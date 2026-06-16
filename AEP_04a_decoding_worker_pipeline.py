"""
Multivariate Decoding - Parallel Modeling Worker Pipeline
========================================================
Normalization, Feature Selection and Classification by ElasticNet, Model Testing.

Dependencies: func_timeout, mne, numpy, pandas, scikit-learn
"""

import mne
import time
import func_timeout
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.linear_model import LogisticRegression

def fit_logistic_regression_timeout(X_train, y_train, kw_model):
    """Wrapper function to fit model for timeout enforcement."""
    clf = LogisticRegression(**kw_model)
    clf.fit(X_train, y_train)
    return clf

def run_subject_decoding_pipeline(params):
    """
    Worker function executed in parallel for a single subject, condition, and session.
    Paths are now dynamically passed via params to ensure modularity.
    """
    # Unpack parameters, including the newly added directory paths
    # decis, stim, day, sid, random_ints, cv_folder_str, mat_folder_str, data_csd_folder_str = params
    decis, stim, sess, sid, random_ints, cv_folder_str, mat_folder_str, data_csd_folder_str = params
    decision = 'Cued' if decis.startswith('Cue') else 'Voluntary'
    stimulus = 'Visual' if stim.startswith('Vis') else 'Auditory'
    session = 'Session1' if sess == 'Sess1' else 'Session2'
    # trialinfo = f'{decis}_{day}_No{stim}PreAction'
    trialinfo = f'{decis}_{sess}_No{stim}PreAction'
    
    # Reconstruct Path objects (strings are safer for multiprocessing serialization)
    cv_folder = Path(cv_folder_str)
    mat_folder = Path(mat_folder_str)
    data_csd_folder = Path(data_csd_folder_str)
    sid_folder = mat_folder / str(sid)

    # Ensure the output directory exists
    if not sid_folder.exists(): sid_folder.mkdir(parents=True, exist_ok=True)
    log_file = sid_folder / f'NoStimPreAction_S{sid}_log.txt'
    
    # Define constant parameters based on protocol requirements
    n_splits = 8
    n_repeat = 20
    n_times_train = 96  # Up to id_time_zero + 1
    n_times_general = n_times_train
    
    # Modeling hyperparameters
    c_default = 0.125
    penalty = 'elasticnet'
    solver = 'saga'
    l1_ratio = 0.5
    tol = 0.001
    max_iter = 1000
    time_limit = 10  # Seconds allowed per time-point execution
    
    t_start_total = time.time()
    
    with open(log_file, "a+") as log:
        log.write(f"--- Starting Processing for S{sid} | {decision} {stimulus} {session} ---\n")
        
        # Load specific CSD epochs data here for the participant
        # fpath_Pre   = data_csd_folder / f'{sid}/s{sid}_Day1_{decision}No{stimulus}StimPre_Dsampled_CSD.set'
        # fpaths_Post = [data_csd_folder / f'{sid}/s{sid}_{day}_{decision}No{stimulus}StimPost_Dsampled_CSD.set',
        #                data_csd_folder / f'{sid}/s{sid}_{day}_{decision}{stimulus}_Dsampled_CSD.set']
        fpath_Pre   = data_csd_folder / f'{sid}/s{sid}_Session1_{decision}No{stimulus}StimPre_CSD_Decode.set'
        fpaths_Post = [data_csd_folder / f'{sid}/s{sid}_{session}_{decision}No{stimulus}StimPost_CSD_Decode.set',
                       data_csd_folder / f'{sid}/s{sid}_{session}_{decision}{stimulus}_CSD_Decode.set']
        data_Pre = mne.io.read_epochs_eeglab(fpath_Pre,verbose='CRITICAL').get_data(copy=True)
        data_Post = np.concatenate([mne.io.read_epochs_eeglab(p, verbose='CRITICAL').get_data(copy=True) for p in fpaths_Post])

        # 为了预分配存储矩阵，先读取第一轮的第一个折叠获取样本总数 len_cv 和通道数 n_chs
        # To pre-allocate the storage matrix, first read the first fold of the first round to obtain len_cv and n_chs
        temp_fname = cv_folder / f'{n_splits}Folds_{trialinfo}_S{sid}_train_Round00_Seed{random_ints[0]}.xlsx'
        if temp_fname.exists():
            len_cv = len(pd.read_excel(temp_fname, sheet_name='Fold1')) + \
                     len(pd.read_excel(str(temp_fname).replace('_train_', '_test_'), sheet_name='Fold1'))
            assert(len_cv == data_Pre.shape[0] * 2)
        else:
            len_cv = data_Pre.shape[0] * 2

        n_chs = data_Pre.shape[1]

        # Initialize storage matrix
        ch_coefs_matrix_cv = np.zeros((n_repeat, n_splits, n_times_train, n_chs))
        prob_matrix_cv     = np.zeros((n_repeat, len_cv, n_times_train, n_times_general))
        pred_matrix_cv     = np.zeros((n_repeat, len_cv, n_times_train, n_times_general))
        
        for id_rand, rand_state in enumerate(random_ints):
            t_round_start = time.time()
            
            # Access dynamic cv_folder instead of global
            fname_train = cv_folder / f'{n_splits}Folds_{trialinfo}_S{sid}_train_Round{id_rand:02d}_Seed{rand_state}.xlsx'
            fname_test  = cv_folder / f'{n_splits}Folds_{trialinfo}_S{sid}_test_Round{id_rand:02d}_Seed{rand_state}.xlsx'
            
            if not (fname_train.exists() and fname_test.exists()):
                log.write(f"Round {id_rand+1} skipped: CV Excel split files missing.\n")
                continue
                
            for kidx in range(n_splits):
                df_train_fold = pd.read_excel(fname_train, sheet_name=f'Fold{kidx+1}')
                df_test_fold = pd.read_excel(fname_test, sheet_name=f'Fold{kidx+1}')
                df_label_fold = pd.concat((df_train_fold,df_test_fold)).sort_values(by=['idx'])
                label_Post = np.array([int(i[-3:]) for i in df_label_fold.loc[df_label_fold['GroundTruth']==1,'TID']])
                X_all = np.concatenate([data_Pre, data_Post[label_Post,:,:]])

                # Z-score Normalization & Modeling Block
                X_train_mean = X_all[df_train_fold.idx,:,:n_times_train].mean()
                X_train_sd   = X_all[df_train_fold.idx,:,:n_times_train].std()
                X_train_normalize = (X_all[df_train_fold.idx,:,:n_times_train]-X_train_mean)/X_train_sd
                X_test_normalize  = (X_all[df_test_fold.idx,:,:n_times_train]-X_train_mean)/X_train_sd

                for tidx_train in range(n_times_train):
                    kw_enet = dict(
                        C=c_default, penalty=penalty, solver=solver,
                        max_iter=max_iter, l1_ratio=l1_ratio, tol=tol
                    )
                    
                    run_flag = True
                    while run_flag:
                        try:
                            # Placeholder logic for training
                            X_dummy_train = X_train_normalize[...,tidx_train] 
                            y_dummy_train = df_train_fold['GroundTruth'].values
                            
                            clf_model = func_timeout.func_timeout(
                                time_limit, fit_logistic_regression_timeout,
                                args=(X_dummy_train, y_dummy_train, kw_enet)
                            )
                            
                            flag_ch_select = clf_model.coef_.flatten() != 0
                            if np.sum(flag_ch_select) > 0:
                                ch_coefs_matrix_cv[id_rand, kidx, tidx_train, :] = clf_model.coef_.flatten()
                                run_flag = False
                            else:
                                kw_enet['C'] += 0.05
                                
                        except func_timeout.FunctionTimedOut:
                            kw_enet['tol'] += 0.001
                            if kw_enet['tol'] > 0.1: 
                                run_flag = False

                    # Predicting of all generalization time points
                    for tidx_general in range(n_times_general):
                        X_test_gen = X_test_normalize[..., tidx_general]
                        idx_test_global = df_test_fold.idx.values # 对应的全局 trial 索引
                        prob_matrix_cv[id_rand, idx_test_global, tidx_train, tidx_general] = clf_model.predict_proba(X_test_gen)[:, 1]
                        pred_matrix_cv[id_rand, idx_test_global, tidx_train, tidx_general] = clf_model.predict(X_test_gen)

                    # avoid using the existed
                    del clf_model 
            
            t_round_end = time.time()
            log.write(f"Round {id_rand+1}/{n_repeat} complete. Duration: {t_round_end - t_round_start:.2f}s\n")
            log.flush()
            
        t_end_total = time.time()
        log.write(f"--- Processing Finished for S{sid}. Total Time: {t_end_total - t_start_total:.2f}s ---\n\n")

        # Store all generated matrices to separate files
        savefile_template = f'{n_repeat}Rounds_{n_splits}Folds_ElasticNet_{trialinfo}_S{sid}'
        np.savez_compressed(sid_folder / f'{savefile_template}_ch_coefs_matrix_cv.npz', ch_coefs_matrix_cv)
        np.savez_compressed(sid_folder / f'{savefile_template}_prob_matrix_cv.npz', prob_matrix_cv)
        np.savez_compressed(sid_folder / f'{savefile_template}_pred_matrix_cv.npz', pred_matrix_cv)
        
    return f"Subject {sid} | {decision} {stimulus} complete."