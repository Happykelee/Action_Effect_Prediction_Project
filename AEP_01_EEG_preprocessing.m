%% EEG Preprocessing and CSD Transformation Pipeline
% =========================================================================
% Platforms and Plugins Used:
% - MATLAB Version:   R2023b
% - EEGLAB Version:   v2024.0
% - ICLabel Plugin:   v1.6
% - PICARD Plugin:    v1.0
% - CSD Toolbox:      v1.1 (Kayser & Tenke, 2006)
%
% Terminology Mapping (Code vs. Manuscript):
% Please note that the raw data files and internal script variables utilize 
% original laboratory naming conventions. These perfectly correspond to the 
% formal terminology reported in the manuscript as follows:
%
%   Code / Filename         Manuscript Terminology
%   -------------------     ----------------------------------
%   "subs","subjects" ==>   "Participants"
%   "noprediction"    ==>   "Pre-association phase"
%   "session1"        ==>   "Post-association Session 1"
%   "session2"        ==>   "Post-association Session 2"
%
% Pipeline Summary:
% This script executes the continuous data preprocessing, ICA-based artifact 
% rejection, epoching, baseline correction, and surface Laplacian (CSD) 
% transformation prior to subsequent univariate and multivariate analyses.
% 
% Pipeline steps:
% 1. Import raw data (.set)
% 2. Bandpass filtered (1 - 40 Hz)
% 3. Re-referenced (TP9, TP10)
% 4. Useless and noised channels rejection
% 5. Edges at start/end trimmed
% 6. Merged datasets (Day 1) & ran ICA
% 7. Rejected artifact components via ICLabel & manual check
% 8. Interpolated missing channels 
% 9. Epoched based on RT
% 10. Rejected bad epochs and subject-level RT outliers
% 11. Shifted baseline and realigned to stimulus onset
% 12. Removed epochs based on group-Level RT outliers
% 13. Current Source Density (CSD) transformation
% =========================================================================

clear; clc; close all;
eeglab; % Initialize EEGLAB

%% ==================== 0. Configuration ====================
% Define all paths and global parameters here to avoid hardcoding below.
% Users reproducing the code should only need to modify this section.

config.root_dir       = '/path/to/your/project/directory';
config.raw_dir        = fullfile(config.root_dir, 'rawdata_set');
config.pre_dir        = fullfile(config.root_dir, 'PipelineData', 'predata');
config.erp_dir        = fullfile(config.root_dir, 'PipelineData', 'erpdata');
config.erp_std_dir    = fullfile(config.root_dir, 'PipelineData', 'erpdata_standard');
config.erp_final_dir  = fullfile(config.root_dir, 'PipelineData', 'erpdata_noout');
config.csd_dir        = fullfile(config.root_dir, 'PipelineData', 'csddata');
config.csd_decode_dir = fullfile(config.root_dir, 'PipelineData', 'csddata_decode');
config.code_dir       = fullfile(config.root_dir, 'code');
config.mat_dir        = fullfile(config.root_dir, 'matdata');
config.stat_dir       = fullfile(config.root_dir, 'statdata');
% config.pre_dir        = '/home/zhengting/BackupData/ParisCite/1_MultivariateProcedure/PipelineData_3rd/prodata';
% config.mat_dir        = '/home/zhengting/Insync/OneDrive/ParisCite/1_MultivariateProcedure/code/Github_pub/matdata';
% config.stat_dir        = '/home/zhengting/Insync/OneDrive/ParisCite/1_MultivariateProcedure/code/Github_pub/statdata';

config.subs           = 1:20; % Subject list
config.sessions       = ["noprediction", "session1", "session2"];
config.conditions     = ["Cued", "Voluntary"];
config.types          = ["single_left_visual", "single_left_auditory", ...
                         "double_left_visual", "double_left_auditory"];
                     
config.filter_band    = [1, 40]; % Hz
config.ref_channels   = {'TP9', 'TP10'};
config.drop_channels  = {'FT9', 'PO9', 'PO10'}; 

addpath(config.code_dir);
if ~exist(config.pre_dir, 'dir'); mkdir(config.pre_dir); end
if ~exist(config.erp_dir, 'dir'); mkdir(config.erp_dir); end
if ~exist(config.erp_std_dir, 'dir'); mkdir(config.erp_std_dir); end
if ~exist(config.erp_final_dir, 'dir'); mkdir(config.erp_final_dir); end
if ~exist(config.csd_dir, 'dir'); mkdir(config.csd_dir); end
if ~exist(config.csd_decode_dir, 'dir'); mkdir(config.csd_decode_dir); end

% --- Artifact Rejection Thresholds (Matches the Methods Section) ---
config.ep_thresh  = [-800, 800]; % Extreme voltage thresholds (µV)
config.ep_trend   = [500, 0.3];  % Trend fitting: [max slope, max R^2]
config.ep_prob    = [5, 8];      % Joint probability: [local SD, global SD limit]
config.ep_kurt    = [5, 8];      % Kurtosis: [local SD, global SD limit]

%% ==================== 1-5. Basic Preprocessing ====================
% Filtered, re-referenced, channels rejection, and edges trimmed.


% Initialize the structure for behavioral accuracy
Cued_Accuracy = struct();

for id_subs = config.subs
    sub_out_dir = fullfile(config.pre_dir, sprintf('s%d', id_subs));
    if ~exist(sub_out_dir, 'dir'); mkdir(sub_out_dir); end

    Cued_Accuracy(id_subs).SID = id_subs;
    
    for id_sess = config.sessions
        fprintf('\n--- Processing Sub %02d | Sess %s ---\n', id_subs, id_sess);
        
        % -----------------------------------------------------
        % 1. Import raw data (.set)
        % -----------------------------------------------------
        filename_raw = sprintf('s%d%s_rawdata.set', id_subs, id_sess);
        EEG = pop_loadset('filename', filename_raw, 'filepath', config.raw_dir);
        
        % -----------------------------------------------------
        % 2. Bandpass filtered (1 - 40 Hz)
        % -----------------------------------------------------
        EEG = pop_eegfiltnew(EEG, 'locutoff', config.filter_band(1), 'hicutoff', config.filter_band(2), 'plotfreqz', 0);
        
        % -----------------------------------------------------
        % 3. Re-referenced (TP9, TP10)
        % -----------------------------------------------------
        EEG = pop_reref(EEG, config.ref_channels);
        
        % -----------------------------------------------------
        % 4. Useless and noised channels rejection
        % -----------------------------------------------------
        % First, automatically remove predefined, globally non-analytical channels.
        % Second, remove subject-specific noised channels identified via manual GUI inspection.
        
        % 4a. Automated exclusion of predefined useless channels
        EEG = pop_select(EEG, 'nochannel', config.drop_channels);
        if ~isfield(config, 'chanlocs')
            config.chanlocs = EEG.chanlocs; 
            config.chanlocs_labels = {EEG_ref.chanlocs.labels}';
            fprintf('  [INFO] Standard channel montage successfully captured into config.\n');
        end
        
        % 4b. Manual intervention for noised channels:
        % Beyond the standard channel drops, channels heavily contaminated by 
        % persistent high-frequency noise, technical drifts, or poor impedance 
        % were meticulously identified via manual visual inspection using the 
        % EEGLAB GUI and subsequently excluded on a subject-by-subject basis.
        
        % (Optional: If reproducing from scratch interactively, uncomment below 
        % to pop up the channel scrolling GUI for manual channel rejection)
        % pop_eegplot(EEG, 1, 1, 1); 
        % keyboard; % Pause script execution to allow manual editing via GUI
        
        % -----------------------------------------------------
        % 5. Edges at start/end trimmed
        % -----------------------------------------------------
        start_time = max(EEG.event(2).latency/EEG.srate - 5, EEG.xmin);
        end_time   = min(EEG.event(end).latency/EEG.srate + 5, EEG.xmax);
        EEG = pop_select(EEG, 'time', [start_time, end_time]);

        % -----------------------------------------------------
        % Calculate Cued Condition Behavioral Accuracy
        % -----------------------------------------------------
        types_EEG = string({EEG.event.type});
        % Count the number of incorrect (S222) and correct (S210) keypress trials
        num_err = sum(types_EEG == "S222");
        num_corr = sum(types_EEG == "S210");
        if (num_err + num_corr) > 0
            acc = 1 - num_err / (num_err + num_corr);
        else
            acc = NaN; % Prevent division by zero if relevant markers are missing
        end
        % Store the calculated accuracy into the structure based on the current session
        if id_sess == "noprediction"
            Cued_Accuracy(id_subs).PreAssn = acc;
        elseif id_sess == "session1"
            Cued_Accuracy(id_subs).PostAssnSess1 = acc;
        elseif id_sess == "session2"
            Cued_Accuracy(id_subs).PostAssnSess2 = acc;
        end

        % Save preprocessed data
        setname = sprintf('s%d%s_1_prep', id_subs, id_sess);
        EEG.setname = setname;
        EEG.comments = pop_comments(EEG.comments, '', ...
            '1-5. Bandpass filtered (1-40Hz), re-referenced, channels rejection, and edges trimmed.', 1);
        pop_saveset(EEG, 'filename', [setname, '.set'], 'filepath', sub_out_dir);

    end
end

acc_table = struct2table(Cued_Accuracy);
csv_file_path = fullfile(config.stat_dir, 'Cued_Accuracy.csv');
writetable(acc_table, csv_file_path);
fprintf('  [INFO] Cued Accuracy successfully calculated and exported to:\n -> -> %s\n', csv_file_path);

%% ==================== 6-7. Merge & ICA ====================
% Merged sessions of day 1 (Session 1 & NoPrediction), ran ICA, 
% and rejected artifact components via ICLabel & manual check


for id_subs = config.subs
    sub_out_dir = fullfile(config.pre_dir, sprintf('s%d', id_subs));
    
    for day = 1:2
        fprintf('\n--- ICA Sub %02d | Day %d ---\n', id_subs, day);
        
        % -----------------------------------------------------
        % 6. Merged datasets (Day 1) & ran ICA
        % -----------------------------------------------------
        if day == 1
            EEG1 = pop_loadset('filename', sprintf('s%d%s_1_prep.set', id_subs, config.sessions(1)), 'filepath', sub_out_dir);
            EEG2 = pop_loadset('filename', sprintf('s%d%s_1_prep.set', id_subs, config.sessions(2)), 'filepath', sub_out_dir);
            OUTEEG = pop_mergeset(EEG1, EEG2, 0);
            setname = sprintf('s%dsession1_2_merge_ica', id_subs);
        else
            OUTEEG = pop_loadset('filename', sprintf('s%d%s_1_prep.set', id_subs, config.sessions(3)), 'filepath', sub_out_dir);
            setname = sprintf('s%dsession2_2_ica', id_subs);
        end
        
        % Run ICA (Picard)
        OUTEEG = pop_runica(OUTEEG, 'icatype', 'picard', 'maxiter', 500, 'mode', 'standard');
        
        % -----------------------------------------------------
        % 7. Rejected artifact components via ICLabel & manual check
        % -----------------------------------------------------
        % Flag artifact components via ICLabel with a conservative probability threshold (>= 80%)
        % Target artifact categories: Muscle (Row 2), Ocular (Row 3), and Channel Noise (Row 6)
        OUTEEG = pop_iclabel(OUTEEG, 'default');
        OUTEEG = pop_icflag(OUTEEG, [NaN NaN; 0.8 1; 0.8 1; NaN NaN; NaN NaN; 0.8 1; NaN NaN]);
        
        OUTEEG.setname = setname;
        pop_saveset(OUTEEG, 'filename', [setname, '.set'], 'filepath', sub_out_dir);

        % Manual intervention & final cleaning:
        % In addition to the automated ICLabel flagging (for the most obvious 
        % muscle, ocular, and channel noises), remaining artifact components 
        % were manually inspected and rejected using the EEGLAB GUI.
        
        % (Optional: Uncomment the next 3 lines if you want the script to pause 
        % here, allowing you to use the GUI before it automatically saves)
        % pop_selectcomps(OUTEEG, 1:size(OUTEEG.icaweights,1));
        % disp('Please manually reject components in the GUI, then type "return" in the command window.');
        % keyboard; 
        
        % Remove the flagged/selected components
        OUTEEG = pop_subcomp(OUTEEG, [], 0); 
        
        % Save the final cleaned continuous dataset
        if day == 1
            setname = sprintf('s%dsession1_2_merge_ica_cleaned', id_subs);
        else
            setname = sprintf('s%dsession2_2_ica_cleaned', id_subs);
        end
        OUTEEG.setname = setname;
        OUTEEG.comments = pop_comments(OUTEEG.comments, '', ...
            '6-7. Merged sessions, ran ICA (Picard), and rejected artifact components via ICLabel & manual check.', 1);
        pop_saveset(OUTEEG, 'filename', [setname, '.set'], 'filepath', sub_out_dir);
    end
end

%% ==================== 8-11. Epoching & Baseline Correction ====================
% Interpolated channels, epoched, rejected bad epochs and subject-level RT outliers,
% shifted baseline and realigned to stimulus onset.


% Ensure RT_subs mapping is generated and loaded prior to this step
rt_file_path = fullfile(config.mat_dir, 'RT_infos.mat');

if ~exist(rt_file_path, 'file')
    warning('Behavioral RT data not found. Automatically executing AEP_01b_Extract_RT.m to generate it...');
    AEP_01b_extract_behavioral_rt(config); 
    fprintf('Behavioral data extraction completed. Resuming main pipeline...\n');
end
load(rt_file_path, 'RT_subs');

for id_subs = config.subs
    sub_out_dir = fullfile(config.pre_dir, sprintf('s%d', id_subs));
    sub_erp_dir = fullfile(config.erp_dir, sprintf('s%d', id_subs));
    if ~exist(sub_erp_dir, 'dir'); mkdir(sub_erp_dir); end
    
    for day = 1:2
        % -----------------------------------------------------
        % 8. Interpolated missing channels
        % -----------------------------------------------------
        filename = dir(fullfile(sub_out_dir, sprintf('s%dsession%d*ica_cleaned.set', id_subs, day)));
        EEG_all = pop_loadset('filename', filename(1).name, 'filepath', sub_out_dir);
        
        if EEG_all.nbchan < 58
            EEG_all = pop_interp(EEG_all, config.chanlocs, 'spherical');
        end

        if day == 1
            sess_idx = 1:2;
        else
            sess_idx = 3;
        end
        
        for id_sess = config.sessions(sess_idx)
            rt_infos = RT_subs(id_subs).(id_sess);
            
            for no = 1:length(rt_infos)
                des = regexprep(rt_infos(no).Description, '(\<[a-z])', '${upper($1)}').erase(" ");
                side = regexprep(rt_infos(no).Side, '(\<[a-z])', '${upper($1)}');
                mark = rt_infos(no).Marker;
                
                if ~contains(des, "Catch")
                    % -----------------------------------------------------
                    % 9. Epoched based on RT
                    % -----------------------------------------------------
                    rt = rt_infos(no).RT;
                    rt_mean = mean(rt(~rt_infos(no).RejSemiManual));
                    rt_std = std(rt(~rt_infos(no).RejSemiManual));
                    
                    rt_thrs_max = min(rt_mean + 2*rt_std, 1.00);
                    rt_thrs_min = max(rt_mean - 2*rt_std, 0.25);
                    rejflag = (rt > rt_thrs_max) | (rt < rt_thrs_min) | rt_infos(no).RejSemiManual;
                    
                    % Define dynamic epoch window based on the maximum valid reaction time.
                    % [NOTE]: May require manual adjustment for atypical samples.
                    epoch_limit = [-(0.81 + mean(rt(~rejflag)) + std(rt(~rejflag))*2), 1.4];

                    % Extract initial large epochs
                    setname = sprintf('s%d_%s_%s_%s', id_subs, id_sess, des, side);
                    EEG = pop_epoch(EEG_all, {mark}, epoch_limit, 'newname', setname, 'epochinfo', 'yes');
                    
                    % -----------------------------------------------------
                    % 10. Rejected bad epochs and subject-level RT outliers
                    % -----------------------------------------------------
                    % Run the 4 mathematical artifact detection algorithms on epochs
                    EEG = pop_eegthresh(EEG, 1, 1:EEG.nbchan, config.ep_thresh(1), config.ep_thresh(2), EEG.xmin, EEG.xmax, 0, 0);
                    EEG = pop_rejtrend(EEG, 1, 1:EEG.nbchan, EEG.pnts, config.ep_trend(1), config.ep_trend(2), 0, 0);
                    EEG = pop_jointprob(EEG, 1, 1:EEG.nbchan, config.ep_prob(1), config.ep_prob(2), 0, 0);
                    EEG = pop_rejkurt(EEG, 1, 1:EEG.nbchan, config.ep_kurt(1), config.ep_kurt(2), 0, 0, 0, [], 0);

                    % Extract automated EEG rejection indices
                    eeg_auto_bad_epochs = [find(EEG.reject.rejthresh == 1), ...
                                           find(EEG.reject.rejconst == 1), ...
                                           find(EEG.reject.rejjp == 1), ...
                                           find(EEG.reject.rejkurt == 1)];
                    
                    % Combine Behavioral (RT) Outliers and EEG Artifact Outliers
                    final_rejflag = rejflag;
                    if ~isempty(eeg_auto_bad_epochs)
                        final_rejflag(eeg_auto_bad_epochs) = true;
                    end
                    
                    % Reject all bad epochs simultaneously
                    if sum(final_rejflag) > 0
                        EEG = pop_rejepoch(EEG, final_rejflag, 0);
                    end

                    % -----------------------------------------------------
                    % 11. Shifted baseline and realigned to stimulus onset
                    % -----------------------------------------------------
                    % Step 11a: Shift zero time point to "Offset of tactile" (S 40) for baseline correction
                    epoch_limit_shift = [-1.0, 1.02];
                    setname_shift = sprintf('%s_shift', EEG.setname);
                    EEG = pop_epoch(EEG, {'S 40'}, epoch_limit_shift, 'newname', setname_shift, 'epochinfo', 'yes');
                    EEG = pop_rmbase(EEG, [-800, -300], []);

                    % Step 11b: Shift back (realign) to "Onset of stimulus" as the actual zero time point
                    epoch_limit_realign = [-0.95, 0.502];
                    EEG = pop_epoch(EEG, {mark}, epoch_limit_realign, 'newname', setname, 'epochinfo', 'yes');

                    % Defensive check: Ensure no trials were lost during the realignment
                    assert(EEG.trials == sum(~final_rejflag), 'Trial count mismatch after realignment!');
                    rt_infos(no).RejSemiManual_shift = final_rejflag;
                    EEG.comments = pop_comments(EEG.comments, '', ...
                        ['8-11. Interpolated channels, epoched, rejected bad epochs and subject-level RT outliers, ' ...
                        'shifted baseline and realigned to stimulus onset.'], 1);
                    pop_saveset(EEG, 'filename', [setname, '.set'], 'filepath', sub_erp_dir);
                end
            end
            
            % Update back into behavioral structure
            RT_subs(id_subs).(id_sess) = rt_infos;
        end
    end
end

% Save updated RT structure containing synchronized rejection flags
save(rt_file_path, 'RT_subs');


%% ==================== 12 Remove epochs based on group-Level RT outliers ====================
% Renames and reorganizes the epoched prior to this step
AEP_01c_standardize_filenames(config);

% -----------------------------------------------------
% 12. Removed epochs based on group-Level RT outliers
% -----------------------------------------------------
AEP_01d_group_rt_cleaning(config);

%% ==================== 13. CSD Transformation (Dual-Pipeline) ====================
% Current Source Density (CSD) transformation using surface Laplacian
% This step splits into two isolated pipelines to prepare data for different endpoints:
%   - Pipeline A: Standard full-resolution CSD for univariate analysis.
%   - Pipeline B: Downsampled CSD designed to optimize multivariate decoding.


disp('Applying Current Source Density (CSD) Transformation...');

% Prepare Montage (Requires CSD Toolbox)
Montage_CSD = ExtractMontage('10-5-System_Mastoids_EGI129.csd', config.chanlocs_labels);
[TransMat_G, TransMat_H] = GetGH(Montage_CSD, 4); 

for id_subs = config.subs
    % Define subject-specific subdirectories for perfect data tracking
    sub_erp_dir   = fullfile(config.erp_final_dir, sprintf('s%d', id_subs));
    sub_csd_std   = fullfile(config.csd_dir, sprintf('s%d', id_subs));
    sub_csd_dec   = fullfile(config.csd_decode_dir, sprintf('s%d', id_subs));
    
    if ~exist(sub_csd_std, 'dir'); mkdir(sub_csd_std); end
    if ~exist(sub_csd_dec, 'dir'); mkdir(sub_csd_dec); end
    
    % Search for all finalized clean epoch files for this subject
    erp_files = dir(fullfile(sub_erp_dir, '*.set'));
    if isempty(erp_files); continue; end
    
    fprintf('\n--- Running Dual-CSD Channels for Subject %02d (%d files found) ---\n', id_subs, length(erp_files));
    
    for i = 1:length(erp_files)
        current_filename = erp_files(i).name;
        
        % =========================================================================
        % PIPELINE A: Standard Full-Resolution CSD
        % =========================================================================
        EEG_std = pop_loadset('filename', current_filename, 'filepath', sub_erp_dir);
        
        % Apply CSD Transformation matrix
        data_CSD_std = CSD(EEG_std.data, TransMat_G, TransMat_H);
        EEG_std.data = data_CSD_std;
        
        % Update Dataset metadata
        EEG_std.setname = [EEG_std.setname, '_CSD'];
        EEG_std.comments = pop_comments(EEG_std.comments, '', ...
            '13a. Standard Pipeline: Surface Laplacian CSD transformation applied to full-resolution clean epochs.', 1);
        
        % Save to standard CSD folder
        pop_saveset(EEG_std, 'filename', [EEG_std.setname, '.set'], 'filepath', sub_csd_std);
        clear EEG_std;
        
        % =========================================================================
        % PIPELINE B: Downsampled CSD for Multivariate Decoding
        % =========================================================================
        EEG_dec = pop_loadset('filename', current_filename, 'filepath', sub_erp_dir);
        
        % 1. Downsample first to reduce computational cost and feature dimensionality
        % (Adjust 100 Hz to your preferred target decoding rate if necessary)
        EEG_dec = pop_resample(EEG_dec, 100); 
        
        % 2. Apply CSD Transformation matrix to the downsampled data
        data_CSD_dec = CSD(EEG_dec.data, TransMat_G, TransMat_H);
        EEG_dec.data = data_CSD_dec;
        
        % Update Dataset metadata with strict tracking info
        EEG_dec.setname = [EEG_dec.setname, '_Dsampled_CSD'];
        EEG_dec.comments = pop_comments(EEG_dec.comments, '', ...
            '13b. Decoding Pipeline: Epochs downsampled to 100 Hz prior to Surface Laplacian CSD transformation.', 1);
        
        % Save to separate decoding CSD folder
        pop_saveset(EEG_dec, 'filename', [EEG_dec.setname, '.set'], 'filepath', sub_csd_dec);
        clear EEG_dec;
    end
end

disp('================ Pipeline Complete ================');