function AEP_01d_group_rt_cleaning(config)
    %% Remove Epochs based on Group-Level RT Outliers (Grand Average)
    % 
    % This script aggregates valid Reaction Times (RTs) across all subjects for 
    % each experimental condition (e.g., Cued vs. Voluntary). It then calculates 
    % the Grand Mean and Grand Standard Deviation (SD). Finally, it drops epochs 
    % falling outside the [GrandMean +/- 2SD] window.
    
    
    %% ==================== 0. Configuration ====================
    if nargin < 1 || isempty(config)
        config.root_dir       = '/path/to/your/project/directory';
        config.erp_std_dir    = fullfile(config.root_dir, 'PipelineData', 'erpdata_standard');
        config.erp_final_dir  = fullfile(config.root_dir, 'PipelineData', 'erpdata_noout');
        config.mat_dir        = fullfile(config.root_dir, 'matdata');
        
        config.subs           = 1:20;
        config.sessions       = ["noprediction", "session1", "session2"];
        config.conditions     = ["Cued", "Voluntary"];
        
        fprintf('[INFO] No config passed. Running EEG_01d with local default settings.\n');
    else
        fprintf('[INFO] EEG_01d executing with configurations passed from the main pipeline.\n');
    end

    load(fullfile(config.mat_dir, 'RT_infos.mat'), 'RT_subs');

    % Outlier Threshold Settings (Matches Methods logic)
    config.grand_sd_limit = 2;   % +/- 2 SD from Grand Mean
    config.abs_rt_min     = 0.25; % Minimum absolute RT (s)
    config.abs_rt_max     = 1.00; % Maximum absolute RT (s)

    %% ==================== 1. Calculate Grand Averages ====================
    fprintf('--- Calculating Grand Averages before cleaning across all subjects ---\n');
    RT_stats_grand = struct();
    
    for condition = config.conditions
        rts_group = [];
        
        % Aggregate RTs from all subjects for the current decision type
        for id_subs = config.subs
            for id_sess_str = config.sessions
                id_sess = char(id_sess_str);
                rt_infos = RT_subs(id_subs).(id_sess);
                
                for no = 1:length(rt_infos)
                    des = regexprep(rt_infos(no).Description, '(\<[a-z])', '${upper($1)}').erase(" ");
                    if contains(des, condition) && ~contains(des, "Catch")
                        % Extract valid RTs (not rejected in earlier steps)
                        rt = rt_infos(no).RT;
                        valid_rts = rt(~rt_infos(no).RejSemiManual_shift);
                        rts_group = [rts_group, valid_rts];
                    end
                end
            end
        end
        
        % Compute Grand Mean and SD
        RT_stats_grand.(condition + "_Mean") = mean(rts_group);
        RT_stats_grand.(condition + "_SD")   = std(rts_group);
        
        fprintf('[%s] Grand Mean: %.3f s | Grand SD: %.3f s (Total trials: %d)\n', ...
                condition, mean(rts_group), std(rts_group), length(rts_group));
    end
    
    %% ==================== 2. Reject Epochs based on Grand Bounds ====================
    fprintf('\n--- Applying Grand Bounds and Saving Final Epochs ---\n');
    
    for id_subs = config.subs
        sub_erp_dir = fullfile(config.erp_std_dir, sprintf('s%d', id_subs));
        sub_out_dir = fullfile(config.erp_final_dir, sprintf('s%d', id_subs));
        if ~exist(sub_out_dir, 'dir'); mkdir(sub_out_dir); end

        aggregated_rt = struct('TrialType', {}, 'RT_original', {}, 'RT_final', {},'RejFlag_final', {});
        agg_idx = 1;
        
        for id_sess_str = config.sessions
            id_sess = char(id_sess_str);
            rt_infos = RT_subs(id_subs).(id_sess);
            
            for no = 1:length(rt_infos)
                des = regexprep(rt_infos(no).Description, '(\<[a-z])', '${upper($1)}').erase(" ");
                % side = regexprep(rt_infos(no).Side, '(\<[a-z])', '${upper($1)}');
                
                if ~contains(des, "Catch")
                    setname = rt_infos(no).TrialType;
                    filename = fullfile(sub_erp_dir, [setname, '.set']);
                    
                    if ~exist(filename, 'file'); continue; end
    
                    
                    % Load dataset
                    EEG = pop_loadset('filename', [setname, '.set'], 'filepath', sub_erp_dir);
                    
                    % Identify condition (Cued vs Voluntary)
                    if contains(des, "Cued")
                        condition = "Cued";
                    else
                        condition = "Voluntary";
                    end
                    
                    % Define thresholds
                    g_mean = RT_stats_grand.(condition + "_Mean");
                    g_sd   = RT_stats_grand.(condition + "_SD");
                    
                    % Apply dynamic limits capped by absolute min/max
                    low_limit = max(config.abs_rt_min, g_mean - config.grand_sd_limit * g_sd);
                    up_limit  = min(config.abs_rt_max, g_mean + config.grand_sd_limit * g_sd);
                    
                    % We must apply this filter ONLY to the remaining valid trials
                    % Since EEG.trials matches the already cleaned RTs
                    rt_remaining = rt_infos(no).RT(~rt_infos(no).RejSemiManual_shift);
                    rm_flag_grand = (rt_remaining < low_limit) | (rt_remaining > up_limit);
    
                    RejSemiManual_final = (rt_infos(no).RT < low_limit) | (rt_infos(no).RT > up_limit);
                    
                    % Reject these specific epochs
                    if sum(rm_flag_grand) > 0
                        EEG = pop_rejepoch(EEG, rm_flag_grand, 0);
                        EEG.comments = pop_comments(EEG.comments, '', ...
                            sprintf('12. Removed epochs exceeding Grand Average +/- %dSD in %s', ...
                            config.grand_sd_limit, condition), 1);
                    end
    
                    % Save the finally cleaned dataset (noout)
                    pop_saveset(EEG, 'filename', [setname, '.set'], 'filepath', sub_out_dir);
                    
                    % Save final record of RTs (Optional but good for stats)
                    rt_infos(no).RT_Noout_Across_2SD = rt_remaining(~rm_flag_grand);
                    rt_infos(no).RejSemiManual_final = RejSemiManual_final;

                    aggregated_rt(agg_idx).TrialType = rt_infos(no).TrialType;
                    aggregated_rt(agg_idx).RT_original = rt_infos(no).RT;
                    aggregated_rt(agg_idx).RT_final = rt_infos(no).RT_Noout_Across_2SD;
                    aggregated_rt(agg_idx).RejFlag_final = rt_infos(no).RejSemiManual_final;
                    agg_idx = agg_idx + 1;
                end
            end
            RT_subs(id_subs).(id_sess) = rt_infos;
        end
        RT_subs(id_subs).RT_infos = aggregated_rt;
    end
    
    % Save the final metadata
    save(fullfile(config.mat_dir, '3_RT_infos.mat'), 'RT_subs');
    fprintf('\n--- Step 12 Complete! Final cleaned datasets saved to erpdata_noout ---\n');
end