function AEP_01c_standardize_filenames(config)
    %% Standardize Epoch Filenames for Group Analysis
    % This script renames and reorganizes the epoched datasets from their 
    % original, condition-specific raw names (e.g., s1_session1_CuedAuditory_Right) 
    % to a standardized format (e.g., s1_Session1_CuedAuditory) suitable for 
    % Group-level averaging and MVPA decoding.
    %
    % Key Mapping Rules:
    % 1. 'noprediction' -> Mapped into 'Session1' with a 'Pre' suffix for NoStim.
    % 2. 'session1/2' NoStim stimulus type -> Appended with 'Post' suffix.
    % 3. 'session1/2' -> Capitalized to 'Session1/Session2'.
    % 4. Left/Right side indicators are dropped.
    
    
    %% ==================== 0. Configuration ====================
    if nargin < 1 || isempty(config)
        config.root_dir    = '/path/to/your/project/directory';
        config.erp_dir     = fullfile(config.root_dir, 'PipelineData', 'erpdata');
        config.erp_std_dir = fullfile(config.root_dir, 'PipelineData', 'erpdata_standard');
        config.mat_dir     = fullfile(config.root_dir, 'matdata');
        
        config.subs        = 1:20;
        config.sessions    = ["noprediction", "session1", "session2"];
        fprintf('[INFO] No config passed. Running EEG_01c with local default settings.\n');
    else
        fprintf('[INFO] EEG_01c executing with configurations passed from the main pipeline.\n');
    end
    
    load(fullfile(config.mat_dir, 'RT_infos.mat'), 'RT_subs');
    
    % Paradigm mapping for Noprediction (Pre) based on Subject ID
    preside_types(1).Left = 'NoVisualStimPre';   preside_types(1).Right = 'NoAuditoryStimPre';
    preside_types(2).Left = 'NoAuditoryStimPre'; preside_types(2).Right = 'NoVisualStimPre';
    preside_types(3).Left = 'NoVisualStimPre';   preside_types(3).Right = 'NoAuditoryStimPre';
    preside_types(4).Left = 'NoAuditoryStimPre'; preside_types(4).Right = 'NoVisualStimPre';
    
    %% ==================== 1. Standardize and Copy Files ====================
    fprintf('--- Standardizing Filenames & Updating RT_infos ---\n');
    
    for id_subs = config.subs
        id_type = mod(id_subs - 1, 4) + 1;
        
        sub_src_dir = fullfile(config.erp_dir, sprintf('s%d', id_subs));
        sub_out_dir = fullfile(config.erp_std_dir, sprintf('s%d', id_subs));
        if ~exist(sub_out_dir, 'dir'); mkdir(sub_out_dir); end
        
        for id_sess_str = config.sessions
            id_sess = char(id_sess_str);
            rt_infos = RT_subs(id_subs).(id_sess);
            
            for no = 1:length(rt_infos)
                des = regexprep(rt_infos(no).Description, '(\<[a-z])', '${upper($1)}').erase(" ");
                side = regexprep(rt_infos(no).Side, '(\<[a-z])', '${upper($1)}');
                
                if ~contains(des, "Catch")
                    % 1. Identify the original raw filename (from Step 11)
                    original_setname = sprintf('s%d_%s_%s_%s', id_subs, id_sess, des, side);
                    original_file = fullfile(sub_src_dir, [original_setname, '.set']);
                    
                    % 2. Map to the new standardized TrialType (matching Image 2)
                    % Capitalize Session
                    if strcmp(id_sess, 'session2')
                        % target_day = 'Session2';
                        target_day = 'Day2';
                    else
                        % target_day = 'Session1'; % noprediction is treated as Session1
                        target_day = 'Day1';
                    end
                    
                    % Base Trial Type
                    if strcmp(id_sess, 'noprediction')
                        % Remove 'Noprediction' redundancy
                        base_des = replace(des, 'Noprediction', '');
                        
                        % Determine if it's NoVis or NoAud based on Subject Type and Side
                        if strcmp(side, 'Left')
                            stim_type = preside_types(id_type).Left;
                        else
                            stim_type = preside_types(id_type).Right;
                        end
                        % Append 'Pre' for noprediction phase
                        standard_des = replace(base_des, side, stim_type);
                        
                    else
                        % For normal sessions, append 'Post' to NoStim conditions
                        if contains(des, 'No')
                            standard_des = strcat(des, 'Post');
                        else
                            standard_des = des;
                        end
                    end
                    
                    % Final Standardized Name
                    standard_setname = sprintf('s%d_%s_%s', id_subs, target_day, standard_des);
                    
                    % 3. Save the mapping to RT_infos
                    rt_infos(no).OriginalFilename = original_setname;
                    rt_infos(no).TrialType = standard_setname;
                    
                    % 4. Load, rename, and save the dataset to the new directory
                    if exist(original_file, 'file')
                        EEG = pop_loadset('filename', [original_setname, '.set'], 'filepath', sub_src_dir);
                        EEG.setname = standard_setname;
                        pop_saveset(EEG, 'filename', [standard_setname, '.set'], 'filepath', sub_out_dir);
                    end
                end
            end
            RT_subs(id_subs).(id_sess) = rt_infos;
        end
    end
    
    % Save the updated metadata with mapping traces
    save(fullfile(config.mat_dir, 'RT_infos.mat'), 'RT_subs');
    fprintf('--- Filename Standardization Complete! ---\n');
end