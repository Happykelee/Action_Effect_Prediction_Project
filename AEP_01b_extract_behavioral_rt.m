function AEP_01b_extract_behavioral_rt(config)
    %% Extract Behavioral Data (Reaction Times) and Identify Outliers
    % This script extracts reaction times (RTs) from the preprocessed continuous 
    % EEG datasets (Step 1 outputs), matches them with experimental trial definitions,
    % and calculates behavioral outliers (e.g., > 1.0s or < 0.25s) to be used 
    % later for epoch rejection.
    %
    % Output: A structured .mat file containing RTs and outlier flags for all subjects.
    
    
    %% ==================== 0. Configuration ====================
    if nargin < 1 || isempty(config)
        config.root_dir       = '/path/to/your/project/directory';
        config.pre_dir        = fullfile(config.root_dir, 'PipelineData', 'predata');
        config.mat_dir        = fullfile(config.root_dir, 'matdata');
        
        config.subs           = 1:20; 
        config.sessions       = ["noprediction", "session1", "session2"];
        fprintf('[INFO] No config passed. Running EEG_01b with local default settings.\n');
    else
        fprintf('[INFO] EEG_01b executing with configurations passed from the main pipeline.\n');
    end

    load(fullfile(config.mat_dir, 'markers_des.mat'));
    
    % Outlier Threshold Settings (Matches Methods logic)
    config.abs_rt_min     = 0.25; % Minimum absolute RT (s)
    config.abs_rt_max     = 1.00; % Maximum absolute RT (s)
    
    % Define the types array as used in your original logic
    paradigm_types = ["single_left_visual", "single_left_auditory",...
                      "double_left_visual", "double_left_auditory"];
    RT_subs = struct();
    
    %% ==================== 1. Extract RTs from EEG Events ====================
    disp('Extracting Reaction Times from EEG Triggers...');
    
    for id_subs = config.subs
        RT_subs(id_subs).ID = id_subs;
        no_type = mod(id_subs - 1, 4) + 1;
        sub_pro_dir = fullfile(config.pre_dir, sprintf('s%d', id_subs));
        
        for id_sess = config.sessions
            fprintf('Processing RTs - Sub %02d | Sess %s\n', id_subs, id_sess);
            
            % Construct the exact variable name dynamically (e.g., 'struct_single_left_visual')
            target_struct_name = sprintf('struct_%s', paradigm_types(no_type));
            
            % Extract the full struct array for this paradigm
            current_paradigm_struct = eval(target_struct_name);
            
            % Slice the struct array based on the session type
            if id_sess == "noprediction"
                Marker_des = current_paradigm_struct(1:4);
            else
                Marker_des = current_paradigm_struct(5:end);
            end
            
            % Load the early preprocessed dataset (before ICA/Epoching shifted events)
            filename = sprintf('s%d%s_1_prep.set', id_subs, id_sess);
            EEG = pop_loadset('filename', filename, 'filepath', sub_pro_dir);
            
            % Extract event attributes
            EEG_types = string({EEG.event.type});
            EEG_latent = [EEG.event.latency];
            
            % Match markers and calculate RT (Stimulus to Response)
            for no_mark = 1:length(Marker_des)
                mark = Marker_des(no_mark).Marker{5};
                index = find(EEG_types == mark);
                
                if ~isempty(index)
                    % RT calculation based on original logic: (Response - Stimulus)/srate
                    Marker_des(no_mark).RT = (EEG_latent(index-1) - EEG_latent(index-3)) / EEG.srate;
                    Marker_des(no_mark).PressLatency = EEG_latent(index-1) / EEG.srate;
                else
                    Marker_des(no_mark).RT = NaN;
                    Marker_des(no_mark).PressLatency = NaN;
                end
                
                % Initialize the manual rejection flag to 0
                Marker_des(no_mark).RejSemiManual = false;
            end
            
            RT_subs(id_subs).(id_sess) = Marker_des;
        end
    end
    
    %% ==================== 2. Global Outlier Identification (Optional but Recommended) ====================
    % To flag extreme RTs (e.g., > 1.0s or < 0.25s).
    
    disp('Flagging Behavioral Outliers...');
    
    for id_subs = config.subs
        for id_sess = config.sessions
            rt_infos = RT_subs(id_subs).(id_sess);
            
            for no = 1:length(rt_infos)
                des = regexprep(rt_infos(no).Description, '(\<[a-z])', '${upper($1)}').erase(" ");
                if ~contains(des, "Catch") && all(~isnan(rt_infos(no).RT))
                    rt = rt_infos(no).RT;
                    
                    % Standardize outlier definition (e.g., > 1.0s or < 0.25s)
                    % You can adjust this to your final paper's exact criteria.
                    is_outlier = (rt > config.abs_rt_max) | (rt < config.abs_rt_min); 
                    
                    rt_infos(no).RejSemiManual = is_outlier;
                end
            end
            RT_subs(id_subs).(id_sess) = rt_infos;
        end
    end
    
    %% ==================== 3. Save Behavioral Data ====================
    out_file = fullfile(config.mat_dir, 'RT_infos.mat');
    save(out_file, 'RT_subs');
    fprintf('\nSuccessfully saved RT structures to: %s\n', out_file);
end