function extract_LFDCS_detections(varargin)
%
%   Saves clips and/or spectrograms of LFDCS detections.
%
%   SYNTAX
%   -----------------------------------------------------------------------
%   extract_LFDCS_detections()
%   extract_LFDCS_detections(Name,Value)
%   -----------------------------------------------------------------------
%
%
%   INPUT ARGUMENTS (optional Name-Value pairs)
%   -----------------------------------------------------------------------
%   "params" - name or path of parameter file to use
%   .......................................................................
%   "input_file" - string specifying path to input file, which should be
%       CSV file exported from LFDCS. If not specified, user will be
%       prompted to select the file.
%   .......................................................................
%   "output_dir" - string specifying path to the folder where clips and/or
%       spectrograms should be saved. If not specified, user will be
%       prompted to choose a folder.
%   .......................................................................
%   "audio_dir" - string specifying path to the folder containing audio
%       files within which calls were detected in LFDCS. If not specified,
%       user will be prompted to choose a folder.
%   .......................................................................
%   "overwrite" - True/false value specifying whether or not existing clips 
%       or spectrograms in the output folder should be overwritten. Default
%       is false.
%   -----------------------------------------------------------------------
%
%
%   DEPENDENCIES
%       MUCA.io.saveFig
%       MUCA.time.readDateTime
%
%
%   NOTES
%   -----------------------------------------------------------------------
%   - If the "SaveClips" parameter is false but "SaveSpecs" is true, the
%   function will search the output folder for existing clips. If clips
%   exist, they will be used to create the spectrograms (in this case no
%   audio folder is needed). If no clips are found, then an input audio
%   folder is required.
%
%   - Upcalls may not always appear where they should in the clips or
%   spectrogram images. There can be two reasons for this, which can occur
%   at the same time:
%       -- Reason 1) The source CSV file may have been saved in Microsoft
%       Excel, which often causes the milliseconds in the detection times
%       to be stripped away, reducing their precision. This problem can be
%       resolved by recreating the CSV file in LFDCS and avoiding saving in
%       Excel. A warning is issued if this problem is detected.
%       -- Reason 2) LFDCS detection times do not always match up with the
%       times in which the detections actually occur in the WAV files. The
%       reason for this is not entirely known - it may be a bug in LFDCS.
%       There is currently no easy solution to this problem.
%   -----------------------------------------------------------------------
%
%   Written by Wilfried Beslin
%   Last updated 2023-12-01 using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    import MUCA.io.saveFig

    % 1) INPUT PARSING ....................................................

    % define paths to resource folders
    rootDir = mfilename('fullpath');
    [rootDir,scriptName,~] = fileparts(rootDir);
    %resDir = fullfile(rootDir,'BrowserResources');
    paramsDir = fullfile(rootDir,'PARAMS',scriptName);

    % parse input
    p = inputParser;
    
    p.addParameter('params', 'default_params.txt', @ischar)
    p.addParameter('input_file', '', @ischar)
    p.addParameter('output_dir', '', @ischar)
    p.addParameter('audio_dir', '', @ischar)
    p.addParameter('overwrite', false)
    
    p.parse(varargin{:})
    paramsFilePath = p.Results.params;
    input_file_path = p.Results.input_file;
    usr_output_dir = p.Results.output_dir;
    audio_dir = p.Results.audio_dir;
    overwrite = p.Results.overwrite;
    
    % import parameters
    [userParamsDir, paramsFilename, paramsExt] = fileparts(paramsFilePath);
    if isempty(paramsExt)
        paramsExt = '.txt';
    end
    if isempty(userParamsDir)
        paramsFilePath = fullfile(paramsDir,[paramsFilename,paramsExt]);
    end
    PARAMS = loadParams(paramsFilePath);
    
    
    % 2) INITIALIZE BASED ON PARAMETERS ...................................
    
    % process file and folder paths based on clip and spectrogram options
    specs_only = PARAMS.SaveSpecs && ~PARAMS.SaveClips;
    specs_from_clips = false;
    
    if ~PARAMS.SaveClips && ~PARAMS.SaveSpecs
        disp('Nothing to do!')
        return
        
    elseif specs_only
        % process output folder first and check if clips exist
        disp('Spectrogram processing only; looking for existing clip files to produce spectrograms')
        [output_folders, existing_clips, existing_specs] = process_output_folders(usr_output_dir, PARAMS.SaveClips, PARAMS.SaveSpecs);
        if isempty(output_folders)
            disp('Cancelling')
            return
        end
        
        if ~isempty(existing_clips)
            specs_from_clips = true;
            data = existing_clips;
            num_detections = numel(data);
        else
            disp('No clips found; will process LFDCS detection output')
        end
    end
    
    if ~specs_from_clips
        % ask user for input file if not specified
        if isempty(input_file_path)
            [input_file_name,input_file_dir] = uigetfile({'*.csv';'*.xlsx'},'Select LFDCS exported results file');
            input_file_path = fullfile(input_file_dir,input_file_name);
            if isnumeric(input_file_name)
                disp('Cancelling')
                return
            end
        end
        
        % ask user for audio folder path if not specified
        if isempty(audio_dir)
            audio_dir = uigetdir(pwd,'Specify audio file folder');
            if isnumeric(audio_dir)
                disp('Cancelling')
                return
            end
        end
        
        % process output folders if not already done
        if ~specs_only
            [output_folders, existing_clips, existing_specs] = process_output_folders(usr_output_dir, PARAMS.SaveClips, PARAMS.SaveSpecs);
            if isempty(output_folders)
                disp('Cancelling')
                return
            end
        end
        
        % read LFDCS file
        disp('Processing LFDCS detections...')
        [data, deployment] = read_LFDCS_file(input_file_path, audio_dir, PARAMS.RecursiveSearch);
        if isempty(data)
            disp('Cancelling')
            return
        end
        num_detections = height(data);
    end
    
    
    % 3) PROCESSING .......................................................
    
    % loop through each detection and save a clip and/or spectrogram
    previous_rec = '';
    error_logs = table((1:num_detections)', repmat({''},num_detections,1), 'VariableNames',{'Detection','Error'});
    for ii = 1:num_detections
        
        try
            if specs_from_clips
                clip_file_path = fullfile(output_folders.clips, existing_clips{ii});
                [~, out_name, ~] = fileparts(clip_file_path);
                spec_file_path = fullfile(output_folders.specs, [out_name,'.png']);

                fprintf('Detection %d/%d (%s)\n', ii, num_detections, out_name);

                [x_det, fs] = audioread(clip_file_path);

            else
                current_rec = data.FilePath{ii};
                current_rec_start = data.FileStart(ii);
                det_call_start = data.DetTime(ii);
                det_call_end = det_call_start + data.DetDur(ii);

                % determine detection start date/time
                det_dt = current_rec_start + seconds(det_call_start);

                % set output file names and paths
                out_name = sprintf('%s_Detection_%s',deployment,char(det_dt,'yyyyMMdd_HHmmss'));
                clip_file_path = fullfile(output_folders.clips, [out_name,'.wav']);
                spec_file_path = fullfile(output_folders.specs, [out_name,'.png']);

                fprintf('Detection %d/%d (%s)\n', ii, num_detections, out_name);

                % process new recording file if necessary
                if ~strcmp(current_rec, previous_rec)
                    fprintf('        Isolating file %s\n', current_rec)
                    current_rec_file = data.FilePath{ii};
                    rec_info = audioinfo(current_rec_file);
                    previous_rec = current_rec;
                end

                % isolate detection
                initial_det_range = [det_call_start, det_call_end];
                [x_det, fs] = isolateDetection(current_rec_file, rec_info, PARAMS.Channel, initial_det_range, PARAMS.SnippetDur);

                % save clip if specified
                if PARAMS.SaveClips
                    if ~(isfile(clip_file_path) && ~overwrite)
                        disp('    Saving clip')
                        audiowrite(clip_file_path, x_det, fs, 'BitsPerSample',rec_info.BitsPerSample);
                    else
                        disp('    Clip already exists')
                    end
                end

            end

            % process spectrogram if specified
            if PARAMS.SaveSpecs
                if ~(isfile(spec_file_path) && ~overwrite)
                    disp('    Saving spectrogram')
                    fig = make_spectrogram(x_det, fs, PARAMS.SpecColorMap, PARAMS.SpecMaxFreq, out_name);
                    saveFig(fig, spec_file_path, PARAMS.SpecFigSize, 'pixels');
                    close(fig)
                else
                    disp('    Spectrogram already exists')
                end
            end
        catch ME
            warning('Failed to process detection %d/%d:\n %s', ii, num_detections, ME.message)
            error_logs.Error{ii} = ME.getReport('extended', 'hyperlinks','off');
        end
        
    end
    
    % check if there are any bad detections and save an error report if so
    has_error = ~cellfun('isempty', error_logs.Error);
    if any(has_error)
        disp('Bad detections encountered, saving report...')
        dt_now = datetime('now');
        dt_now.Format = 'yyyyMMddHHmmss';
        error_file_name = ['bad_detections_',char(dt_now),'.xlsx'];
        error_file_path = fullfile(output_folders.root, error_file_name);
        writetable(error_logs(has_error,:), error_file_path)
    end
    disp('Done')
end


% loadParams --------------------------------------------------------------
function PARAMS = loadParams(paramFile)
% Reads in program parameters from file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    import BWAV_code.readParam
    import BWAV_code.buildColormaps

    % read parameter file as a block of text
    params_text = fileread(paramFile);

    % initialize output
    PARAMS = struct;
    
    % set parameters
    PARAMS.Channel = readParam(params_text, 'Channel', {@(var)validateattributes(var,{'numeric'},{'scalar','positive','integer'})});
    PARAMS.RecursiveSearch = readParam(params_text, 'RecursiveSearch', {@(var)validateattributes(var,{'logical'},{'scalar'})});
    PARAMS.SaveClips = readParam(params_text, 'SaveClips', {@(var)validateattributes(var,{'logical'},{'scalar'})});
    PARAMS.SaveSpecs = readParam(params_text, 'SaveSpecs', {@(var)validateattributes(var,{'logical'},{'scalar'})});
    PARAMS.SnippetDur = readParam(params_text, 'SnippetDur', {@(var)validateattributes(var,{'numeric'},{'scalar','positive'})});
    PARAMS.SpecMaxFreq = readParam(params_text, 'SpecMaxFreq', {@(var)validateattributes(var,{'numeric'},{'scalar','positive'})});
    PARAMS.SpecColorMap = readParam(params_text, 'SpecColorMap', {@(var)validateattributes(var,{'char'},{'row'})});
    PARAMS.SpecFigSize = readParam(params_text, 'SpecFigSize', {@(var)validateattributes(var,{'numeric'},{'numel',2,'integer','positive'})});
    
    % assign colormap matrix
    cmaps = buildColormaps();
    try
        PARAMS.SpecColorMap = feval(cmaps.(PARAMS.SpecColorMap));
    catch
        error('Invalid colormap "%s"', PARAMS.SpecColorMap)
    end
end


% read_LFDCS_file ---------------------------------------------------------
function [data, deployment] = read_LFDCS_file(LFDCS_file_path, audio_dir, recursive_search)
% Extract relevant info from CSV file exported by LFDCS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    import MUCA.time.readDateTime

    dt_ref = datetime(1970,1,1,0,0,0);
    LFDCS_header_rows = 23;
    start_time_col = 2;
    duration_col = 4;
    deployment_expr = '[A-Z]{3,4}_\d{4}_\d{2}';
    
    % get table
    import_opts = detectImportOptions(LFDCS_file_path,'NumHeaderLines',LFDCS_header_rows);
    LFDCS_table = readtable(LFDCS_file_path,import_opts);
    num_detections = height(LFDCS_table);
    
    % check detection times to see if Excel has dropped the milliseconds. 
    % Issue a warning if so.
    if mean(LFDCS_table.Var2 - round(LFDCS_table.Var2) == 0) > 0.8
        time_prompt_cell = {...
            'It appears that the LFDCS detection times have been rounded, likely because the CSV file was opened and saved in Microsoft Excel. This will cause the position of the detections within the clips or spectrograms to be slightly inaccurate.';...
            'It is STRONGLY RECOMMENDED to use a CSV file that contains the true detection times. If an unaltered backup of the original file is not available, it will have to be recreated using the "export_autodetections" command in LFDCS.';...
            'Are you sure you wish to extract and/or plot detections using the current (inaccurate) file anyway?'};
        time_prompt = sprintf('%s\n\n%s\n\n%s\n', time_prompt_cell{:});
        time_btn1str = 'Yes, proceed anyway (NOT RECOMMENDED)';
        time_btn2str = 'No, I will use a better file';
        usr_opt = questdlg(time_prompt, 'WARNING: Imprecise Detection Times', time_btn1str, time_btn2str, time_btn2str);
        
        if strcmp(usr_opt,time_btn2str) || isempty(usr_opt)
            data = [];
            deployment = [];
            return
        end
    end
    
    % get WAV file for each detection
    disp('Getting WAV file times...')
    [rec_file_names,rec_file_paths] = Utilities.getFileNames(audio_dir, 'wav', recursive_search);
    rec_times = readDateTime(rec_file_names);
    %rec_times_secs = seconds(rec_times - dt_ref);
    
    %%% sort files in case they are out of order
    [rec_times_sorted, idx_sort] = sort(rec_times);
    rec_file_paths_sorted = rec_file_paths(idx_sort);
    
    det_times = dt_ref + seconds(LFDCS_table{:, start_time_col});
    det_rec_file_idx = NaN(num_detections,1);
    for ii = 1:num_detections
        det_time_ii = det_times(ii);
        rec_idx_ii = find(rec_times_sorted <= det_time_ii,1,'last');
        det_rec_file_idx(ii) = rec_idx_ii;
    end
    
    % reorganize and keep useful data
    data_table_headers = {'FilePath', 'FileStart', 'DetTime', 'DetDur'};
    data_FilePath = rec_file_paths_sorted(det_rec_file_idx);
    data_FileStart = rec_times_sorted(det_rec_file_idx);
    data_DetTime = seconds(det_times - rec_times_sorted(det_rec_file_idx));
    data_DetDur = LFDCS_table{:, duration_col};
    
    data = table(data_FilePath, data_FileStart, data_DetTime, data_DetDur, 'VariableNames',data_table_headers);
    
    % determine deployment
    [~, LFDCS_file_name, ~] = fileparts(LFDCS_file_path);
    deployment = regexp(LFDCS_file_name, deployment_expr, 'match');
    deployment = deployment{:};
end


% process_output_folders --------------------------------------------------
function [output_folders, existing_clips, existing_specs] = process_output_folders(usr_output_dir, do_clips, do_specs)
% Process user input relating to output to set and/or create output folders
% and query existing files.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ask user for output folder path if not specified
    if isempty(usr_output_dir)
        output_root_dir = uigetdir(pwd, 'Specify root output folder');
        if isnumeric(output_root_dir)
            output_folders = [];
            existing_clips = [];
            existing_specs = [];
            return
        end
    else
        output_root_dir = usr_output_dir;
    end
    
    % if folders exist, query them; otherwise, create them
    
    %%% clips
    output_clip_dir = fullfile(output_root_dir, 'clips');
    existing_clips = [];
    if isfolder(output_clip_dir)
        existing_clips = Utilities.getFileNames(output_clip_dir, 'wav');
    elseif do_clips
        mkdir(output_clip_dir);
    end
    
    %%% spectrograms
    output_spec_dir = fullfile(output_root_dir, 'spectrograms');
    existing_specs = [];
    if isfolder(output_spec_dir)
        existing_specs = Utilities.getFileNames(output_spec_dir, 'png');
    elseif do_specs
        mkdir(output_spec_dir);
    end

    % save folder paths to struct
    output_folders = struct(...
        'root',output_root_dir,...
        'clips',output_clip_dir,...
        'specs',output_spec_dir);
end


% isolateDetection --------------------------------------------------------
function [x, fs] = isolateDetection(rec_file, rec_info, channel, detection_call_range, snippet_dur)
% Isolate call using centre of detection window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fs = rec_info.SampleRate;
    N = rec_info.TotalSamples;
    
    % determine centre of call
    t_centre = diff(detection_call_range)/2 + detection_call_range(1);
    
    % determine final call range
    final_call_range = [t_centre - snippet_dur/2, t_centre + snippet_dur/2];
    
    % convert range to samples
    final_call_range_samples = round(final_call_range.*fs);
    
    % correct range in case it is out of bounds
    final_call_range_samples_adjusted = [max([1,final_call_range_samples(1)]), min([N,final_call_range_samples(2)])];
    if final_call_range_samples_adjusted(1) > final_call_range_samples_adjusted(2)
        error('Adjusted sample range was [%d, %d].\nNumber of samples in audio file is %d.\nThis detection may occur outside of its assigned audio file.\nSource file path: %s', final_call_range_samples_adjusted(1), final_call_range_samples_adjusted(2), N, rec_file)
    end
    start_pad_length = final_call_range_samples_adjusted(1) - final_call_range_samples(1);
    end_pad_length = final_call_range_samples(2) - final_call_range_samples_adjusted(2);
    
    % read WAV file
    x = audioread(rec_file, final_call_range_samples_adjusted);
    
    % adjust snippet with zero padding if needed
    x = [zeros(start_pad_length,1); x(:,channel); zeros(end_pad_length,1)];
end


% make_spectrogram --------------------------------------------------------
function fig = make_spectrogram(x, fs, colmap, spec_f_max, name)
% Does just that. Settings are based on the ones used by MERIDIAN.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    specparam_win = round(fs*0.256);
    specparam_ovl = specparam_win - round(fs*0.032);
    specparam_nfft = 2^nextpow2(fs*(2048/8000));
    
    % calculate spectrogram
    [~, f, t, P] = spectrogram(x, specparam_win, specparam_ovl, specparam_nfft, fs);
    
    % figure
    fig = figure('Visible', 'off');
    ax = axes(fig);
    surf(ax, t, f, 10*log10(P), 'EdgeColor','none');
    view(ax, 2);
    ax.YLim = [0, spec_f_max];
    colormap(ax, colmap)
    
    title(ax, name, 'Interpreter','none');
    xlabel(ax, 'Time [s]')
    ylabel(ax, 'Frequency [Hz]')
    
end