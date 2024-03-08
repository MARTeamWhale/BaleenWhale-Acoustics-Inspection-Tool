function BAIT_convertLFDCSTable(varargin)
%
%   Convert LFDCS autodetection results into xlsx file compatible with 
%   the BAIT validator program.
%
%   NOTE: this function makes the assumption that the source audio file of 
%   each detection in the input spreadsheet exists. If there are missing
%   files in the middle or end of the recording sequence, results may be
%   incorrect. Files missing at the start are not as bad - in that case, 
%   the affected detections will simply be removed.
%
%
%   SYNTAX
%   -----------------------------------------------------------------------
%   BAIT_convertLFDCSTable
%   BAIT_convertLFDCSTable(Name,Value)
%   -----------------------------------------------------------------------
%
%
%   INPUT ARGUMENTS (optional Name-Value pairs)
%   -----------------------------------------------------------------------
%   params -> name or path of parameter to use. The parameter file for this
%       function defines some filtering variables that decide which
%       detections will be included in the converted spreadsheet.
%
%   input_file -> path of the input LFDCS autodetetions CSV file. If not
%       specified, user is prompted to select file.
%
%   wav_dir -> path of folder containing raw audio files from which the
%       detections originate. If not specified, user will be prompted to
%       choose the folder.
%
%   output_file -> path in which to write the output BAIT-style XLSX 
%       spreadsheet. If not specified, user is prompted to save the file.
%   
%   wav_subfolders -> True/False value specifying whether to search for
%       audio files within any subfolders that may exist in the root audio
%       folder. Default is true.
%
%   read_stop_times -> True/False value specifying whether or not audio
%       file stop times should be read or not. Setting this to true results
%       in more accurate and robust determination of which audio file each
%       detection is contained in, but at the cost of significantly greater
%       processing time. If false, it is assumed that the stop time of an
%       audio file corresponds to the start time of the next audio file, 
%       which is much faster and works in most cases, but may result in 
%       errors if the timeseries is not monotonic or contains missing 
%       files. Default is false.
% -------------------------------------------------------------------------
%
% DEPENDENCIES:
%   BAIT.readLFDCSTable
%   BAIT.getRecTimes
%   MUCA.filepaths.listFiles
%   MUCA.time.readDateTime
%
%
%   Written by Wilfried Beslin
%   Last updated 2024-03-08 using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    import BAIT.readLFDCSTable
    import BAIT.getRecTimes
    import MUCA.filepaths.listFiles
    import MUCA.time.readDateTime

    % 1) INITIALIZATION AND INPUT ARGUMENT PARSING ........................
    disp('Initializing...')
    
    % define common variables
    dtRef = datetime(1970,1,1,0,0,0);
    
    % define paths to resource folders and files
    scriptPath = mfilename('fullpath');
    [rootDir,~,~] = fileparts(scriptPath);
    outputTemplatePath = fullfile(rootDir,'+BAIT','OutputTemplate.xlsx');
    
    % parse input args
    p = inputParser;
    p.addParameter('params','', @ischar)
    p.addParameter('input_file', '', @ischar)
    p.addParameter('wav_dir', '', @ischar)
    p.addParameter('output_file', '', @ischar)
    p.addParameter('wav_subfolders', true, @islogical)
    p.addParameter('read_stop_times', false, @islogical)
    p.parse(varargin{:})

    paramsFileInput = p.Results.params;
    wavDir = p.Results.wav_dir;
    outputFilePath = p.Results.output_file;
    search_wav_subfolders = p.Results.wav_subfolders;
    getRecStopTimes = p.Results.read_stop_times;
    
    % get and validate input file paths
    
    %%% template BAIT output xlsx
    if ~isfile(outputTemplatePath)
        error('Could not find Detection Browser output template file')
    end
    
    %%% LFDCS spreadsheet
    inFilePath = p.Results.input_file;
    if isempty(inFilePath)
        [inFileName,inFileDir] = uigetfile('*.csv','Select LFDCS autodetections CSV file');
        inFilePath = fullfile(inFileDir,inFileName);
        if isnumeric(inFileName)
            return
        end
    end
    
    %%% WAV folder
    if isempty(wavDir)
        wavDir = uigetdir(rootDir,'Select WAV folder');
        if isnumeric(wavDir)
            return
        end
    end
    
    
    % 2) READ FILTERING PARAMETERS ........................................
    disp('Reading parameter file...')
    PARAMS = loadParams(paramsFileInput);
    

    % 3) EXTRACT LFDCS DATA ...............................................
    disp('Extracting LFDCS data...')
    
    % read the LFDCS autodetections table (filtered as needed)
    tableFiltParams = rmfield(PARAMS, setdiff(fieldnames(PARAMS),{'ManualSpeciesCodes','AutoCallTypes','StartDateTime','StopDateTime'}));
    [tableLFDCS, detTimes, ~, precisionLoss] = readLFDCSTable(inFilePath, tableFiltParams);
    
    % Issue a warning if Excel has dropped milliseconds and prompt user for
    % action
    if precisionLoss
        time_prompt_cell = {...
            'It appears that the LFDCS detection times have been rounded, likely because the CSV file was opened and saved in Microsoft Excel. This will result in inaccurate bounding boxes when viewing the detections, and may cause further issues for other people using the data.';...
            'It is STRONGLY RECOMMENDED to use a CSV file that contains the true detection times. If an unaltered backup of the original file is not available, it will have to be recreated using the "export_autodetections" command in LFDCS.';...
            'Are you sure you wish to convert the current (inaccurate) file anyway?'};
        time_prompt = sprintf('%s\n\n%s\n\n%s\n', time_prompt_cell{:});
        time_btn1str = 'Yes, proceed anyway (NOT RECOMMENDED)';
        time_btn2str = 'No, I will use a better file';
        usr_opt = questdlg(time_prompt, 'WARNING: Imprecise Detection Times', time_btn1str, time_btn2str, time_btn2str);
        
        if strcmp(usr_opt,time_btn2str) || isempty(usr_opt)
            return
        end
    end
    
    % get number of filtered detections
    n = height(tableLFDCS);
    
    
    % 4) GET WAV FILE LIST AND RECORDING TIMES ............................
    disp('Getting WAV file times...')
    [wavFilePaths, ~] = listFiles(wavDir, 'wav', 'Recursive',search_wav_subfolders);
    
    % extract datetime from WAV files
    %dtWav = readDateTime(wavFileNames);
    if getRecStopTimes
        [dtWavStart, dtWavStop] = getRecTimes(wavFilePaths);
    else
        dtWavStart = getRecTimes(wavFilePaths);
    end
    
    % sort WAV files based on start time
    [dtWavStart, iSort] = sort(dtWavStart);
    if getRecStopTimes
        dtWavStop = dtWavStop(iSort);
    else
        dtWavStop = [dtWavStart(2:end); datetime(Inf, Inf, Inf)];
    end
    %wavFileNames = wavFileNames(iSort);
    wavFilePaths = wavFilePaths(iSort);
    
    % get WAV file paths relative to root
    wavFileRelPaths = erase(wavFilePaths, [wavDir,filesep]);
    
    
    % 5) ASSIGN WAV FILES TO EACH AUTODETECTION ...........................
    disp('Finding origin WAV file for each detection...')
    
    % get WAV file for each detection
    iDetWav = NaN(n,1);
    for ii = 1:n
        detStartii = detTimes(ii,1);
        %iWavii = find(dtWavStart <= detStartii, 1, 'last');
        iWavii = find(dtWavStart <= detStartii & dtWavStop > detStartii);
        if numel(iWavii) == 1
            iDetWav(ii) = iWavii;
        elseif numel(iWavii) > 1
            warning('Detection %d/%d: multiple audio files occur simultaneously at the time of this detection. Detection will be ignored.', ii, n)
        end
    end
    
    % remove entries that have no or overlapping WAV files
    good_files = ~isnan(iDetWav);
    n_good = sum(good_files);
    
    
    % 6) CREATE OUTPUT ....................................................
    disp('Writing output...')
    
    % create output table
    outTableHeader = {'FileName','FileStart','SigStart','SigEnd','SigStartDateTime','Class_LFDCS','Class_MATLAB','ReasonForUNK','Comments'};
    %FileName = wavFileNames(iDetWav(good_files));
    FileName = wavFileRelPaths(iDetWav(good_files)); % relative paths are needed if files are spread across subfolders
    FileStart = seconds(dtWavStart(iDetWav(good_files)) - dtRef);
    SigStart = seconds(detTimes(good_files,1) - dtRef) - FileStart;
    SigEnd = seconds(detTimes(good_files,2) - dtRef) - FileStart;
    SigStartDateTime = detTimes(good_files,1);
    SigStartDateTime.Format = 'dd-MMM-yyyy HH:mm:ss';
    Class_LFDCS = repmat({''},n_good,1);
    Class_LFDCS(tableLFDCS.ManualSpeciesCode(good_files)==9999) = {'Correct'};
    Class_LFDCS(tableLFDCS.ManualSpeciesCode(good_files)==0) = {'Unknown'};
    Class_LFDCS(tableLFDCS.ManualSpeciesCode(good_files)==-9999) = {'Incorrect'};
    Class_MATLAB = Class_LFDCS;
    ReasonForUNK = repmat({''},n_good,1);
    Comments = repmat({''},n_good,1);
    
    tableOut = table(...
        FileName,FileStart,SigStart,SigEnd,SigStartDateTime,Class_LFDCS,Class_MATLAB,ReasonForUNK,Comments,...
        'VariableNames',outTableHeader);
    
    % initialize output file
    if isempty(outputFilePath)
        [outputFileName,outputFileDir] = uiputfile('*.xlsx','Save output file');
        outputFilePath = fullfile(outputFileDir,outputFileName);
    end
    [copyOK,copyMsg] = copyfile(outputTemplatePath,outputFilePath);
    if ~copyOK
        error('Could not create output file:\n%s',copyMsg)
    end
    
    % update output table
    writetable(tableOut,outputFilePath,'Sheet','Detected');
    
    disp('Done')
end


% loadParams --------------------------------------------------------------
function PARAMS = loadParams(paramFileInput)
% Reads in program parameters from file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    import BAIT.processParamFile
    import BAIT.readParam

    % read parameter file as a block of text
    scriptPath = mfilename('fullpath');
    [rootDir,scriptName,~] = fileparts(scriptPath);
    paramsText = processParamFile(paramFileInput, rootDir, scriptName);

    % initialize output
    PARAMS = struct;
    
    % set parameters
    PARAMS.ManualSpeciesCodes = readParam(paramsText, 'ManualSpeciesCodes', {@(var)validateattributes(var,{'numeric'},{'integer'}), @(var)assert(isnan(var))});
    PARAMS.AutoCallTypes = readParam(paramsText, 'AutoCallTypes', {@(var)validateattributes(var,{'numeric'},{'integer'}), @(var)assert(isnan(var))});
    PARAMS.StartDateTime = readParam(paramsText, 'StartDateTime', {@(var)validateattributes(var,{'numeric'},{'numel',6}), @(var)assert(isnan(var))});
    PARAMS.StopDateTime = readParam(paramsText, 'StopDateTime', {@(var)validateattributes(var,{'numeric'},{'numel',6}), @(var)assert(isnan(var))});
    
    % change time parameters to datetime Infs if they are NaNs
    if isnan(PARAMS.StartDateTime)
        PARAMS.StartDateTime = datetime('-Inf');
    else
        PARAMS.StartDateTime = datetime(PARAMS.StartDateTime);
    end
    if isnan(PARAMS.StopDateTime)
        PARAMS.StopDateTime = datetime('Inf');
    else
        PARAMS.StopDateTime = datetime(PARAMS.StopDateTime);
    end
end