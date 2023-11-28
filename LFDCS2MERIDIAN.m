function LFDCS2MERIDIAN(varargin)
%
%   Convert LFDCS autodetection results into xlsx file compatible with 
%   MERIDIAN detection validator program.
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
%   LFDCS2MERIDIAN
%   LFDCS2MERIDIAN(Name,Value)
%   -----------------------------------------------------------------------
%
%
%   INPUT ARGUMENTS (optional Name-Value pairs)
%   -----------------------------------------------------------------------
%   input_file -> path of the input LFDCS autodetetions spreadsheet. May be 
%       CSV or XLSX. If not specified, user is prompted to select file.
%
%   wav_dir -> path of folder containing raw audio files from which the
%       detections originate. If not specified, user will be prompted to
%       choose the folder.
%
%   output_file -> path in which to write the output MERIDIAN-style XLSX 
%       spreadsheet. If not specified, user is prompted to save the file.
%
%   sp_codes -> vector of numeric LFDCS validated species codes to keep in 
%       the converted spreadsheet. These should be limited to the 
%       following codes:
%           9999    (correct)
%           -9999   (incorrect)
%           0       (unknown)
%           -32767  (unclassified)
%       The converted spreadsheet will only retain detections with a manual 
%       validation code that corresponds to one of the entries in the list. 
%       For example, set this to [9999, 0] to only keep correct and unknown
%       calls. The converted spreadsheet will also translate these numeric 
%       codes into strings. If this list is empty, then all detections will
%       be retained (default behaviour)
%   
%   wav_subfolders -> True/False value specifying whether to search for
%       audio files within any subfolders that may exist in the root audio
%       folder. Default is true.
% -------------------------------------------------------------------------
%
%   Written by Wilfried Beslin
%   Last updated Nov 28, 2023 using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES
% 2022-06-17
% ----------
% - WARNING! Currently, this script may fail if there are subfolders that
% contain earlier recordings, as the recordings may be out of order

    import BWAV_code.utilities.getFileNames
    import BWAV_code.utilities.readDateTime

    % 1) INITIALIZATION AND INPUT ARGUMENT PARSING ........................
    disp('Initializing...')
    
    % define common variables
    dtRef = datetime(1970,1,1,0,0,0);
    
    % define paths to resource folders
    rootDir = mfilename('fullpath');
    [rootDir,~,~] = fileparts(rootDir);
    resDir = fullfile(rootDir,'BrowserResources');
    paramsDir = fullfile(rootDir,'PARAMS','LFDCS2MERIDIAN');
    
    % parse input args
    p = inputParser;
    p.addParameter('params','default_params.txt', @(v)ischar(v))
    p.addParameter('input_file', '', @ischar)
    p.addParameter('wav_dir', '', @ischar)
    p.addParameter('output_file', '', @ischar)
    p.addParameter('wav_subfolders', true, @islogical)
    p.parse(varargin{:})

    paramsFilePath = p.Results.params;
    wavDir = p.Results.wav_dir;
    outFilePath = p.Results.output_file;
    search_wav_subfolders = p.Results.wav_subfolders;
    
    % get and validate input file paths
    
    %%% parameter file
    [userParamsDir, paramsFilename, paramsExt] = fileparts(paramsFilePath);
    if isempty(paramsExt)
        paramsExt = '.txt';
    end
    if isempty(userParamsDir)
        paramsFilePath = fullfile(paramsDir,[paramsFilename,paramsExt]);
    end
    
    %%% template MERIDIAN output xlsx
    outTemplatePath = fullfile(resDir,'OutputTemplate.xlsx');
    if ~isfile(outTemplatePath)
        error('Could not find Detection Browser output template file')
    end
    
    %%% LFDCS spreadsheet
    inFilePath = p.Results.input_file;
    if isempty(inFilePath)
        [inFileName,inFileDir] = uigetfile({'*.csv';'*.xlsx'},'Select LFDCS spreadsheet');
        inFilePath = fullfile(inFileDir,inFileName);
        if isnumeric(inFileName)
            return
        end
    end
    
    %%% WAV folder
    if isempty(wavDir)
        wavDirRoot = '\\142.2.83.111\teamwhalenas2\MOORED_PAM_DATA';
        wavDir = uigetdir(wavDirRoot,'Select WAV folder');
        if isnumeric(wavDir)
            return
        end
    end
    
    
    % 2) READ FILTERING PARAMETERS ........................................
    disp('Reading parameter file...')
    PARAMS = loadParams(paramsFilePath);
    

    % 3) EXTRACT LFDCS DATA ...............................................
    disp('Extracting LFDCS data...')
    
    % set number of LFDCS header lines
    LFDCSHeader = 23;
    
    % set columns containing manual species codes and auto call type codes
    iColSpecies = 10; %15 for commented XLSX files...
    iColCallType = 1;
    
    % read LFDCS autodetections file as table
    importOpts = detectImportOptions(inFilePath, 'NumHeaderLines',LFDCSHeader, 'DatetimeType','text');
    tableLFDCS = readtable(inFilePath, importOpts);
    
    % truncate table to include only the species and call type codes of
    % interest
    %%% manual species codes
    if ~isnan(PARAMS.ManualSpeciesCodes)
        rowsInclude = ismember(tableLFDCS{:,iColSpecies},PARAMS.ManualSpeciesCodes);
        tableLFDCS = tableLFDCS(rowsInclude,:);
    end
    %%% auto call type codes
    if ~isnan(PARAMS.AutoCallTypes)
        rowsInclude = ismember(tableLFDCS{:,iColCallType},PARAMS.AutoCallTypes);
        tableLFDCS = tableLFDCS(rowsInclude,:);
    end
    
    % get absolute detection times based on time type
    iColStartTime = 2;
    if isnumeric(tableLFDCS{:,iColStartTime})
        % process for case where start-end times are provided relative
        % to 1970/01/01
        iColEndTime = 3;

        det_start_absolute = dtRef + seconds(tableLFDCS{:,iColStartTime});
        det_end_absolute = dtRef + seconds(tableLFDCS{:,iColEndTime});
    else
        % process for case where start times are provided using
        % absolute datetimes and fractional seconds
        iColStartFracSec = 3;
        iColDuration = 4;

        det_start_rounded = datetime(tableLFDCS{:,iColStartTime}, 'InputFormat','MM/dd/yy HH:mm:ss', 'PivotYear',year(datetime('now'))-99);
        det_start_absolute = det_start_rounded + seconds(tableLFDCS{:,iColStartFracSec});
        det_end_absolute = det_start_absolute + seconds(tableLFDCS{:,iColDuration});
    end
    
    % truncate table to include only detections within time period of
    % interest
    rowsInclude = det_start_absolute >= PARAMS.StartDateTime & det_end_absolute <= PARAMS.StopDateTime;
    det_start_absolute = det_start_absolute(rowsInclude);
    det_end_absolute = det_end_absolute(rowsInclude);
    tableLFDCS = tableLFDCS(rowsInclude,:);
    
    n = height(tableLFDCS);
    
    % check detection times to see if Excel has dropped the milliseconds. 
    % Issue a warning if so.
    if mean(second(det_start_absolute) - round(second(det_start_absolute)) == 0) > 0.8
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
    
    
    % 4) GET WAV FILE LIST AND RECORDING TIMES ............................
    disp('Getting WAV file times...')
    [wavFileNames,~] = getFileNames(wavDir, 'wav', search_wav_subfolders);
    
    % extract datetime from WAV files
    dtWav = readDateTime(wavFileNames);
    
    
    % 5) ASSIGN WAV FILES TO EACH AUTODETECTION ...........................
    disp('Finding origin WAV file for each detection...')
    
    % get WAV file for each detection
    iDetWav = NaN(n,1);
    for ii = 1:n
        detStartii = det_start_absolute(ii);
        iWavii = find(dtWav <= detStartii, 1, 'last');
        if ~isempty(iWavii)
            iDetWav(ii) = iWavii;
        end
    end
    
    % remove entries that have no WAV files
    good_files = ~isnan(iDetWav);
    n_good = sum(good_files);
    
    
    % 6) CREATE OUTPUT ....................................................
    disp('Writing output...')
    
    % create output table
    outTableHeader = {'FileName','FileStart','SigStart','SigEnd','SigStartDateTime','Class_LFDCS','Class_MATLAB','ReasonForUNK','Comments'};
    FileName = wavFileNames(iDetWav(good_files));
    FileStart = seconds(dtWav(iDetWav(good_files)) - dtRef);
    SigStart = seconds(det_start_absolute(good_files) - dtRef) - FileStart;
    SigEnd = seconds(det_end_absolute(good_files) - dtRef) - FileStart;
    SigStartDateTime = det_start_absolute(good_files);
    SigStartDateTime.Format = 'dd-MMM-yyyy HH:mm:ss';
    Class_LFDCS = repmat({''},n_good,1);
    Class_LFDCS(tableLFDCS{good_files,iColSpecies}==9999) = {'Correct'};
    Class_LFDCS(tableLFDCS{good_files,iColSpecies}==0) = {'Unknown'};
    Class_LFDCS(tableLFDCS{good_files,iColSpecies}==-9999) = {'Incorrect'};
    Class_MATLAB = Class_LFDCS;
    ReasonForUNK = repmat({''},n_good,1);
    Comments = repmat({''},n_good,1);
    
    tableOut = table(...
        FileName,FileStart,SigStart,SigEnd,SigStartDateTime,Class_LFDCS,Class_MATLAB,ReasonForUNK,Comments,...
        'VariableNames',outTableHeader);
    
    % initialize output file
    if isempty(outFilePath)
        [outFileName,outFileDir] = uiputfile('*.xlsx','Save output file');
        outFilePath = fullfile(outFileDir,outFileName);
    end
    [copyOK,copyMsg] = copyfile(outTemplatePath,outFilePath);
    if ~copyOK
        error('Could not create output file:\n%s',copyMsg)
    end
    
    % update output table
    writetable(tableOut,outFilePath,'Sheet','Detected');
    
    disp('Done')
end


% loadParams --------------------------------------------------------------
function PARAMS = loadParams(paramFile)
% Reads in program parameters from file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    import BWAV_code.readParam

    % read parameter file as a block of text
    params_text = fileread(paramFile);

    % initialize output
    PARAMS = struct;
    
    % set parameters
    PARAMS.ManualSpeciesCodes = readParam(params_text, 'ManualSpeciesCodes', {@(var)validateattributes(var,{'numeric'},{'integer'}), @(var)assert(isnan(var))});
    PARAMS.AutoCallTypes = readParam(params_text, 'AutoCallTypes', {@(var)validateattributes(var,{'numeric'},{'integer'}), @(var)assert(isnan(var))});
    PARAMS.StartDateTime = readParam(params_text, 'StartDateTime', {@(var)validateattributes(var,{'numeric'},{'numel',6}), @(var)assert(isnan(var))});
    PARAMS.StopDateTime = readParam(params_text, 'StopDateTime', {@(var)validateattributes(var,{'numeric'},{'numel',6}), @(var)assert(isnan(var))});
    
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