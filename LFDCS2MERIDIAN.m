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
%   Last updated May 30, 2022 using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES
% 2022-06-17
% ----------
% - WARNING! Currently, this script may fail if there are subfolders that
% contain earlier recordings, as the recordings may be out of order

    % parse input
    p = inputParser;
    p.addParameter('input_file', '', @ischar)
    p.addParameter('wav_dir', '', @ischar)
    p.addParameter('output_file', '', @ischar)
    p.addParameter('sp_codes', [], @isnumeric)  % original default was [9999,0,-32767]; empty selects all
    p.addParameter('wav_subfolders', true, @islogical)
    
    p.parse(varargin{:})
    inFilePath = p.Results.input_file;
    wavDir = p.Results.wav_dir;
    outFilePath = p.Results.output_file;
    sp_codes = p.Results.sp_codes;
    search_wav_subfolders = p.Results.wav_subfolders;

    % define common variables
    dtRef = datetime(1970,1,1,0,0,0);

    
    % 1) GET INPUT ........................................................
    
    % template MERIDIAN output xlsx
    browserPath = which('validateMeridianDetections');
    if isempty(browserPath)
        error('Could not find Detection Browser output template file')
    end
    [browserDir,~,~] = fileparts(browserPath);
    outTemplatePath = fullfile(browserDir,'BrowserResources','OutputTemplate.xlsx');
    
    % LFDCS spreadsheet
    if isempty(inFilePath)
        [inFileName,inFileDir] = uigetfile({'*.csv';'*.xlsx'},'Select LFDCS spreadsheet');
        inFilePath = fullfile(inFileDir,inFileName);
        if isnumeric(inFileName)
            return
        end
    end
    
    % WAV folder
    if isempty(wavDir)
        wavDirRoot = '\\142.2.83.111\teamwhalenas2\MOORED_PAM_DATA';
        wavDir = uigetdir(wavDirRoot,'Select WAV folder');
        if isnumeric(wavDir)
            return
        end
    end
    
    
    % 2) GET WAV FILE LIST ................................................
    disp('Getting WAV file times...')
    [wavFileNames,~] = Utilities.getFileNames(wavDir, 'wav', search_wav_subfolders);
    
    % extract datetime from WAV files
    dtWav = Utilities.readDateTime(wavFileNames);
 

    % 3) EXTRACT LFDCS DATA ...............................................
    disp('Extracting LFDCS data...')
    
    % set number of LFDCS header lines
    LFDCSHeader = 23;
    
    % set column containing manual species codes
    iColSpecies = 10; %15 for commented XLSX files...
    
    % read LFDCS autodetections file as table
    importOpts = detectImportOptions(inFilePath, 'NumHeaderLines',LFDCSHeader, 'DatetimeType','text');
    tableLFDCS = readtable(inFilePath, importOpts);
    
    % truncate table to include only entries of interest
    if ~isempty(sp_codes)
        rowsInclude = ismember(tableLFDCS{:,iColSpecies},sp_codes);
        tableLFDCS = tableLFDCS(rowsInclude,:);
    end
    n = height(tableLFDCS);
    
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

    
    % 4) CREATE OUTPUT ....................................................
    disp('Writing output...')
    
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