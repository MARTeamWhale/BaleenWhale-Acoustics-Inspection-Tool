%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function "validateMeridianDetections"
%   Written by Wilfried Beslin
%   Last updated Nov 27, 2023 using MATLAB R2018b
%
%   Description:
%   Allows a user to browse through and validate spectrograms of potential 
%   right whale upcalls (or other whale calls of interest) detected by the 
%   MERIDIAN neural network. It is also possible to use this for LFDCS
%   detections by converting the LFDCS autodetections outputs into
%   MERIDIAN-style spreadsheets using the LFDCS2MERIDIAN script.
%   Also supports marking of undetected (missed) calls and general 
%   annotations.
%
%   See reference manual for details.
%
%   Syntax
%   -----------------------------------------------------------------------
%   validateMeridianDetections - launches the browser with default 
%       settings, and prompts the user to specify the input spreadsheet 
%       file and folder of WAV files
%
%   validateMeridianDetections(Name,Value) - launches the browser with
%       specified Name-Value pair input arguments. Multiple Name-Value
%       pairs can be specified at once. These named arguments are the
%       following:
%           'params' - string specifying name or path of parameter file to
%               use. If not specified, the program will load the default
%               parameters sheet.
%           'data_file' - string specifying the input CSV or XLSX 
%               spreadsheet. If not specified, user will be prompted to 
%               choose one.
%           'audio_folder' - string specifying the path of the folder
%               containing the raw WAV files. If not specified, user will
%               be prompted to select the folder. 
%
%   Usage
%   -----------------------------------------------------------------------
%   Most navigation is done via the keyboard, as specified below. Note that
%   the main figure must be in focus for the keyboard commands to execute.
%   
%       W = zoom in
%       S = zoom out
%       A = previous detection
%       D = next detection
%
%       UpArrow/Shift+W = increase frequency range
%       DownArrow/Shift+S = decrease frequency range
%       LeftArrow/Shift+A = pan left
%       RightArrow/Shift+D = pan right
%
%       R = reset zoom and translation
%       T = toggle detection windows on/off
%       C = change colormap
%       Q = quit
%       Esc = Cancel operations:
%           - marking/unmarking of undetected upcalls
%           - audio playback
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES
%{
% - Some elements are a bit buggy. Unfortunately, there's nothing that can 
%   be done (easily). Fortunately, they're all very minor:
%       -- The bar shown during audio playback doesn't go all the way to
%       the end
%           --- BUT, checkout the animatedline() function - it may or may
%           not be a solution
%       -- Buttons don't all resize properly when expanding figure
%       -- Audio playback does not work for extreme sampling rates and/or
%       playback speeds. This can be a limitation with both MATLAB and the
%       computer's audio hardware.
%           --- For now, the program is designed to issue a simple warning 
%           if playback fails. But it may be possible to mitigate the
%           problem somewhat by resampling.
% - Figure updating can be very slow especially when dense spectrograms
%   have to be displayed. Profiling has shown that the drawnow() function 
%   is what consumes the most time by far. Unfortunately, all my attempts 
%   to overcome the limitations of drawnow() have met with limited success.
%   It is costly to render dense data in MATLAB and there is only so 
%   much optimization that can be done.
%}


function validateMeridianDetections(varargin)

    % Clear profiler (DEBUG ONLY)
    %profile clear

    % define global variables
    % NOTE: don't forget to clear these in the figure deletion routine!
    global OUTPUT PARAMS DATA UI PLOT
    
    disp('Initializing...')
    
    % define paths to resource folders
    rootDir = mfilename('fullpath');
    [rootDir,~,~] = fileparts(rootDir);
    resDir = fullfile(rootDir,'BrowserResources');
    paramsDir = fullfile(rootDir,'PARAMS');
    
    % parse input
    p = inputParser;
    p.addParameter('params','default_params.txt', @(v)ischar(v))
    p.addParameter('data_file', '', @(v)ischar(v));
    p.addParameter('audio_folder', '', @(v)ischar(v));
    p.parse(varargin{:})
    
    %%% process params file
    paramsFilePath = p.Results.params;
    [userParamsDir, paramsFilename, paramsExt] = fileparts(paramsFilePath);
    if isempty(paramsExt)
        paramsExt = '.txt';
    end
    if isempty(userParamsDir)
        paramsFilePath = fullfile(paramsDir,[paramsFilename,paramsExt]);
    end
    
    %%% process data file
    inFilePath = p.Results.data_file;
    if isempty(inFilePath)
        % prompt user for input spreadsheet file
        [inFile,inDir] = uigetfile('*.csv;*.xlsx','Select initial MERIDIAN or working spreadsheet file');
        inFilePath = fullfile(inDir,inFile);
        if isnumeric(inFile)
            return
        end
    end
    
    %%% process audio folder
    wavDir = p.Results.audio_folder;
    if isempty(wavDir)
        % prompt user for WAV folder
        wavDir = uigetdir(pwd,'Specify WAV folder');
        if isnumeric(wavDir)
            return
        end
    end
    
    % read input file
    [OUTPUT,cancel] = readInput(inFilePath,resDir);
    if cancel
        return
    end
    
    % set static parameters
    PARAMS = loadSetParams(paramsFilePath);
    
    % initialize program data
    DATA = initializeData(OUTPUT,PARAMS,wavDir);

    % build UI components
    UI = buildUI(PARAMS);
    
    % initialize plot elements
    PLOT = initializePlotElements(UI.axesMain,UI.axesHidden);

    % assign callbacks
    assignCallbacks(UI);
    
    % set state
    changeState('normal',false)

    % initialize first detection
    jumpStart = do_GotoUnclassified();
    if ~jumpStart
        DATA.currentDet.idx = DATA.detInfo.n;
        processNewDet();
    end
    
    disp('Initialization complete')
end


%% INITIALIZATION FUNCTIONS ===============================================
% readInput ---------------------------------------------------------------
function [OUTPUT,cancel] = readInput(inFilePath,resDir)
% Reads input file, determines its type, and setup output as appropriate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % determine type of input spreadsheet (initial MERIDIAN CSV or 
    % processed XLSX) based on file extension
    [~,~,inFileExt] = fileparts(inFilePath);
    switch inFileExt
        case '.csv'
            [OUTPUT,cancel] = parseInitialCSV(inFilePath,resDir);
        case '.xlsx'
            OUTPUT = parseWorkingXLSX(inFilePath);
            cancel = false;
        otherwise
            error('Expected input file to be a CSV or XLSX file')
    end
    
    % Make sure undetected table has correct field types if it's empty
    if isempty(OUTPUT.tableMissedCalls)
        OUTPUT.tableMissedCalls.FileName = cell.empty(0,1);
        OUTPUT.tableMissedCalls.FileStart = double.empty(0,1);
        OUTPUT.tableMissedCalls.SigStart = double.empty(0,1);
        OUTPUT.tableMissedCalls.SigEnd = double.empty(0,1);
        OUTPUT.tableMissedCalls.SigStartDateTime = datetime.empty(0,1);
        OUTPUT.tableMissedCalls.Type = cell.empty(0,1);
        OUTPUT.tableMissedCalls.Comments = cell.empty(0,1);
    end
    
    % Make sure annotations table has correct field types if it's empty
    if isempty(OUTPUT.tableAnnotations)
        OUTPUT.tableAnnotations.FileName = cell.empty(0,1);
        OUTPUT.tableAnnotations.FileStart = double.empty(0,1);
        OUTPUT.tableAnnotations.SigStart = double.empty(0,1);
        OUTPUT.tableAnnotations.SigEnd = double.empty(0,1);
        OUTPUT.tableAnnotations.SigMinFreq = double.empty(0,1);
        OUTPUT.tableAnnotations.SigMaxFreq = double.empty(0,1);
        OUTPUT.tableAnnotations.SigStartDateTime = datetime.empty(0,1);
        OUTPUT.tableAnnotations.Comments = cell.empty(0,1);
    end
    
    
    % NESTED FUNCTIONS 
    % parse Initial CSV ...................................................
    function [OUTPUT,cancel] = parseInitialCSV(inFilePath,resDir)
    % Processes initial MERIDIAN CSV file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % initialize output struct
        OUTPUT = struct();
        
        % read input file
        inTable = readtable(inFilePath);
        
        % check columns to make sure file is as expected
        cols = {'filename','det_start','det_stop','threshold','file_start','label','start','stop','fmin','fmax','val'};
        colsInput = inTable.Properties.VariableNames;
        assert(all(ismember(colsInput,cols)),'File columns not recognized')
        
        % create output file and variable
        [outFileName,outFileDir] = uiputfile('*.xlsx','Specify output file');
        outFilePath = fullfile(outFileDir,outFileName);
        if ischar(outFileName)
            % delete existing file if it exists
            if isfile(outFilePath)
                delete(outFilePath);
            end
            
            % create output file and variable
            [tableDetections,tableMissedCalls,tableAnnotations] = createOutput(inTable,outFilePath,resDir);
            OUTPUT.tableDetections = tableDetections;
            OUTPUT.tableMissedCalls = tableMissedCalls;
            OUTPUT.tableAnnotations = tableAnnotations;
            OUTPUT.path = outFilePath;
            cancel = false;
        else
            % cancel if user doesn't specify path
            cancel = true;
        end
    end

    % parseWorkingXLSX ....................................................
    function OUTPUT = parseWorkingXLSX(inFilePath)
    % Processes existing XLSX file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % initialize output struct
        OUTPUT = struct();
        
        % read input file
        tableDetections = readtable(inFilePath,'Sheet','Detected','DateTimeType','datetime');
        tableMissedCalls = readtable(inFilePath,'Sheet','Undetected','DateTimeType','datetime');
        tableAnnotations = readtable(inFilePath,'Sheet','Annotations','DateTimeType','datetime');
        
        % check columns to make sure file is as expected
        colsDetectedExpected = {'FileName','FileStart','SigStart','SigEnd','SigStartDateTime','Class_LFDCS','Class_MATLAB','ReasonForUNK','Comments'};
        colsUndetectedExpected = {'FileName','FileStart','SigStart','SigEnd','SigStartDateTime','Type','Comments'};
        colsAnnotationsExpected = {'FileName','FileStart','SigStart','SigEnd','SigMinFreq','SigMaxFreq','SigStartDateTime','Comments'};
        colsDetected = tableDetections.Properties.VariableNames;
        colsUndetected = tableMissedCalls.Properties.VariableNames;
        colsAnnotations = tableAnnotations.Properties.VariableNames;
        assert(all(ismember(colsDetected,colsDetectedExpected)),'Columns of Detected sheet not recognized')
        assert(all(ismember(colsUndetected,colsUndetectedExpected)),'Columns of Undetected sheet not recognized')
        assert(all(ismember(colsAnnotations,colsAnnotationsExpected)),'Columns of Annotations sheet not recognized')
        
        % make sure fields expected to be cellstrs are indeed cellstrs
        tableDetections = forceTableCellStr(tableDetections,{'Class_LFDCS','Class_MATLAB','ReasonForUNK','Comments'});
        tableMissedCalls = forceTableCellStr(tableMissedCalls,{'Type','Comments'});
        tableAnnotations = forceTableCellStr(tableAnnotations,{'Comments'});
        
        % convert datetime format
        tableDetections.SigStartDateTime.Format = 'dd-MMM-yyyy HH:mm:ss';
        if height(tableMissedCalls) > 0
            tableMissedCalls.SigStartDateTime.Format = 'dd-MMM-yyyy HH:mm:ss';
        else
            tableMissedCalls.SigStartDateTime = datetime.empty(0,1);
        end
        if height(tableAnnotations) > 0
            tableAnnotations.SigStartDateTime.Format = 'dd-MMM-yyyy HH:mm:ss';
        else
            tableAnnotations.SigStartDateTime = datetime.empty(0,1);
        end
        
        % set output struct
        OUTPUT.tableDetections = tableDetections;
        OUTPUT.tableMissedCalls = tableMissedCalls;
        OUTPUT.tableAnnotations = tableAnnotations;
        OUTPUT.path = inFilePath;
    end
    
    % createOutput ........................................................
    function [tableDetections,tableMissedCalls,tableAnnotations] = createOutput(inTable,outFilePath,resDir)
    % initializes output variables and creates the file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % define variables
        cols = {'FileName','FileStart','SigStart','SigEnd','SigStartDateTime','Class_LFDCS','Class_MATLAB','ReasonForUNK','Comments'};
        blankCellstr = repmat({''},height(inTable),1);
        dtRef = datetime(1970,1,1,0,0,0);
    
        % copy template output file to specified output path
        templateFilePath = fullfile(resDir,'OutputTemplate.xlsx');
        [copyOK,copyMsg] = copyfile(templateFilePath,outFilePath);
        if ~copyOK
            error('Could not create output file:\n%s',copyMsg)
        end
        
        % get detection start datetimes
        dtStart = dtRef + seconds(inTable.file_start) + seconds(inTable.det_start);

        % get LFDCS classification
        classLFDCS_num = inTable.val;
        classLFDCS_str = blankCellstr;
        if isnumeric(classLFDCS_num)
            classLFDCS_str(classLFDCS_num == 9999) = {'Correct'};
            classLFDCS_str(classLFDCS_num == 0) = {'Unknown'};
            classLFDCS_str(classLFDCS_num == -9999) = {'Incorrect'};
            classLFDCS_str(classLFDCS_num == -32767) = {'Unclassified'};
        end

        % get data for each column
        FileName = inTable.filename;
        FileStart = inTable.file_start;
        SigStart = inTable.det_start;
        SigEnd = inTable.det_stop;
        SigStartDateTime = dtStart;
        Class_LFDCS = classLFDCS_str;
        Class_MATLAB = blankCellstr;
        ReasonForUNK = blankCellstr;
        Comments = blankCellstr;
        
        % create Detected table
        tableDetections = table(FileName,FileStart,SigStart,SigEnd,SigStartDateTime,Class_LFDCS,Class_MATLAB,ReasonForUNK,Comments,...
            'VariableNames',cols);
        writetable(tableDetections,outFilePath,'Sheet','Detected');
        
        % initialize Undetected table
        tableMissedCalls = readtable(outFilePath,'Sheet','Undetected');
        
        % initialize Annotations table
        tableAnnotations = readtable(outFilePath,'Sheet','Annotations');
    end

    % forceTableCellStr ...................................................
    function tbl = forceTableCellStr(tbl,fields)
    % Ensures that specified fields in a table are cell arrays of strings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        nFields = numel(fields);
        for ii = 1:nFields
            fieldii = fields{ii};
            if isnumeric(tbl.(fieldii))
                tbl.(fieldii) = num2cell(tbl.(fieldii));
                for jj = 1:height(tbl)
                    if isnan(tbl.(fieldii){jj})
                        tbl.(fieldii){jj} = '';
                    else
                        tbl.(fieldii){jj} = num2str(tbl.(fieldii){jj});
                    end
                end
            end
        end
    end
end

% loadSetParams -----------------------------------------------------------
function PARAMS = loadSetParams(paramFile)
% Reads in program parameters from file (some may also be hard-coded)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    import BWAV_code.readParam

    % read parameter file as a block of text
    params_text = fileread(paramFile);

    % initialize output
    PARAMS = struct;
    
    % set spectrogram settings
    PARAMS.spec = struct;
    PARAMS.spec.WinSize = readParam(params_text, 'WinSize', {@(var)validateattributes(var,{'numeric'},{'scalar','positive'})});
    PARAMS.spec.StepSize = readParam(params_text, 'StepSize', {@(var)validateattributes(var,{'numeric'},{'scalar','positive'})});
    PARAMS.spec.NFFT_8kHz = readParam(params_text, 'NFFT_8kHz', {@(var)validateattributes(var,{'numeric'},{'scalar','positive','integer'})});
    
    % set plot parameters
    PARAMS.plot = struct();
    PARAMS.plot.TSpanDefault = readParam(params_text, 'TSpanDefault', {@(var)validateattributes(var,{'numeric'},{'scalar','positive'})});
    PARAMS.plot.FMaxDefault = readParam(params_text, 'FMaxDefault', {@(var)validateattributes(var,{'numeric'},{'scalar','positive'})});
    PARAMS.plot.TPanFactor = readParam(params_text, 'TPanFactor', {@(var)validateattributes(var,{'numeric'},{'scalar','positive'})});
    PARAMS.plot.TZoomFactor = readParam(params_text, 'TZoomFactor', {@(var)validateattributes(var,{'numeric'},{'scalar','positive'})});
    PARAMS.plot.FZoomFactor = readParam(params_text, 'FZoomFactor', {@(var)validateattributes(var,{'numeric'},{'scalar','positive'})});
    PARAMS.plot.Colormaps = struct(...
        'parula', @parula,...
        'parula_inverted', @()flipud(parula),...
        'jet', @jet,...
        'jet_inverted', @()flipud(jet),...
        'hsv', @hsv,...
        'hsv_inverted', @()flipud(hsv),...
        'hot', @hot,...
        'hot_inverted', @()flipud(hot),...
        'cool', @cool,...
        'cool_inverted', @()flipud(cool),...
        'spring', @spring,...
        'spring_inverted', @()flipud(spring),...
        'summer', @summer,...
        'summer_inverted', @()flipud(summer),...
        'autumn', @autumn,...
        'autumn_inverted', @()flipud(autumn),...
        'winter', @winter,...
        'winter_inverted', @()flipud(winter),...    
        'gray', @gray,...
        'gray_inverted', @()flipud(gray),...
        'bone', @bone,...
        'bone_inverted', @()flipud(bone),...    
        'copper', @copper,...
        'copper_inverted', @()flipud(copper),...    
        'pink', @pink,...
        'pink_inverted', @()flipud(pink));
    PARAMS.plot.DefaultColormap = lower(readParam(params_text, 'DefaultColormap', {@(var)validateattributes(var,{'char'},{'vector'})}));
    % ensure colormap is valid
    if ~ismember(PARAMS.plot.DefaultColormap,fieldnames(PARAMS.plot.Colormaps))
        bad_colmap = PARAMS.plot.DefaultColormap;
        PARAMS.plot.DefaultColormap = 'bone';
        warning('"%s" is not a valid colormap; defaulting to "%s"', bad_colmap, PARAMS.plot.DefaultColormap)
    end
    
    % set colour data
    PARAMS.col = struct();
    PARAMS.col.Unclassified = [1 1 1];
    PARAMS.col.Correct = [0.3922 0.8314 0.0745];
    PARAMS.col.Unknown = [0.0745 0.6235 1.0000];
    PARAMS.col.Incorrect = [0.9020 0.2000 0.2000];
    PARAMS.col.Annotations = [0.93 0.69 0.13];  % [1, 0.8, 0.2]
    
    % set marker parameters
    PARAMS.markers = struct();
    PARAMS.markers.LineWidth = readParam(params_text, 'LineWidth', {@(var)validateattributes(var,{'numeric'},{'scalar','positive'})});
    PARAMS.markers.FaceAlpha = readParam(params_text, 'FaceAlpha', {@(var)validateattributes(var,{'numeric'},{'scalar','nonnegative','<=',1})});
    PARAMS.markers.StandardFRange = readParam(params_text, 'StandardFRange', {@(var)validateattributes(var,{'numeric'},{'numel',2,'nonnegative','increasing'})});
    %%% Colors
    PARAMS.markers.Color = struct();
    PARAMS.markers.Color.Detections = zeros(2,3);
    PARAMS.markers.Color.Detections(1,:) = readParam(params_text, 'Color_DetectedOther', {@(var)validateattributes(var,{'numeric'},{'numel',3,'nonnegative','<=',1})});
    PARAMS.markers.Color.Detections(2,:) = readParam(params_text, 'Color_DetectedFocal', {@(var)validateattributes(var,{'numeric'},{'numel',3,'nonnegative','<=',1})});
    PARAMS.markers.Color.MissedCalls = zeros(2,3);
    PARAMS.markers.Color.MissedCalls(1,:) = setAutoColor(readParam(params_text, 'Color_MissedDefinite', {@(var)validateattributes(var,{'numeric'},{'numel',3,'nonnegative','<=',1}),@(var)assert(strcmpi(var,'auto'))}), PARAMS.col.Correct);
    PARAMS.markers.Color.MissedCalls(2,:) = setAutoColor(readParam(params_text, 'Color_MissedPotential', {@(var)validateattributes(var,{'numeric'},{'numel',3,'nonnegative','<=',1}),@(var)assert(strcmpi(var,'auto'))}), PARAMS.col.Unknown);
    PARAMS.markers.Color.Annotations = setAutoColor(readParam(params_text, 'Color_Annotations', {@(var)validateattributes(var,{'numeric'},{'numel',3,'nonnegative','<=',1}),@(var)assert(strcmpi(var,'auto'))}), PARAMS.col.Annotations);
    %%% Frequency range type
    PARAMS.markers.FRangeType = struct();
    PARAMS.markers.FRangeType.Detections = 'standard';
    PARAMS.markers.FRangeType.MissedCalls = 'standard';
    PARAMS.markers.FRangeType.Annotations = 'variable';
    
    % NESTED FUNCTIONS 
    % setAutoColor ........................................................
    function newcol = setAutoColor(colval, autocol)
    % Decide colour if colour value is "auto"
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmpi(colval, 'auto')
            newcol = autocol;
        else
            newcol = colval;
        end
    end
end

% initializeData ----------------------------------------------------------
function DATA = initializeData(OUTPUT,PARAMS,wavDir)
% Initialize dependent and/or variable program variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get/define common variables
    outTable = OUTPUT.tableDetections;
    dtRef = datetime(1970,1,1,0,0,0);
    [wavList,iUniqueWav,iDetWav] = unique(outTable.FileName);

    % initialize output
    DATA = struct;
    
    % initialize plotting data
    DATA.plot = struct;
    DATA.plot.tSpan = PARAMS.plot.TSpanDefault;
    DATA.plot.tOffset = 0;
    DATA.plot.fMax = PARAMS.plot.FMaxDefault;
    DATA.plot.colmapActive = PARAMS.plot.DefaultColormap;
    DATA.plot.tRange = [NaT, NaT];
    DATA.plot.fRange = [NaN, NaN];

    % get info on all detections
    DATA.detInfo = struct;
    DATA.detInfo.n = height(outTable);
    DATA.detInfo.startTimes = outTable.SigStart;
    DATA.detInfo.endTimes = outTable.SigEnd;
    DATA.detInfo.iWav = iDetWav;

    % info on all WAV files
    DATA.wavInfo = struct;
    DATA.wavInfo.dir = wavDir;
    DATA.wavInfo.Names = wavList;
    DATA.wavInfo.dtStarts = dtRef + seconds(outTable.FileStart(iUniqueWav));
    DATA.wavInfo.dtStartRef = dtRef;

    % info on current detection (including WAV file data)
    DATA.currentDet = struct;
    DATA.currentDet.idx = 1;
    DATA.currentDet.tStart = NaN;
    DATA.currentDet.tEnd = NaN;
    DATA.currentDet.wav = struct();
    DATA.currentDet.wav.Name = '';
    DATA.currentDet.wav.dtStart = NaT;
    DATA.currentDet.wav.idx = NaN;
    DATA.currentDet.wav.x = [];
    DATA.currentDet.wav.Fs = NaN;
    DATA.currentDet.wav.BitsPerSample = NaN;
    
    % initialize marker data
    DATA.markers = struct();
    DATA.markers.show = true;
    %%% Screen indices
    DATA.markers.ScreenIndices = struct();
    DATA.markers.ScreenIndices.Detections = NaN(height(OUTPUT.tableDetections),1);
    DATA.markers.ScreenIndices.MissedCalls = NaN(height(OUTPUT.tableMissedCalls),1);
    DATA.markers.ScreenIndices.Annotations = NaN(height(OUTPUT.tableAnnotations),1);
    %%% Colour indices
    DATA.markers.ColorIndices = struct();
    DATA.markers.ColorIndices.Detections = NaN(height(OUTPUT.tableDetections),1);
    DATA.markers.ColorIndices.MissedCalls = NaN(height(OUTPUT.tableMissedCalls),1);
    DATA.markers.ColorIndices.Annotations = NaN(height(OUTPUT.tableAnnotations),1);
    
    % playback data
    DATA.playback = struct();
    DATA.playback.player = audioplayer.empty(0);
    DATA.playback.speed = 1;
    DATA.playback.vol = 1;
    
    % program state
    DATA.program = struct();
    DATA.program.mode = 'normal'; %%% May be 'normal', 'draw', 'edit', or 'playback'
    DATA.program.drawState = ''; %%% May be 'start', 'drag', 'end', or ''
    DATA.program.drawPurpose = ''; %%% May be 'play_audio', 'save_audio', 'mark_undetected', or 'annotate'
    DATA.program.drawType = ''; %%% May be 'x' or 'xy'
    DATA.program.undetectedType = ''; %%% May be 'Definite' or 'Potential'
    DATA.program.editState = ''; %%% May be 'select' or 'resize'
    DATA.program.editType = ''; %%% May be 'MissedCalls' or 'Annotations'
    DATA.program.sheetChanged = struct();
    DATA.program.sheetChanged.Detections = false;  %%% indicates if "Detected" sheet has been updated or not
    DATA.program.sheetChanged.MissedCalls = false;  %%% indicates if "Undetected" sheet has been updated or not
    DATA.program.sheetChanged.Annotations = false;  %%% indicates if "Annotations" sheet has been updated or not 
end

% buildUI -----------------------------------------------------------------
function UI = buildUI(PARAMS)
% Creates and initializes the figure and UI components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: UI figure resizing is bugged for reasons I don't understand.
% This is particularly an issue with UI Buttons - they don't seem to resize
% except when a text area is present within their panels. And even when a
% text area is present, they may not resize properly.
% The only sensible solution may to to implement a custom resize function
% and resize all UI components programatically. The auto system is broken.

    % define common properties
    axesPosCorrection = [50 50 -60 -60];
    axesFontSize = 12;
    axesFontWeight = 'b';
    
    % get monitor positions (use last monitor) and figure position vector
    mpAll = get(0,'MonitorPosition');
    mp = mpAll(end,:);
    fHeight = 720;
    fWidth = 1280;
    xPad = 20;
    yPad = 40;
    fX = mp(1)+xPad;
    fY = mp(4)-yPad-fHeight + (mp(2)-1);   
    fPos = [fX,fY,fWidth,fHeight];
    
    % initialize output
    UI = struct;

    % MAIN ................................................................
    % Create figure
    UI.fig = uifigure;
    UI.fig.Position = fPos;
    UI.fig.Name = 'Whale Call Validator';

    % Create axesMain
    UI.axesMain = axes(UI.fig);
    xlabel(UI.axesMain, 'Time',...
        'FontSize', axesFontSize,...
        'FontWeight', axesFontWeight)
    ylabel(UI.axesMain, 'Frequency [Hz]',...
        'FontSize', axesFontSize,...
        'FontWeight', axesFontWeight)
    UI.axesMain.Color = [0 0 0];
    UI.axesMain.Units = 'pixels';
    UI.axesMain.Position = [191 21 880 680] + axesPosCorrection;
    UI.axesMain.Units = 'normalized';
    UI.axesMain.NextPlot = 'add';
    UI.axesMain.PickableParts = 'none';

    % Create axesHidden
    UI.axesHidden = axes(UI.fig);
    UI.axesHidden.Units = 'normalized';
    UI.axesHidden.Position = UI.axesMain.Position;
    UI.axesHidden.NextPlot = 'add';
    UI.axesHidden.PickableParts = 'none';


    % SPECTROGRAM PANEL ...................................................
    % Create panelSpectrogram
    UI.panelSpectrogram = uipanel(UI.fig);
    UI.panelSpectrogram.Title = 'SPECTROGRAM';
    UI.panelSpectrogram.FontWeight = 'bold';
    UI.panelSpectrogram.Position = [21 446 150 255];
    
    % Create buttonResetTimeSpan
    UI.buttonResetTimeSpan = uibutton(UI.panelSpectrogram, 'push');
    UI.buttonResetTimeSpan.Position = [10 190 130 35];
    UI.buttonResetTimeSpan.Text = 'Reset Time Span';

    % Create buttonResetTimeOffset
    UI.buttonResetTimeOffset = uibutton(UI.panelSpectrogram, 'push');
    UI.buttonResetTimeOffset.Position = [10 145 130 35];
    UI.buttonResetTimeOffset.Text = 'Reset Time Offset';

    % Create buttonResetFreqRange
    UI.buttonResetFreqRange = uibutton(UI.panelSpectrogram, 'push');
    UI.buttonResetFreqRange.Position = [10 100 130 35];
    UI.buttonResetFreqRange.Text = 'Reset Freq. Range';

    % Create buttonChangeColormap
    UI.buttonChangeColormap = uibutton(UI.panelSpectrogram, 'push');
    UI.buttonChangeColormap.Position = [10 55 130 35];
    UI.buttonChangeColormap.Text = 'Change Colormap';

    % Create buttonToggleMarkers
    UI.buttonToggleMarkers = uibutton(UI.panelSpectrogram, 'push');
    UI.buttonToggleMarkers.Position = [10 10 130 35];
    UI.buttonToggleMarkers.Text = 'Toggle Markers';
    
    
    % DETECTION PANEL .....................................................
    % Create panelDetection
    UI.panelDetection = uipanel(UI.fig);
    UI.panelDetection.Title = 'DETECTION';
    UI.panelDetection.FontWeight = 'bold';
    UI.panelDetection.Position = [21 301 150 120];
    
    % Create buttonJump
    UI.buttonJump = uibutton(UI.panelDetection, 'push');
    UI.buttonJump.Position = [10 55 130 35];
    UI.buttonJump.Text = 'Jump';
    
    % Create buttonFirstUnclassified
    UI.buttonFirstUnclassified = uibutton(UI.panelDetection, 'push');
    UI.buttonFirstUnclassified.Position = [10 10 130 35];
    UI.buttonFirstUnclassified.Text = 'First Unclassified';


    % PLAYBACK PANEL ......................................................
    % Create panelPlayback
    UI.panelPlayback = uipanel(UI.fig);
    UI.panelPlayback.Title = 'PLAYBACK';
    UI.panelPlayback.FontWeight = 'bold';
    UI.panelPlayback.Position = [21 21 150 255];

    % Create buttonPlayDetection
    UI.buttonPlayDetection = uibutton(UI.panelPlayback, 'push');
    UI.buttonPlayDetection.Position = [10 190 130 35];
    UI.buttonPlayDetection.Text = 'Play Detection';
    
    % Create buttonPlayRange
    UI.buttonPlayRange = uibutton(UI.panelPlayback, 'push');
    UI.buttonPlayRange.Position = [10 145 130 35];
    UI.buttonPlayRange.Text = 'Play Range';
    
    % Create buttonChangeSpeedVol
    UI.buttonChangeSpeedVol = uibutton(UI.panelPlayback, 'push');
    UI.buttonChangeSpeedVol.Position = [10 100 130 35];
    UI.buttonChangeSpeedVol.Text = 'Change Speed/Vol';
    
    % Create buttonSaveClip
    UI.buttonSaveClip = uibutton(UI.panelPlayback, 'push');
    UI.buttonSaveClip.Position = [10 55 130 35];
    UI.buttonSaveClip.Text = 'Save Clip';

    % Create buttonStopPlayback
    UI.buttonStopPlayback = uibutton(UI.panelPlayback, 'push');
    UI.buttonStopPlayback.BackgroundColor = [0.902 0.2 0.2];
    UI.buttonStopPlayback.Position = [10 10 130 35];
    UI.buttonStopPlayback.Text = 'Stop Playback';
    UI.buttonStopPlayback.Enable = 'off';


    % VALIDATION PANEL ....................................................
    % Create panelValidation
    UI.panelValidation = uipanel(UI.fig);
    UI.panelValidation.Title = 'VALIDATION';
    UI.panelValidation.FontWeight = 'bold';
    UI.panelValidation.Position = [1091 386 170 315];

    % Create buttonResetValidation
    UI.buttonResetValidation = uibutton(UI.panelValidation, 'push');
    UI.buttonResetValidation.BackgroundColor = PARAMS.col.Unclassified;
    UI.buttonResetValidation.FontWeight = 'bold';
    UI.buttonResetValidation.Position = [10 245 150 40];
    UI.buttonResetValidation.Text = 'Reset';

    % Create buttonMarkCorrect
    UI.buttonMarkCorrect = uibutton(UI.panelValidation, 'push');
    UI.buttonMarkCorrect.BackgroundColor = PARAMS.col.Correct;
    UI.buttonMarkCorrect.FontWeight = 'bold';
    UI.buttonMarkCorrect.Position = [10 195 150 40];
    UI.buttonMarkCorrect.Text = 'Correct';

    % Create buttonMarkIncorrect
    UI.buttonMarkIncorrect = uibutton(UI.panelValidation, 'push');
    UI.buttonMarkIncorrect.BackgroundColor = PARAMS.col.Incorrect;
    UI.buttonMarkIncorrect.FontWeight = 'bold';
    UI.buttonMarkIncorrect.Position = [10 145 150 40];
    UI.buttonMarkIncorrect.Text = 'Incorrect';

    % Create buttonMarkUnknown
    UI.buttonMarkUnknown = uibutton(UI.panelValidation, 'push');
    UI.buttonMarkUnknown.BackgroundColor = PARAMS.col.Unknown;
    UI.buttonMarkUnknown.FontWeight = 'bold';
    UI.buttonMarkUnknown.Position = [10 95 150 40];
    UI.buttonMarkUnknown.Text = 'Unknown';

    % Create labelReasonForUnknown
    UI.labelReasonForUnknown = uilabel(UI.panelValidation);
    UI.labelReasonForUnknown.Position = [10 60 150 25];
    UI.labelReasonForUnknown.Text = 'Reason for Unknown:';
    UI.labelReasonForUnknown.Enable = 'off';

    % Create textareaReasonForUnknown
    UI.textareaReasonForUnknown = uitextarea(UI.panelValidation);
    UI.textareaReasonForUnknown.Position = [10 10 150 50];
    UI.textareaReasonForUnknown.Editable = 'off';
    UI.textareaReasonForUnknown.Enable = 'off';


    % COMMENT PANEL .......................................................
    % Create panelComment
    UI.panelComment = uipanel(UI.fig);
    UI.panelComment.Title = 'COMMENT';
    UI.panelComment.FontWeight = 'bold';
    UI.panelComment.Position = [1091 246 170 125];

    % Create textareaComment
    UI.textareaComment = uitextarea(UI.panelComment);
    UI.textareaComment.Position = [10 35 150 60];
    UI.textareaComment.Editable = 'off';
    
    % Create buttonEditComment
    UI.buttonEditComment = uibutton(UI.panelComment, 'push');
    UI.buttonEditComment.Position = [10 10 150 20];
    UI.buttonEditComment.Text = 'Edit';


    % MISSED CALLS PANEL ..................................................
    % Create panelMissedCalls
    UI.panelMissedCalls = uipanel(UI.fig);
    UI.panelMissedCalls.Title = 'MISSED CALLS';
    UI.panelMissedCalls.FontWeight = 'bold';
    UI.panelMissedCalls.Position = [1091 111 170 120];  

    % Create buttonAddDefiniteMissed
    UI.buttonAddDefiniteMissed = uibutton(UI.panelMissedCalls, 'push');
    UI.buttonAddDefiniteMissed.BackgroundColor = [0 0 0];
    UI.buttonAddDefiniteMissed.FontWeight = 'bold';
    UI.buttonAddDefiniteMissed.FontColor = PARAMS.col.Correct;
    UI.buttonAddDefiniteMissed.Position = [10 55 90 35];
    UI.buttonAddDefiniteMissed.Text = 'Add Definite';
    
    % Create buttonAddPotentialMissed
    UI.buttonAddPotentialMissed = uibutton(UI.panelMissedCalls, 'push');
    UI.buttonAddPotentialMissed.BackgroundColor = [0 0 0];
    UI.buttonAddPotentialMissed.FontWeight = 'bold';
    UI.buttonAddPotentialMissed.FontColor = PARAMS.col.Unknown;
    UI.buttonAddPotentialMissed.Position = [10 10 90 35];
    UI.buttonAddPotentialMissed.Text = 'Add Potential';

    % Create statebuttonEditMissedCalls
    UI.statebuttonEditMissedCalls = uibutton(UI.panelMissedCalls, 'state');
    UI.statebuttonEditMissedCalls.Position = [110 10 50 80];
    UI.statebuttonEditMissedCalls.Text = 'Edit';
    
    
    % ANNOTATIONS PANEL ...................................................
    % Create panelAnnotations
    UI.panelAnnotations = uipanel(UI.fig);
    UI.panelAnnotations.Title = 'ANNOTATIONS';
    UI.panelAnnotations.FontWeight = 'bold';
    UI.panelAnnotations.Position = [1091 21 170 75];  
    
    % Create buttonAddAnnotation
    UI.buttonAddAnnotation = uibutton(UI.panelAnnotations, 'push');
    UI.buttonAddAnnotation.BackgroundColor = PARAMS.col.Annotations;
    UI.buttonAddAnnotation.FontWeight = 'bold';
    UI.buttonAddAnnotation.FontColor = [0 0 0];
    UI.buttonAddAnnotation.Position = [10 10 90 35];
    UI.buttonAddAnnotation.Text = 'Add New';
    
    % Create statebuttonEditAnnotations
    UI.statebuttonEditAnnotations = uibutton(UI.panelAnnotations, 'state');
    UI.statebuttonEditAnnotations.Position = [110 10 50 35];
    UI.statebuttonEditAnnotations.Text = 'Edit';
end

% initializePlotElements --------------------------------------------------
function PLOT = initializePlotElements(haMain,haHidden)
% Creates and initializes elements that appear in plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define common properties
    textInfoFontSize = 12;
    textInfoFontColor = [1 1 1];
    textLoadingFontSize = 20;
    patchLoadingFaceColor = [190 200 215]/255;
    patchLoadingEdgeColor = [165 180 200]/255;
    patchLoadingAlpha = 0.3;
    
    % initialize output
    PLOT = struct();
    
    % MAIN ................................................................
    % Initialize spectrogram surface
    PLOT.surfaceSpectrogram = matlab.graphics.GraphicsPlaceholder;
    
    % initialize marker patche arrays
    PLOT.patchesDetections = matlab.graphics.GraphicsPlaceholder;
    PLOT.patchesMissedCalls = matlab.graphics.GraphicsPlaceholder;
    PLOT.patchesAnnotations = matlab.graphics.GraphicsPlaceholder;
    
    % Initialize selection patch
    PLOT.patchSelection = matlab.graphics.GraphicsPlaceholder;
    
    % Initialize playback line
    PLOT.linePlayback = matlab.graphics.GraphicsPlaceholder;
    
    % Initialize highlight lines
    PLOT.highlight = BWAV_code.MarkerHighlight(haMain);
 
    
    % HIDDEN ..............................................................
    % Create textDetInfo
    PLOT.textDetInfo = text(haHidden, 0, 1, '',...
        'HorizontalAlignment', 'left',...
        'VerticalAlignment', 'top',...
        'FontSize', textInfoFontSize,...
        'Color', textInfoFontColor,...
        'Interpreter','none');

    % Create textWindowInfo
    PLOT.textWindowInfo = text(haHidden, 1, 1, '',...
        'HorizontalAlignment', 'right',...
        'VerticalAlignment', 'top',...
        'FontSize', textInfoFontSize,...
        'Color', textInfoFontColor,...
        'Interpreter','none');

    % Create textValidation
    PLOT.textValidation = text(haHidden, 0.5, 1, '',...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'top',...
        'FontSize', textInfoFontSize,...
        'FontWeight','bold',...
        'Interpreter','none');

    % Create patchLoading
    PLOT.patchLoading = fill(haHidden, [0 0 1 1], [0 1 1 0], patchLoadingFaceColor,...
        'FaceAlpha', patchLoadingAlpha,...
        'EdgeColor', patchLoadingEdgeColor);

    % Create textLoading
    PLOT.textLoading = text(haHidden, 0.5, 0.5, 'Initializing...',...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle',...
        'FontSize',textLoadingFontSize); 
    
    % Make hidden axes invisible
    haHidden.Visible = 'off';
end

% assignCallbacks ---------------------------------------------------------
function assignCallbacks(UI)
% Assigns callback functions to UI components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    UI.fig.DeleteFcn = @FigDeleteFcn;
    UI.buttonResetTimeSpan.ButtonPushedFcn = @buttonResetTimeSpan_Callback;
    UI.buttonResetTimeOffset.ButtonPushedFcn = @buttonResetTimeOffset_Callback;
    UI.buttonResetFreqRange.ButtonPushedFcn = @buttonResetFreqRange_Callback;
    UI.buttonChangeColormap.ButtonPushedFcn = @buttonChangeColormap_Callback;
    UI.buttonToggleMarkers.ButtonPushedFcn = @buttonToggleMarkers_Callback;
    UI.buttonJump.ButtonPushedFcn = @buttonJump_Callback;
    UI.buttonFirstUnclassified.ButtonPushedFcn = @buttonFirstUnclassified_Callback;
    UI.buttonPlayDetection.ButtonPushedFcn = @buttonPlayDetection_Callback;
    UI.buttonPlayRange.ButtonPushedFcn = @buttonPlayRange_Callback;
    UI.buttonChangeSpeedVol.ButtonPushedFcn = @buttonChangeSpeedVol_Callback;
    UI.buttonSaveClip.ButtonPushedFcn = @buttonSaveClip_Callback;
    UI.buttonStopPlayback.ButtonPushedFcn = @buttonStopPlayback_Callback;
    UI.buttonResetValidation.ButtonPushedFcn = @buttonResetValidation_Callback;
    UI.buttonMarkCorrect.ButtonPushedFcn = @buttonValidation_Callback;
    UI.buttonMarkIncorrect.ButtonPushedFcn = @buttonValidation_Callback;
    UI.buttonMarkUnknown.ButtonPushedFcn = @buttonValidation_Callback;
    UI.buttonEditComment.ButtonPushedFcn = @buttonEditComment_Callback;
    UI.buttonAddDefiniteMissed.ButtonPushedFcn = @buttonAddMissedCall_Callback;
    UI.buttonAddPotentialMissed.ButtonPushedFcn = @buttonAddMissedCall_Callback;
    UI.statebuttonEditMissedCalls.ValueChangedFcn = @statebuttonEditMissedCalls_Callback;
    UI.buttonAddAnnotation.ButtonPushedFcn = @buttonAddAnnotation_Callback;
    UI.statebuttonEditAnnotations.ValueChangedFcn = @statebuttonEditAnnotations_Callback;
end


%% PROCESS FUNCTIONS ======================================================
% processNewDet -----------------------------------------------------------
function processNewDet()
% processes a new detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global OUTPUT DATA UI PLOT

    % Get index of current detection
    iDet = DATA.currentDet.idx;

    % check if WAV file needs to be updated
    iWavPrevious = DATA.currentDet.wav.idx;
    iWavCurrent = DATA.detInfo.iWav(iDet);
    if iWavCurrent ~= iWavPrevious
        PLOT.textLoading.String = 'Loading WAV file...';
        showLoadScreen();

        % get WAV file name and path
        wavName = DATA.wavInfo.Names{iWavCurrent};
        wavPath = fullfile(DATA.wavInfo.dir,wavName);

        % load WAV 
        w = audioinfo(wavPath);
        [x,Fs] = audioread(wavPath);

        % remove DC offset
        x = x - mean(x);
        
        %**TESTING - downsample to 1000 Hz
        %x = decimate(x,8,'fir');
        %Fs = 1000;

        % update WAV file data for current detection
        DATA.currentDet.wav.idx = iWavCurrent;
        DATA.currentDet.wav.Name = wavName;
        DATA.currentDet.wav.dtStart = DATA.wavInfo.dtStarts(iWavCurrent);
        DATA.currentDet.wav.x = x;
        DATA.currentDet.wav.Fs = Fs;
        DATA.currentDet.wav.BitsPerSample = w.BitsPerSample;

        % clean up
        clear w x Fs
    else
        showLoadScreen();
    end

    % update other info in currentDet (start/end times)
    DATA.currentDet.tStart = DATA.detInfo.startTimes(iDet);
    DATA.currentDet.tEnd = DATA.detInfo.endTimes(iDet);

    % reset tOffset (will also update plot)
    do_ResetPlot('offset');
    
    % update comment boxes
    UI.textareaComment.Value = OUTPUT.tableDetections.Comments(iDet);
    UI.textareaReasonForUnknown.Value = OUTPUT.tableDetections.ReasonForUNK(iDet);
    if strcmp(PLOT.textValidation.String,'Unknown')
        UI.textareaReasonForUnknown.Enable = 'on';
        UI.labelReasonForUnknown.Enable = 'on';
    else
        UI.textareaReasonForUnknown.Enable = 'off';
        UI.labelReasonForUnknown.Enable = 'off';
    end
    drawnow
    
    % save detection spreadsheet
    saveSpreadsheet();
    
    % reset loading string
    PLOT.textLoading.String = 'Updating...';
end

% updatePlot --------------------------------------------------------------
function updatePlot()
% Updates the spectrogram according to current settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % resume profiler (DEBUG ONLY)
    %profile resume

    % checkout global variables
    global OUTPUT PARAMS DATA UI PLOT

    % get relevant stuff
    ha = UI.axesMain;
    tDetStart = DATA.currentDet.tStart;
    tDetEnd = DATA.currentDet.tEnd;
    tDetMid = tDetStart + (tDetEnd - tDetStart)/2;
    dtStart = DATA.currentDet.wav.dtStart;
    dtDetStart = dtStart + seconds(tDetStart);
    Fs = DATA.currentDet.wav.Fs;
    tSpan = DATA.plot.tSpan;
    tOffset = DATA.plot.tOffset;
    tBufferOffscreen = PARAMS.spec.WinSize;
    fRange = [0, min([DATA.plot.fMax,Fs/2])];
    specWinSize = round(PARAMS.spec.WinSize*Fs);
    ovlSize = round((PARAMS.spec.WinSize - PARAMS.spec.StepSize)*Fs);
    nfft = Fs*(PARAMS.spec.NFFT_8kHz/8000);
    FColmap = PARAMS.plot.Colormaps.(DATA.plot.colmapActive);

    % get t-range of plot (relative to WAV start)
    tPlotStart = tDetMid + tOffset - tSpan/2;
    tPlotEnd = tDetMid + tOffset + tSpan/2;
    tRange = [tPlotStart,tPlotEnd];
    dtRange = dtStart + seconds(tRange);
    rangeDiff = diff(dtRange);
    rangeDiff.Format = 's';

    % extract target waveform samples (including off-screen buffer)
    startSample = round((tPlotStart - tBufferOffscreen)*Fs);
    endSample = round((tPlotEnd + tBufferOffscreen)*Fs);
    sampleVec = (startSample:endSample)';
    outOfRange = sampleVec < 1 | sampleVec > numel(DATA.currentDet.wav.x);
    xTarget = zeros(size(sampleVec));
    xTarget(~outOfRange) = DATA.currentDet.wav.x(sampleVec(~outOfRange));
    %%% areas that are out of range will be zero-padded. Fortunately, it
    %%% seems the "spectrogram" function ignores zeros at the start and end

    % get and replot spectrogram, if needed
    if ~isequal(tRange,DATA.plot.tRange) || ~isequal(fRange,DATA.plot.FRange)
        [~,F,T,P] = spectrogram(xTarget,specWinSize,ovlSize,nfft,Fs);

        %%% remove excess frequencies
        FExcess = F > fRange(2) + 2;
        F = F(~FExcess);
        P = P(~FExcess,:);

        %%% log amplitude
        S = 10*log10(abs(P));

        %%% adjust T vector (convert to datetime)
        T = dtStart + seconds(T + tPlotStart - tBufferOffscreen);
        
        delete(ha.Children);
        PLOT.surfaceSpectrogram = surf(ha, T, F, S, 'EdgeColor','none');
        axis(ha,'xy'); view(ha,0,90);
    else
        % Only reset non-spectrogram elements
        nonSurfKids = findobj(ha.Children, '-not','Type','Surface');
        delete(nonSurfKids)
    end
    colormap(ha,feval(FColmap));
    
    % update figure ranges in DATA
    DATA.plot.tRange = tRange;
    DATA.plot.FRange = fRange;
    
    % reset screen indices
    DATA.markers.ScreenIndices.Detections = NaN(height(OUTPUT.tableDetections),1);
    DATA.markers.ScreenIndices.MissedCalls = NaN(height(OUTPUT.tableMissedCalls),1);
    DATA.markers.ScreenIndices.Annotations = NaN(height(OUTPUT.tableAnnotations),1);

    % set color indices
    %%% Detections
    isTarget = false(DATA.detInfo.n,1);
    isTarget(DATA.currentDet.idx) = true;
    DATA.markers.ColorIndices.Detections = isTarget + 1;
    %%% MissedCalls
    isPotentialMissed = strcmp(OUTPUT.tableMissedCalls.Type,'Potential');
    DATA.markers.ColorIndices.MissedCalls = isPotentialMissed + 1;
    %%% Annotations
    DATA.markers.ColorIndices.Annotations = ones(height(OUTPUT.tableAnnotations),1);
    
    if DATA.markers.show
        % process markers in reverse order of importance:
        
        %%% Annotations
        drawMarkers(ha, 'Annotations', [], dtRange, fRange);
        
        %%% MissedCalls
        isPotentialMissed = strcmp(OUTPUT.tableMissedCalls.Type,'Potential');
        drawMarkers(ha, 'MissedCalls', isPotentialMissed, dtRange, fRange);
        
        %%% Detections
        isTarget = false(DATA.detInfo.n,1);
        isTarget(DATA.currentDet.idx) = true;
        drawMarkers(ha, 'Detections', isTarget, dtRange, fRange);
        
        % edit appearance depending on mode
        if strcmp(DATA.program.mode,'edit')
            switch DATA.program.editType
                case 'Annotations'
                    hpFocus = PLOT.patchesAnnotations;
                    hpNonFocus = [PLOT.patchesDetections;PLOT.patchesMissedCalls];
                case 'MissedCalls'
                    hpFocus = PLOT.patchesMissedCalls;
                    hpNonFocus = [PLOT.patchesDetections;PLOT.patchesAnnotations];
                otherwise
                    error('Bad case')
            end
            set(hpNonFocus,'LineStyle','--');
            set(hpNonFocus,'LineWidth',1); % EdgeAlpha still doesn't work, but LineWidth does
            
            % update status of highlight object
            if PLOT.highlight.hasSelection
                iHighlight = PLOT.highlight.signalIdx;
                iHighlightScreen = DATA.markers.ScreenIndices.(DATA.program.editType)(iHighlight);
                % if selected marker is on screen, update highlight object
                % with that info
                if ~isnan(iHighlightScreen)
                    refreshArgs = {hpFocus(iHighlightScreen),iHighlightScreen};
                else
                    refreshArgs = {};
                end
                PLOT.highlight.refresh(refreshArgs{:});
            end
        end
    end

    % set axes ranges and title
    ylim(ha,fRange)
    xlim(ha,dtRange)

    % set text for message strings
    %%% detection info string
    strDetNum = sprintf('Detection %d/%d',DATA.currentDet.idx,DATA.detInfo.n);
    strWav = DATA.currentDet.wav.Name;
    strTDet = char(dtDetStart);
    strDur = sprintf('duration = %s',char(seconds(tDetEnd - tDetStart)));
    strDetInfo = strcat({' '},{strDetNum;strWav;strTDet;strDur});
    %%% window info string
    strWinSize = sprintf('Time Span = %s',char(rangeDiff));
    strOffset = sprintf('Time Offset = %s',char(seconds(tOffset)));
    strWinInfo = strcat({strWinSize;strOffset},{' '});
    %%% validation string
    strVal = OUTPUT.tableDetections.Class_MATLAB{DATA.currentDet.idx};
    if isempty(strVal)
        strVal = 'Unclassified';
    end

    % set info strings in hidden axes
    PLOT.textDetInfo.String = strDetInfo;
    PLOT.textWindowInfo.String = strWinInfo;
    PLOT.textValidation.String = strVal;
    
    % edit colour of validation string
    PLOT.textValidation.Color = PARAMS.col.(strVal);

    % remove load screen if present
    hideLoadScreen();
    
    % stop profiler (DEBUG ONLY)
    %profile off
    
    % NESTED FUNCTIONS ....................................................
    function drawMarkers(ha,markerType,col2Idx,plotDTRange,plotFRange)
    % Draws the marker patches for a particular signal type.
    % The function reads color data from "PARAMS.markers" automatically 
    % based on "markerType". Color matrices in "PARAMS.markers" should 
    % be either 1-by-3 (1 color) or 2-by-3 (2 colors) - if the current 
    % "markerType" has 2 colors associated with it, then "col2Idx" should 
    % be a N-by-1 logical vector controlling assignment of the 2nd color 
    % (N corresponds to the number of signals). If there is only one color,
    % "col2Idx" should be empty.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOTE: nested functions read global variables from their parent
    % functions, and so need not (and cannot) declare them explicitly.
    
        % initialize
        outTableName = ['table',markerType];
        patchObjName = ['patches',markerType];
        col = PARAMS.markers.Color.(markerType);
    
        % get signal datetime and frequency ranges
        [sigDTRange,sigFRange] = getSignalRanges(OUTPUT.(outTableName));
    
        % get logical indices of markers that appear on-screen
        markerOnScreen = ...
            any((sigDTRange >= plotDTRange(1) & sigDTRange <= plotDTRange(2)),2) &...
            any((sigFRange >= plotFRange(1) & sigFRange <= plotFRange(2)),2);
        nMarkersOnScreen = sum(markerOnScreen);
        
        % set up X and Y data, and initial color
        patchXData = sigDTRange(markerOnScreen,[1,1,2,2])';
        patchYData = sigFRange(markerOnScreen,[1,2,2,1])';
        col1 = col(1,:);
        
        % draw patches (note "fill" returns a patch array)
        hp = fill(ha,patchXData,patchYData,col1,...
            'EdgeColor',col1,...
            'FaceAlpha',PARAMS.markers.FaceAlpha,...
            'LineWidth',PARAMS.markers.LineWidth);
        
        % edit secondary color if specified
        if ~isempty(col2Idx)
            col2 = col(2,:);
            set(hp(col2Idx(markerOnScreen)),{'FaceColor','EdgeColor'},{col2,col2});
        end
        
        % assign handles to PLOT struct
        PLOT.(patchObjName) = hp;
        
        % update screen index
        DATA.markers.ScreenIndices.(markerType)(markerOnScreen) = 1:nMarkersOnScreen;
    end
end

% showLoadScreen ----------------------------------------------------------
function showLoadScreen()
% Turns the objects on the load screen visible
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global PLOT

    % make text and patch visible
    PLOT.patchLoading.Visible = 'on';
    PLOT.textLoading.Visible = 'on';

    % draw
    drawnow
end

% hideLoadScreen ----------------------------------------------------------
function hideLoadScreen()
% Makes the objects on the load screen invisible
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global PLOT

    % hide stuff
    PLOT.patchLoading.Visible = 'off';
    PLOT.textLoading.Visible = 'off';

    % draw
    drawnow
end

% changeState -------------------------------------------------------------
function changeState(newState,doUpdate,varargin)
% Changes program state and initializes the new state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global DATA UI PLOT
    
    hf = UI.fig;
    
    % assign new state
    DATA.program.mode = newState;
    
    % process new state
    switch newState
        case 'normal'
            narginchk(2,2);
            
            % update keyboard and mouse callbacks
            hf.WindowKeyPressFcn = @FigKeyPressFcn_normal;
            hf.WindowButtonDownFcn = '';
            hf.WindowButtonMotionFcn = '';
            hf.WindowButtonUpFcn = '';
            
            % activate buttons (except "Playback"), and set edit mode state
            % buttons to false
            changeButtonState('on');
            UI.buttonStopPlayback.Enable = 'off';
            UI.statebuttonEditMissedCalls.Value = false;
            UI.statebuttonEditAnnotations.Value = false;
            
            % reset draw and edit state data
            DATA.program.drawState = '';
            DATA.program.drawPurpose = '';
            DATA.program.drawType = '';
            DATA.program.undetectedType = '';
            DATA.program.editState = '';
            DATA.program.editType = '';
            
            % clear marker highlight if there is one
            PLOT.highlight.reset();
            
        case 'draw'
            narginchk(4,4);
            disp('Highlight area on plot')
            
            % update keyboard and mouse callbacks
            hf.WindowKeyPressFcn = @FigKeyPressFcn_draw;
            hf.WindowButtonDownFcn = @FigButtonDownFcn_draw;
            hf.WindowButtonMotionFcn = @FigButtonMotionFcn_draw;
            hf.WindowButtonUpFcn = @FigButtonUpFcn_draw;
            
            % disable buttons
            changeButtonState('off');
            %** Note that disabling the source button causes a loss of
            %focus. Consider changing it to state buttons...
            
            % initialize draw state
            DATA.program.drawState = 'start';
            DATA.program.drawPurpose = varargin{1};
            DATA.program.drawType = varargin{2};
            
        case 'playback'
            narginchk(2,2);
            
            % update keyboard and mouse callbacks
            hf.WindowKeyPressFcn = @FigKeyPressFcn_playback;
            hf.WindowButtonDownFcn = '';
            hf.WindowButtonMotionFcn = '';
            hf.WindowButtonUpFcn = '';
            
            % disable buttons except "Stop Playback"
            changeButtonState('off');
            UI.buttonStopPlayback.Enable = 'on';
            
        case 'edit'
            narginchk(3,3);
            
            % set edit state and type
            DATA.program.editState = 'select';
            DATA.program.editType = varargin{1};
            
            % update keyboard and mouse callbacks
            hf.WindowKeyPressFcn = @FigKeyPressFcn_edit;
            hf.WindowButtonDownFcn = @FigButtonDownFcn_edit;
            hf.WindowButtonMotionFcn = @FigButtonMotionFcn_edit;
            hf.WindowButtonUpFcn = @FigButtonUpFcn_edit;
            
            % disable all buttons except appropriate edit button and
            % (most) spectrogram panel buttons
            changeButtonState('off');
            editButtonName = ['statebuttonEdit',DATA.program.editType];
            UI.(editButtonName).Enable = 'on';
            UI.buttonResetTimeSpan.Enable = 'on';
            UI.buttonResetTimeOffset.Enable = 'on';
            UI.buttonResetFreqRange.Enable = 'on';
            UI.buttonChangeColormap.Enable = 'on';
            
            % show markers if they're hidden
            if ~DATA.markers.show
                DATA.markers.show = true;
                doUpdate = true;
            end
                
        otherwise
            error('Bad case')
    end
    
    % update plot if specified
    if doUpdate
        showLoadScreen();
        updatePlot
    end
    
    % NESTED FUNCTIONS 
    % changeButtonState ...................................................
    function changeButtonState(state)
    % Change the state of all UIButtons and UIStateButtons
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hButtons = findobj(UI.fig,'Type','uibutton','-or','Type','uistatebutton');
        set(hButtons,'Enable',state);
    end
    
end

% saveSpreadsheet ---------------------------------------------------------
function saveStatus = saveSpreadsheet()
% Updates the contents of the output spreadsheet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global OUTPUT DATA

    % initialize save status
    saveStatus = zeros(3,1); %%% one value for each sheet: 1 = success, 0 = no save attempt, -1 = failed
    
    % save first sheet if it's been updated
    if DATA.program.sheetChanged.Detections
        outTable1 = OUTPUT.tableDetections;
        try
            writetable(outTable1,OUTPUT.path,'Sheet','Detected');
            DATA.program.sheetChanged.Detections = false;
            disp('"Detected" sheet saved successfully')
            saveStatus(1) = 1;
        catch ME1
            warning('Could not save "Detected" sheet: %s',ME1.message) %%% ignore CodeAnalyzer warning, it's incorrect
            saveStatus(1) = -1;
        end
    end
    
    % save second sheet if it's been updated
    if DATA.program.sheetChanged.MissedCalls
        saveStatus(2) = processAuxillarySheet('Undetected',OUTPUT.tableMissedCalls,OUTPUT.path);
        if saveStatus(2) == 1
            DATA.program.sheetChanged.MissedCalls = false;
        end
    end
    
    % save third sheet if it's been updated
    if DATA.program.sheetChanged.Annotations
        saveStatus(3) = processAuxillarySheet('Annotations',OUTPUT.tableAnnotations,OUTPUT.path);
        if saveStatus(3) == 1
            DATA.program.sheetChanged.Annotations = false;
        end
    end
 
    
    % NESTED FUNCTIONS
    % processAuxillarySheet ...............................................
    function saveStatusSub = processAuxillarySheet(sheetName,sheetDataNew,outPath)
    % Try saving Undetected or Annotations sheet
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        sheetDataOld = readtable(outPath,'Sheet',sheetName);
        nRowsNew = height(sheetDataNew);
        nRowsOld = height(sheetDataOld);

        if nRowsNew >= nRowsOld
            outTable = sheetDataNew;
        else
            %workaround to ensure old excel rows are deleted
            sheetDataCellNew = cell(size(sheetDataOld));
            sheetDataCellNew(1:nRowsNew,:) = table2cell(sheetDataNew);
            outTable = cell2table(sheetDataCellNew,'VariableNames',sheetDataNew.Properties.VariableNames);
        end
        try
            writetable(outTable,outPath,'Sheet',sheetName);
            fprintf('"%s" sheet saved successfully\n',sheetName)
            saveStatusSub = 1;
        catch MESub
            warning('Could not save "%s" sheet: %s',sheetName,MESub.message)
            saveStatusSub = -1;
        end
    end
end

% do_GotoNext -------------------------------------------------------------
function do_GotoNext()
% Go to next detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global DATA

    if DATA.currentDet.idx < DATA.detInfo.n
        DATA.currentDet.idx = DATA.currentDet.idx + 1;
        processNewDet();
    else
        disp('Last detection reached!')
    end
end

% do_GotoPrev -------------------------------------------------------------
function do_GotoPrev()
% Go to previous detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global DATA

    if DATA.currentDet.idx > 1
        DATA.currentDet.idx = DATA.currentDet.idx - 1;
        processNewDet();
    else
        disp('First detection reached!')
    end
end

% do_Jump -----------------------------------------------------------------
function do_Jump()
% Select a detection to go to
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global OUTPUT DATA UI

    % generate detection list strings
    numStr = strtrim(cellstr(num2str((1:DATA.detInfo.n)')));
    dateStr = cellstr(DATA.wavInfo.dtStarts(DATA.detInfo.iWav) + seconds(DATA.detInfo.startTimes));
    classStr = OUTPUT.tableDetections.Class_MATLAB;
    detList = strcat('(', numStr, {')     '}, dateStr, {'     '},classStr);

    % create list dialogue
    [iDetNew,detSelected] = listdlg(...
        'ListString',detList,...
        'PromptString','Choose detection',...
        'SelectionMode','single',...
        'ListSize',[250, 300],...
        'InitialValue',DATA.currentDet.idx);
    figure(UI.fig);

    if detSelected
        DATA.currentDet.idx = iDetNew;
        processNewDet();
    end
end

% do_GotoUnclassified -----------------------------------------------------
function doJump = do_GotoUnclassified()
% Jump to the first unvalidated detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global DATA OUTPUT
    
    % get first unvalidated detection
    unclassified = strcmp(OUTPUT.tableDetections.Class_MATLAB,'');
    iFirstUnclassified = find(unclassified,1,'first');
    
    % change index
    if ~isempty(iFirstUnclassified)
        doJump = true;
        DATA.currentDet.idx = iFirstUnclassified;
        processNewDet();
    else
        doJump = false;
        disp('All detections classified!')
    end
end

% do_ZoomIn ---------------------------------------------------------------
function do_ZoomIn()
% Reduces the span of the t axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global PARAMS DATA

    % reduce tSpan
    DATA.plot.tSpan = DATA.plot.tSpan/PARAMS.plot.TZoomFactor;

    % update plot
    showLoadScreen();
    updatePlot();
end

% do_ZoomOut --------------------------------------------------------------
function do_ZoomOut()
% Increases the span of the t axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global PARAMS DATA

    % increase tSpan
    DATA.plot.tSpan = DATA.plot.tSpan*PARAMS.plot.TZoomFactor;

    % update plot
    showLoadScreen();
    updatePlot();
end

% do_PanLeft --------------------------------------------------------------
function do_PanLeft()
% Translates the plot to the left
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global PARAMS DATA UI

    % get span of t axis and edit offset based on that
    tRangeDur = diff(UI.axesMain.XLim);
    DATA.plot.tOffset = DATA.plot.tOffset - seconds(tRangeDur)*PARAMS.plot.TPanFactor;

    % update plot
    showLoadScreen();
    updatePlot();
end

% do_PanRight -------------------------------------------------------------
function do_PanRight()
% Translates the plot to the right
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global PARAMS DATA UI

    % get span of t axis and edit offset based on that
    tRangeDur = diff(UI.axesMain.XLim);
    DATA.plot.tOffset = DATA.plot.tOffset + seconds(tRangeDur)*PARAMS.plot.TPanFactor;

    % update plot
    showLoadScreen();
    updatePlot();
end

% do_RaiseFMax ------------------------------------------------------------
function do_RaiseFMax()
% Increases fMax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global PARAMS DATA

    % change fRange
    fMaxOld = DATA.plot.fMax;
    if fMaxOld < DATA.currentDet.wav.Fs/2
        DATA.plot.fMax = fMaxOld*PARAMS.plot.FZoomFactor;
        
        % update plot
        showLoadScreen();
        updatePlot();
    else
        disp('Cannot increase frequency range any further')
    end
end

% do_LowerFMax ------------------------------------------------------------
function do_LowerFMax()
% Lowers fMax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global PARAMS DATA

    % change fRange
    DATA.plot.fMax = DATA.plot.fMax/PARAMS.plot.FZoomFactor;
        
    % update plot
    showLoadScreen();
    updatePlot();
end

% do_ResetPlot ------------------------------------------------------------
function do_ResetPlot(type)
% Resets the zoom, offset, and/or fRange
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global PARAMS DATA

    % reset parameters according to reset type
    switch type
        case 'span'
            DATA.plot.tSpan = PARAMS.plot.TSpanDefault;
        case 'offset'
            DATA.plot.tOffset = 0;
        case 'frange'
            DATA.plot.fMax = PARAMS.plot.FMaxDefault;
        case 'all'
            DATA.plot.tSpan = PARAMS.plot.TSpanDefault;
            DATA.plot.tOffset = 0;
            DATA.plot.fMax = PARAMS.plot.FMaxDefault;
        otherwise
            error('Bad case')
    end

    % update plot
    showLoadScreen();
    updatePlot();
end

% do_ToggleMarkers --------------------------------------------------------
function do_ToggleMarkers()
% Changes the flag to show or hide detection patches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global DATA

    % change flag
    DATA.markers.show = ~DATA.markers.show;

    % update plot
    showLoadScreen();
    updatePlot();
end

% do_ChangeColormap  ------------------------------------------------------
function do_ChangeColormap()
% Changes the colormap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global PARAMS DATA UI

    % define list of colormaps
    colmapList = fieldnames(PARAMS.plot.Colormaps);

    % get current colormap
    iMapCurrent = find(strcmp(DATA.plot.colmapActive,colmapList));

    % prompt user to choose map
    opt = listdlg(...
        'ListString',colmapList,...
        'PromptString','Choose colormap',...
        'Name','colormap',...
        'SelectionMode','single',...
        'InitialValue',iMapCurrent);
    figure(UI.fig);

    % change colormap
    if ~isempty(opt) && opt ~= iMapCurrent
        DATA.plot.colmapActive = colmapList{opt};

        % update plot
        showLoadScreen();
        updatePlot();
    end
end

% do_ResetValidation ------------------------------------------------------
function do_ResetValidation()
% Resets the classification of a detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global OUTPUT PARAMS DATA UI PLOT
    
    % update relevant fields in table
    OUTPUT.tableDetections.Class_MATLAB{DATA.currentDet.idx} = '';
    OUTPUT.tableDetections.ReasonForUNK{DATA.currentDet.idx} = '';
    
    % clear and deactivate reason for unknown text area
    UI.textareaReasonForUnknown.Value = {''};
    UI.textareaReasonForUnknown.Enable = 'off';
    UI.labelReasonForUnknown.Enable = 'off';
    
    % Change validation text
    PLOT.textValidation.String = 'Unclassified';
    PLOT.textValidation.Color = PARAMS.col.Unclassified;
    drawnow
    
    % mark table change status
    DATA.program.sheetChanged.Detections = true;
end

% do_Validate -------------------------------------------------------------
function do_Validate(type)
% Sets the classification of a detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global OUTPUT PARAMS DATA UI PLOT
    
    % get index of current detection
    iDet = DATA.currentDet.idx;
    
    % update relevant field in table
    OUTPUT.tableDetections.Class_MATLAB{iDet} = type;
    
    % add reason for unknown
    if strcmp(type,'Unknown')
        % get current reason for unknown
        currentRFU = OUTPUT.tableDetections.ReasonForUNK(iDet);

        % prompt user for reason and process
        RFU = inputdlg('Enter reason for unknown:','Reason for Unknown',[1 80],currentRFU);
        figure(UI.fig);
        if isempty(RFU)
            RFU = {''};
        end
        
        % activate text area
        UI.textareaReasonForUnknown.Enable = 'on';
        UI.labelReasonForUnknown.Enable = 'on';
        
    else
        RFU = {'N/A'};
        UI.textareaReasonForUnknown.Enable = 'off';
        UI.labelReasonForUnknown.Enable = 'off';
    end
    
    % update reason for unknown
    UI.textareaReasonForUnknown.Value = RFU;
    OUTPUT.tableDetections.ReasonForUNK(iDet) = RFU;
    
    % Change validation text
    col = PARAMS.col.(type);
    PLOT.textValidation.String = type;
    PLOT.textValidation.Color = col;
    drawnow
    
    % mark table change status
    DATA.program.sheetChanged.Detections = true;
end

% do_AddComment -----------------------------------------------------------
function do_AddComment()
% Enters a comment in the spreadsheet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global OUTPUT DATA UI
    
    iDet = DATA.currentDet.idx;
    
    % get current comment
    currentComment = OUTPUT.tableDetections.Comments(iDet);

    % prompt user for comment
    comment = inputdlg('Enter comment:','Comment',[1 80],currentComment);
    figure(UI.fig);
    
    % process user input
    if ~isempty(comment)
        UI.textareaComment.Value = comment;
        OUTPUT.tableDetections.Comments(iDet) = comment;
        
        % mark table change status
        DATA.program.sheetChanged.Detections = true;
    end
end

% do_EditComment_NonDetection ---------------------------------------------
function do_EditComment_NonDetection()
% Add or edit comment of a missed call or annotation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global PLOT OUTPUT DATA UI
    
    iSig = PLOT.highlight.signalIdx;
    sigType = PLOT.highlight.signalType;
    outTableName = ['table',sigType];
    sheetStatVarName = [lower(sigType(1)),sigType(2:end),'SheetChanged'];
    
    % get current comment
    currentComment = OUTPUT.(outTableName).Comments(iSig);

    % prompt user for comment
    comment = inputdlg('Enter comment:','Comment',[1 80],currentComment);
    figure(UI.fig);
    
    % process user input
    if ~isempty(comment)
        OUTPUT.(outTableName).Comments(iSig) = comment;
        
        % mark table change status
        DATA.program.(sheetStatVarName) = true;
    end
end

% do_AddMissedCall --------------------------------------------------------
function do_AddMissedCall(dtRange)
% Add missed call to Undetected spreadsheet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % checkout global variables
    global OUTPUT DATA UI
    
    % prompt user for comment
    comment = inputdlg(sprintf('Enter comment for new %s missed call\n(Or press Cancel to stop marking this area)',lower(DATA.program.undetectedType)),'Missed Call',[1 80]);
    figure(UI.fig)
    
    % process user input
    if ~isempty(comment)
    
        % show load screen
        showLoadScreen();

        % round dtRange
        %dtRange = dateshift(dtRange,'start','second','nearest');

        % create table entry for current call
        FileName = {DATA.currentDet.wav.Name};
        FileStart = seconds(DATA.currentDet.wav.dtStart - DATA.wavInfo.dtStartRef);
        SigStart = seconds(dtRange(1) - DATA.currentDet.wav.dtStart);
        SigEnd = seconds(dtRange(2) - DATA.currentDet.wav.dtStart);
        SigStartDateTime = dtRange(1);
        Type = {DATA.program.undetectedType};
        Comments = comment;
        tableAddition = table(FileName,FileStart,SigStart,SigEnd,SigStartDateTime,Type,Comments,...
            'VariableNames',{'FileName','FileStart','SigStart','SigEnd','SigStartDateTime','Type','Comments'});

        % create new table
        tableNewUnsorted = [OUTPUT.tableMissedCalls;tableAddition];

        % sort table
        sigStarts = tableNewUnsorted.FileStart + tableNewUnsorted.SigStart;
        [~,iSort] = sort(sigStarts,'ascend');
        tableNew = tableNewUnsorted(iSort,:);
        OUTPUT.tableMissedCalls = tableNew;

        % save spreadsheet
        DATA.program.sheetChanged.MissedCalls = true;
        saveSpreadsheet();

        fprintf('%s call added to "Undetected" sheet\n',DATA.program.undetectedType)
    end
    
    % change program state
    changeState('normal',true)
end

% do_AddAnnotation --------------------------------------------------------
function do_AddAnnotation(dtRange,fRange)
% Add an annotation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % checkout global variables
    global OUTPUT DATA UI
    
    % prompt user to add comment or cancel
    comment = inputdlg(sprintf('Enter comment for new annotation\n(Or press Cancel to stop marking this area)'),'New Annotation',[1 80]);
    figure(UI.fig);
    
    % process user input
    if ~isempty(comment)
        
        % show load screen
        showLoadScreen();
        
        % create table entry for current annotation
        FileName = {DATA.currentDet.wav.Name};
        FileStart = seconds(DATA.currentDet.wav.dtStart - DATA.wavInfo.dtStartRef);
        SigStart = seconds(dtRange(1) - DATA.currentDet.wav.dtStart);
        SigEnd = seconds(dtRange(2) - DATA.currentDet.wav.dtStart);
        SigMinFreq = fRange(1);
        SigMaxFreq = fRange(2);
        SigStartDateTime = dtRange(1);
        Comments = comment;
        tableAddition = table(FileName,FileStart,SigStart,SigEnd,SigMinFreq,SigMaxFreq,SigStartDateTime,Comments,...
            'VariableNames',{'FileName','FileStart','SigStart','SigEnd','SigMinFreq','SigMaxFreq','SigStartDateTime','Comments'});
        
        % create new table
        tableNewUnsorted = [OUTPUT.tableAnnotations;tableAddition];
        
        % sort table
        SigStarts = tableNewUnsorted.FileStart + tableNewUnsorted.SigStart;
        [~,iSort] = sort(SigStarts,'ascend');
        tableNew = tableNewUnsorted(iSort,:);
        OUTPUT.tableAnnotations = tableNew;
        
        % save spreadsheet
        DATA.program.sheetChanged.Annotations = true;
        saveSpreadsheet(); 
        
        disp('Annotation added to "Annotations" sheet')
    end
    
    % change state
    changeState('normal',true)
end

% do_UnmarkNonDetection ---------------------------------------------------
function do_UnmarkNonDetection()
% Removes the currently highlighted MissedCall or Annotation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global DATA PLOT OUTPUT
    
    % extract useful data
    hh = PLOT.highlight;
    
    if hh.isDisplayed
        % get type and index of currently selected signal
        sigType = hh.signalType;
        sigIdx = hh.signalIdx;
        
        % set vars according to signal type
        switch sigType
            case 'MissedCalls'
                promptStr = 'Unmark this area as a missed call?';
                outTableName = 'tableMissedCalls';
            case 'Annotations'
                promptStr = 'Remove this annotation?';
                outTableName = 'tableAnnotations';
            otherwise
                error('Bad case')
        end

        % prompt user for action
        opt = questdlg(promptStr,'Unmark Non-Detection',...
            'Yes','No','Yes');
        if strcmp(opt,'Yes')
            showLoadScreen();

            % remove selection
            OUTPUT.(outTableName)(sigIdx,:) = [];
            DATA.program.sheetChanged.(sigType) = true;
            
            % clear highlight
            hh.reset();

            % update plot
            updatePlot();

            fprintf('Removed entry from "%s" sheet\n',sigType)
        end
    end
end

% do_HighlightMarker ------------------------------------------------------
function do_HighlightMarker(x,y,markerType)
% Change the appearance of a marker located at a certain position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global PLOT DATA
    
    % look for markers that overlap with cursor
    patchObjName = ['patches',markerType];
    hp = PLOT.(patchObjName);
    if ~isempty(hp)
        markerXData = [hp.XData];
        markerYData = [hp.YData];
        markerXRanges = markerXData([1,3],:)';
        markerYRanges = markerYData([1,2],:)';
        markerHit_X = x <= markerXRanges(:,2) & x >= markerXRanges(:,1);
        markerHit_Y = y <= markerYRanges(:,2) & y >= markerYRanges(:,1);
        markerHit = markerHit_X & markerHit_Y;
    else
        markerHit = [];
    end
    %nMarkers = numel(markerHit);
    nHits = sum(markerHit);
    
    % get screen index of previous highlighted marker, if any, and clear 
    % the highlight
    iHitScreenOld = PLOT.highlight.screenIdx;
    PLOT.highlight.reset();
    
    % highlight a marker if there are any hits
    if any(markerHit)
        % check the hit stack to see if previously highlighted marker (if
        % any) is still within. If so, move to the next hit in the stack;
        % otherwise, select the first hit in the stack.
        iHitScreenAll = find(markerHit);
        iHitScreenOldIdx = find(iHitScreenAll == iHitScreenOld);
        if isempty(iHitScreenOldIdx)
            % first hit
            iHitScreenNewIdx = 1;
        else
            % next in stack
            iHitScreenNewIdx = iHitScreenOldIdx + 1;
            if iHitScreenNewIdx > nHits
                iHitScreenNewIdx = 1;
            end
        end
        
        % highlight new marker
        iHitScreen = iHitScreenAll(iHitScreenNewIdx);
        iHitList = find(DATA.markers.ScreenIndices.(markerType) == iHitScreen,1,'first'); % expect only one hit
        hp = PLOT.(patchObjName)(iHitScreen);
        PLOT.highlight.setHighlight(hp, markerType, iHitList, iHitScreen);
    end
end

% do_PlaySaveAudio --------------------------------------------------------
function do_PlaySaveAudio(dtRange,action)
% Play or save audio clip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global DATA
    
    % get Fs and bits per sample
    Fs = DATA.currentDet.wav.Fs;
    bps = DATA.currentDet.wav.BitsPerSample;

    % isolate area to play/save
    tRange = seconds(dtRange - DATA.currentDet.wav.dtStart);
    sampleRangeRaw = round(tRange*Fs);
    sampleRange = zeros(1,2);
    sampleRange(1) = max([1,sampleRangeRaw(1)]);
    sampleRange(2) = min([numel(DATA.currentDet.wav.x),sampleRangeRaw(2)]);
    sampleStartOffset = -(sampleRangeRaw(1) - sampleRange(1));
    x = DATA.currentDet.wav.x(sampleRange(1):sampleRange(2));
    
    % do playback or save
    if ~isempty(x)
        switch action
            case 'play' 
            
                % set playback variables
                FsPlay = Fs*DATA.playback.speed;
                xPlay = x*DATA.playback.vol;
                
                try
                    % create audio player object
                    DATA.playback.player = audioplayer(xPlay,FsPlay,bps);

                    % assign player start time to user data
                    DATA.playback.player.UserData = struct();
                    DATA.playback.player.UserData.dtStart = dtRange(1) + seconds(sampleStartOffset);
                    DATA.playback.player.UserData.FsActual = Fs;

                    % assign audio player callbacks
                    DATA.playback.player.StartFcn = @PlaybackStartFcn;
                    DATA.playback.player.TimerFcn = @PlaybackTimerFcn;
                    DATA.playback.player.StopFcn = @PlaybackStopFcn;

                    % play audio 
                    play(DATA.playback.player);
                    % Note player will take care of state change after
                    % completing playback
                catch ME
                    warning('Audio playback failed.\nThe sampling rate and/or playback speed may be too high or too low.\nOriginal error message:\n\t%s', ME.message)
                    
                    % manually trigger post-playback actions
                    PlaybackStopFcn(DATA.playback.player,[])
                end
                
            case 'save'
            % Save audio clip
            
                % prompt user for output file path
                [clipName,clipDir] = uiputfile('*.wav','Save audio clip');
                clipPath = fullfile(clipDir,clipName);
                if ischar(clipName)
                    % save clip
                    audiowrite(clipPath,x,Fs,'BitsPerSample',bps)
                    fprintf('Audio clip saved to file:\n%s\n',clipPath);
                end
                
                % revert to normal state
                changeState('normal',false)
            
            otherwise
                error('Bad case')
        end
    else
        warning('No audio data to %s!',action)
        changeState('normal',false)
    end
end

% do_StopPlayback ---------------------------------------------------------
function do_StopPlayback()
% stop audio playback
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global DATA
    
    % stop playback
    stop(DATA.playback.player);
end


%% CALLBACK FUNCTIONS =====================================================
% FigDeleteFcn ------------------------------------------------------------
function FigDeleteFcn(src,edata)
% Saves spreadsheet and destroys global variables at same time as figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % save spreadsheet
    saving = true;
    while saving
        saveStatus = saveSpreadsheet();
        
        if any(saveStatus == -1)
            OKPressed = false;
            warnStr = sprintf('Spreadsheet file did not save successfully.\n\nIf the file is open, close it and press OK to try again.\n\nOtherwise, press X to cancel (unsaved data will be lost).');
            hfWarn = warndlg(warnStr,'Warning: save failed','modal');
            
            % overwrite callback function for OK button
            hpbOK = findobj(hfWarn,'Tag','OKButton');
            hpbOK.Callback = @warnOK_Callback;
            
            % wait for user action
            waitfor(hfWarn)
            
            % check if OK button was pressed. If not, cancel the save.
            if ~OKPressed
                saving = false;
                disp('Save canceled')
            end
        else
            saving = false;
        end
    end

    % clear global variables
    clear global OUTPUT PARAMS DATA UI PLOT
    
    % view profiler results (DEBUG ONLY)
    %profile viewer
    
    % nested function: OK button callback .................................
    function warnOK_Callback(src,edata)
        OKPressed = true;
        delete(gcbf)
    end
end

% FigKeyPressFcn_normal ---------------------------------------------------
function FigKeyPressFcn_normal(src,edata)
% Executes actions based on keyboard key presses for the 'normal' mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES: edata in this case will be a KeyData object, which 
% contains info on which key was pressed.

    % checkout global variables
    global UI

    % get name of key pressed
    key = edata.Key;

    % check if shift was pressed too
    noShift = ~any(contains(edata.Modifier,'shift'));

    % take action depending on key
    switch key
        case 'w'
            if noShift
                do_ZoomIn();
            else
                do_RaiseFMax();
            end
        case 's'
            if noShift
                do_ZoomOut();
            else
                do_LowerFMax();
            end
        case 'a'
            if noShift
                do_GotoPrev();
            else
                do_PanLeft();
            end
        case 'd'
            if noShift
                do_GotoNext();
            else
                do_PanRight();
            end
        case 'uparrow'
            do_RaiseFMax();
        case 'downarrow'
            do_LowerFMax();
        case 'rightarrow'
            do_PanRight();
        case 'leftarrow'
            do_PanLeft();
        case 'space'
            %do_Jump(); %%% doesn't work well with UIFigure
        case 'r'
            do_ResetPlot('all');
        case 't'
            do_ToggleMarkers();
        case 'c'
            do_ChangeColormap();
        case 'q'
            close(UI.fig)
    end
end

% FigKeyPressFcn_draw -----------------------------------------------------
function FigKeyPressFcn_draw(src,edata)
% Executes actions based on keyboard key presses for the 'draw' mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global DATA
    
    % get name of key pressed and draw state
    key = edata.Key;
    drawState = DATA.program.drawState;
    
    % take action depending on key and draw state
    if strcmp(key,'escape') && strcmp(drawState,'start')
        changeState('normal',false)
    end
end

% FigKeyPressFcn_playback -------------------------------------------------
function FigKeyPressFcn_playback(src,edata)
% Executes actions based on keyboard key presses for the 'playback' mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get name of key pressed
    key = edata.Key;
    
    % take action depending on key and draw state
    if strcmp(key,'escape')
        do_StopPlayback();
    end
end

% FigKeyPressFcn_edit -----------------------------------------------------
function FigKeyPressFcn_edit(src,edata)

    % get name of key pressed
    key = edata.Key;
    
    % check if shift was pressed too
    noShift = ~any(contains(edata.Modifier,'shift'));
    
    % take action depeding on key press
    switch key
        case {'escape','return'}
            saveSpreadsheet();
            changeState('normal',true)
        case 'w'
            if noShift
                do_ZoomIn();
            else
                do_RaiseFMax();
            end
        case 's'
            if noShift
                do_ZoomOut();
            else
                do_LowerFMax();
            end
        case 'a'
            if noShift
                %** consider allowing selection cycling
            else
                do_PanLeft();
            end
        case 'd'
            if noShift
                %** consider allowing selection cycling
            else
                do_PanRight();
            end
        case 'uparrow'
            do_RaiseFMax(); 
        case 'downarrow'
            do_LowerFMax();
        case 'rightarrow'
            do_PanRight();
        case 'leftarrow'
            do_PanLeft();
        case 'space'
            do_EditComment_NonDetection();
        case {'x','delete'}
            do_UnmarkNonDetection();
        case 'r'
            do_ResetPlot('all');
        case 'c'
            do_ChangeColormap();
        case 'q'
            close(UI.fig)
    end
end

% FigButtonDownFcn_draw ---------------------------------------------------
function FigButtonDownFcn_draw(src,edata)
% Creates selection patch object if mouse was clicked over plot area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global DATA UI PLOT
    
    % check mouse clicks if app is in correct state
    if strcmp(DATA.program.drawState,'start')
        ha = UI.axesMain;
        
        % get datetime-numeric ruler
        ruler = ha.XAxis;
        
        % get position of mouse cursor relative to axes
        cp = ha.CurrentPoint;
        x = cp(1,1);
        y = cp(1,2);
        
        % determine if cursor position is within axes
        xMin = ruler2num(ha.XLim(1),ruler);
        xMax = ruler2num(ha.XLim(2),ruler);
        yMin = ha.YLim(1);
        yMax = ha.YLim(2);

        withinX = x <= xMax && x >= xMin;
        withinY = y <= yMax && y >= yMin;
        if withinX && withinY
            % create patch object
            xDT = num2ruler(x,ruler);
            xPatch = [xDT;xDT;xDT;xDT];
            switch DATA.program.drawType
                case 'x'
                    yPatch = [yMin;yMax;yMax;yMin]; 
                case 'xy'
                    yPatch = [y;y;y;y]; 
                otherwise
                    error ('Bad case')
            end
            
            %{
            PLOT.patchSelection = fill(ha,xPatch,yPatch,[0.65 0.8 0.9],...
                'EdgeColor',[0.3 0.4 0.45],...
                'FaceAlpha',0.25,...
                'EdgeAlpha',0.25);
            %}
            
            PLOT.patchSelection = fill(ha,xPatch,yPatch,[0.65 0.8 0.9],...
                'EdgeColor',[0.65 0.8 0.9],...
                'FaceAlpha',0.3,...
                'EdgeAlpha',0.5);

            % change draw state
            DATA.program.drawState = 'drag';
        end
    end
end

% FigButtonDownFcn_edit ---------------------------------------------------
function FigButtonDownFcn_edit(src,edata)
% Pass position of cursor upon mouse click for marker highlighting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global UI DATA PLOT
    
    % check if there exists a highlighted marker with edge hits.
    % If YES, then entre resize mode.
    % If NO, then check if a new marker can be highlighted.
    if any(PLOT.highlight.lineData.Hit)
        DATA.program.editState = 'resize';
    else
        % get axis and datetime-numeric ruler
        ha = UI.axesMain;
        ruler = ha.XAxis;

        % get position of mouse cursor relative to axes
        cp = ha.CurrentPoint;
        x = cp(1,1);
        y = cp(1,2);

        % determine if cursor position is within axes
        xMin = ruler2num(ha.XLim(1),ruler);
        xMax = ruler2num(ha.XLim(2),ruler);
        yMin = ha.YLim(1);
        yMax = ha.YLim(2);

        withinX = x <= xMax && x >= xMin;
        withinY = y <= yMax && y >= yMin;
        if withinX && withinY
            % get marker type
            markerType = DATA.program.editType;

            % highlight marker if possible
            do_HighlightMarker(x,y,markerType);
        end
    end
end

% FigButtonMotionFcn_draw -------------------------------------------------
function FigButtonMotionFcn_draw(src,edata)
% Updates the area of the selection patch as the mouse moves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global DATA UI PLOT
    
    % update patch if app is in correct state
    if strcmp(DATA.program.drawState,'drag')
        ha = UI.axesMain;
        hp = PLOT.patchSelection;
        
        % get position of mouse cursor relative to axes
        cp = ha.CurrentPoint;
        x = cp(1,1);
        y = cp(1,2);
        
        % edit patch object
        switch DATA.program.drawType
            case 'x'
                hp.XData = [hp.XData(1:2);x;x];
            case 'xy'
                hp.XData = [hp.XData(1:2);x;x];
                hp.YData = [hp.YData(1);y;y;hp.YData(4)];
            otherwise
                error('Bad case')
        end
    end
end

% FigButtonMotionFcn_edit -------------------------------------------------
function FigButtonMotionFcn_edit(src,edata)
% Depending on sub-state, checks if the mouse cursor is over the edge of a 
% marker, or enables resizing of marker.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global DATA UI PLOT
    
    % check if something is highlighted
    if PLOT.highlight.isDisplayed
        
        % get position of mouse cursor relative to axes
        ha = UI.axesMain;
        cp = ha.CurrentPoint;
        x = cp(1,1);
        y = cp(1,2);
    
        % If mouse is down over edge, resize marker. 
        % Otherwise, check for edge hits.
        switch DATA.program.editState
            case 'resize'
                PLOT.highlight.resize(x,y);
            case 'select'
                PLOT.highlight.checkMouseHits(x,y);
            otherwise
                error('Bad case')
        end  
    end
end

% FigButtonUpFcn_draw -----------------------------------------------------
function FigButtonUpFcn_draw(src,edata)
% Finalizes selection patch object position and ends draw state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global DATA UI PLOT
    
    % end patch dragging if app is in correct state
    if strcmp(DATA.program.drawState,'drag')
        ha = UI.axesMain;
        hp = PLOT.patchSelection;
        DATA.program.drawState = 'end';
        
        % get patch x-coordinates (in datetime)
        ruler = ha.XAxis;
        xPatchRange = sort(num2ruler(hp.XData([1,3]),ruler))';
        yPatchRange = sort(hp.YData([1,2]))';
        
        % process according to purpose
        switch DATA.program.drawPurpose
            
            case 'play_audio'
                do_PlaySaveAudio(xPatchRange,'play');
                
            case 'save_audio'
                do_PlaySaveAudio(xPatchRange,'save');
                delete(hp);
                
            case 'mark_undetected'
                do_AddMissedCall(xPatchRange);
                delete(hp);
                
            case 'annotate'
                do_AddAnnotation(xPatchRange,yPatchRange);
                delete(hp);
                
            otherwise
                error('Bad case')
        end
    end
end

% FigButtonUpFcn_edit -----------------------------------------------------
function FigButtonUpFcn_edit(src,edata)
% Finalizes marker resize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % checkout global variables
    global DATA OUTPUT PLOT
    
    % If in resize state, end resize and update dimensions in spreadsheet
    if strcmp(DATA.program.editState,'resize')
        DATA.program.editState = 'select';
        
        iHighlight = PLOT.highlight.signalIdx;
        SigStartDateTimeNew = PLOT.highlight.DTMin;
        SigStartNew = seconds(PLOT.highlight.DTMin - DATA.currentDet.wav.dtStart);
        SigEndNew = seconds(PLOT.highlight.DTMax - DATA.currentDet.wav.dtStart);
        switch DATA.program.editType
            case 'MissedCalls'
                OUTPUT.tableMissedCalls.SigStartDateTime(iHighlight) = SigStartDateTimeNew;
                OUTPUT.tableMissedCalls.SigStart(iHighlight) = SigStartNew;
                OUTPUT.tableMissedCalls.SigEnd(iHighlight) = SigEndNew;
                DATA.program.sheetChanged.MissedCalls = true;
            case 'Annotations'
                SigMinFreqNew = PLOT.highlight.FMin;
                SigMaxFreqNew = PLOT.highlight.FMax;
                OUTPUT.tableAnnotations.SigStartDateTime(iHighlight) = SigStartDateTimeNew;
                OUTPUT.tableAnnotations.SigStart(iHighlight) = SigStartNew;
                OUTPUT.tableAnnotations.SigEnd(iHighlight) = SigEndNew;
                OUTPUT.tableAnnotations.SigMinFreq(iHighlight) = SigMinFreqNew;
                OUTPUT.tableAnnotations.SigMaxFreq(iHighlight) = SigMaxFreqNew;
                DATA.program.sheetChanged.Annotations = true;
            otherwise
                error('Bad case')
        end
    end
end

% buttonResetTimeSpan_Callback --------------------------------------------
function buttonResetTimeSpan_Callback(src,edata)
    do_ResetPlot('span');
end

% buttonResetTimeOffset_Callback ------------------------------------------
function buttonResetTimeOffset_Callback(src,edata) 
    do_ResetPlot('offset');
end

% buttonResetFreqRange_Callback -------------------------------------------
function buttonResetFreqRange_Callback(src,edata)
    do_ResetPlot('frange');
end

% buttonChangeColormap_Callback -------------------------------------------
function buttonChangeColormap_Callback(src,edata)
    do_ChangeColormap();
end

% buttonToggleMarkers_Callback --------------------------------------------
function buttonToggleMarkers_Callback(src,edata)
    do_ToggleMarkers();
end

% buttonJump_Callback -----------------------------------------------------
function buttonJump_Callback(src,edata)
    do_Jump();
end

% buttonFirstUnclassified_Callback ----------------------------------------
function buttonFirstUnclassified_Callback(src,edata)
    do_GotoUnclassified();
end

% buttonPlayDetection_Callback --------------------------------------------
function buttonPlayDetection_Callback(src,edata)

    % checkout global variables
    global DATA OUTPUT
    
    % get time range of current detection
    iDet = DATA.currentDet.idx;
    [dtRangeAll,~] = getSignalRanges(OUTPUT.tableDetections);
    dtRange = dtRangeAll(iDet,:);
    
    % do playback
    do_PlaySaveAudio(dtRange,'play');
end

% buttonPlayRange_Callback ------------------------------------------------
function buttonPlayRange_Callback(src,edata)
    changeState('draw',false,'play_audio','x')
end

% buttonChangeSpeedVol_Callback -------------------------------------------
function buttonChangeSpeedVol_Callback(src,edata)

    % checkout global variables
    global DATA UI
    
    % prompt user for input
    playSettingsCellCurrent = {num2str(DATA.playback.speed);num2str(DATA.playback.vol)};
    playSettingsCell = inputdlg({'Enter playback speed:';'Enter playback volume:'},'Speed/Volume',1,playSettingsCellCurrent);
    figure(UI.fig);
    
    % process input
    playSettings = str2double(playSettingsCell);
    if ~isempty(playSettings)
        % validate play settings
        validSettings = ~any(isnan(playSettings)) && all(playSettings > 0);
        if validSettings
            DATA.playback.speed = playSettings(1);
            DATA.playback.vol = playSettings(2);
        else
            warning('Invalid speed and/or volume')
        end
    end
end

% buttonSaveClip_Callback -------------------------------------------------
function buttonSaveClip_Callback(src,edata)
    changeState('draw',false,'save_audio','x')
end

% buttonStopPlayback_Callback ---------------------------------------------
function buttonStopPlayback_Callback(src,edata)
    do_StopPlayback()
end

% buttonResetValidation_Callback ------------------------------------------
function buttonResetValidation_Callback(src,edata)
    do_ResetValidation();
end

% buttonValidation_Callback -----------------------------------------------
function buttonValidation_Callback(src,edata)
    do_Validate(src.Text);
end

% buttonEditComment_Callback ----------------------------------------------
function buttonEditComment_Callback(src,edata)
    do_AddComment();
end

% buttonAddMissedCall_Callback --------------------------------------------
function buttonAddMissedCall_Callback(src,edata)

    % checkout global variables
    global DATA
    
    % activate patch windows if they're not active
    if ~DATA.markers.show
        showLoadScreen();
        DATA.markers.show = true;
        updatePlot();
    end
    
    % set type of undetected call
    typeCell = regexp(src.Text,'(?<=Add )\w*','match');
    DATA.program.undetectedType = typeCell{:};
    
    % change state
    changeState('draw',false,'mark_undetected','x')
end

% statebuttonEditMissedCalls_Callback -------------------------------------
function statebuttonEditMissedCalls_Callback(src,edata)

    % activate or deactivate edit mode
    isOn = edata.Value;
    if isOn
        changeState('edit',true,'MissedCalls');
    else
        saveSpreadsheet();
        changeState('normal',true);
    end
end

% buttonAddAnnotation_Callback --------------------------------------------
function buttonAddAnnotation_Callback(src,edata)

    % checkout global variables
    global DATA
    
    % activate markers if they're not active
    if ~DATA.markers.show
        showLoadScreen();
        DATA.markers.show = true;
        updatePlot();
    end

    % change state
    changeState('draw',false,'annotate','xy');
end

% statebuttonEditAnnotations_Callback -------------------------------------
function statebuttonEditAnnotations_Callback(src,edata)
    
    % activate or deactivate edit mode
    isOn = edata.Value;
    if isOn
        changeState('edit',true,'Annotations');
    else
        saveSpreadsheet();
        changeState('normal',true);
    end
end

% PlaybackStartFcn --------------------------------------------------------
function PlaybackStartFcn(src,edata)

    % checkout global variables
    global UI PLOT

    % create line
    ha = UI.axesMain;
    xLine = repelem(src.UserData.dtStart,1,2);
    yLine = ha.YLim;
    PLOT.linePlayback = line(UI.axesMain,xLine,yLine,...
        'LineWidth',2,...
        'LineStyle','-',...
        'Color',[0 0 0]);

    % change state
    changeState('playback',false)
end

% PlaybackTimerFcn --------------------------------------------------------
function PlaybackTimerFcn(src,edata)

    % checkout global variables
    global DATA PLOT
    
    % get current sample in datetime
    lag = -0.5*DATA.playback.speed;
    dtStart = src.UserData.dtStart;
    dtNow = dtStart + seconds(src.CurrentSample/src.UserData.FsActual + lag);
    dtNow = max([dtStart,dtNow]);

    % update line
    PLOT.linePlayback.XData = [dtNow,dtNow];
    drawnow
end

% PlaybackStopFcn ---------------------------------------------------------
function PlaybackStopFcn(src,edata)

    % checkout global variables
    global PLOT

    % change state
    changeState('normal',false)
    
    % remove cursor line and selection patch (in case it exists)
    delete(PLOT.linePlayback);
    delete(PLOT.patchSelection);
    
    % delete audioplayer
    delete(src)
end


%% UTILITY FUNCTIONS ======================================================
% getSignalRanges ---------------------------------------------------------
function [DTRange,FRange] = getSignalRanges(outTable)
% Returns the start and end datetimes and min/max frequencies of all 
% signals of a particular type based on its output spreadsheet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: Do not use the SigStartDateTime field directly, it results in loss
% of accuracy

    % checkout global variables
    global PARAMS DATA
    
    % get reference datetime
    dtRef = DATA.wavInfo.dtStartRef;
    
    % get datetime range of signals
    %DTRange = [outTable.SigStartDateTime, outTable.SigStartDateTime + seconds(outTable.SigEnd - outTable.SigStart)];
    DTRange = dtRef + seconds(outTable.FileStart + [outTable.SigStart, outTable.SigEnd]);
    
    % get frequency range
    try
        FRange = [outTable.SigMinFreq,outTable.SigMaxFreq];
    catch
        FRange = repmat(PARAMS.markers.StandardFRange,size(DTRange,1),1);
    end
end