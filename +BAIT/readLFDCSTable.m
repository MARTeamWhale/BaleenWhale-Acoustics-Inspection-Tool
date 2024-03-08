function [tableLFDCS, absDetTimes, usingRelativeTimeFormat, hasPrecisionLoss] = readLFDCSTable(csvFilePath, varargin)
%
% Read an LFDCS autodetections CSV file as a MATLAB table. May also isolate
% a subset of the file by applying filters.
%
%   INPUT ARGUMENTS
%   -----------------------------------------------------------------------
%   "csvFilePath" - char string specifying path to LFDCS autodetections CSV
%       file
%   .......................................................................
%   "filtParams" [OPTIONAL] - struct specifying filtering parameters to 
%       isolate a subset of the CSV file. Fields must include:
%           'ManualSpeciesCodes' -> list of manual species codes to include
%           'AutoCallTypes' -> list of auto call types to include
%           'StartDateTime' -> earliest occurence time a detection can have
%           'StopDateTime' -> latest occurrence time a detection can have
%       To ignore filters, set them to NaN (for species and call type 
%       codes) or datetime -InF/+Inf (for start and stop date times,
%       respectively).
%   -----------------------------------------------------------------------
%
%   OUTPUT ARGUMENTS
%   -----------------------------------------------------------------------
%   "tableLFDCS" - MATLAB table corresponding to the original data in the
%       LFDCS CSV file (does not include the header). Variable names have
%       been added to this table. The names and contents of columns 2 and 3
%       vary depending on the time format used in the file.
%   .......................................................................
%   "absDateTime" - N-by-2 matrix of detection datetimes, where column 1 
%       corresponds to start times and column 2 corresponds to stop times.
%   .......................................................................
%   "usingRelativeTimeFormat" - true/false value that indicates if the CSV
%       file describes detection times as seconds relative to 1970/01/01
%       (true), or as absolute date-times with fractional seconds (false).
%   .......................................................................
%   "hasPrecisionLoss" - true/false value indicating if detection times in
%       this CSV file were rounded by MS Excel or not.
%   -----------------------------------------------------------------------
%
%
%   Written by Wilfried Beslin
%   Last updated 2024-03-08 using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    import MUCA.io.importTextFile
    import MUCA.time.isoFormat
    
    narginchk(1,2)

    % determine number of header lines in LFDCS CSV file
    fileText = importTextFile(csvFilePath);
    numHeaderLines = find(cellfun('isempty',fileText), 1, 'last');
    if isempty(numHeaderLines)
        headerExpr = '^,+$'; % expression for getting a row of commas
        headerMatch = regexp(fileText, headerExpr);
        numHeaderLines = find(~cellfun('isempty',headerMatch), 1, 'last');
    end
    
    % read LFDCS autodetections csv as a table
    importOpts = detectImportOptions(csvFilePath, 'NumHeaderLines',numHeaderLines, 'DatetimeType','text');
    tableLFDCS = readtable(csvFilePath, importOpts);
    
    % Determine which time format is used in the CSV file.
    %%% LFDCS has two different ways of specifying detection times,
    %%% depending on which settings were used when running
    %%% "export_autodetections":
    %%% Option 1) specifies detection start times in column 2 and end times 
    %%%     in column 3 as the number of seconds relative to 
    %%%     1970/01/01 00:00:00.
    %%% Option 2) specifies absolute detection start date-times in column 2
    %%%     using the format MM/dd/yy HH:mm:ss, with added fractional
    %%%     seconds in column 3.
    %%% Both options include duration in column 4.
    iColStartTime = 2;
    usingRelativeTimeFormat = isnumeric(tableLFDCS{:,iColStartTime});
    
    % determine table variable names based on format
    varnames = {'CallType', '', '', 'Duration', 'MinFreq', 'MaxFreq', 'Bandwidth', 'Amplitude', 'MDist', 'ManualSpeciesCode', 'ManualCallTypeCode'};
    if usingRelativeTimeFormat
        varnames(2:3) = {'RelStartTime','RelStopTime'};
    else
        varnames(2:3) = {'AbsStartTime','FracSec'};
    end
    tableLFDCS.Properties.VariableNames = varnames;
    
    % get absolute detection times based on time type
    if usingRelativeTimeFormat
        % process for case where start-end times are provided relative
        % to 1970/01/01
        dtRef = datetime(1970, 1, 1, 0, 0, 0);
        absDetTimes = dtRef + seconds([tableLFDCS.RelStartTime, tableLFDCS.RelStopTime]);
    else
        % process for case where start times are provided using
        % absolute datetimes and fractional seconds
        absDetStartRounded = datetime(tableLFDCS.AbsStartTime, 'InputFormat','MM/dd/yy HH:mm:ss', 'PivotYear',year(datetime('now'))-99);
        absDetTimes = NaT(height(tableLFDCS),2);
        absDetTimes(:,1) = absDetStartRounded + seconds(tableLFDCS.FracSec);
        absDetTimes(:,2) = absDetTimes(:,1) + seconds(tableLFDCS.Duration);
    end
    absDetTimes.Format = isoFormat('long','simplified','milliseconds');
    
    % check detection times to see if Excel has dropped the milliseconds 
    hasPrecisionLoss = mean(second(absDetTimes(:)) - round(second(absDetTimes(:))) == 0) > 0.8;
    
    % identify rows to keep based on filters, if any
    if nargin > 1
        filtParams = varargin{1};
        rowsInclude = true(height(tableLFDCS),1);
        
        %%% manual species codes
        if ~isnan(filtParams.ManualSpeciesCodes)
            rowsInclude = rowsInclude & ismember(tableLFDCS.ManualSpeciesCode,filtParams.ManualSpeciesCodes);
        end
        
        %%% auto call type codes
        if ~isnan(filtParams.AutoCallTypes)
            rowsInclude = rowsInclude & ismember(tableLFDCS.CallType,filtParams.AutoCallTypes);
        end
        
        %%% date-time range
        rowsInclude = rowsInclude & absDetTimes(:,1) >= filtParams.StartDateTime & absDetTimes(:,1) <= filtParams.StopDateTime;
    
        % truncate table
        tableLFDCS = tableLFDCS(rowsInclude,:);
        absDetTimes = absDetTimes(rowsInclude,:);
    end
end