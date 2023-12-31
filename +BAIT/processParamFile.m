function [paramFileText, usingDefault] = processParamFile(paramFileInputRaw, rootDir, scriptNameFull)
% Parse an input parameter file name or path and try to locate the file. If
% the file is found, return the file contents.
%
%   Written by Wilfried Beslin
%   Last updated 2023-12-05 using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % remove package name from script name
    scriptNameSlim = erase(scriptNameFull, 'BAIT_');

    % set path to PARAMS folder
    paramsDir = fullfile(rootDir, 'PARAMS', scriptNameSlim);

    % check if the input is empty, in which case default parameters should
    % be used
    usingDefault = isempty(paramFileInputRaw);
    if usingDefault
        % get default params
        defaultFileFullName = ['DefaultParams.',scriptNameSlim,'.txt'];
        paramFileInitPath = fullfile(paramsDir, defaultFileFullName);
        assert(isfile(paramFileInitPath), sprintf('Could not find default parameter file "%s"',defaultFileFullName))
    else
        % parse user input

        % break up input string
        [paramFileInitDir, paramFileName, paramFileExt] = fileparts(paramFileInputRaw);

        % add TXT extension if it doesn't exist
        if isempty(paramFileExt)
            paramFileExt = '.txt';
            paramFileInput = [paramFileInputRaw, paramFileExt];
        else
            paramFileInput = paramFileInputRaw;
        end

        % set initial file path
        paramFileInitPath = fullfile(paramFileInitDir, [paramFileName,paramFileExt]);

        % check if file exists - if it doesn't, begin loop to look for it
        while ~isfile(paramFileInitPath)
            % if only a name was specified, look for the file in the PARAMS
            % folder
            if isempty(paramFileInitDir)
                paramFileInitDir = paramsDir;
                paramFileInitPath = fullfile(paramFileInitDir, [paramFileName,paramFileExt]);
            else
                % if file cannot be found, then throw an error
                error('Could not find parameter file "%s"', paramFileInput)
            end
        end
    end
    
    
    % get the full path of the paramfile
    paramFileInfo = dir(paramFileInitPath);
    paramFileDir = paramFileInfo.folder;
    paramFileFullName = paramFileInfo.name;
    paramFilePath = fullfile(paramFileDir, paramFileFullName);

    % load the parameters
    if usingDefault
        disp('Loading default parameters');
    else
        fprintf('Loading parameter file "%s"\n', paramFilePath);
    end
    paramFileText = fileread(paramFilePath);
end