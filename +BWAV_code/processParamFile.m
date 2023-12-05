function [paramFileText, usingDefault] = processParamFile(paramFileInputRaw, rootDir, scriptName)
% Parse an input parameter file name or path and try to locate the file. If
% the file is found, return the file contents.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % set path to PARAMS folder
    paramsDir = fullfile(rootDir, 'PARAMS', scriptName);

    % check if the input is empty, in which case default parameters should
    % be used
    usingDefault = isempty(paramFileInputRaw);
    if usingDefault
        % get default params
        defaultFileFullName = ['DefaultParams.',scriptName,'.txt'];
        paramFileInitPath = fullfile(paramsDir, defaultFileFullName);
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