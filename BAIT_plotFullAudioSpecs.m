function BAIT_plotFullAudioSpecs(varargin)
%
% Create image files of spectrograms for entire WAV files within a folder.
%
% SYNTAX
% -------------------------------------------------------------------------
% BAIT_plotFullAudioSpecs()
% BAIT_plotFullAudioSpecs(Name,Value)
% -------------------------------------------------------------------------
%
%
%   INPUT ARGUMENTS (optional Name-Value pairs)
%   -----------------------------------------------------------------------
%   "params" - name or path of parameter file to use
%   .......................................................................
%   "input_dir" - string specifying path to the folder containing audio
%       files to be plotted. If not specified, user will be prompted to 
%       choose a folder.
%   .......................................................................
%   "output_dir" - string specifying path to the folder where spectrogram 
%       images should be saved. If not specified, user will be prompted to 
%       choose a folder.
%   .......................................................................
%   "overwrite" - True/false value specifying whether or not existing clips 
%       or spectrograms in the output folder should be overwritten. Default
%       is false.
%   -----------------------------------------------------------------------
%
%
%   DEPENDENCIES
%       BAIT.processParamFile
%       BAIT.readParam
%       BAIT.buildColormaps
%       MUCA.filepaths.listFiles
%       MUCA.audio.plotSpectrogram
%       MUCA.io.saveFig
%       
%
%   Written by Wilfried Beslin
%   Last updated 2025-12-02 using MATLAB R2024a
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    import MUCA.filepaths.listFiles
    import MUCA.audio.plotSpectrogram
    import MUCA.io.saveFig

    % 1) INPUT PARSING ....................................................

    p = inputParser;
    
    p.addParameter('params', '', @ischar)
    p.addParameter('input_dir', '', @ischar)
    p.addParameter('output_dir', '', @ischar)
    p.addParameter('overwrite', false)
    
    p.parse(varargin{:})
    paramsFileInput = p.Results.params;
    input_dir = p.Results.input_dir;
    output_dir = p.Results.output_dir;
    overwrite = p.Results.overwrite;
    
    % import parameters
    PARAMS = loadParams(paramsFileInput);

    
    % 2) PROCESSING .......................................................

    % get audio files from input folder
    [wav_filepaths, wav_filenames] = listFiles(input_dir, 'wav', 'Recursive',PARAMS.RecursiveSearch);
    num_files = numel(wav_filenames);

    % abort if no files are found
    if num_files == 0
        warning('No WAV files detected')
        return
    end

    % initialize spectrogram parameters
    stft_params = struct(...
        'winfcn', {'hann'},...
        'winsize', [],...
        'novl', [],...
        'nfft', []...
        );

    % initialize hidden figure
    fig = figure('Visible','off');
    ax = axes('Parent',fig);

    % loop through each file
    for ii = 1:num_files
        in_filename_ii = wav_filenames{ii};
        in_filepath_ii = wav_filepaths{ii};
        fprintf('Processing file %d of %d: ''%s''\n', ii, num_files, in_filename_ii);

        % check if output file already exists
        out_filepath_ii = replace(in_filepath_ii, input_dir, output_dir);
        out_filepath_ii = replace(out_filepath_ii,'.wav','.png');
        if ~overwrite && isfile(out_filepath_ii)
            disp('    Output file already exists. To overwrite this file, set ''overwrite'' = true')
            continue
        end

        try
            % get file sampling rate
            finfo_ii = audioinfo(in_filepath_ii);
            Fs_ii = finfo_ii.SampleRate;
    
            % determine spectrogram parameters for this file
            stft_params.winsize = round(Fs_ii*PARAMS.SpecFrameDur);
            stft_params.novl = stft_params.winsize - round(Fs_ii*PARAMS.SpecStepDur);
            stft_params.nfft = 2.^nextpow2(stft_params.winsize);
    
            % clear current figure
            cla(ax);
    
            % create spectrogram
            plotSpectrogram(ax, in_filepath_ii,...
                'Channel', PARAMS.Channel,...
                'FreqRange', [0,PARAMS.SpecMaxFreq],...
                'SpecParams', stft_params,...
                'Smooth', PARAMS.SmoothSpec,...
                'LogFreqs', PARAMS.LogFreqs,...
                'ColMap', PARAMS.SpecColorMap);

            % add title
            title(ax, in_filename_ii, 'Interpreter','none')

            % create output folder if it doesn't exist
            % (if the folder is within one or more folders that also don't
            % exist, all the required folders will be created at once)
            [out_dir_ii, ~, ~] = fileparts(out_filepath_ii);
            if ~isfolder(out_dir_ii)
                mkdir(out_dir_ii)
            end

            % save spectrogram as an image
            if isfile(out_filepath_ii) && overwrite
                fprintf('    Overwriting image ''%s''\n', out_filepath_ii)
                delete(out_filepath_ii)
            else
                fprintf('    Saving spectrogram image ''%s''\n', out_filepath_ii)
            end
            saveFig(fig, out_filepath_ii, PARAMS.SpecFigSize, 'pixels');

        catch ME
            warnmsg = sprintf('Failed to process file ''%s'':\n%s', in_filepath_ii, ME.message);
            warning(warnmsg)
        end
    end

    % wrap up
    close(fig);
    disp('Done')
end


% loadParams --------------------------------------------------------------
function PARAMS = loadParams(paramFileInput)
% Reads in program parameters from file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    import BAIT.processParamFile
    import BAIT.readParam
    import BAIT.buildColormaps

    % read parameter file as a block of text
    scriptPath = mfilename('fullpath');
    [rootDir,scriptName,~] = fileparts(scriptPath);
    paramsText = processParamFile(paramFileInput, rootDir, scriptName);

    % initialize output
    PARAMS = struct;
    
    % set parameters
    PARAMS.Channel = readParam(paramsText, 'Channel', {@(var)validateattributes(var,{'numeric'},{'scalar','positive','integer'})});
    PARAMS.RecursiveSearch = readParam(paramsText, 'RecursiveSearch', {@(var)validateattributes(var,{'logical'},{'scalar'})});
    PARAMS.SpecFrameDur = readParam(paramsText, 'SpecFrameDur', {@(var)validateattributes(var,{'numeric'},{'scalar','positive'})});
    PARAMS.SpecStepDur = readParam(paramsText, 'SpecStepDur', {@(var)validateattributes(var,{'numeric'},{'scalar','positive','<=',PARAMS.SpecFrameDur})});
    PARAMS.SmoothSpec = readParam(paramsText, 'SmoothSpec', {@(var)validateattributes(var,{'logical'},{'scalar'})});
    PARAMS.LogFreqs = readParam(paramsText, 'LogFreqs', {@(var)validateattributes(var,{'logical'},{'scalar'})});
    PARAMS.SpecMaxFreq = readParam(paramsText, 'SpecMaxFreq', {@(var)validateattributes(var,{'numeric'},{'scalar','positive'})});
    PARAMS.SpecColorMap = readParam(paramsText, 'SpecColorMap', {@(var)validateattributes(var,{'char'},{'row'})});
    PARAMS.SpecFigSize = readParam(paramsText, 'SpecFigSize', {@(var)validateattributes(var,{'numeric'},{'numel',2,'integer','positive'})});
    
    % assign colormap matrix
    cmaps = buildColormaps();
    try
        PARAMS.SpecColorMap = feval(cmaps.(PARAMS.SpecColorMap));
    catch
        error('Invalid colormap "%s"', PARAMS.SpecColorMap)
    end
end