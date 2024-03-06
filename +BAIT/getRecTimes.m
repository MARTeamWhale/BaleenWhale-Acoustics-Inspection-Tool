function varargout = getRecTimes(recFilePaths)
%
% Read the start and optionally end times of WAV files. If end times are
% requested, each WAV file will be inspected for its duration, which can
% take a long time if there are many files.
% Includes a waitbar to estimate processing time remaining.
%
%   Written by Wilfried Beslin
%   Last updated 2024-03-06 using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    import MUCA.time.readDateTime
    
    nargoutchk(1,2)
    getStopTimes = nargout > 1;
    
    % get number of WAV files
    numRecs = numel(recFilePaths);
    
    % initialize output
    dtStart = NaT(numRecs,1);
    dtStop = dtStart;
    
    % initialize waitbar
    if getStopTimes
        waitmsg = 'Reading audio file start and stop times...';
    else
        waitmsg = 'Reading audio file start times...';
    end
    f = waitbar(0, waitmsg);
    tic
    
    % loop through each detection
    for ii = 1:numRecs
        % update waitbar
        t_elapsed = toc;
        t_rem = ceil(t_elapsed.*((numRecs-(ii-1))/(ii-1)));
        waitbar(ii/numRecs, f, sprintf('%s\nEstimated Time remaining: %s', waitmsg, duration(0,0,t_rem)))
        
        % isolate WAV file
        recPathii = recFilePaths{ii};
        [~,recNameii,~] = fileparts(recPathii);
        
        % get start time
        dtStart(ii) = readDateTime(recNameii);
        
        % get stop time if requested
        if getStopTimes
            try
                recInfoii = audioinfo(recPathii);
                recDurii = duration(0,0,recInfoii.Duration);
                dtStop(ii) = dtStart(ii) + recDurii;
            catch
                warning('Could not read stop time for file "%s"', recNameii)
            end
        end
    end
    
    % close waitbar
    close(f);
    
    % process output
    varargout{1} = dtStart;
    if getStopTimes
        varargout{2} = dtStop;
    end
end