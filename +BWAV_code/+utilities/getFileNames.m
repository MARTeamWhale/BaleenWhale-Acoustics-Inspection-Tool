% Utility function for obtaining filenames of a particular type of file
% within a directory.
% Written by Wilfried Beslin
% Last updated Nov 2020
%
% Input: (dirpath,ext)
%   dirpath [required] = string specifying full or relative path of 
%       directory of interest
%   ext [required] = string specifying the file extension of interest (don't include
%       the period)
%   recursive [optional] = logical specifying if subdirectories should be
%       included or not. Default is false.
%   MustContain [parameter] = regular expression specifying elements that a
%       path must have for it to be included
%   MustNotContain [parameter] = regular expression specifying elements
%       that a path should not have for it to be included
%
% Output: (filenames)
%   filenames = cell array of strings containg the names of each file
%       containing the extension. The extension itself is included too.
%   filepaths = cell array of strings containg the path of each file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [filenames,filepaths] = getFileNames(dirpath,ext,varargin)

    % input validation
    p = inputParser;
    
    % dirpath
    addRequired(p,'dirpath',@isdir)
    % ext
    validExt = @(arg) validateattributes(arg,{'char'},{'row'});
    addRequired(p,'ext',validExt)
    % recursive [OPTIONAL]
    defaultRecursive = false;
    validRecursive = @(arg) validateattributes(arg,{'logical'},{'scalar'});
    addOptional(p,'recursive',defaultRecursive,validRecursive)
    % MustContain [PARAMETER]
    defaultLimitExpr = '.*';
    validLimitExpr = @(arg) validateattributes(arg,{'char'},{'row'});
    addParameter(p,'MustContain',defaultLimitExpr,validLimitExpr);
    % MustNotContain [Parameter]
    defaultExclusionExpr = '';
    validExclusionExpr = @(arg) validateattributes(arg,{'char'},{'row'});
    addParameter(p,'MustNotContain',defaultExclusionExpr,validExclusionExpr);
    
    parse(p,dirpath,ext,varargin{:})
    recursive = p.Results.recursive;
    limitExpr = p.Results.MustContain;
    exclusionExpr = p.Results.MustNotContain;
    % end input parsing
    
    % initialize variables
    filenames = {};
    filepaths = {};
    
    % get file names and paths
    [filenames,filepaths] = getFilesFromDir(dirpath,ext,recursive,limitExpr,exclusionExpr,filenames,filepaths);
end

%%-------------------------------------------------------------------------
function [allFilenames,allFilepaths] = getFilesFromDir(dirPath,ext,recursive,limitExpr,exclusionExpr,allFilenames,allFilepaths)

    % get info on all files from specified directory
    dirFileInfo = dir([dirPath,filesep,'*.',ext]);
    wantedFilenames = {dirFileInfo.name}';
    wantedFilepaths = fullfile(dirPath,wantedFilenames);
    
    % only keep files that meet limit expression and don't meet exclusion expression
    keep = meetsRegEx(wantedFilepaths,limitExpr) & ~meetsRegEx(wantedFilepaths,exclusionExpr);
    wantedFilenames = wantedFilenames(keep);
    wantedFilepaths = wantedFilepaths(keep);
    
    % concatenate results with existing file names/paths
    allFilenames = [allFilenames;wantedFilenames];
    allFilepaths = [allFilepaths;wantedFilepaths];
    
    % recursion, if specified
    if recursive
        % get subdirectories
        dirInfoFull = dir(dirPath);
        subdirs = {dirInfoFull.name}';
        subdirs = subdirs([dirInfoFull.isdir]);
        dotExpr = '^(\.+)$'; % regular expression to find '.' or '..'
        keepDir = ~meetsRegEx(subdirs,dotExpr) & ~meetsRegEx(subdirs,exclusionExpr); % check for exclusions ONLY (because we're dealing with partial paths)
        subdirs = subdirs(keepDir);
        subdirPaths = fullfile(dirPath,subdirs);
        
        % get files in subdirectories
        for ii = 1:numel(subdirPaths)
            [allFilenames,allFilepaths] = getFilesFromDir(subdirPaths{ii},ext,recursive,limitExpr,exclusionExpr,allFilenames,allFilepaths);
        end 
    end
end

%% ------------------------------------------------------------------------
function m = meetsRegEx(str,expr)
    exprOut = regexp(str,expr);
    m = ~cellfun('isempty',exprOut);
end