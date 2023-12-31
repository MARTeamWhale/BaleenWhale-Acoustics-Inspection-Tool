function var = readParam(filestr, varname, validation_fcns)
% Read a parameter value from a parameter file.
% The validation_fcns argument is a cell array of function handles for
% ensuring the extracted variable is valid. It should be populated with
% functions that throw errors when conditions are not met, such as 'assert'
% and 'validateattributes'. The variable is considered valid if it passes
% at least one condition.
%
%   Written by Wilfried Beslin
%   Last updated 2023-12-06 using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define standard regex for reading parameters
    param_prefix_regex = 'VAR[ \t]*=[ \t]*';  % regular expression part to isolate parameter name and the equal sign that comes after it
    param_std_value_regex = '[\w.-]+';  % regular expression part to isolate standard (scalar) parameter values
    param_array_value_regex = '\[[\d.-, ]*\]';  % regular expression part to isolate multivalue parameters enclosed in braces
    param_regex = ['(?<=',param_prefix_regex,')((',param_std_value_regex,')|(',param_array_value_regex,'))'];

    % read variable string
    varstr_raw = char(regexp(filestr, strrep(param_regex,'VAR',varname), 'match'));
    varstr_num = regexpi(varstr_raw, '[\d.-]+|NaN', 'match');
    varstr_bool = regexpi(varstr_raw, 'true|false', 'match');

    % check what kind of variable it is
    if ~isempty(varstr_num)
        % numeric inputs
        var = str2double(varstr_num);
    elseif ~isempty(varstr_bool)
        % logical inputs
        var = strcmpi(varstr_bool, 'true');
    elseif ~isempty(varstr_raw) && ~strcmp(varstr_raw,'[]')
        % char input
        var = varstr_raw;
    else
        % empty
        var = [];
    end

    % validate input
    valid_input = true;
    for ii = 1:numel(validation_fcns)
        try
            feval(validation_fcns{ii}, var)
            valid_input = true;
            break
        catch
            valid_input = false;
        end
    end
    if ~valid_input
        error('The value of parameter "%s" is invalid or could not be read', varname)
    end
end