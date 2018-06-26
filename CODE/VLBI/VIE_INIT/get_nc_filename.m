% ************************************************************************
%   Description:
%   This function matches the string(s) inputted via the cell-array "field_name_pattern" with the entries (stings) in "wrapper_data_files".
%   It returns the matchig entry.
%
%   Input:	
%   - field_name_pattern (Cell-array):      Cell-array with 1 or more strings, which are matched with entries in "wrapper_data_files"
%   - wrapper_data_files (Cell-array):      Containing name(s) of .nc files in vgosDB (sub-)folder
%   - flag_generate_error_msg (optional):   If this flag is set to "1" (default = 0), an error will be generated if there is no match. Otherwise just a warnign appears.
%
%   Output:
%   - fieldname_str (string):               Name of entry in "wrapper_data_files" that matches the string(s) in "field_name_pattern"
% 
%   External calls: 	
%       
%   Coded for VieVS: 
%   2018-06-26, A. Hellerschmied
%
%   Revision: 
%   yyyy-mm-dd,FIRSTNAME SECONDNAME:
%

function fieldname_str = get_nc_filename(field_name_pattern, wrapper_data_files, varargin)

    % Init.:
    flag_generate_error_msg = 0;
    switch(nargin)
        case 3
            flag_generate_error_msg = varargin{1};
    end

    % Check input:
    if ~iscell(field_name_pattern)
        field_name_pattern = {field_name_pattern};
    end
    if ~iscell(wrapper_data_files)
        wrapper_data_files = {wrapper_data_files};
    end
    
    % Init.:
    flaglist_match_wrapper_data_files = true(length(wrapper_data_files), length(field_name_pattern));
    
    % Select file from list according to "field_name_default"  
    for i_1 = 1 : length(field_name_pattern)
        for i_2 = 1 : length(wrapper_data_files)
            flaglist_match_wrapper_data_files(i_2, i_1) = flaglist_match_wrapper_data_files(i_2, i_1) && ~isempty(strfind(wrapper_data_files{i_2},field_name_pattern{i_1}));
        end
    end
    
    % Check, wich file in "wrapper_data_files" matches all strings in "field_name_default":
    ind = true(length(wrapper_data_files), 1);
    for i_1 = 1 : length(field_name_pattern)
        ind = ind & flaglist_match_wrapper_data_files(:, i_1);
    end
    
    % Check, if there is a match:
    if sum(ind) == 1 % one match
        tmp_str = wrapper_data_files{ind};
        % Remove the ".nc" file ending:
        fieldname_str = tmp_str(1 : strfind(tmp_str, '.nc')-1);
    elseif sum(ind) >= 1 % more than 1 match
        error('More than 1 match was found!');
    else % sum == 0, no match
        fieldname_str = '';
        switch(flag_generate_error_msg)
            case 0
                warning('No match for .nc file (from wrapper) found in out_struct!')
            case 1
                error('No match for .nc file (from wrapper) found in out_struct!')
        end
    end
    
return