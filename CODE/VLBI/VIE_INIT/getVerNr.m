% #########################################################################
% #     getVerNr
% #########################################################################
%
% DESCRIPTION
%   Takes a string and looks for version number. Output of the function is the version number in double. If the a version nr does not exist in the string the output is -1
%   
%
% CREATED  
%   2017-11-14     Jakob Gruber
%
% REFERENCES
%
%
% COUPLING
%
%
% INPUT:
%	- str (string)
%       - id_ver (identifier of version in string)
% OUTPUT:
%	- num_curr_ver (version number)
%
% CHANGES:
%
function [ num_curr_ver ] = getVerNr( str,id_ver )

i_str_ver = strfind(str,id_ver);
if isempty(i_str_ver)
    num_curr_ver = -1;
else
    str_curr_ver = str(i_str_ver+length(id_ver):end);
    i_str_end = strfind(str_curr_ver,'_'); % find end indices in string
    if isempty(i_str_end) % if string doesn't consist of further identifier cut at "."
        i_str_end = strfind(str_curr_ver,'.');
    end
    if ~isempty(i_str_end) % if file doesn't consist of ending ".nc"
        str_curr_ver = str_curr_ver(1:i_str_end-1);
    else
        str_curr_ver = str_curr_ver;
    end

    % Correct version number with letter
    if ~isempty(find(double(str_curr_ver)>96))
        i_no_double_ver = find(double(str_curr_ver)>96);
        i_yes_double_ver = find(double(str_curr_ver)>47 & double(str_curr_ver)<58);
        num_curr_ver = str2double(str_curr_ver(i_yes_double_ver))+double(str_curr_ver(i_no_double_ver))*10^-3;
    else
        num_curr_ver = str2double(str_curr_ver);
    end
end

end

