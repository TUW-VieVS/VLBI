% #########################################################################
% #     getFreqBand
% #########################################################################
%
% DESCRIPTION
%   Takes a string and looks for frequency band and store it into str_inst
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
%       - id_fband (identifier of frequency band in string)
% OUTPUT:
%	- str_curr_fband (frequency band name)
%
% CHANGES:
%
function [ str_curr_fband ] = getFreqBand( str,id_fband )

i_str_fband = strfind(str,id_fband);

if isempty(i_str_fband)
    str_curr_fband = [];
else
    str_curr_fband = str(i_str_fband+length(id_fband):end);
    i_str_end = strfind(str_curr_fband,'_'); % find end indices in string
    if isempty(i_str_end) % if string doesn't consist of further identifier cut at "."
        i_str_end = strfind(str_curr_fband,'.');
    end
    if ~isempty(i_str_end) % if file doesn't consist of ending ".nc"
        str_curr_fband = str_curr_fband(1:i_str_end-1);
    end
end

end

