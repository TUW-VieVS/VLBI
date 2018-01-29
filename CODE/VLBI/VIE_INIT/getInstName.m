% #########################################################################
% #     getinstName
% #########################################################################
%
% DESCRIPTION
%   Takes a string and looks for institute name and store it into str_inst
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
%       - id_inst (identifier of institute in string)
% OUTPUT:
%	- str_inst (institute name)
%
% CHANGES:
%
function [ str_inst ] = getInstName( str, id_inst )
 
i_str_inst = strfind(str,id_inst); % find start indices in string

if isempty(i_str_inst)
    str_inst = [];
else
    str_inst = str(i_str_inst+length(id_inst):end);
    i_str_end = strfind(str_inst,'_'); % find end indices in string
    if isempty(i_str_end) % if string doesn't consist of further identifier cut at "."
        i_str_end = strfind(str_inst,'.');
    end
    if ~isempty(i_str_end) % if file doesn't consist of ending ".nc"
        str_inst = str_inst(1:i_str_end-1);
    else
        str_inst = str_inst;
    end
end

end

