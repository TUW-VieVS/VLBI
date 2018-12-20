% #########################################################################
% #     saveProcessList
% #########################################################################
%
% DESCRITPION
% This function saves a process list to the file-path/name defined in teh input variable "fullOutName".
%
% AUTHOR 
%   ???
%
% INPUT
% - hObject
% - handles     : GUI handles structure
% - fullOutName : Filepath and filename of the process list
%
% OUTPUT
%
% COUPLING
% - eop_out
% - saveProcessList
%
% CHANGES
% - 2015-08-20, A. Hellerschmied: Function extraced from vievs2_3.m

function saveProcessList(hObject, handles, fullOutName)

% get sessions
allSelectedSessions=get(handles.listbox_setInput_processList, 'String');

if ~isempty(allSelectedSessions)
    process_list=char(allSelectedSessions);

    % save the parameter file to file
    save(fullOutName, 'process_list');
end