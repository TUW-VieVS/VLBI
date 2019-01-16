% #########################################################################
% #     handles2OptfilePath
% #########################################################################
%
% DESCRITPION
% This function takes the information in the handles structure and creates 
% filename and filepath of the OPT file
%
% AUTHOR 
%       A. Hellerschmied (2016-10-11)
%
% INPUT
%   handles     structure from the GUI (also containing data, eg residuals)
%
% OUTPUT
%   OPTfolder           Filepath of the OPt file (relatice to the VieVS WORK dir.)  - string
%   OPTFileName         Filename of the OPT file                                    - string
%   selectedOPTdir      Subdirectory of the OPT file (selected in the GUI)          - string
%
% CHANGES
%   

function [OPTfolder, OPTFileName, selectedOPTdir] = handles2OptfilePath(handles)

% get selected file
allSessionsInList   = get(handles.popupmenu_plot_residuals_session, 'String');
session             = allSessionsInList{get(handles.popupmenu_plot_residuals_session, 'Value')};

% get currently selected OPT directory
OPTdirs         = get(handles.popupmenu_setInput_optDir, 'String');
selectedOPTdir  = OPTdirs{get(handles.popupmenu_setInput_optDir, 'value')};


% ##### Check, if the input dataset file followed the standard naming convention #####
% => YYMMMDDcc_Nnnn (c...alphabetic character, n...number)
flag_std_naming_convention = true;

% Total length = 14 char
if length(session) ~= 14
    flag_std_naming_convention = false; 
end

if flag_std_naming_convention
    % "_" at car. 10
    if ~strcmp(session(10), '_')
        flag_std_naming_convention = false;
    end
    % first two characters are numbers:
    [num, status_1] = str2num(session(1:2));
    % char. 6+7 are numbers:
    [num, status_2] = str2num(session(6:7));
    % char. 12-14 are numbers:
    [num, status_3] = str2num(session(12:14));
    if ~(status_1 && status_2 && status_3)
        flag_std_naming_convention = false;
    end
end

if flag_std_naming_convention
    % ##### Standard naming convention is used for this session: #####
    
    % Get the year from the session name:
    if str2double(session(1:2)) > 75
        yearStr = ['19', session(1:2)];
    else
        yearStr = ['20', session(1:2)];
    end
    
    OPTFileName = [session(1:9), '.OPT'];
    
    
    
else
    % ##### Non-Standard naming convention is used for this session: #####
    % e.g. when using .vso input files
    % => Get yearStr from the opt_ file!
    
    % #### Load opt_ file from LEVEL 3 (sub-)directory: ####
    % Get sub-dir.:
    allSubfolders   = get(handles.popupmenu_plot_residuals_folder, 'string');
    curSubfolder    = allSubfolders{get(handles.popupmenu_plot_residuals_folder, 'Value')};
    
    % load opt_ file and get the year:
    load(['../DATA/LEVEL3/', curSubfolder, '/opt_', session, '.mat']);
    yearStr = opt_.year;
    
    OPTFileName = [session, '.OPT'];
end

OPTfolder = ['../../VLBI_OPT/', selectedOPTdir, '/', yearStr, '/'];

