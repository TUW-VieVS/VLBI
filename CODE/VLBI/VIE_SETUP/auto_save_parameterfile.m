% #########################################################################
% #     auto_save_parameterfile
% #########################################################################
%
% DESCRITPION
% This function saves a parameterfile (containing all GUI options) to a
% file (name eg. auto_save_fromGUI_yyyymmdd) which can be loaded by a menu
% entry 'load current'.
% This function is called whenever an option in the GUI is changed
%
% AUTHOR 
%   Matthias Madzak
%
% INPUT
%   handles      structure from the GUI (also containing data, eg residuals)
%   hObject      ...unused so far.
%
% OUTPUT
%
% CHANGES
%   2014-11-05, A. Hellerschmied: Exctacted code from "vievs2_3.m" file and
%     created this m-file. Header added.

function auto_save_parameterfile(hObject, handles)

% define base for filename
baseFileName='auto_save_fromGUI_';

% get current date and time
c=clock;
autoFilename=sprintf('%s%04.0f%02.0f%02.0f', baseFileName, c(1), c(2), c(3));
%[baseFileName, num2str(c(1)), num2str(c(2)), num2str(c(3)), '.mat'];

% delete all older entries ( if they have a different name, they would not
% be overwritten)
filesInParamFolder=dir(['../WORK/PARAMETERS/', baseFileName, '*.mat']);

% remove the filename which would be overwritten anyway
if ~isempty(filesInParamFolder)
    filesInParamFolder(strcmp({filesInParamFolder.name}, autoFilename))=[];
end

% delete entries
if ~isempty(filesInParamFolder)
    for k=1:length(filesInParamFolder)
        delete(['../WORK/PARAMETERS/', filesInParamFolder(k).name])
    end
end
    
% save new parameter file (and overwrite if already exists)
saveParamFile(hObject, handles, ['../WORK/PARAMETERS/', autoFilename]);