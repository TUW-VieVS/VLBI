% #########################################################################
% #     loadSupersourceFile
% #########################################################################
%
% DESCRITPION
% This function loads the supersource file (containing crf catalogs) and
% update the popupmenu to be able to select those.
%
% AUTHOR 
%   Matthias Madzak
%
% INPUT
%   handles      structure from the GUI (also containing data, eg residuals)
%   hObject      ...unused so far.
%   filename     filename of loaded Supersource file
%
% OUTPUT
%   handles      structure from the GUI (also containing data, eg residuals)
%
% CHANGES
%   2014-11-05, A. Hellerschmied: Exctacted code from "vievs2_3.m" file and
%     created this m-file
%   2016-10-28, H. Krasna: further catalogue fields added 


function handles=loadSupersourceFile(hObject, handles, filename)

load(filename);

supersourceFields=fieldnames(source);

% define all crf frames which should be written to the popupmenu
allCrfFrames={'icrf2','icrf2nonVCS', 'icrf3sx', 'icrf3k', 'VieCRF13','GSFC2015b', 'vievsCrf', 'userCrf'};

newPopupmenuEntries=supersourceFields(ismember(supersourceFields, allCrfFrames));
set(handles.popupmenu_parameters_refFrames_supersourceCRF, 'Value', 1)
set(handles.popupmenu_parameters_refFrames_supersourceCRF, 'String', newPopupmenuEntries)