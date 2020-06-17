% #########################################################################
% #     loadSupersourceFile
% #########################################################################
%
% DESCRITPION
% This function loads the "superstation" file (contatining all station
% dependant information, such as loading, coordinates, velocities,...) and
% updates the popupmenu to be able to select e.g. different trfs.
%
% AUTHOR 
%   Matthias Madzak
%
% INPUT
%   handles      structure from the GUI (also containing data, eg residuals)
%   hObject      ...unused so far.
%   filename     filename of loaded Superstation file
%
% OUTPUT
%   handles      structure from the GUI (also containing data, eg residuals)
%
% CHANGES
%   2014-11-05, A. Hellerschmied: Extracted code from "vievs2_3.m" file and
%     created this m-file
%   2015-10-27, A. Hellerschmied: vrtf2014 and VieTRF13 added, vievsTrf is set as default whenever a superstation file is loaded.
%   2016-02-02, Y. Kwak: VieTRF10a is removed since it is not available in superstation file anymore
%                        userOwnTrf can be selected if there is a field for it in a superstation file   
%   2016-05-10, M. Madzak: ITRF2014, VTRF2014, ivsTrf2014b added
%   2017-02-19, H. Krasna: popup menu for APLrg and GIA
%   2017-10-20, H. Krasna: DTRF2014 added

function handles=loadSuperstationFile(hObject, handles, filename)


load(filename);

superstatFields=fieldnames(superstations);

% define all coordinate frames which exist in superstations file (those
% should be used for the popupmenu in vievs)
allCoordinateFrames={'vievsTrf', 'itrf2005', 'itrf2008', 'itrf2014', 'dtrf2014',...
    'vtrf2008', 'vtrf2014', 'ivsTrf2014b', 'VieTRF13', 'userOwnTrf'};

% update popupmenu
newPopupmenuEntries=superstatFields(ismember(superstatFields, allCoordinateFrames));
% set(handles.popupmenu_parameters_refFrames_superstationTRF, 'Value', 1);
set(handles.popupmenu_parameters_refFrames_superstationTRF, 'String', newPopupmenuEntries);

% make vievsTrf default if available
logTRFFound=~cellfun(@isempty, strfind(get(handles.popupmenu_parameters_refFrames_superstationTRF, 'String'), 'itrf2014'));
if sum(logTRFFound)>0
    set(handles.popupmenu_parameters_refFrames_superstationTRF, 'Value', find(logTRFFound));
else
    set(handles.popupmenu_parameters_refFrames_superstationTRF, 'Value', 1);
end

% update ocean/atmo tide popupmenu
ocLoadFields=fieldnames(superstations(1).ocean_loading);
atLoadFields=fieldnames(superstations(1).atmosphere_tidal_loading);
APLrgLoadFields=fieldnames(superstations(1).aplrg);
GIALoadFields=fieldnames(superstations(1).gia);
set(handles.popupmenu_parameters_statCorr_tidalOceanLoad, 'String',...
    ocLoadFields)
set(handles.popupmenu_parameters_statCorr_tidalAtmoOceanLoad, 'String',...
    atLoadFields)
set(handles.popupmenu_parameters_statCorr_APLrg, 'String',...
    APLrgLoadFields)
set(handles.popupmenu_parameters_statCorr_GIA, 'String',...
    GIALoadFields(2:end)) % 1st field is the reference epoch
% make gsfc, fes2004
otldefaultLog=~cellfun(@isempty, strfind(ocLoadFields, 'TPXO72'));
if sum(otldefaultLog)>0
    set(handles.popupmenu_parameters_statCorr_tidalOceanLoad, 'Value', ...
        find(otldefaultLog))
end
gsfcLog=~cellfun(@isempty, strfind(atLoadFields, 'GSFC'));
if sum(gsfcLog)>0
    set(handles.popupmenu_parameters_statCorr_tidalAtmoOceanLoad, 'Value',...
        find(gsfcLog))
end
