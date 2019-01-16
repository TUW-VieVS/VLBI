% #########################################################################
% #     updateSeveralPopupmenus
% #########################################################################
%
% DESCRITPION
% This function updates some of the popupmenus where e.g. models can be
% chosen. This might be nice when you created a new file/folder and want to
% use it in VieVS
%
% AUTHOR
%   Matthias Madzak ???
%
% INPUT
%   handles      structure from the GUI (also containing data, e.g. residuals)
%   hObject
%
% OUTPUT
%
% COUPLING
% - getTLEFileNames
%
% CHANGES
% 2015-09-30 by A. Hellerschmied: Extracted from vievs2_3.m
% 2015-10-01 by A. Hellerschmied: Updated (now all popupmenues are updated!)
% 2015-10-02 by A. Hellerschmied: Added possibility to add filenames to the popupmenu string by defining them as input argument for the subroutine "update_popupmenu_files_in_dir". And Bug-fix!
% 2017-01-27 by A. Hellerschmied: bug fix with TLE directory
% 2017-02-19 by H. Krasna: ../DATA/GLOB/TRF/AO/*.txt, ../DATA/GLOB/TRF/STSEASON/*.txt, ../DATA/GLOB/TRF/APLRG/*.txt added
% 2017-02-28 by A. Hellerschmied: Fixed a problem with updating OUTLIER directories
% 2017-06-08 by A. Hellerschmied: Fixed a bug which caused problems, if no own TRF or won CRF was available in tghe TRF or CRF folder.
% 2017-08-31 by A. Hellerschmied: Fixed a problem when the dir /ION/FILES/ was empty.
% 2018-02-11 by D. Landskron: external troposphere files section removed
% 

function updateSeveralPopupmenus(hObject, handles)
% This function updates some of the popupmenus where e.g. models can be
% chosen. This might be nice when you created a new file/folder and want to
% use it in VieVS


% #########################################################################
% # File
% #########################################################################

% ### OPT-directory ###
path_dir               = '../../VLBI_OPT/';
popupmenu_tag          = 'popupmenu_setInput_optDir';
folder_description_str = 'OPT file directory';
update_popupmenu_folder_in_dir(path_dir, popupmenu_tag, folder_description_str, handles)

% ### Outlier-directory ###
curContent=get(handles.popupmenu_setInput_outDir, 'String');
if ~iscell(curContent)
   curSelected=curContent;
   if strcmp(curSelected, ' ')
       curSelected = '';
   end
else
    curSelected=curContent{get(handles.popupmenu_setInput_outDir, 'Value')};
end
dirsInOutlierFolder=dir('../DATA/OUTLIER/');
yrsToDelete=~cellfun(@isnan, cellfun(@str2double, {dirsInOutlierFolder.name}, 'UniformOutput', false));
dirsInOutlierFolder(strcmp({dirsInOutlierFolder.name}, '.')|strcmp({dirsInOutlierFolder.name}, '..')|~[dirsInOutlierFolder.isdir]|yrsToDelete)=[];
if isempty({dirsInOutlierFolder.name})
    set(handles.popupmenu_setInput_outDir, 'Value', 1)
    set(handles.popupmenu_setInput_outDir, 'string', ' ');
else
    set(handles.popupmenu_setInput_outDir, 'string', {'', dirsInOutlierFolder.name});
end
oldFoundInNew=~cellfun(@isempty, strfind({'', dirsInOutlierFolder.name}, curSelected));
if sum(oldFoundInNew)>0
    set(handles.popupmenu_setInput_outDir, 'Value', find(oldFoundInNew))
else
    % maybe it was '':
    if isempty(curSelected)
        set(handles.popupmenu_setInput_outDir, 'Value', 1)
    else
        msgbox('Previously selected Outlier directory was not found!', 'Warning', 'warn')
    end
end



% #########################################################################
% # Models
% #########################################################################

% ##### Reference Frames #####

% ### TRF file (ASCII) ###
path_dir             = '../TRF/*.txt';
popupmenu_tag        = 'popupmenu_parameters_refFrames_otherTRF';
file_description_str = 'TRF file';
additonal_filenames = {};
remove_filenames = {'SavedGuiData_superstations.txt'};
flag_no_files = update_popupmenu_files_in_dir(path_dir, popupmenu_tag, file_description_str, handles, additonal_filenames, remove_filenames);
if flag_no_files
    set(handles.popupmenu_parameters_refFrames_otherTRF, 'Enable', 'off')
    set(handles.radiobutton_parameters_refFrames_otherTRF, 'Enable', 'off')
else
    set(handles.popupmenu_parameters_refFrames_otherTRF, 'Enable', 'on')
    set(handles.radiobutton_parameters_refFrames_otherTRF, 'Enable', 'on')
end

% ### CRF file (ASCII) ###
path_dir             = '../CRF/*.txt';
popupmenu_tag        = 'popupmenu_parameters_refFrames_otherCRF';
file_description_str = 'CRF file';
flag_no_files = update_popupmenu_files_in_dir(path_dir, popupmenu_tag, file_description_str, handles);
if flag_no_files
    set(handles.popupmenu_parameters_refFrames_otherCRF, 'Enable', 'off')
    set(handles.radiobutton_parameters_refFrames_otherCRF, 'Enable', 'off')
else
    set(handles.popupmenu_parameters_refFrames_otherCRF, 'Enable', 'on')
    set(handles.radiobutton_parameters_refFrames_otherCRF, 'Enable', 'on')
end


% ##### Ionosphere #####

% ### folder of external iono files ###
path_dir = '../ION/';
popupmenu_tag          = 'popupmenu_parameters_iono_ext';
folder_description_str = 'Ionospheric delay folder';
update_popupmenu_folder_in_dir(path_dir, popupmenu_tag, folder_description_str, handles)


% ##### Station models #####

% ### non-tidal atmo loading ###
curContent=get(handles.popupmenu_parameters_statCorr_nonTidalAtmoOceanLoad, 'String');
curSelected=curContent{get(handles.popupmenu_parameters_statCorr_nonTidalAtmoOceanLoad, 'Value')};
dirsInAtmFolder=dir('../ATM/');
dirsInAtmFolder(strcmp({dirsInAtmFolder.name}, '.')|strcmp({dirsInAtmFolder.name}, '..')|strcmp({dirsInAtmFolder.name}, 'temp')|~[dirsInAtmFolder.isdir])=[]; % exclude the "temp" folder!
set(handles.popupmenu_parameters_statCorr_nonTidalAtmoOceanLoad, 'String', {dirsInAtmFolder.name})
oldFoundInNew=~cellfun(@isempty, strfind({dirsInAtmFolder.name}, curSelected));
if sum(oldFoundInNew)>0
    set(handles.popupmenu_parameters_statCorr_nonTidalAtmoOceanLoad, 'Value', find(oldFoundInNew))
else
    msgbox('Previously selected non-tidal atmo loading directory was not found!\n', 'Warning', 'warn')
end

% ### hydrology loading ###
path_dir               = '../HYDLO/';
popupmenu_tag          = 'popupmenu_parameters_statCorr_hydroLoading';
folder_description_str = 'Hydrology loading data';
update_popupmenu_folder_in_dir(path_dir, popupmenu_tag, folder_description_str, handles)


% ##### EOP #####

% ### EOP a priori time series loading ###
path_dir             = '../EOP/*.txt';
popupmenu_tag        = 'popupmenu_parameters_eop_aPriori_other';
file_description_str = 'EOP a priori time series loading';
update_popupmenu_files_in_dir(path_dir, popupmenu_tag, file_description_str, handles);


% ### High frequency ocean tides ###
% Here, the following strinsg have to be added to the popupmenu.
% - 'interpf (Conventions)'
% - 'Combi_IGG_Bonn'
% ...The according data for those corrections is not taken from dedicated files in the /EOP/eophf/ folder, but theyare hardcoded in the related functions.
path_dir                 = '../EOP/eophf/*.dat';
popupmenu_tag        = 'popupmenu_parameters_eop_oceanTideModel';
file_description_str = 'High frequency ocean tides file';
add_folder_names     = {'interpf (Conventions)', 'Combi_IGG_Bonn'};
update_popupmenu_files_in_dir(path_dir, popupmenu_tag, file_description_str, handles, add_folder_names);

%set(handles.popupmenu_parameters_eop_oceanTideModel, 'String', ['interpf (Conventions)', {dirsInEophfFolder.name}, 'Combi_IGG_Bonn'])


% ##### Source Structure #####

% ### source structure catalog ###
path_dir                 = '../CRF/SOURCE_STRUCTURE_CAT/*.cat';
popupmenu_tag        = 'popupmenu_parameters_ss_catalog';
file_description_str = 'Source structure catalog';
update_popupmenu_files_in_dir(path_dir, popupmenu_tag, file_description_str, handles);


% #########################################################################
% # Global Solution
% #########################################################################

% ##### CRF/TRF parameterization #####

% TRF:

% ### Station list for datum definition (NNR/NNT) ###
path_dir             = '../DATA/GLOB/TRF/DATUM/*.txt';
popupmenu_tag        = 'popupmenu_vie_glob_trfCrf_trf_stations4Nnt';
file_description_str = 'Station list for datum definition (NNR/NNT)';
update_popupmenu_files_in_dir(path_dir, popupmenu_tag, file_description_str, handles);

% ### List of stations which are reduced sessionwise ###
path_dir             = '../DATA/GLOB/TRF/REDUCE/*.txt';
popupmenu_tag        = 'popupmenu_vie_glob_trfCrf_trf_sessionWiseRedStat';
file_description_str = 'List of stations which are reduced sessionwise';
update_popupmenu_files_in_dir(path_dir, popupmenu_tag, file_description_str, handles);

% ### List of stations with breaks, but constant velocities ###
path_dir             = '../DATA/GLOB/TRF/VELOC/*.txt';
popupmenu_tag        = 'popupmenu_vie_glob_trfCrf_trf_keepConstVel';
file_description_str = 'List of stations with breaks, but constant velocities';
update_popupmenu_files_in_dir(path_dir, popupmenu_tag, file_description_str, handles);

% ### List of velocity ties between stations ###
path_dir             = '../DATA/GLOB/TRF/VELOC/TIES/*.txt';
popupmenu_tag        = 'popupmenu_vie_glob_trfCrf_trf_velTies';
file_description_str = 'List of velocity ties between stations';
update_popupmenu_files_in_dir(path_dir, popupmenu_tag, file_description_str, handles);

% ### List of discontinuities at stations (.mat) ###
path_dir             = '../DATA/GLOB/TRF/DISCONT/*.mat';
popupmenu_tag        = 'popupmenu_vie_glob_trfCrf_trf_positDisc';
file_description_str = 'List of discontinuities at stations';
update_popupmenu_files_in_dir(path_dir, popupmenu_tag, file_description_str, handles);

% CRF:

% ### List of sources for datum definition(NNR) ###
path_dir             = '../DATA/GLOB/CRF/DATUM/*.txt';
popupmenu_tag        = 'popupmenu_vie_glob_trfCrf_crf_sourcesForNnr';
file_description_str = 'List of sources for datum definition(NNR)';
update_popupmenu_files_in_dir(path_dir, popupmenu_tag, file_description_str, handles);

% ### List of sources which are fixed to a priori coordinates ###
path_dir             = '../DATA/GLOB/CRF/FIXED_SOURCES/*.txt';
popupmenu_tag        = 'popupmenu_vie_glob_trfCrf_crf_fixedSources';
file_description_str = 'List of sources which are fixed to a priori coordinates';
update_popupmenu_files_in_dir(path_dir, popupmenu_tag, file_description_str, handles);

% ### List of sources which are sessionwise reduced ###
path_dir             = '../DATA/GLOB/CRF/REDUCE/*.txt';
popupmenu_tag        = 'popupmenu_vie_glob_trfCrf_crf_sessionWiseRedSources';
file_description_str = 'List of sources which are sessionwise reduced';
update_popupmenu_files_in_dir(path_dir, popupmenu_tag, file_description_str, handles);

% ##### special parameters #####
% ### List of stations for antenna offset estimation###
path_dir             = '../DATA/GLOB/TRF/AO/*.txt';
popupmenu_tag        = 'listbox_vie_glob_axisOffsets';
file_description_str = 'List of stations for antenna offset estimation';
update_popupmenu_files_in_dir(path_dir, popupmenu_tag, file_description_str, handles);

% ### List of stations for seasonal harmonic signal estimation###
path_dir             = '../DATA/GLOB/TRF/STSEASON/*.txt';
popupmenu_tag        = 'listbox_vie_glob_stseaspos';
file_description_str = 'List of stations for seasonal harmonic signal estimation';
update_popupmenu_files_in_dir(path_dir, popupmenu_tag, file_description_str, handles);

% ### List of stations for APL regression coefficients estimation###
path_dir             = '../DATA/GLOB/TRF/APLRG/*.txt';
popupmenu_tag        = 'listbox_vie_glob_APLrg';
file_description_str = 'List of stations for APL regression coefficients estimation';
update_popupmenu_files_in_dir(path_dir, popupmenu_tag, file_description_str, handles);



% #########################################################################
% # Simulation
% #########################################################################

% ### turbulence parameters (listbox) ###
curContent=get(handles.listbox_vie_sim_paramFile, 'String');
curSelected=curContent{get(handles.listbox_vie_sim_paramFile, 'Value')};
dirsInTurbFolder=dir('../DATA/TURB/*.dat');
dirsInTurbFolder([dirsInTurbFolder.isdir])=[]; % remove all folders (i just want .dat files)
set(handles.listbox_vie_sim_paramFile, 'String', {dirsInTurbFolder.name})
oldFoundInNew=~cellfun(@isempty, strfind({dirsInTurbFolder.name}, curSelected));
if sum(oldFoundInNew)>0
    set(handles.listbox_vie_sim_paramFile, 'Value', find(oldFoundInNew))
else
    msgbox('Previously selected simulation turbulance file was not found!\n', 'Warning', 'warn')
end


% #########################################################################
% # Scheduling
% #########################################################################

% ### Local TLE files ###
path_tle_dir = '../ORBIT/TLE/';    % TLE data folder
[tle_str] = getTLEFileNames(path_tle_dir);
set(handles.popupmenu_vie_sched_select_tle,'String',tle_str);



% #########################################################################
% # Plotting
% #########################################################################

% ##### Residuals #####

% ### LEVEL3 subfolder ###
path_dir               = '../DATA/LEVEL3/';
popupmenu_tag          = 'popupmenu_plot_residuals_folder';
folder_description_str = 'LEVEL3 subfolder for plotting residuals';
update_popupmenu_folder_in_dir(path_dir, popupmenu_tag, folder_description_str, handles)


% ##### Parameters #####

% ### LEVEL3 subfolder 1 ###
path_dir               = '../DATA/LEVEL3/';
popupmenu_tag          = 'popupmenu_plot_folder1_subfolder';
folder_description_str = 'LEVEL3 subfolder 1 for plotting parameters';
update_popupmenu_folder_in_dir(path_dir, popupmenu_tag, folder_description_str, handles)

% ### LEVEL3 subfolder 2 ###
path_dir               = '../DATA/LEVEL3/';
popupmenu_tag          = 'popupmenu_plot_folder2_subfolder';
folder_description_str = 'LEVEL3 subfolder 2 for plotting parameters';
update_popupmenu_folder_in_dir(path_dir, popupmenu_tag, folder_description_str, handles)

% ### LEVEL3 subfolder 3 ###
path_dir               = '../DATA/LEVEL3/';
popupmenu_tag          = 'popupmenu_plot_folder3_subfolder';
folder_description_str = 'LEVEL3 subfolder 3 for plotting parameters';
update_popupmenu_folder_in_dir(path_dir, popupmenu_tag, folder_description_str, handles)


% ##### Session analysis #####

% ### LEVEL3 subfolder 1 (black) ###
path_dir               = '../DATA/LEVEL3/';
popupmenu_tag          = 'popupmenu_plot_sessionAnalysis_subfolder';
folder_description_str = 'LEVEL3 subfolder 1 for session analysis';
update_popupmenu_folder_in_dir(path_dir, popupmenu_tag, folder_description_str, handles)

% ### LEVEL3 subfolder 2 (red) ###
path_dir               = '../DATA/LEVEL3/';
popupmenu_tag          = 'popupmenu_plot_sessionAnalysis_subfolder2';
folder_description_str = 'LEVEL3 subfolder 2 for session analysis';
update_popupmenu_folder_in_dir(path_dir, popupmenu_tag, folder_description_str, handles)

% ### LEVEL3 subfolder 3 (blue) ###
path_dir               = '../DATA/LEVEL3/';
popupmenu_tag          = 'popupmenu_plot_sessionAnalysis_subfolder3';
folder_description_str = 'LEVEL3 subfolder 3 for session analysis';
update_popupmenu_folder_in_dir(path_dir, popupmenu_tag, folder_description_str, handles)

% ### LEVEL3 subfolder 4 (green) ###
path_dir               = '../DATA/LEVEL3/';
popupmenu_tag          = 'popupmenu_plot_sessionAnalysis_subfolder4';
folder_description_str = 'LEVEL3 subfolder 4 for session analysis';
update_popupmenu_folder_in_dir(path_dir, popupmenu_tag, folder_description_str, handles)




% ##### EOP/BAS output #####

% ### process list ###
path_dir             = '../WORK/PROCESSLIST/*.mat';
popupmenu_tag        = 'popupmenu_plot_eopOut_pl';
file_description_str = 'Process list for EOP/BAS output.';
update_popupmenu_files_in_dir(path_dir, popupmenu_tag, file_description_str, handles);

% ### LEVEL3 subfolder ###
path_dir               = '../DATA/LEVEL3/';
popupmenu_tag          = 'popupmenu_plot_eopOut_subfolder';
folder_description_str = 'LEVEL3 subfolder for EOP/BAS output';
update_popupmenu_folder_in_dir(path_dir, popupmenu_tag, folder_description_str, handles)



msgbox('Folders have been updated!', 'DONE!', 'help');

guidata(hObject, handles);

end

%% Subroutines

% ##### Function to update popupmenus containing names of files found in a defined deirectory #####
function flag_no_files = update_popupmenu_files_in_dir(path_dir, popupmenu_tag, file_description_str, handles, varargin)

flag_no_files = false;

switch nargin
    case 4 % Just add the names of all files found in the derfined directory to the popupmenu
        additonal_filenames = {};
        remove_filenames = {};
    case 5 % additional strings (which are provided as input argument) have to be added to the popupmenu
        additonal_filenames = varargin{1};
        remove_filenames = {};
    case 6
        additonal_filenames = varargin{1};
        remove_filenames = varargin{2};
    otherwise
        additonal_filenames = {}; % Error init.
end
curContent = eval(['get(handles.', popupmenu_tag, ', ''String'')']);
if ischar(curContent)
    curContent = {curContent};
end
curSelected = eval( ['curContent{get(handles.', popupmenu_tag, ', ''Value'')}'] );
filesInDir = dir(path_dir);
filesInDir(strcmp({filesInDir.name}, '.')|strcmp({filesInDir.name}, '..'))=[]; % Remove "." and ".."

if ~isempty(remove_filenames) % Remove filenames defined in input argument "remove_filenames"
    for i = 1 : length(remove_filenames)
        filesInDir(strcmp({filesInDir.name}, remove_filenames))=[]; % Remove "." and ".."
    end
end

filenames = {filesInDir.name};
filenames = [filenames, additonal_filenames]; % Add additional filesnames
if isempty(filenames)
    filenames = ' '; % It is not allowed to set an empty string!
    flag_no_files = true;
end
eval(['set(handles.', popupmenu_tag, ', ''String'', filenames)']);
oldFoundInNew = strcmp(filenames, curSelected); % ~cellfun(@isempty, strfind(filenames, curSelected));

if sum(oldFoundInNew) == 1
    eval(['set(handles.', popupmenu_tag, ', ''Value'', find(oldFoundInNew))']);
elseif sum(oldFoundInNew) == 0
    msgbox(['Previously selected file was not found: ', file_description_str, '(filename: ',curSelected ,')'], 'Warning', 'warn')
else
    error('ERROR in function: updateSeveralPopupmenus: The "Value" argument has to be a scalar!');
end

end % function


% ##### Function to update popupmenus containing names of folders found in a defined deirectory #####
function update_popupmenu_folder_in_dir(path_dir, popupmenu_tag, folder_description_str, handles)
curContent = eval(['get(handles.', popupmenu_tag, ', ''String'')']);
if ~iscell(curContent)
    curContent = {curContent};
end
% try
curSelected = eval( ['curContent{get(handles.', popupmenu_tag, ', ''Value'')}'] );
% catch
%     disp('iii');
% end
foldersInDir = dir(path_dir);
foldersInDir(strcmp({foldersInDir.name}, '.')|strcmp({foldersInDir.name}, '..')|~[foldersInDir.isdir])=[]; % Only get the names of the sub-directories!
if strcmp(curSelected, '/') % ...because "/" is set as default in the LEVEL3-subfolder popupmenus on program start in Vie_SETUP (plotting tools and EOP/BAS output)!
    foldersInDir(end + 1).name = '/';
end
eval(['set(handles.', popupmenu_tag, ', ''String'', {foldersInDir.name})']);
oldFoundInNew = strcmp({foldersInDir.name}, curSelected);% ~cellfun(@isempty, strfind({foldersInDir.name}, curSelected));
if sum(oldFoundInNew) == 1
    eval(['set(handles.', popupmenu_tag, ', ''Value'', find(oldFoundInNew))']);
elseif sum(oldFoundInNew) == 0
    msgbox(['Previously selected file was not found: ', folder_description_str, '(folder name: ',curSelected ,')'], 'Warning', 'warn')
else
    error('ERROR in function: updateSeveralPopupmenus: The "Value" argument has to be a scalar!');
end

end % function











