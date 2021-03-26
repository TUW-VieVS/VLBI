function varargout = vie_setup(varargin)
% VIE_SETUP MATLAB code for vie_setup.fig
%      VIE_SETUP, by itself, creates a new VIE_SETUP or raises the existing
%      singleton*.
%
%      H = VIE_SETUP returns the handle to a new VIE_SETUP or the handle to
%      the existing singleton*.
%
%      VIE_SETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIE_SETUP.M with the given input arguments.
%
%      VIE_SETUP('Property','Value',...) creates a new VIE_SETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vie_setup_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vie_setup_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vie_setup

% Last Modified by GUIDE v2.5 30-Jun-2017 14:44:55


% 07 Jan 2014 by Matthias Madzak: LEVEL2 bug corrected
% 27 Jan 2014 by Hana Krasna: icrf2nonVCS set to default
% 31 Jan 2014 by Lucia Plank: time added for excluded station in OPT file
% 12 Feb 2014 by Lucia Plank: source structure simulation added
% 25 Mar 2014 by Lucia Plank: typo & bug corrected (load parameter file)
% 06 Jun 2014 by A. Hellerschmied: Time scale for residual plot
% 18 Jun 2014 by Hana Krasna: sources can be estimated in vie_lsm with
%   NNR condition
% 29 Jul 2014 by Hana Krasna: button added for baseline dependent weighting
% (function calc_bas_weights.m in VIE_LSM from Minttu Uunila)
% 18 Aug 2014 by A. Hellerschmied: option added to plot source coordinates with plotting tool
% 28 Aug 2014 by Hana Krasna: Version 2.3
% 09 Sep 2014 by David Mayer: scheduling update
% 26 Sep 2014 by Daniel Landskron: TRP update
% 05 Nov 2014 by A. Hellerschmied: default parameter file is deleted and
%   re-created during VieVS startup
% 02 Dez 2014 by M. Madzak: Show Warning, if subdir strings are not equal and vie_lsm is ticked
% 02 Dez 2014 by M. Madzak: Bugfix - fields of handles.data.plot.res structure.
% 08 Jan 2015 by M. Madzak & Anastasiia: Minor bug-fix
% 13 Jan 2015 by M. Madzak: Contact infos for sinex now added
% 20 Jan 2015 by A. Hellerschmied: Move function saveSchedParameters to a
%   separate file.
% 09 Feb 2015 by A. Hellerschmied: General GUI update for VIE_SCHED and the Run-Options menu;
%    vievsTrf and vievsCrf set as default for Vie_MOD
% 06 Mar 2015 by A. Hellerschmied: Minor bug fix
% 23 Apr 2015 by A. Hellerschmied: Error handling routine for process lists added 
%   (button in Run/Run Options added to the GUI)
% 19 Aug 2015 by A. Girdiuk: New Output options were added for writing EOP files (eop_out.m).
% 30 Sep 2015 by A. Hellerschmied: function "updateSeveralPopupmenus" extracted to a separate file and updated.
% 02 Oct 2015 by A. Hellerschmied: A sub-folder named "DEFAULT" is created automatically, if the "../../../VLBI_OPT/" directory is empty!
% 27 Oct 2015 by A. Hellerschmied: Minor change: default TRF is now defined/set in function "loadSuperstationFile".
% 06 Nov 2015 by S. Boehm: option to estimate tidal ERP variations added
% 16 Nov 2015 by A. Hellerschmied: Bug-fix: Problem with finding "superstation.mat" in the TRF folder while starting vievs fixed.
% 03 Dez 2015 by A. Hellerschmied: Web address of the VieVS webpage changed (to 'http://vievs.geo.tuwien.ac.at/'); Help/Documentation menu removed in the GUI
% 15 Dez 2015 by A. Hellerschmied:  - OPT and OUTLIER files will now be opened with the default editor, defined for these file types in Windows 
%                                   - Function OLDplotResidualsToAxes removed
%                                   - Bug-fix: Path for opening the outlier file via the GUI corrected
%                                   - Moved function plotSessionAnalysisToAxes to a separate file and corrected a minor bug in there
% 18 Dez 2015 by A. Hellerschmied:  - Adapted for VieVS 3.0
%                                   - Renamed from vievs2_3.m to vie_setup.m
% 07 Jan 2016 by M. Madzak: net CDF support and other minor changes
% 21 Jun 2016 by A. Hellerschmied: Changes to start netCDF analysis tool (analyseNetcdf.m) from VieVS GUI
% 07 Jul 2016 by A. Hellerschmied: Added checkbox "checkbox_estimation_leastSquares_sources_ICRF2_def" 
% 11 Jul 2016 by D. Mayer: Added edit field "edit_estimation_leastSquares_sources_abs_constr"
% 14 Jul 2016 by A. Hellerschmied: - Added possibility to select VSO files as input data (changes in: updateInputFilesBox, pushbutton_setInput_browseForSessions_Callback, openOPTfile) 
%                                  - Open/Create OPT and Outlier files also for vgosdb and vso files
% 25 Jul 2016 by M. Schartner: Sched analyser tool added to GUI menu
% 02 Aug 2016 by A. Hellerschmied: - Menu "Models/Space Crafts" added.
%                                  - Possibility to select SP3 orbit file added
% 08 Aug 2016 by A. Hellerschmied: - DATA/SIM/.. is now a valid directory for NGS files.
% 30 Aug 2016 by A. Hellerschmied: - DATA/SCHED/.. is now a valid directory for NGS files.
% 31 Aug 2016 by A. Girdiuk: - several function are moved in separate files and bug-fixed 
% 06 Sep 2016 by M. Schartner: changes for sched_manually and sched_analyser
% 10 Sep 2016 by A. Girdiuk: residuals plot: new panel is added: select data: different units can be shown 
%                            also: selected data can be written to OPT-file for stations and sources only
%                            WriteSelectDatatobeExcluded to write in OPT-file is added
% 19 Sep 2016 by A. Girdiuk: eop/bas out panel: select sessions in the subfolder intensive sessions output review
% 06 Sep 2016 by M. Schartner: changes for sched_manually and sched_analyser
% 05 Oct 2016 by M. Schartner: changes to select Bands and to select CATALOGS folder
% 17 Oct 2016 by A. Hellerschmied: Correct folder for ext. ion. files selected (/ION/FILES/)% 06 Sep 2016 by M. Schartner: changes for sched_manually and sched_analyser
% 19 Oct 2016 by A. Hellerschmied: - Moved function saveSimParameters to separte m-file
%                                  - Added additional items to the VIE_SIM GUI (rng setting, wn for satellite obs.)
% 05 Nov 2016 by M. Schartner: - no more warning if process list is empty and vie_sched is selected
%							   - different warnings if you select process list and run vie_sched
% 15 Dec 2016 by H. Krasna: field for elevation dependent noise added in GUI
% 22 Dec 2016 by A. Hellerschmied: - Edited some text fields in VIE_SCHED GUI
%                                  - Changed TLE directory from "/CATALOGS/TLE/" to "/ORBIT/TLE/" in the VieVS root dir.
%                                  - Changed SP3 directory from "/CATALOGS/SP3/" to "/ORBIT/SP3/" in the VieVS root dir.
% 17 Jan 2017 by H. Krasna: station seasonal harmonic variation and pole
%                           tide Love and Shida numbers added to GUI (as extra parameters in the N-global matrix and in the Global solution GUI )
% 19 Jan 2017 by D. Landskron: choice of temperature source transferred from "Models - Troposphere" to "Models - Station models", and bug fixed
% 23 Jan 2017 by D. Landskron: general shape of Troposphere changed
% 02 Feb 2017 by A. Hellerschmied: Minor change in VIE_SCHED GUI
% 19 Feb 2017 by H. Krasna: APL regression coefficient added (modelling in vie_mod and
%                           estimation in vie_glob); GIA uplift rates for
%                           modelling added; vie_glob GUI split to main and
%                           special parameter window
% 22 Feb 2017 by A. Hellerschmied: Error handling in case no input file was selected
% 23 Feb 2017 by A. Hellerschmied: - "edit_run_add_noise" added
%                                  - The filenames and -paths of the superstation- and supersource files are now displayed in the GUI
% 27.Feb 2017 by M. Schartner: changes for vie sched conditions
% 01 Mar 2017 by A. Hellerschmied: "checkbox_run_runOptions_manuallyFindBreaks" removed
% 10 Mar 2017 by D. Landskron: hard-coded part added for troposphere settings
% 04 Apr 2017 by A. Girdiuk: Plotting tool was supplied with warning message that you need to change default 'hh' before push the button 'Display'
% 08 May 2017 by A. Hellerschmied: Added unit for satellite pos. estimates in parameter plotting tool
% 29 May 2017 by A. Hellerschmied: Revised code to open Outlier und OPT files in text editor
% 31 May 2017 by A. Hellerschmied: VieVS logo removed from GUI title
% 08 Jun 2017 by A. Hellerschmied: Fixed a bug which caused problems, if no own TRF or won CRF was available in tghe TRF or CRF folder.
% 14 Jun 2017 by A. Hellerschmied: wrong vie_batch3_1.m is now called.
% 30 Jun 2017 by A. Hellerschmied: Added option to show errorbars when plotting parameters (checkbox_plot_show_errorbars added)
% 24 Aug 2017 by A. Girdiuk: output in the text file modified to supply the case when the data were removed from the plot window by clicking Reset
% 31 Aug 2017 by A. Hellerschmied: Problem with empty /ION/Files/ dir solved (popupmenu in GUI crashed)
% 11 Jan 2018 by D. Landskron: external troposphere modeling removed
% 18 Jan 2018 by A. Hellerschmied: Changes for transition to GIT (call of vie_batch.m)
% 06 Jul 2018 by D. Landskron: VMF3 added to the troposphere models
% 25 Sep 2018 by D. Landskron: specific warning message suppressed
% 18 Dec 2018 by D. Landskron: VieVS now starts with File - Set input files
% 06 Mar 2019 by D. Landskron: suffix checkbox added to the sinex files
% 15 Jan 2020 by M. Mikschi: added gravitational deformation checkbox
%
%*************************************************************************



% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vie_setup_OpeningFcn, ...
                   'gui_OutputFcn',  @vie_setup_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
   gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



% --- Executes just before vie_setup is made visible.
function vie_setup_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vie_setup (see VARARGIN)

%profile on
%set(handles.figure_vievs2, 'WindowButtonMotionFcn', {@mouseMove,handles});
% change logo icon
% if exist('vievs_logo.gif', 'file')
%     warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
%     jframe=get(handles.figure_vievs2,'javaframe');
%     jIcon=javax.swing.ImageIcon('vievs_logo.gif');
%     jframe.setFigureIcon(jIcon);
% end

% suppress the warning 'Setting the "WindowButtonUpFcn" property is not permitted while this mode is active.' which appears in Plotting - Residuals when scrolling through stations while "Zoom in" is active
warning('off','MATLAB:modes:mode:InvalidPropertySet')


allMainUiPanels=[handles.uipanel_file_setInputFiles, ...
    handles.uipanel_parameters_referenceFrames, ...
    handles.uipanel_parameters_ephemerides, ...
    handles.uipanel_parameters_troposphere, ...
    handles.uipanel_parameters_ionosphere, ...
    handles.uipanel_parameters_stationCorrections, ...
    handles.uipanel_parameters_eop, ...
    handles.uipanel_parameters_observationRestrictions,...
    handles.uipanel_parameters_sourcestructure,...
    handles.uipanel_estimation_leastSquares_tropo,...
    handles.uipanel_estimation_leastSquares_clock,...
    handles.uipanel_estimation_leastSquares_eop,...
    handles.uipanel_estimation_leastSquares_stationCoordinates,...
    handles.uipanel_estimation_leastSquares,...
    handles.uipanel_run_sinexOutput,...
    handles.uipanel_run_runOptions,...
    handles.uipanel_run_vievsProcSettings,...
    handles.uipanel_vie_glob,...
    handles.uipanel_vie_glob_TRFCRFparam,...
    handles.uipanel_vie_glob_special_param,...
    handles.uipanel_plot_parameters,...
    handles.uipanel_plot_residuals,...
    handles.uipanel_vie_sched,...
	handles.uipanel_vie_sched_minor_parameters,...
    handles.uipanel_vie_sim,...
    handles.uipanel_estimation_globalSolution,...
    handles.uipanel_plot_sessionAnalysis,...
    handles.uipanel_plot_scheduling,...
	handles.uipanel_plot_eopOut,...
    handles.uipanel_models_space_crafts,...
    handles.uipanel_welcome];


% save all ui panels to handles struct
handles.allMainUiPanels=allMainUiPanels;

% get directories' content
dirsInOptFolder=dir('../../VLBI_OPT/');
dirsInOutlierFolder=dir('../DATA/OUTLIER/');
dirsInTrpFolder=dir('../TRP/OUTPUT_DATA/');
dirsInIonFolder=dir('../ION/FILES/');
dirsInAtmFolder=dir('../ATM/');
dirsInWorkPlFolder=dir('../WORK/PROCESSLIST/*.mat');
% dirsInAtideFolder=dir('../ATIDE/*.mat');
% dirsInOtideFolder=dir('../OTIDE/*.mat');
dirsInHydloFolder=dir('../HYDLO/');
dirsInTrfFolder=dir('../TRF/*');
dirsInCrfFolder=dir('../CRF/*.txt');
dirsInCrfFolderMat=dir('../CRF/*.mat');
dirsInEopFolder=dir('../EOP/*.txt');
dirsInEophfFolder=dir('../EOP/eophf/*.dat');
dirsInDataFolder=dir('../DATA/LEVEL3/');
dirsInTurbFolder=dir('../DATA/TURB/*.dat');
dirsInGlobTrfDatum=dir('../DATA/GLOB/TRF/DATUM/*.txt');
dirsInGlobTrfReduce=dir('../DATA/GLOB/TRF/REDUCE/*.txt');
dirsInGlobTrfDiscont=dir('../DATA/GLOB/TRF/DISCONT/*.mat');
dirsInGlobTrfAo=dir('../DATA/GLOB/TRF/AO/*.txt');
dirsInGlobTrfSeason=dir('../DATA/GLOB/TRF/STSEASON/*.txt');
dirsInGlobTrfAPLrg=dir('../DATA/GLOB/TRF/APLRG/*.txt');
dirsInGlobTrfVel=dir('../DATA/GLOB/TRF/VELOC/*.txt');
dirsInGlobTrfVelTies=dir('../DATA/GLOB/TRF/VELOC/TIES/*.txt');
dirsInGlobCrfDatum=dir('../DATA/GLOB/CRF/DATUM/*.txt');
dirsInGlobCrfFixed=dir('../DATA/GLOB/CRF/FIXED_SOURCES/*.txt');
dirsInGlobCrfReduce=dir('../DATA/GLOB/CRF/REDUCE/*.txt');
dirsInSoucatFolder=dir('../CRF/SOURCE_STRUCTURE_CAT/*.cat');
statlistFolder='../WORK/STATIONLIST/';
if ~exist(statlistFolder, 'dir')
    mkdir(statlistFolder);
end
dirsInStatlistFolder=dir([statlistFolder, '*.mat']);

% try to find superstations file in TRF folder
superstatfile=find(strcmp({dirsInTrfFolder.name}, 'superstation.mat'));
if ~isempty(superstatfile)
    % get coordinate frames of that file to popupmenu
    handles=loadSuperstationFile(hObject, handles, ['../TRF/', dirsInTrfFolder(superstatfile).name]);
    handles.data.superstationFile=['../TRF/', dirsInTrfFolder(superstatfile).name];
    set(handles.text_parameters_refFrames_selected_superstation_file, 'String', handles.data.superstationFile);
    
%     % make vievsTrf default if available
%     logTRFFound=~cellfun(@isempty, strfind(get(handles.popupmenu_parameters_refFrames_superstationTRF, 'String'), 'vievsTrf'));
%     if sum(logTRFFound)>0
%         set(handles.popupmenu_parameters_refFrames_superstationTRF, 'Value', find(logTRFFound));
%     end
else
    msgbox('There was no superstation (.mat) file found in the TRF folder. Make sure you get one!', 'No superstation.mat file found in /TRF/', 'warn');
end

% try to find source file in CRF folder
supersourcefile=find(~cellfun(@isempty, strfind({dirsInCrfFolderMat.name}, 'supersource.mat')));
if ~isempty(supersourcefile)
    % get CRF frames of that file to popupmenu
    handles=loadSupersourceFile(hObject, handles, ['../CRF/', dirsInCrfFolderMat(supersourcefile).name]);
    handles.data.supersourceFile = ['../CRF/', dirsInCrfFolderMat(supersourcefile).name];
    set(handles.text_parameters_refFrames_selected_supersource_file, 'String', handles.data.supersourceFile);
    
    % make vievsCRF default if available
    logCRFFound=~cellfun(@isempty, strfind(get(handles.popupmenu_parameters_refFrames_supersourceCRF, 'String'), 'icrf3sx'));
    if sum(logCRFFound)>0
        set(handles.popupmenu_parameters_refFrames_supersourceCRF, 'Value', find(logCRFFound));
    end
else
    msgbox('There was no supersource (.mat) file found in the CRF folder. Make sure you get one!', 'No supersource.mat file found in /CRF/', 'warn');
end

% remove '.', '..' and entries which are not directories (eg files)
dirsInOptFolder(strcmp({dirsInOptFolder.name}, '.')|strcmp({dirsInOptFolder.name}, '..')|strcmp({dirsInOptFolder.name}, '.git')|~[dirsInOptFolder.isdir])=[];
yrsToDelete=~cellfun(@isnan, cellfun(@str2double, {dirsInOutlierFolder.name}, 'UniformOutput', false));
dirsInOutlierFolder(strcmp({dirsInOutlierFolder.name}, '.')|strcmp({dirsInOutlierFolder.name}, '..')|~[dirsInOutlierFolder.isdir]|yrsToDelete)=[];
dirsInTrpFolder(strcmp({dirsInTrpFolder.name}, '.')|strcmp({dirsInTrpFolder.name}, '..')|~[dirsInTrpFolder.isdir])=[];
dirsInIonFolder(strcmp({dirsInIonFolder.name}, '.')|strcmp({dirsInIonFolder.name}, '..')|~[dirsInIonFolder.isdir])=[];
dirsInAtmFolder(strcmp({dirsInAtmFolder.name}, '.')|strcmp({dirsInAtmFolder.name}, '..')|strcmp({dirsInAtmFolder.name}, 'temp')|~[dirsInAtmFolder.isdir])=[];
dirsInHydloFolder(strcmp({dirsInHydloFolder.name}, '.')|strcmp({dirsInHydloFolder.name}, '..')|~[dirsInHydloFolder.isdir])=[];
dirsInTrfFolder( strcmp({dirsInTrfFolder.name}, '.') | strcmp({dirsInTrfFolder.name}, '..') | strcmp({dirsInTrfFolder.name}, 'SavedGuiData_superstations.txt') | cellfun(@isempty, strfind({dirsInTrfFolder.name}, '.txt')) )=[]; % Exception for superstation GUI settup savings file ("SavedGuiData_superstations.txt")
dirsInCrfFolder(strcmp({dirsInCrfFolder.name}, '.')|strcmp({dirsInCrfFolder.name}, '..'))=[];
dirsInEopFolder(strcmp({dirsInEopFolder.name}, '.')|strcmp({dirsInEopFolder.name}, '..'))=[];
dirsInEophfFolder(strcmp({dirsInEophfFolder.name}, '.')|strcmp({dirsInEophfFolder.name}, '..'))=[];
dirsInDataFolder(strcmp({dirsInDataFolder.name}, '.')|strcmp({dirsInDataFolder.name}, '..')|~[dirsInDataFolder.isdir])=[];
dirsInTurbFolder([dirsInTurbFolder.isdir])=[]; % remove all folders (i just want .dat files)
dirsInWorkPlFolder([dirsInWorkPlFolder.isdir])=[];
dirsInGlobTrfDatum([dirsInGlobTrfDatum.isdir])=[];
dirsInGlobTrfReduce([dirsInGlobTrfReduce.isdir])=[];
dirsInGlobTrfDiscont([dirsInGlobTrfDiscont.isdir])=[];
dirsInGlobTrfVel([dirsInGlobTrfVel.isdir])=[];
dirsInGlobTrfVelTies([dirsInGlobTrfVelTies.isdir])=[];
dirsInGlobTrfAo([dirsInGlobTrfAo.isdir])=[];
dirsInGlobTrfSeason([dirsInGlobTrfSeason.isdir])=[];
dirsInGlobTrfAPLrg([dirsInGlobTrfAPLrg.isdir])=[];
dirsInGlobCrfDatum([dirsInGlobCrfDatum.isdir])=[];
dirsInGlobCrfFixed([dirsInGlobCrfFixed.isdir])=[];
dirsInGlobCrfReduce([dirsInGlobCrfReduce.isdir])=[];
dirsInSoucatFolder([dirsInSoucatFolder.isdir])=[];

% set new entries for popup menu
set(handles.popupmenu_setInput_optDir, 'String', {dirsInOptFolder.name})
if isempty({dirsInOutlierFolder.name})
    set(handles.popupmenu_setInput_outDir, 'string', ' ');
else
    set(handles.popupmenu_setInput_outDir, 'string', {'', dirsInOutlierFolder.name});
end

if isempty(dirsInIonFolder)
    set(handles.popupmenu_parameters_iono_ext, 'String', ' ')
else
    set(handles.popupmenu_parameters_iono_ext, 'String', {dirsInIonFolder.name})
end

if isempty(dirsInAtmFolder)
    set(handles.popupmenu_parameters_statCorr_nonTidalAtmoOceanLoad, 'String', ' ')
else
    set(handles.popupmenu_parameters_statCorr_nonTidalAtmoOceanLoad, 'String', {dirsInAtmFolder.name})
end
if isempty(dirsInTrfFolder)
    set(handles.popupmenu_parameters_refFrames_otherTRF, 'String', ' ')
    set(handles.popupmenu_parameters_refFrames_otherTRF, 'Enable', 'off')
    set(handles.radiobutton_parameters_refFrames_otherTRF, 'Enable', 'off')
else
    set(handles.popupmenu_parameters_refFrames_otherTRF, 'String', {dirsInTrfFolder.name})
    set(handles.popupmenu_parameters_refFrames_otherTRF, 'Enable', 'on')
    set(handles.radiobutton_parameters_refFrames_otherTRF, 'Enable', 'on')
end
if isempty(dirsInCrfFolder)
    set(handles.popupmenu_parameters_refFrames_otherCRF, 'String', ' ')
    set(handles.popupmenu_parameters_refFrames_otherCRF, 'Enable', 'off')
    set(handles.radiobutton_parameters_refFrames_otherCRF, 'Enable', 'off')
else
    set(handles.popupmenu_parameters_refFrames_otherCRF, 'String', {dirsInCrfFolder.name})
    set(handles.popupmenu_parameters_refFrames_otherCRF, 'Enable', 'on')
    set(handles.radiobutton_parameters_refFrames_otherCRF, 'Enable', 'on')
end
% set(handles.popupmenu_parameters_statCorr_tidalOceanLoad, 'String', {dirsInOtideFolder.name})
% set(handles.popupmenu_parameters_statCorr_tidalAtmoOceanLoad, 'String', {dirsInAtideFolder.name})
if isempty(dirsInHydloFolder)
    set(handles.popupmenu_parameters_statCorr_hydroLoading, 'String', ' ');
else
    set(handles.popupmenu_parameters_statCorr_hydroLoading, 'String', {dirsInHydloFolder.name});
end
set(handles.popupmenu_parameters_eop_aPriori_other, 'String', {dirsInEopFolder.name})
set(handles.popupmenu_parameters_eop_oceanTideModel, 'String', [{dirsInEophfFolder.name},'interpf (Conventions)','Combi_IGG_Bonn'])
set(handles.popupmenu_plot_folder1_subfolder, 'String', ['/', {dirsInDataFolder.name}])
set(handles.popupmenu_plot_folder2_subfolder, 'String', ['/', {dirsInDataFolder.name}])
set(handles.popupmenu_plot_folder3_subfolder, 'String', ['/', {dirsInDataFolder.name}])
set(handles.popupmenu_plot_residuals_folder,  'String', ['/', {dirsInDataFolder.name}])
set(handles.popupmenu_plot_sessionAnalysis_subfolder, 'String', ['/', {dirsInDataFolder.name}])
set(handles.popupmenu_plot_sessionAnalysis_subfolder2, 'String', ['/', {dirsInDataFolder.name}])
set(handles.popupmenu_plot_sessionAnalysis_subfolder3, 'String', ['/', {dirsInDataFolder.name}])
set(handles.popupmenu_plot_sessionAnalysis_subfolder4, 'String', ['/', {dirsInDataFolder.name}])
set(handles.popupmenu_plot_eopOut_subfolder, 'String', ['/', {dirsInDataFolder.name}])

if isempty(dirsInWorkPlFolder)
    set(handles.popupmenu_plot_eopOut_pl, 'String', ' ')
else
    set(handles.popupmenu_plot_eopOut_pl, 'String', {dirsInWorkPlFolder.name})
end

set(handles.listbox_vie_sim_paramFile, 'String', {dirsInTurbFolder.name})
if isempty(dirsInGlobTrfDatum)
    set(handles.popupmenu_vie_glob_trfCrf_trf_stations4Nnt, 'String', ' ')
else
    set(handles.popupmenu_vie_glob_trfCrf_trf_stations4Nnt, 'String', {dirsInGlobTrfDatum.name})
end
if isempty(dirsInGlobTrfReduce)
    set(handles.popupmenu_vie_glob_trfCrf_trf_sessionWiseRedStat, 'String', ' ')
else
    set(handles.popupmenu_vie_glob_trfCrf_trf_sessionWiseRedStat, 'String', {dirsInGlobTrfReduce.name})
end
if isempty(dirsInGlobTrfDiscont)
    set(handles.popupmenu_vie_glob_trfCrf_trf_positDisc, 'String', ' ')
else
    set(handles.popupmenu_vie_glob_trfCrf_trf_positDisc, 'String', {dirsInGlobTrfDiscont.name})
end
if isempty(dirsInGlobTrfVel)
    set(handles.popupmenu_vie_glob_trfCrf_trf_keepConstVel, 'String', ' ')
else
    set(handles.popupmenu_vie_glob_trfCrf_trf_keepConstVel, 'String', {dirsInGlobTrfVel.name})
end
if isempty(dirsInGlobTrfVelTies)
    set(handles.popupmenu_vie_glob_trfCrf_trf_velTies, 'String', ' ')
else
    set(handles.popupmenu_vie_glob_trfCrf_trf_velTies, 'String', {dirsInGlobTrfVelTies.name})
end
if isempty(dirsInGlobCrfDatum)
    set(handles.popupmenu_vie_glob_trfCrf_crf_sourcesForNnr, 'String', ' ')
else
    set(handles.popupmenu_vie_glob_trfCrf_crf_sourcesForNnr, 'String', {dirsInGlobCrfDatum.name})
end
if isempty(dirsInGlobCrfFixed)
    set(handles.popupmenu_vie_glob_trfCrf_crf_fixedSources, 'String', ' ')
else
    set(handles.popupmenu_vie_glob_trfCrf_crf_fixedSources, 'String', {dirsInGlobCrfFixed.name})
end
if isempty(dirsInGlobCrfReduce)
    set(handles.popupmenu_vie_glob_trfCrf_crf_sessionWiseRedSources, 'String', ' ')
else
    set(handles.popupmenu_vie_glob_trfCrf_crf_sessionWiseRedSources, 'String', {dirsInGlobCrfReduce.name})
end
if isempty(dirsInGlobTrfAo)
    set(handles.listbox_vie_glob_axisOffsets, 'String', ' ')
else
    set(handles.listbox_vie_glob_axisOffsets, 'String', {dirsInGlobTrfAo.name})
end
if isempty(dirsInGlobTrfSeason)
    set(handles.listbox_vie_glob_stseaspos, 'String', ' ')
else
    set(handles.listbox_vie_glob_stseaspos, 'String', {dirsInGlobTrfSeason.name})
end
if isempty(dirsInGlobTrfAPLrg)
    set(handles.listbox_vie_glob_APLrg, 'String', ' ')
else
    set(handles.listbox_vie_glob_APLrg, 'String', {dirsInGlobTrfAPLrg.name})
end
if isempty(dirsInSoucatFolder)
    set(handles.listbox_vie_sim_ss, 'String', ' ')
    set(handles.popupmenu_parameters_ss_catalog, 'String', ' ')
else
    set(handles.listbox_vie_sim_ss, 'String', {dirsInSoucatFolder.name});
    set(handles.popupmenu_parameters_ss_catalog, 'String', {dirsInSoucatFolder.name});
end
if isempty(dirsInStatlistFolder)
    set(handles.listbox_vie_sched_prenet, 'String', ' ')
else
    set(handles.listbox_vie_sched_prenet, 'String', {dirsInStatlistFolder.name})
end


% set default files for otide atide, a-non-tide:

% OTIDE
%defaultOtideFile='FES2004.mat';

% if we have found it -> set this one as selected
%if sum(strcmp({dirsInOtideFolder.name}, defaultOtideFile))==1
%   set(handles.popupmenu_parameters_statCorr_tidalOceanLoad, 'Value', find(strcmp({dirsInOtideFolder.name}, defaultOtideFile)))
%end
% if we have found 'FES2004.mat'
% if sum(strcmp({dirsInOtideFolder.name}, 'FES2004.mat'))==1
%     set(handles.popupmenu_parameters_statCorr_tidalOceanLoad, 'String', [dirsInOtideFolder(strcmp({dirsInOtideFolder.name}, 'FES2004.mat')).name])
% else
%     fprintf('Default tidal ocean loading file (FES2004.mat) not found!\nCorrection might be not used!\n--> Be sure to select a proper file!\n\n');
% end

% ATIDE
%defaultAtideFile='s12_cm_noib_leonid.mat';

% if we have found it -> set this one as selected
%if sum(strcmp({dirsInAtideFolder.name}, defaultAtideFile))==1
%    set(handles.popupmenu_parameters_statCorr_tidalAtmoOceanLoad, 'Value', find(strcmp({dirsInAtideFolder.name}, defaultAtideFile)))
%end


% % if we have found 's12_cm_noib_leonid.mat'
% if sum(strcmp({dirsInAtideFolder.name}, 's12_cm_noib_leonid.mat'))==1
%     set(handles.popupmenu_parameters_statCorr_tidalAtmoOceanLoad, 'String', [dirsInAtideFolder(strcmp({dirsInAtideFolder.name}, 's12_cm_noib_leonid.mat')).name])
% else
%     msgbox('Default tidal atmosphere loading file (s12_cm_noib_leonid.mat) not found!\nCorrection might be not used!\n--> Be sure to select a proper file!\n\n','file not found in default folder!','warn');
% end

% if default parameter file does not exist - create it
defaultParamFilename='PARAMETERS/default_fromGUI.mat';
if ~exist(defaultParamFilename, 'file')
    saveParamFile(hObject, handles, defaultParamFilename)
else % if it exists - delete it and create a new one
    delete(defaultParamFilename)
    saveParamFile(hObject, handles, defaultParamFilename)
end

% get vie_sched GUI startup options ++++++++++++++++++++++++
filename = '../CATALOGS/antenna.cat';
if exist(filename, 'file')
    fid = fopen(filename, 'r');
    stanum = 0;
    staname = '';
    while ~feof(fid)
        line = fgetl(fid);
        linelength = length(line);
        if ((linelength < 11) || (~strcmp(line(1), ' ')))   %%%
            continue;
        end
        stanum = stanum + 1;
        staname(stanum,1:8) = line(4:11);
    end
    fclose(fid);
    %
    selectstastr='';
    for i = 1 : stanum-1
        selectstastr = strcat(selectstastr,staname(i,1:8),'|');
    end
    selectstastr = strcat(selectstastr,staname(stanum,1:8));
    set(handles.listbox_vie_sched_selectsta,'String',selectstastr);
	set(handles.listbox_vie_sched_sat_obs_network_available,'String',selectstastr);
end
% get vie_sched GUI startup options ------------------------------

setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_welcome, 'Visible', 'On');
bc = get(handles.uipanel_welcome, 'BackgroundColor');

axes(handles.axes_logo);
set(gca, 'Color', 'none');

RGB = imread('../CODE/VIE_SETUP/vievs.png');
mask = sum(RGB,3)==765;
R = RGB(:,:,1);
R(mask)=bc(1)*255;

G = RGB(:,:,2);
G(mask)=bc(2)*255;

B = RGB(:,:,3);
B(mask)=bc(3)*255;
RGB = cat(3,R,G,B);

imshow(RGB)

try
    [status_hash,hash] = system('git rev-parse --short HEAD');
    hash = strtrim(hash);
    [status_tag,tag] = system('git describe --abbrev=0');
    tag = strtrim(tag);
    
    if status_tag == 0 && length(tag)>1
        f = gcf;
        f.Name = sprintf('Vienna VLBI and Satellite Software %s',tag);
    end
    
    if status_hash == 0 && status_tag== 0 && length(tag) > 1 && length(hash) == 7
        handles.text_version.String = sprintf('%s (%s)',tag, hash);
    elseif status_hash == 0 && length(hash) == 7
        handles.text_version.String = sprintf('%s', hash);
    elseif status_tag == 0 && length(tag) > 1
        handles.text_version.String = sprintf('%s',tag);
    else
        handles.text_version.String = 'no git';
    end
catch
    handles.text_version.String = 'no git';
end

% Choose default command line output for vie_setup
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vie_setup wait for user response (see UIRESUME)
% uiwait(handles.figure_vievs2);

%profile viewer


% --- Outputs from this function are returned to the command line.
function varargout = vie_setup_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_parameterFiles_loadCurrent_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_parameterFiles_loadCurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get automatically loaded file
allAutoSavedFiles=dir('../WORK/PARAMETERS/auto_save_fromGUI_*.mat');

% get sorting index
[temp, sortInd]=sort({allAutoSavedFiles.name});

% sort all found files according to index
allAutoSavedFiles=allAutoSavedFiles(sortInd);

% load the latest file
if isempty(allAutoSavedFiles)
    msgbox('No current parameter file exists. Use the GUI to create it automatically.', 'Current not available', 'warn');
else
    loadParamFile(hObject, handles, ['../WORK/PARAMETERS/', allAutoSavedFiles(end).name]);
    
    % message box for user
    msgbox('Current parameters have been loaded', 'Load current', 'help');
end




% --------------------------------------------------------------------
function menu_file_parameterFiles_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_parameterFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_vie_glob_special_param, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_file_reloadFolders_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_reloadFolders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

updateSeveralPopupmenus(hObject, handles)

% --------------------------------------------------------------------
function menu_help_Callback(hObject, eventdata, handles)
% hObject    handle to menu_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_outlier_simple.
function checkbox_outlier_simple_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_outlier_simple (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_outlier_simple


% --------------------------------------------------------------------
function Untitled_8_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in checkbox_run_sinex_write.
function checkbox_run_sinex_write_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_sinex_write (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_sinex_write

if get(hObject, 'Value')
    newState='On';
else
    newState='Off';
end
set(handles.text_run_sinex_header, 'Enable', newState)
set(handles.text_run_sinex_clockParam, 'Enable', newState)
set(handles.text_run_sinex_zwd, 'Enable', newState)
set(handles.text_run_sinex_tropoGradients, 'Enable', newState)
set(handles.text_run_sinex_sources, 'Enable', newState)
set(handles.text_run_sinex_stationCoords, 'Enable', newState)
set(handles.text_run_sinex_eop, 'Enable', newState)
set(handles.checkbox_run_sinex_sources, 'Enable', newState)
set(handles.text300, 'Enable', newState)

set(handles.radiobutton_run_sinex_clockParam_incl, 'Enable', 'off')
set(handles.radiobutton_run_sinex_zwd_incl, 'Enable', newState)
set(handles.radiobutton_run_sinex_tropoParam_incl, 'Enable', newState)
set(handles.radiobutton_run_sinex_sources_incl, 'Enable', newState)
set(handles.radiobutton_run_sinex_stationCoords_incl, 'Enable', newState)
set(handles.radiobutton_run_sinex_eop_incl, 'Enable', newState)

set(handles.radiobutton_run_sinex_clockParam_excl, 'Enable', newState)
set(handles.radiobutton_run_sinex_zwd_excl, 'Enable', newState)
set(handles.radiobutton_run_sinex_tropoParam_excl, 'Enable', newState)
set(handles.radiobutton_run_sinex_eop_excl, 'Enable', newState)

% if neither eop is estimtated -> disable EOP-options in sinex output
if (get(handles.checkbox_estimation_leastSquares_eop_xpEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_ypEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_dut1Est, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_nutdxEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_nutdyEst, 'Value'))==0
    set(handles.radiobutton_run_sinex_eop_incl, 'Enable', 'off');
    set(handles.radiobutton_run_sinex_eop_excl, 'Enable', 'off');
    set(handles.text_run_sinex_eop, 'Enable', 'off')
end

% if stations are not estimated -> disable option
if get(handles.checkbox_estimation_leastSquares_coordinates_estimate, 'Value')==0
    set(handles.text_run_sinex_stationCoords, 'Enable', 'off')
    set(handles.radiobutton_run_sinex_stationCoords_incl, 'Enable', 'off')
end

% if gradients are not estimated -> disable sinex option
if get(handles.checkbox_estimation_leastSquares_tropo_ngr, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_tropo_egr, 'Value')==0
    set(handles.text_run_sinex_tropoGradients, 'Enable', 'off')
    set(handles.radiobutton_run_sinex_tropoParam_incl, 'Enable', 'off')
    set(handles.radiobutton_run_sinex_tropoParam_excl, 'Enable', 'off')
end

% if zwd is not esimtated -> disable sinex option
if get(handles.checkbox_estimation_leastSquares_tropo_zwd, 'value')==0
    set(handles.text_run_sinex_zwd, 'Enable', 'off')
    set(handles.radiobutton_run_sinex_zwd_incl, 'Enable', 'off')
    set(handles.radiobutton_run_sinex_zwd_excl, 'Enable', 'off')
end

% if source checkbox is unticked
if get(handles.checkbox_run_sinex_sources, 'Value')==0
    set(handles.text_run_sinex_sources, 'Enable', 'off')
    set(handles.radiobutton_run_sinex_sources_incl, 'Enable', 'off')
end

set(handles.checkbox_run_sinex_changeAnalystsName, 'Enable', newState)
set(handles.edit_run_sinex_firstname, 'Enable', newState)
set(handles.edit_run_sinex_lastname, 'Enable', newState)
set(handles.edit_run_sinex_email, 'Enable', newState)
if get(handles.checkbox_run_sinex_changeAnalystsName, 'Value')==0
    set(handles.edit_run_sinex_firstname, 'Enable', 'Off')
    set(handles.edit_run_sinex_lastname, 'Enable', 'Off')
    set(handles.edit_run_sinex_email, 'Enable', 'Off')
end

set(handles.checkbox_run_sinex_addSuffix, 'Enable', newState)
set(handles.edit_run_sinex_suffix, 'Enable', newState)
if get(handles.checkbox_run_sinex_addSuffix, 'Value')==0
    set(handles.edit_run_sinex_suffix, 'Enable', 'Off')
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --------------------------------------------------------------------
function menu_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to menu_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_parameters_referenceFrames_Callback(hObject, eventdata, handles)
% hObject    handle to menu_parameters_referenceFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_parameters_referenceFrames, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_parameters_troposphere_Callback(hObject, eventdata, handles)
% hObject    handle to menu_parameters_troposphere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_parameters_troposphere, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_parameters_ionosphere_Callback(hObject, eventdata, handles)
% hObject    handle to menu_parameters_ionosphere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_parameters_ionosphere, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_help_vievsWebsite_Callback(hObject, eventdata, handles)
% hObject    handle to menu_help_vievsWebsite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

web 'http://vievs.geo.tuwien.ac.at/'


% --------------------------------------------------------------------
function menu_estimation_Callback(hObject, eventdata, handles)
% hObject    handle to menu_estimation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_15_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_globalSolution_Callback(hObject, eventdata, handles)
% hObject    handle to menu_globalSolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_19_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_20_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_22_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function menu_plotting_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to menu_plotting_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% shift the whole window a little bit down because the menu bar appears
curPos=get(handles.figure_vievs2, 'Position');
curPos(2)=curPos(2)-27;
set(handles.figure_vievs2, 'Position', curPos);

% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_plot_parameters, 'Visible', 'On');

% set bar for plotting options to visible
set(handles.uitoolbar_plotOptions, 'Visible', 'On');

% update the folders (since there could have been created another subfolder)
dirsInDataFolder=dir('../DATA/LEVEL3/');
dirsInDataFolder(strcmp({dirsInDataFolder.name}, '.')|strcmp({dirsInDataFolder.name}, '..')|~[dirsInDataFolder.isdir])=[];
set(handles.popupmenu_plot_folder1_subfolder, 'String', ['/', {dirsInDataFolder.name}])
set(handles.popupmenu_plot_folder2_subfolder, 'String', ['/', {dirsInDataFolder.name}])
set(handles.popupmenu_plot_folder3_subfolder, 'String', ['/', {dirsInDataFolder.name}])

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_run_sinexOutput_Callback(hObject, eventdata, handles)
% hObject    handle to menu_run_sinexOutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_run_sinexOutput, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_run_outputDirectories_Callback(hObject, eventdata, handles)
% hObject    handle to menu_run_outputDirectories (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_run_runOptions, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_run_runOptions_Callback(hObject, eventdata, handles)
% hObject    handle to menu_run_runOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_run_vievsProcSettings, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_estimation_leastSquares_Callback(hObject, eventdata, handles)
% hObject    handle to menu_estimation_leastSquares (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_estimation_kalmanFilter_Callback(hObject, eventdata, handles)
% hObject    handle to menu_estimation_kalmanFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_parameters_stationModels_Callback(hObject, eventdata, handles)
% hObject    handle to menu_parameters_stationModels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_parameters_stationCorrections, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function menu_parameters_eop_Callback(hObject, eventdata, handles)
% hObject    handle to menu_parameters_eop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_parameters_eop, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function Untitled_33_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_34_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_35_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_36_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_setInputFiles_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_setInputFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_file_setInputFiles, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_file_exit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(handles.figure_vievs2);

% --------------------------------------------------------------------
function menu_file_parameterFiles_loadDefaults_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_parameterFiles_loadDefaults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% define default parameter filename
defaultParamFilename='PARAMETERS/default_fromGUI.mat';
if exist(defaultParamFilename, 'file')
    loadParamFile(hObject, handles, defaultParamFilename)
    msgbox('Default parameters have been successfully loaded', 'Done', 'help');
else
    msgbox('No default parameters file found (must have been deleted!)', 'Error', 'warn');
end


% --------------------------------------------------------------------
function menu_file_parameterFiles_loadParameters_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_parameterFiles_loadParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get file from explorer
[FileName, PathName] = uigetfile({'*.mat', 'Matlab Binary Format (*.mat)'}, ...
    'Select a parameter file',...
    '../WORK/PARAMETERS/');

if ~isempty(FileName)
    
    if ~strcmp('0', num2str(FileName))
        loadParamFile(hObject, handles, [PathName, FileName])

        % write message box for information for user
        msgbox('Parameter file has been successfully loaded.', 'Load parameter file', 'help')
    end
end
        

% --------------------------------------------------------------------
function menu_file_parameterFiles_saveParameters_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_parameterFiles_saveParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get a file
[filename, pathname, filterindex] = uiputfile( ...
{'*.mat', 'Matlab Binary Format (*.mat)'},...
 'Save as',...
 '../WORK/PARAMETERS/');

% see if something was chosen
if ~isempty(filename)
    if ~strcmp('0', num2str(filename))
        % call function for saving a parameter file
        saveParamFile(hObject, handles, [pathname, filename])
        
        % write message box
        msgbox('Parameter file successfully saved!', 'Save parameter file', 'help');
    end
end



% --------------------------------------------------------------------
function menu_estimation_leastSquares_troposphere_Callback(hObject, eventdata, handles)
% hObject    handle to menu_estimation_leastSquares_troposphere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_estimation_leastSquares_tropo, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_estimation_leastSquares_clock_Callback(hObject, eventdata, handles)
% hObject    handle to menu_estimation_leastSquares_clock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_estimation_leastSquares_clock, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_estimation_leastSquares_eop_Callback(hObject, eventdata, handles)
% hObject    handle to menu_estimation_leastSquares_eop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_estimation_leastSquares_eop, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_estimation_leastSquares_stationCoordinates_Callback(hObject, eventdata, handles)
% hObject    handle to menu_estimation_leastSquares_stationCoordinates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_estimation_leastSquares_stationCoordinates, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_estimation_leastSquares_sourceCoordinates_Callback(hObject, eventdata, handles)
% hObject    handle to menu_estimation_leastSquares_sourceCoordinates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_estimation_leastSquares, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function Untitled_50_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_parameters_observationRestrictions_Callback(hObject, eventdata, handles)
% hObject    handle to menu_parameters_observationRestrictions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_parameters_observationRestrictions, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function menu_parameters_ephemerides_Callback(hObject, eventdata, handles)
% hObject    handle to menu_parameters_ephemerides (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_parameters_ephemerides, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);


% --- Executes on selection change in listbox_sessions.
function listbox_sessions_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_sessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_sessions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_sessions


% --- Executes during object creation, after setting all properties.
function listbox_sessions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_sessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listbox6.
function listbox6_Callback(hObject, eventdata, handles)
% hObject    handle to listbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox6


% --- Executes during object creation, after setting all properties.
function listbox6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox7.
function listbox7_Callback(hObject, eventdata, handles)
% hObject    handle to listbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox7


% --- Executes during object creation, after setting all properties.
function listbox7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox49.
function checkbox49_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox49


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox11.
function listbox11_Callback(hObject, eventdata, handles)
% hObject    handle to listbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox11 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox11


% --- Executes during object creation, after setting all properties.
function listbox11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox12.
function listbox12_Callback(hObject, eventdata, handles)
% hObject    handle to listbox12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox12 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox12


% --- Executes during object creation, after setting all properties.
function listbox12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox13.
function listbox13_Callback(hObject, eventdata, handles)
% hObject    handle to listbox13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox13 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox13


% --- Executes during object creation, after setting all properties.
function listbox13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

keyboard;


% --- Executes on selection change in popupmenu9.
function popupmenu9_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu9 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu9


% --- Executes during object creation, after setting all properties.
function popupmenu9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu12.
function popupmenu12_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu12 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu12


% --- Executes during object creation, after setting all properties.
function popupmenu12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu13.
function popupmenu13_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu13 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu13


% --- Executes during object creation, after setting all properties.
function popupmenu13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu14.
function popupmenu14_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu14 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu14


% --- Executes during object creation, after setting all properties.
function popupmenu14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox51.
function checkbox51_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox51


% --- Executes on button press in checkbox52.
function checkbox52_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox52


% --- Executes on button press in checkbox53.
function checkbox53_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox53


% --- Executes on button press in checkbox54.
function checkbox54_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox54


% --- Executes on button press in checkbox55.
function checkbox55_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox55


% --- Executes on button press in checkbox56.
function checkbox56_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox56


% --- Executes during object creation, after setting all properties.
function uipanel_setInputFiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_setInputFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in checkbox58.
function checkbox58_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox58


% --- Executes on button press in checkbox59.
function checkbox59_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox59


% --- Executes on button press in checkbox60.
function checkbox60_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox60 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox60


% --- Executes on button press in checkbox61.
function checkbox61_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox61 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox61


% --- Executes on button press in checkbox62.
function checkbox62_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox62 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox62


% --- Executes on button press in checkbox63.
function checkbox63_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox63 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox63


% --- Executes on button press in radiobutton_parameters_statCorr_poleModel_lin.
function radiobutton_parameters_statCorr_poleModel_lin_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_parameters_statCorr_poleModel_lin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_parameters_statCorr_poleModel_lin


% --- Executes on button press in radiobutton_parameters_statCorr_poleModel_cub.
function radiobutton_parameters_statCorr_poleModel_cub_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_parameters_statCorr_poleModel_cub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_parameters_statCorr_poleModel_cub


% --- Executes on button press in checkbox_parameters_statCorr_solidEarthTides.
function checkbox_parameters_statCorr_solidEarthTides_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameters_statCorr_solidEarthTides (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameters_statCorr_solidEarthTides

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)

% --- Executes on button press in checkbox_parameters_statCorr_gravitationalDef_Callback.
function checkbox_parameters_statCorr_gravitationalDef_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameters_statCorr_solidEarthTides (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameters_statCorr_gravitationalDef_Callback

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)

% --- Executes on button press in checkbox_parameters_statCorr_tidalOceanLoad.
function checkbox_parameters_statCorr_tidalOceanLoad_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameters_statCorr_tidalOceanLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameters_statCorr_tidalOceanLoad

if get(hObject, 'Value')
    newState='on';
else
    newState='off';
end

set(handles.popupmenu_parameters_statCorr_tidalOceanLoad, 'Enable', newState)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox_parameters_statCorr_tidalAtmoLoad.
function checkbox_parameters_statCorr_tidalAtmoLoad_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameters_statCorr_tidalAtmoLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameters_statCorr_tidalAtmoLoad

if get(hObject, 'Value')
    newState='on';
else
    newState='off';
end

set(handles.popupmenu_parameters_statCorr_tidalAtmoOceanLoad, 'Enable', newState)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox_parameters_statCorr_nonTidalAtmoLoad.
function checkbox_parameters_statCorr_nonTidalAtmoLoad_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameters_statCorr_nonTidalAtmoLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameters_statCorr_nonTidalAtmoLoad

if get(hObject, 'Value')
    newState='on';
else
    newState='off';
end

set(handles.popupmenu_parameters_statCorr_nonTidalAtmoOceanLoad, 'Enable', newState)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox_parameters_statCorr_poleTides.
function checkbox_parameters_statCorr_poleTides_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameters_statCorr_poleTides (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameters_statCorr_poleTides

if get(hObject, 'Value')
    set(handles.radiobutton_parameters_statCorr_poleModel_lin, 'Enable', 'on')
    set(handles.radiobutton_parameters_statCorr_poleModel_cub, 'Enable', 'on')
    set(handles.radiobutton_parameters_statCorr_poleModel_iers2015, 'Enable', 'on')
else
    % only set to disable when also other checkbox is 0
    if get(handles.checkbox_parameters_statCorr_oceanPoleTides, 'Value')==0
        set(handles.radiobutton_parameters_statCorr_poleModel_lin, 'Enable', 'off')
        set(handles.radiobutton_parameters_statCorr_poleModel_cub, 'Enable', 'off')
        set(handles.radiobutton_parameters_statCorr_poleModel_iers2015, 'Enable', 'off')
    end
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox_parameters_statCorr_thermalDef.
function checkbox_parameters_statCorr_thermalDef_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameters_statCorr_thermalDef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameters_statCorr_thermalDef

if get(hObject, 'Value')
    set(handles.checkbox_parameters_statCorr_temp_fromInSitu, 'Enable', 'on')
    set(handles.checkbox_parameters_statCorr_temp_GPT3, 'Enable', 'on')
else
    % only set to disable when also other checkbox is 0
    if get(handles.checkbox_parameters_statCorr_thermalDef, 'Value')==0
        set(handles.checkbox_parameters_statCorr_temp_fromInSitu, 'Enable', 'off')
        set(handles.checkbox_parameters_statCorr_temp_GPT3, 'Enable', 'off')
    end
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on selection change in listbox_setInput_processList.
function listbox_setInput_processList_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_setInput_processList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_setInput_processList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_setInput_processList


% --- Executes during object creation, after setting all properties.
function listbox_setInput_processList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_setInput_processList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_setInput_browseForSessions.
function pushbutton_setInput_browseForSessions_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_setInput_browseForSessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName, PathName] = uigetfile('*.*','Select VLBI sessions', '../DATA/NGS', 'multiselect', 'on');

if ischar(FileName) || iscell(FileName)

    FileName = cellstr(FileName);

    fileWasChosen=1;

    if isempty(FileName)
        fileWasChosen=0;
    elseif ~iscell(FileName)
        if FileName == 0
            fileWasChosen=0;
        end
    end

    if fileWasChosen
        % Distinguish between NGS and VSO (VieVS simple observation files) files:
        if ~isempty(strfind(PathName, 'NGS'))
            filetype_str = 'ngs';
        elseif ~isempty(strfind(PathName, 'SIM')) || ~isempty(strfind(PathName, 'SCHED'))
            % Check file name:
            switch(FileName{1}(end-3:end))
                case '.vso'
                    filetype_str = 'vso';
                otherwise % => NGS
                    filetype_str = 'ngs';
            end
        elseif ~isempty(strfind(PathName, 'VSO'))
            filetype_str = 'vso';
        else
            fileWasChosen = 0;
            fprintf('ERROR: Invalid filepath for input data! Valid fileapths: <VieVS root>/DATA/NGS/, <VieVS root>/DATA/VSO/, <VieVS root>/DATA/SIM/ or <VieVS root>/DATA/SCHED/\n');
        end
    end % if fileWasChosen

    if fileWasChosen

        % make relative paths from PathName
        switch(filetype_str)
            case 'ngs'
                subfolder_name_str = 'NGS';
            case 'vso'
                subfolder_name_str = 'VSO';
                % Add ' [VSO]' to filename
                for i = 1 : length(FileName)
                    FileName{i} = [FileName{i}, ' [VSO]'];
                end
        end % switch(filetype_str)

        % remove part until /DATA/<subfolder_name_str>/ from PathName
        if ispc
            indOfNGSfolder=strfind(PathName, ['DATA\', subfolder_name_str]) + 6 + length(subfolder_name_str);
        else
            indOfNGSfolder=strfind(PathName, ['DATA/', subfolder_name_str]) + 6 + length(subfolder_name_str);
        end
        if ~isempty(indOfNGSfolder)
            PathName=PathName(indOfNGSfolder:end);
        end

        % get newly selected file(s) to one cellstr
        newFiles=cellstr([(repmat(PathName, size(cellstr(FileName),2),1)), char(FileName)]);

        % update listbox
        updateInputFilesBox(hObject, eventdata,handles,newFiles)

        % save parameter file automatically 
        auto_save_parameterfile(hObject, handles)

    end
    
end

% Function to update input-files-listbox
function updateInputFilesBox(hObject, eventdata,handles,newFiles)

% if something was given
if ~isempty(newFiles) && sum(strcmpi(newFiles,''))==0 % second: if cell, returned is [0 0]
    newFiles=cellstr(newFiles);
    
    % \ -> /
    newFiles=strrep(newFiles, '\', '/');
    
    % if the popupmenu is empty, just take the new sessions
    if isempty(get(handles.listbox_setInput_processList, 'String'))
        files4listbox=newFiles;
        %allSelectedFilenames=FileName;
    else
        files4listbox=unique([handles.allSelectedFiles; newFiles]);
    end

    % \ -> /
    files4listbox=strrep(files4listbox, '\', '/');

    % set only filename (for clear view) to listbox string
    set(handles.listbox_setInput_processList, 'String', files4listbox);

    % save all selected sessions to handles struct
    handles.allSelectedFiles=files4listbox;

    % save changes to handles struct
    guidata(hObject, handles);
end
    


% --- Executes on selection change in listbox19.
function listbox19_Callback(hObject, eventdata, handles)
% hObject    handle to listbox19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox19 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox19


% --- Executes during object creation, after setting all properties.
function listbox19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox20.
function listbox20_Callback(hObject, eventdata, handles)
% hObject    handle to listbox20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox20 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox20


% --- Executes during object creation, after setting all properties.
function listbox20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox21.
function listbox21_Callback(hObject, eventdata, handles)
% hObject    handle to listbox21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox21 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox21


% --- Executes during object creation, after setting all properties.
function listbox21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_setInput_optDir.
function popupmenu_setInput_optDir_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_setInput_optDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_setInput_optDir contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_setInput_optDir

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu_setInput_optDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_setInput_optDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_setInput_outDir.
function popupmenu_setInput_outDir_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_setInput_outDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_setInput_outDir contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_setInput_outDir

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu_setInput_outDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_setInput_outDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_setInput_eliminOutliers.
function checkbox_setInput_eliminOutliers_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_setInput_eliminOutliers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_setInput_eliminOutliers

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)



function edit_parameter_obsRestr_qualityCode_Callback(hObject, eventdata, handles)
% hObject    handle to edit_parameter_obsRestr_qualityCode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_parameter_obsRestr_qualityCode as text
%        str2double(get(hObject,'String')) returns contents of edit_parameter_obsRestr_qualityCode as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_parameter_obsRestr_qualityCode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_parameter_obsRestr_qualityCode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_parameter_obsRestr_cutOff_Callback(hObject, eventdata, handles)
% hObject    handle to edit_parameter_obsRestr_cutOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_parameter_obsRestr_cutOff as text
%        str2double(get(hObject,'String')) returns contents of edit_parameter_obsRestr_cutOff as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_parameter_obsRestr_cutOff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_parameter_obsRestr_cutOff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_parameters_eop_models_inclAPrioriNutOffs.
function checkbox_parameters_eop_models_inclAPrioriNutOffs_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameters_eop_models_inclAPrioriNutOffs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameters_eop_models_inclAPrioriNutOffs

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)

% --- Executes on button press in checkbox_parameters_eop_inclHf_LibrationXpYp.
function checkbox_parameters_eop_inclHf_LibrationXpYp_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameters_eop_inclHf_LibrationXpYp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameters_eop_inclHf_LibrationXpYp

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox76.
function checkbox76_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox76 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox76


% --- Executes on button press in checkbox77.
function checkbox77_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox77 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox77


% --- Executes on button press in checkbox78.
function checkbox78_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox78 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox78


% --- Executes on button press in checkbox79.
function checkbox79_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox79 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox79


% --- Executes on selection change in popupmenu19.
function popupmenu19_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu19 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu19


% --- Executes during object creation, after setting all properties.
function popupmenu19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_parameters_eop_tidalUtVariations.
function checkbox_parameters_eop_tidalUtVariations_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameters_eop_tidalUtVariations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameters_eop_tidalUtVariations
% if lagrange is used -> enable "tidal UT1 variations"
if get(handles.checkbox_parameters_eop_tidalUtVariations,'Value')==1
    set(handles.rb_eop_inter_UT1R, 'Enable', 'on')
    set(handles.rb_eop_inter_UT1S, 'Enable', 'on')
else
    set(handles.rb_eop_inter_UT1R, 'Enable', 'off')
    set(handles.rb_eop_inter_UT1S, 'Enable', 'off')
end
% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on selection change in popupmenu_parameters_eop_oceanTideModel.
function popupmenu_parameters_eop_oceanTideModel_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_eop_oceanTideModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_parameters_eop_oceanTideModel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_parameters_eop_oceanTideModel

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu_parameters_eop_oceanTideModel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_eop_oceanTideModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_parameters_eop_inclHf_LibrationUt1.
function checkbox_parameters_eop_inclHf_LibrationUt1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameters_eop_inclHf_LibrationUt1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameters_eop_inclHf_LibrationUt1

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox86.
function checkbox86_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox86 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox86


% --- Executes on button press in checkbox87.
function checkbox87_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox87 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox87


% --- Executes on selection change in popupmenu22.
function popupmenu22_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu22 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu22


% --- Executes during object creation, after setting all properties.
function popupmenu22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox91.
function checkbox91_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox91 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox91


% --- Executes on button press in checkbox_run_runMainSolution.
function checkbox_run_runMainSolution_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_runMainSolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_runMainSolution

if get(hObject, 'Value')
    set(handles.checkbox_run_simpleOutlierTest, 'Enable', 'on')
    set(handles.checkbox_run_normalOutlierTest, 'Enable', 'on')
    if get(handles.checkbox_run_simpleOutlierTest, 'Value') || get(handles.checkbox_run_normalOutlierTest, 'Value') 
        set(handles.text_run_runOptions_c, 'Enable', 'on')
        set(handles.edit_run_outlierTestC, 'Enable', 'on')
    end
else
    set(handles.checkbox_run_simpleOutlierTest, 'Enable', 'off')
    set(handles.checkbox_run_normalOutlierTest, 'Enable', 'off')
    set(handles.text_run_runOptions_c, 'Enable', 'off')
    set(handles.edit_run_outlierTestC, 'Enable', 'off')
end
    
% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox_run_simpleOutlierTest.
function checkbox_run_simpleOutlierTest_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_simpleOutlierTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_simpleOutlierTest
if get(hObject, 'Value')
    set(handles.text_run_runOptions_c, 'Enable', 'on')
    set(handles.edit_run_outlierTestC, 'Enable', 'on')
    set(handles.checkbox_run_normalOutlierTest, 'Value', 0)
else
    % only if also other checkbox is off
    if get(handles.checkbox_run_normalOutlierTest, 'Value')==0
        set(handles.text_run_runOptions_c, 'Enable', 'off')
        set(handles.edit_run_outlierTestC, 'Enable', 'off')
    end
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox_run_normalOutlierTest.
function checkbox_run_normalOutlierTest_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_normalOutlierTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_normalOutlierTest

if get(hObject, 'Value')
    set(handles.text_run_runOptions_c, 'Enable', 'on')
    set(handles.edit_run_outlierTestC, 'Enable', 'on')
    set(handles.checkbox_run_simpleOutlierTest, 'Value', 0)
else
    % only if also other checkbox is off
    if get(handles.checkbox_run_simpleOutlierTest, 'Value') == 0
        set(handles.text_run_runOptions_c, 'Enable', 'off')
        set(handles.edit_run_outlierTestC, 'Enable', 'off')
    end
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)



function edit_run_outlierTestC_Callback(hObject, eventdata, handles)
% hObject    handle to edit_run_outlierTestC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_run_outlierTestC as text
%        str2double(get(hObject,'String')) returns contents of edit_run_outlierTestC as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_run_outlierTestC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_run_outlierTestC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_run_runFirstSolution.
function checkbox_run_runFirstSolution_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_runFirstSolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_runFirstSolution

if get(hObject, 'Value')
    newState='On';
else
    newState='off';
end

set(handles.radiobutton_run_oneOffsPerClock, 'Enable', newState)
set(handles.radiobutton_run_oneOffsAndRatePerClock, 'Enable', newState)
set(handles.radiobutton_run_oneOffsAndRateAndQuPerClock, 'Enable', newState)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox88.
function checkbox88_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox88 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox88


% --- Executes on button press in checkbox89.
function checkbox89_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox89 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox89


% --- Executes on selection change in popupmenu23.
function popupmenu23_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu23 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu23


% --- Executes during object creation, after setting all properties.
function popupmenu23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox90.
function checkbox90_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox90


% --- Executes on button press in checkbox99.
function checkbox99_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox99 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox99


% --- Executes on button press in checkbox_estimation_leastSquares_clocks.
function checkbox_estimation_leastSquares_clocks_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_clocks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_clocks

% set enabling
if get(hObject, 'Value')
    set(handles.radiobutton_estimation_leastSquares_pwlClock, 'Enable', 'on');
    set(handles.radiobutton_estimation_leastSquares_pwlAndRateClock, 'Enable', 'on');
    set(handles.radiobutton_estimation_leastSquares_pwlRateAndQuClock, 'Enable', 'on');
    set(handles.text_estimation_leastSquares_clockInt, 'Enable', 'on');
    set(handles.edit_estimation_leastSquares_clockInterval, 'Enable', 'on');
    set(handles.checkbox_estimation_leastSquares_clocksRelConstr, 'Enable', 'on');
    if get(handles.checkbox_estimation_leastSquares_clocksRelConstr, 'Value')
        set(handles.text_estimation_leastSquares_clockConstr, 'Enable', 'on');
        set(handles.edit_estimation_leastSquares_clockConstr, 'Enable', 'on');
        set(handles.text308, 'Enable', 'on')
    end
else
    set(handles.radiobutton_estimation_leastSquares_pwlClock, 'Enable', 'off');
    set(handles.radiobutton_estimation_leastSquares_pwlAndRateClock, 'Enable', 'off');
    set(handles.radiobutton_estimation_leastSquares_pwlRateAndQuClock, 'Enable', 'off');
    set(handles.text_estimation_leastSquares_clockInt, 'Enable', 'off');
    set(handles.edit_estimation_leastSquares_clockInterval, 'Enable', 'off');
    set(handles.checkbox_estimation_leastSquares_clocksRelConstr, 'Enable', 'off');
    set(handles.text_estimation_leastSquares_clockConstr, 'Enable', 'off');
    set(handles.edit_estimation_leastSquares_clockConstr, 'Enable', 'off');
    set(handles.text308, 'Enable', 'off')
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox_estimation_leastSquares_clocksBasDepOffset.
function checkbox_estimation_leastSquares_clocksBasDepOffset_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_clocksBasDepOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_clocksBasDepOffset

% set enabling
if get(hObject, 'Value')
    if get(handles.checkbox_setInput_useOptFiles, 'Value')
        set(handles.radiobutton_estimation_leastSquares_basdepClockoff_OPT, 'Enable', 'on');
    end
    set(handles.radiobutton_estimation_leastSquares_basdepClockoff_automatic, 'Enable', 'on');
    set(handles.edit_estimation_leastSquares_clockBasDepO_minNobs, 'Enable', 'on');
    set(handles.text403, 'Enable', 'on');
    set(handles.text404, 'Enable', 'on');
else
    set(handles.radiobutton_estimation_leastSquares_basdepClockoff_OPT, 'Enable', 'off');
    set(handles.radiobutton_estimation_leastSquares_basdepClockoff_automatic, 'Enable', 'off');
    set(handles.edit_estimation_leastSquares_clockBasDepO_minNobs, 'Enable', 'off');
    set(handles.text403, 'Enable', 'off');
    set(handles.text404, 'Enable', 'off');

end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)



% --- Executes on button press in checkbox_estimation_leastSquares_clocksRelConstr.
function checkbox_estimation_leastSquares_clocksRelConstr_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_clocksRelConstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_clocksRelConstr

if get(hObject, 'Value')
    set(handles.edit_estimation_leastSquares_clockConstr, 'Enable', 'on');
    set(handles.text_estimation_leastSquares_clockConstr, 'Enable', 'on');
    set(handles.text308, 'Enable', 'on')
    
else
    set(handles.edit_estimation_leastSquares_clockConstr, 'Enable', 'off');
    set(handles.text_estimation_leastSquares_clockConstr, 'Enable', 'off');
    set(handles.text308, 'Enable', 'off')
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


function edit_estimation_leastSquares_clockConstr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_clockConstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_clockConstr as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_clockConstr as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_clockConstr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_clockConstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_estimation_leastSquares_clockInterval_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_clockInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_clockInterval as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_clockInterval as a double

% save parameter file automatically 
% set the minutes value also to the static text
set(handles.text308, 'String', sprintf('after %s minutes', get(hObject, 'String')))

auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_clockInterval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_clockInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_estimation_leastSquares_tropo_ngr.
function checkbox_estimation_leastSquares_tropo_ngr_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_tropo_ngr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_tropo_ngr

if get(hObject, 'Value')
    set(handles.text_estimation_leastSquares_tropo_ngrInt, 'Enable', 'on')
    set(handles.edit_estimation_leastSquares_tropo_ngrInterval, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_tropo_ngrRelConstr, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_tropo_ngrAbsConstr, 'Enable', 'on')
    if get(handles.checkbox_estimation_leastSquares_tropo_ngrRelConstr, 'Value')
        set(handles.text_estimation_leastSquares_tropo_ngrRelConstr, 'Enable', 'on')
        set(handles.edit_estimation_leastSquares_tropo_ngrRelConstr, 'Enable', 'on')
    end
    if get(handles.checkbox_estimation_leastSquares_tropo_ngrAbsConstr, 'Value')
        set(handles.text_estimation_leastSquares_tropo_ngrAbsConstr, 'Enable', 'on')
        set(handles.edit_estimation_leastSquares_tropo_ngrAbsConstr, 'Enable', 'on')
    end
    % if sinex is written -> update option
    if get(handles.checkbox_run_sinex_write, 'Value')
        set(handles.text_run_sinex_tropoGradients, 'Enable', 'on')
        set(handles.radiobutton_run_sinex_tropoParam_incl, 'Enable', 'on')
        set(handles.radiobutton_run_sinex_tropoParam_excl, 'Enable', 'on')
    end
else
    set(handles.text_estimation_leastSquares_tropo_ngrInt, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_tropo_ngrInterval, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_tropo_ngrRelConstr, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_tropo_ngrAbsConstr, 'Enable', 'off')
    set(handles.text_estimation_leastSquares_tropo_ngrRelConstr, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_tropo_ngrRelConstr, 'Enable', 'off')
    set(handles.text_estimation_leastSquares_tropo_ngrAbsConstr, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_tropo_ngrAbsConstr, 'Enable', 'off')
    
    % if also egr is not estimated -> disable sinex option
    if get(handles.checkbox_estimation_leastSquares_tropo_egr, 'Value')==0
        set(handles.text_run_sinex_tropoGradients, 'Enable', 'off')
        set(handles.radiobutton_run_sinex_tropoParam_incl, 'Enable', 'off')
        set(handles.radiobutton_run_sinex_tropoParam_excl, 'Enable', 'off')
    end

end

% save parameter file automatically
auto_save_parameterfile(hObject, handles)   


% --- Executes on button press in checkbox_estimation_leastSquares_tropo_ngrRelConstr.
function checkbox_estimation_leastSquares_tropo_ngrRelConstr_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_tropo_ngrRelConstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_tropo_ngrRelConstr

if get(hObject, 'Value')
    newState='on';
else
    newState='off';
end

set(handles.text_estimation_leastSquares_tropo_ngrRelConstr, 'Enable', newState)
set(handles.edit_estimation_leastSquares_tropo_ngrRelConstr, 'Enable', newState)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


function edit_estimation_leastSquares_tropo_ngrRelConstr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_tropo_ngrRelConstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_tropo_ngrRelConstr as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_tropo_ngrRelConstr as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_tropo_ngrRelConstr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_tropo_ngrRelConstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_estimation_leastSquares_tropo_ngrInterval_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_tropo_ngrInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_tropo_ngrInterval as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_tropo_ngrInterval as a double

% set two static texts
set(handles.text310, 'String', sprintf('after %s minutes', get(hObject, 'String')))

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_tropo_ngrInterval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_tropo_ngrInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_estimation_leastSquares_tropo_ngrAbsConstr.
function checkbox_estimation_leastSquares_tropo_ngrAbsConstr_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_tropo_ngrAbsConstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_tropo_ngrAbsConstr

if get(hObject, 'Value')
    newState='on';
else
    newState='off';
end
set(handles.text_estimation_leastSquares_tropo_ngrAbsConstr, 'Enable', newState)
set(handles.edit_estimation_leastSquares_tropo_ngrAbsConstr, 'Enable', newState)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


function edit_estimation_leastSquares_tropo_ngrAbsConstr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_tropo_ngrAbsConstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_tropo_ngrAbsConstr as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_tropo_ngrAbsConstr as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_tropo_ngrAbsConstr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_tropo_ngrAbsConstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_estimation_leastSquares_tropo_zwd.
function checkbox_estimation_leastSquares_tropo_zwd_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_tropo_zwd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_tropo_zwd

if get(hObject, 'Value')
    set(handles.text_estimation_leastSquares_tropo_zwdInt, 'Enable', 'on')
    set(handles.edit_estimation_leastSquares_tropo_zwdInterval, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_tropo_zwdRelConstr, 'Enable', 'on')
    if get(handles.checkbox_estimation_leastSquares_tropo_zwdRelConstr, 'Value')
        set(handles.text_estimation_leastSquares_tropo_zwdRelConstr, 'Enable', 'on')
        set(handles.edit_estimation_leastSquares_tropo_zwdRelConstr, 'Enable', 'on')
    end
    if get(handles.checkbox_run_sinex_write, 'Value')
        newSinexState='on';
    else
        newSinexState='off';
    end
else
    set(handles.text_estimation_leastSquares_tropo_zwdInt, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_tropo_zwdInterval, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_tropo_zwdRelConstr, 'Enable', 'off')
    set(handles.text_estimation_leastSquares_tropo_zwdRelConstr, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_tropo_zwdRelConstr, 'Enable', 'off')
    
    newSinexState='off';
end

set(handles.text_run_sinex_zwd, 'Enable', newSinexState)
set(handles.radiobutton_run_sinex_zwd_incl, 'Enable', newSinexState)
set(handles.radiobutton_run_sinex_zwd_excl, 'Enable', newSinexState)
        
% save parameter file automatically 
auto_save_parameterfile(hObject, handles)   


% --- Executes on button press in checkbox_estimation_leastSquares_tropo_zwdRelConstr.
function checkbox_estimation_leastSquares_tropo_zwdRelConstr_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_tropo_zwdRelConstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_tropo_zwdRelConstr

if get(hObject, 'Value')
    newState='on';
else
    newState='off';
end
set(handles.text_estimation_leastSquares_tropo_zwdRelConstr, 'Enable', newState);
set(handles.edit_estimation_leastSquares_tropo_zwdRelConstr, 'Enable', newState);

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


function edit_estimation_leastSquares_tropo_zwdRelConstr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_tropo_zwdRelConstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_tropo_zwdRelConstr as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_tropo_zwdRelConstr as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_tropo_zwdRelConstr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_tropo_zwdRelConstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_estimation_leastSquares_tropo_zwdInterval_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_tropo_zwdInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_tropo_zwdInterval as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_tropo_zwdInterval as a double

% change static text below
set(handles.text309, 'String', sprintf('after %s minutes', get(hObject, 'String')))
% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_tropo_zwdInterval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_tropo_zwdInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_estimation_leastSquares_tropo_egr.
function checkbox_estimation_leastSquares_tropo_egr_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_tropo_egr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_tropo_egr

if get(hObject, 'Value')
    set(handles.text_estimation_leastSquares_tropo_egrInt, 'Enable', 'on')
    set(handles.edit_estimation_leastSquares_tropo_egrInterval, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_tropo_egrRelConstr, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_tropo_egrAbsConstr, 'Enable', 'on')
    if get(handles.checkbox_estimation_leastSquares_tropo_egrRelConstr, 'Value')
        set(handles.text_estimation_leastSquares_tropo_egrRelConstr, 'Enable', 'on')
        set(handles.edit_estimation_leastSquares_tropo_egrRelConstr, 'Enable', 'on')
    end
    if get(handles.checkbox_estimation_leastSquares_tropo_egrAbsConstr, 'Value')
        set(handles.text_estimation_leastSquares_tropo_egrAbsConstr, 'Enable', 'on')
        set(handles.edit_estimation_leastSquares_tropo_egrAbsConstr, 'Enable', 'on')
    end
    
    % if sinex is written -> enable this option
    if get(handles.checkbox_run_sinex_write, 'Value')
        set(handles.text_run_sinex_tropoGradients, 'Enable', 'on')
        set(handles.radiobutton_run_sinex_tropoParam_incl, 'Enable', 'on')
        set(handles.radiobutton_run_sinex_tropoParam_excl, 'Enable', 'on')
    end
else
    set(handles.text_estimation_leastSquares_tropo_egrInt, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_tropo_egrInterval, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_tropo_egrRelConstr, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_tropo_egrAbsConstr, 'Enable', 'off')
    set(handles.text_estimation_leastSquares_tropo_egrRelConstr, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_tropo_egrRelConstr, 'Enable', 'off')
    set(handles.text_estimation_leastSquares_tropo_egrAbsConstr, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_tropo_egrAbsConstr, 'Enable', 'off')
    
    % if also ngr is not estimated -> disable sinex option
    if get(handles.checkbox_estimation_leastSquares_tropo_ngr, 'Value')==0
        set(handles.text_run_sinex_tropoGradients, 'Enable', 'off')
        set(handles.radiobutton_run_sinex_tropoParam_incl, 'Enable', 'off')
        set(handles.radiobutton_run_sinex_tropoParam_excl, 'Enable', 'off')
    end
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox_estimation_leastSquares_tropo_egrRelConstr.
function checkbox_estimation_leastSquares_tropo_egrRelConstr_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_tropo_egrRelConstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_tropo_egrRelConstr

if get(hObject, 'Value')
    newState='on';
else
    newState='off';
end

set(handles.text_estimation_leastSquares_tropo_egrRelConstr, 'Enable', newState)
set(handles.edit_estimation_leastSquares_tropo_egrRelConstr, 'Enable', newState)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


function edit_estimation_leastSquares_tropo_egrRelConstr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_tropo_egrRelConstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_tropo_egrRelConstr as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_tropo_egrRelConstr as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_tropo_egrRelConstr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_tropo_egrRelConstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_estimation_leastSquares_tropo_egrInterval_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_tropo_egrInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_tropo_egrInterval as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_tropo_egrInterval as a double

% set two static texts
set(handles.text313, 'String', sprintf('after %s minutes', get(hObject, 'String')))

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_tropo_egrInterval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_tropo_egrInterval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_estimation_leastSquares_tropo_egrAbsConstr.
function checkbox_estimation_leastSquares_tropo_egrAbsConstr_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_tropo_egrAbsConstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_tropo_egrAbsConstr

if get(hObject, 'Value')
    newState='on';
else
    newState='off';
end
set(handles.text_estimation_leastSquares_tropo_egrAbsConstr, 'Enable', newState)
set(handles.edit_estimation_leastSquares_tropo_egrAbsConstr, 'Enable', newState)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


function edit_estimation_leastSquares_tropo_egrAbsConstr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_tropo_egrAbsConstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_tropo_egrAbsConstr as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_tropo_egrAbsConstr as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_tropo_egrAbsConstr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_tropo_egrAbsConstr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_estimation_leastSquares_coordinates_estimate.
function checkbox_estimation_leastSquares_coordinates_estimate_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_coordinates_estimate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_coordinates_estimate

% define wether sinex options (write station coords) should be enabled
newSinexState='off';

if get(hObject, 'Value')
    set(handles.checkbox_estimation_leastSquares_coordinates_NNT, 'Enable', 'on');
    set(handles.checkbox_estimation_leastSquares_coordinates_NNR, 'Enable', 'on');
    set(handles.checkbox_estimation_leastSquares_coordinates_NNS, 'Enable', 'on');
    % disable option for sinex output
    if get(handles.checkbox_run_sinex_write, 'Value')==1
        newSinexState='on';
    end
    
else
    set(handles.checkbox_estimation_leastSquares_coordinates_NNT, 'Enable', 'off');
    set(handles.checkbox_estimation_leastSquares_coordinates_NNR, 'Enable', 'off');
    set(handles.checkbox_estimation_leastSquares_coordinates_NNS, 'Enable', 'off');
end

set(handles.radiobutton_run_sinex_stationCoords_incl, 'Enable', newSinexState)
set(handles.text_run_sinex_stationCoords, 'Enable', newSinexState)
    
% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox_estimation_leastSquares_coordinates_NNT.
function checkbox_estimation_leastSquares_coordinates_NNT_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_coordinates_NNT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_coordinates_NNT

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox_estimation_leastSquares_coordinates_NNR.
function checkbox_estimation_leastSquares_coordinates_NNR_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_coordinates_NNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_coordinates_NNR

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox_estimation_leastSquares_coordinates_NNS.
function checkbox_estimation_leastSquares_coordinates_NNS_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_coordinates_NNS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_coordinates_NNS

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox_estimation_leastSquares_eop_xpEst.
function checkbox_estimation_leastSquares_eop_xpEst_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_eop_xpEst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_eop_xpEst

if get(hObject, 'Value')
    set(handles.edit_estimation_leastSquares_eop_interval_xp, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_eop_constr_xp, 'Enable', 'on')
    if get(handles.checkbox_estimation_leastSquares_eop_constr_xp, 'Value')
        set(handles.edit_estimation_leastSquares_constr_xp, 'Enable', 'on')
    end
else
    set(handles.edit_estimation_leastSquares_eop_interval_xp, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_eop_constr_xp, 'Enable', 'off')    
    set(handles.edit_estimation_leastSquares_constr_xp, 'Enable', 'off')
    
end

% if neither eop is estimtated -> disable options in sinex output
if (get(handles.checkbox_estimation_leastSquares_eop_xpEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_ypEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_dut1Est, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_nutdxEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_nutdyEst, 'Value'))==0 || ...
        get(handles.checkbox_run_sinex_write, 'Value')==0      
    newState='off';
    
else
    newState='on';
end

set(handles.radiobutton_run_sinex_eop_incl, 'Enable', newState);
set(handles.radiobutton_run_sinex_eop_excl, 'Enable', newState);
set(handles.text_run_sinex_eop, 'Enable', newState)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


function edit_estimation_leastSquares_eop_interval_xp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_eop_interval_xp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_eop_interval_xp as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_eop_interval_xp as a double

% set static text
set(handles.text316, 'String', sprintf('after %s minutes', get(hObject, 'String')))
% save parameter file automatically 
auto_save_parameterfile(hObject, handles)

% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_eop_interval_xp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_eop_interval_xp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_estimation_leastSquares_eop_constr_xp.
function checkbox_estimation_leastSquares_eop_constr_xp_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_eop_constr_xp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_eop_constr_xp

if get(hObject, 'Value')
    set(handles.edit_estimation_leastSquares_constr_xp, 'Enable', 'on')
else
    set(handles.edit_estimation_leastSquares_constr_xp, 'Enable', 'off')
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


function edit_estimation_leastSquares_constr_xp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_constr_xp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_constr_xp as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_constr_xp as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_constr_xp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_constr_xp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_estimation_leastSquares_eop_ypEst.
function checkbox_estimation_leastSquares_eop_ypEst_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_eop_ypEst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_eop_ypEst

if get(hObject, 'Value')
    set(handles.edit_estimation_leastSquares_eop_interval_yp, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_eop_constr_yp, 'Enable', 'on')
    if get(handles.checkbox_estimation_leastSquares_eop_constr_yp, 'Value')
        set(handles.edit_estimation_leastSquares_constr_yp, 'Enable', 'on')
    end
else
    set(handles.edit_estimation_leastSquares_eop_interval_yp, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_eop_constr_yp, 'Enable', 'off')    
    set(handles.edit_estimation_leastSquares_constr_yp, 'Enable', 'off')
end

% if neither eop is estimtated -> disable options in sinex output
if (get(handles.checkbox_estimation_leastSquares_eop_xpEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_ypEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_dut1Est, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_nutdxEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_nutdyEst, 'Value'))==0 || ...
        get(handles.checkbox_run_sinex_write, 'Value')==0      
    newState='off';
    
else
    newState='on';
end

set(handles.radiobutton_run_sinex_eop_incl, 'Enable', newState);
set(handles.radiobutton_run_sinex_eop_excl, 'Enable', newState);
set(handles.text_run_sinex_eop, 'Enable', newState)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)



function edit_estimation_leastSquares_eop_interval_yp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_eop_interval_yp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_eop_interval_yp as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_eop_interval_yp as a double

set(handles.text317, 'String', sprintf('after %s minutes', get(hObject, 'String')))

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_eop_interval_yp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_eop_interval_yp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_estimation_leastSquares_eop_constr_yp.
function checkbox_estimation_leastSquares_eop_constr_yp_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_eop_constr_yp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_eop_constr_yp

if get(hObject, 'Value')
    set(handles.edit_estimation_leastSquares_constr_yp, 'Enable', 'on')
else
    set(handles.edit_estimation_leastSquares_constr_yp, 'Enable', 'off')
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


function edit_estimation_leastSquares_constr_yp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_constr_yp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_constr_yp as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_constr_yp as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_constr_yp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_constr_yp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_estimation_leastSquares_eop_dut1Est.
function checkbox_estimation_leastSquares_eop_dut1Est_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_eop_dut1Est (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_eop_dut1Est

if get(hObject, 'Value')
    set(handles.edit_estimation_leastSquares_eop_interval_dut1, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_eop_constr_dut1, 'Enable', 'on')
    if get(handles.checkbox_estimation_leastSquares_eop_constr_dut1, 'Value')
        set(handles.edit_estimation_leastSquares_constr_dut1, 'Enable', 'on')
    end
else
    set(handles.edit_estimation_leastSquares_eop_interval_dut1, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_eop_constr_dut1, 'Enable', 'off')    
    set(handles.edit_estimation_leastSquares_constr_dut1, 'Enable', 'off')
end

% if neither eop is estimtated -> disable options in sinex output
if (get(handles.checkbox_estimation_leastSquares_eop_xpEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_ypEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_dut1Est, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_nutdxEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_nutdyEst, 'Value'))==0 || ...
        get(handles.checkbox_run_sinex_write, 'Value')==0      
    newState='off';
    
else
    newState='on';
end

set(handles.radiobutton_run_sinex_eop_incl, 'Enable', newState);
set(handles.radiobutton_run_sinex_eop_excl, 'Enable', newState);
set(handles.text_run_sinex_eop, 'Enable', newState)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)



function edit_estimation_leastSquares_eop_interval_dut1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_eop_interval_dut1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_eop_interval_dut1 as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_eop_interval_dut1 as a double

set(handles.text318, 'String', sprintf('after %s minutes', get(hObject, 'String')))

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_eop_interval_dut1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_eop_interval_dut1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_estimation_leastSquares_eop_constr_dut1.
function checkbox_estimation_leastSquares_eop_constr_dut1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_eop_constr_dut1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_eop_constr_dut1

if get(hObject, 'Value')
    set(handles.edit_estimation_leastSquares_constr_dut1, 'Enable', 'on')
else
    set(handles.edit_estimation_leastSquares_constr_dut1, 'Enable', 'off')
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


function edit_estimation_leastSquares_constr_dut1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_constr_dut1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_constr_dut1 as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_constr_dut1 as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_constr_dut1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_constr_dut1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_estimation_leastSquares_eop_nutdxEst.
function checkbox_estimation_leastSquares_eop_nutdxEst_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_eop_nutdxEst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_eop_nutdxEst

if get(hObject, 'Value')
    set(handles.edit_estimation_leastSquares_eop_interval_nutdx, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_eop_constr_nutdx, 'Enable', 'on')
    if get(handles.checkbox_estimation_leastSquares_eop_constr_nutdx, 'Value')
        set(handles.edit_estimation_leastSquares_constr_nutdx, 'Enable', 'on')
    end
else
    set(handles.edit_estimation_leastSquares_eop_interval_nutdx, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_eop_constr_nutdx, 'Enable', 'off')    
    set(handles.edit_estimation_leastSquares_constr_nutdx, 'Enable', 'off')
end

% if neither eop is estimtated -> disable options in sinex output
if (get(handles.checkbox_estimation_leastSquares_eop_xpEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_ypEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_dut1Est, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_nutdxEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_nutdyEst, 'Value'))==0 || ...
        get(handles.checkbox_run_sinex_write, 'Value')==0      
    newState='off';
    
else
    newState='on';
end

set(handles.radiobutton_run_sinex_eop_incl, 'Enable', newState);
set(handles.radiobutton_run_sinex_eop_excl, 'Enable', newState);
set(handles.text_run_sinex_eop, 'Enable', newState)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)



function edit_estimation_leastSquares_eop_interval_nutdx_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_eop_interval_nutdx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_eop_interval_nutdx as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_eop_interval_nutdx as a double

set(handles.text319, 'String', sprintf('after %s minutes', get(hObject, 'String')))

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_eop_interval_nutdx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_eop_interval_nutdx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_estimation_leastSquares_eop_constr_nutdx.
function checkbox_estimation_leastSquares_eop_constr_nutdx_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_eop_constr_nutdx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_eop_constr_nutdx

if get(hObject, 'Value')
    set(handles.edit_estimation_leastSquares_constr_nutdx, 'Enable', 'on')
else
    set(handles.edit_estimation_leastSquares_constr_nutdx, 'Enable', 'off')
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


function edit_estimation_leastSquares_constr_nutdx_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_constr_nutdx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_constr_nutdx as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_constr_nutdx as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_constr_nutdx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_constr_nutdx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_estimation_leastSquares_eop_nutdyEst.
function checkbox_estimation_leastSquares_eop_nutdyEst_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_eop_nutdyEst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_eop_nutdyEst

if get(hObject, 'Value')
    set(handles.edit_estimation_leastSquares_eop_interval_nutdy, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_eop_constr_nutdy, 'Enable', 'on')
    if get(handles.checkbox_estimation_leastSquares_eop_constr_nutdy, 'Value')
        set(handles.edit_estimation_leastSquares_constr_nutdy, 'Enable', 'on')
    end
else
    set(handles.edit_estimation_leastSquares_eop_interval_nutdy, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_eop_constr_nutdy, 'Enable', 'off')    
    set(handles.edit_estimation_leastSquares_constr_nutdy, 'Enable', 'off')
end

% if neither eop is estimtated -> disable options in sinex output
if (get(handles.checkbox_estimation_leastSquares_eop_xpEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_ypEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_dut1Est, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_nutdxEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_nutdyEst, 'Value'))==0 || ...
        get(handles.checkbox_run_sinex_write, 'Value')==0      
    newState='off';
    
else
    newState='on';
end

set(handles.radiobutton_run_sinex_eop_incl, 'Enable', newState);
set(handles.radiobutton_run_sinex_eop_excl, 'Enable', newState);
set(handles.text_run_sinex_eop, 'Enable', newState)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)



function edit_estimation_leastSquares_eop_interval_nutdy_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_eop_interval_nutdy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_eop_interval_nutdy as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_eop_interval_nutdy as a double

set(handles.text320, 'String', sprintf('after %s minutes', get(hObject, 'String')))

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_eop_interval_nutdy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_eop_interval_nutdy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_estimation_leastSquares_eop_constr_nutdy.
function checkbox_estimation_leastSquares_eop_constr_nutdy_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_eop_constr_nutdy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_eop_constr_nutdy

if get(hObject, 'Value')
    set(handles.edit_estimation_leastSquares_constr_nutdy, 'Enable', 'on')
else
    set(handles.edit_estimation_leastSquares_constr_nutdy, 'Enable', 'off')
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


function edit_estimation_leastSquares_constr_nutdy_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_constr_nutdy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_constr_nutdy as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_constr_nutdy as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_constr_nutdy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_constr_nutdy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_estimation_leastSquares_sources_est.
function checkbox_estimation_leastSquares_sources_est_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_sources_est (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_sources_est

% set enabling
if get(hObject, 'Value')==0 % disabled
    set(handles.edit_estimation_leastSquares_sources_interval, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_sources_constr, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_sources_constr, 'Enable', 'off')
    % disable option for sinex output
    set(handles.radiobutton_run_sinex_sources_incl, 'Enable', 'off')
    set(handles.text_run_sinex_sources, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_sources_ICRF2_def, 'Enable', 'off')
	set(handles.edit_estimation_leastSquares_sources_abs_constr, 'Enable', 'off')
else % enabled
    set(handles.edit_estimation_leastSquares_sources_interval, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_sources_constr, 'Enable', 'on')
    if get(handles.checkbox_estimation_leastSquares_sources_constr, 'Value')
        set(handles.edit_estimation_leastSquares_sources_constr, 'Enable', 'on')
    else
        set(handles.edit_estimation_leastSquares_sources_constr, 'Enable', 'off')
    end
    % disable option for sinex output
    if get(handles.checkbox_run_sinex_write, 'Value')
        if get(handles.checkbox_run_sinex_sources, 'Value')
            set(handles.radiobutton_run_sinex_sources_incl, 'Enable', 'on')
            set(handles.text_run_sinex_sources, 'Enable', 'on')
        end
    end
    
    % set(handles.checkbox_estimation_leastSquares_sources_NNR, 'Enable', 'off')
	set(handles.checkbox_estimation_leastSquares_sources_NNR, 'Value', 0)
    set(handles.checkbox_estimation_leastSquares_sources_ICRF2_def, 'Enable', 'off')
	set(handles.edit_estimation_leastSquares_sources_abs_constr, 'Enable', 'off')
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


function edit_estimation_leastSquares_sources_interval_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_sources_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_sources_interval as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_sources_interval as a double

% set static text
set(handles.text321, 'String', sprintf('after %s minutes', get(hObject, 'String')))

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_sources_interval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_sources_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_estimation_leastSquares_sources_constr.
function checkbox_estimation_leastSquares_sources_constr_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_sources_constr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_sources_constr

if get(hObject, 'Value')
    set(handles.edit_estimation_leastSquares_sources_constr, 'Enable', 'on')
else
    set(handles.edit_estimation_leastSquares_sources_constr, 'Enable', 'off')
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


function edit_estimation_leastSquares_sources_constr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_sources_constr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_sources_constr as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_sources_constr as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_sources_constr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_sources_constr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on button press in checkbox_estimation_leastSquares_sources_NNR.
function checkbox_estimation_leastSquares_sources_NNR_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_sources_NNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_sources_NNR

% set enabling
if get(hObject, 'Value')==0 % Checkbox disabled
    % set(handles.checkbox_estimation_leastSquares_sources_est, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_sources_ICRF2_def, 'Enable', 'off')
	set(handles.edit_estimation_leastSquares_sources_abs_constr, 'Enable', 'off')
	set(handles.checkbox_estimation_leastSquares_sources_abs_constr, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_sources_obs_per_source, 'Enable', 'off')    
    set(handles.edit_estimation_leastSquares_sources_obs_per_source, 'Enable', 'off')

else % == 1; Checkbox enabled
    set(handles.checkbox_estimation_leastSquares_sources_ICRF2_def, 'Enable', 'on')
    
	set(handles.checkbox_estimation_leastSquares_sources_abs_constr, 'Enable', 'on')
    if get(handles.checkbox_estimation_leastSquares_sources_abs_constr,'Value')
        set(handles.edit_estimation_leastSquares_sources_abs_constr, 'Enable', 'on')
    else
        set(handles.edit_estimation_leastSquares_sources_abs_constr, 'Enable', 'off')
    end

    set(handles.checkbox_estimation_leastSquares_sources_obs_per_source, 'Enable', 'on')    
    if get(handles.checkbox_estimation_leastSquares_sources_obs_per_source, 'Value')
        set(handles.edit_estimation_leastSquares_sources_obs_per_source, 'Enable', 'on')
    else
        set(handles.edit_estimation_leastSquares_sources_obs_per_source, 'Enable', 'off')
    end
    
    set(handles.checkbox_estimation_leastSquares_sources_est, 'Value',0)
	set(handles.edit_estimation_leastSquares_sources_interval, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_sources_constr, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_sources_constr, 'Enable', 'off')
    % disable option for sinex output
    set(handles.radiobutton_run_sinex_sources_incl, 'Enable', 'off')
    set(handles.text_run_sinex_sources, 'Enable', 'off')
	
end
% save parameter file automatically 
auto_save_parameterfile(hObject, handles)

% --- Executes on button press in checkbox_estimation_leastSquares_sources_NNR.
function checkbox_estimation_leastSquares_sources_abs_constr_Callb(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_sources_NNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_sources_NNR

% set enabling
if get(hObject, 'Value')==0 % Checkbox disabled
	set(handles.edit_estimation_leastSquares_sources_abs_constr, 'Enable', 'off')
else % == 1; Checkbox enabled
	set(handles.edit_estimation_leastSquares_sources_abs_constr, 'Enable', 'on')

end
% save parameter file automatically 
auto_save_parameterfile(hObject, handles)

% --- Executes on button press in checkbox_estimation_leastSquares_sources_NNR.
function checkbox_estimation_leastSquares_sources_obs_per_source_Callb(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_sources_NNR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_sources_NNR

% set enabling
if get(hObject, 'Value')==0 % Checkbox disabled
	set(handles.edit_estimation_leastSquares_sources_obs_per_source, 'Enable', 'off')
else % == 1; Checkbox enabled
	set(handles.edit_estimation_leastSquares_sources_obs_per_source, 'Enable', 'on')

end
% save parameter file automatically 
auto_save_parameterfile(hObject, handles)

    






% --- Executes on button press in checkbox_run_allowStationwise.
function checkbox_run_allowStationwise_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_allowStationwise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_allowStationwise

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox_run_prepareGlobParam.
function checkbox_run_prepareGlobParam_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_prepareGlobParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_prepareGlobParam

if get(hObject, 'Value')
    newState='on';
    if get(handles.checkbox_run_globalPram_statVel, 'Value')
        set(handles.edit_run_globalPram_statVelEpoch, 'Enable', 'on')
        set(handles.text_run_runOptions_refEpochInYrs, 'Enable', 'on')
    end
else
    newState='off';
    set(handles.edit_run_globalPram_statVelEpoch, 'Enable', 'off')
    set(handles.text_run_runOptions_refEpochInYrs, 'Enable', 'off')
end

set(handles.checkbox_run_globalPram_source, 'Enable', newState)
set(handles.text_run_runOptions_attentionSources, 'Enable', newState)
set(handles.checkbox_run_globalPram_statVel, 'Enable', newState)
set(handles.checkbox_run_globalPram_axisOffset, 'Enable', newState)
set(handles.checkbox_run_globalPram_stseaspos, 'Enable', newState)
set(handles.checkbox_run_globalPram_tidERPvar, 'Enable', newState)
set(handles.checkbox_run_globalPram_hlpole, 'Enable', newState)
set(handles.checkbox_run_globalPram_APLrg, 'Enable', newState)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox_run_globalPram_source.
function checkbox_run_globalPram_source_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_globalPram_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_globalPram_source

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox_run_globalPram_statVel.
function checkbox_run_globalPram_statVel_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_globalPram_statVel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_globalPram_statVel
if get(hObject, 'Value')
    newState='on';
else
    newState='off';
end

set(handles.edit_run_globalPram_statVelEpoch, 'Enable', newState)
set(handles.text_run_runOptions_refEpochInYrs, 'Enable', newState)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


function edit_run_globalPram_statVelEpoch_Callback(hObject, eventdata, handles)
% hObject    handle to edit_run_globalPram_statVelEpoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_run_globalPram_statVelEpoch as text
%        str2double(get(hObject,'String')) returns contents of edit_run_globalPram_statVelEpoch as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_run_globalPram_statVelEpoch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_run_globalPram_statVelEpoch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_estimation_leastSquares_clocks_useClockBreaks.
function checkbox_estimation_leastSquares_clocks_useClockBreaks_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_clocks_useClockBreaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_clocks_useClockBreaks

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox173.
function checkbox173_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox173 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox173


% --- Executes on button press in checkbox174.
function checkbox174_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox174 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox174



function edit_run_outDirs_oneSub_Callback(hObject, eventdata, handles)
% hObject    handle to edit_run_outDirs_oneSub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_run_outDirs_oneSub as text
%        str2double(get(hObject,'String')) returns contents of edit_run_outDirs_oneSub as a double

% add the subfolder also to the eop_out panel stuff
set(handles.radiobutton_plot_eopOut_subfolder_current, 'String',...
    ['Current subfolder (''', get(hObject, 'String'), ''')']);

set(handles.radiobutton_plot_eopOut_outFile_default, 'String',...
    ['Default: ''/VieVS/OUT/eop_', get(hObject, 'String'), '.txt or ''bas_', ...
    get(hObject, 'String'),'.txt''']);
	
% also to the baseline length (basel) output
set(handles.edit_plot_eopOut_basRepOptions_writeBasOutFname, 'String',...
	['basel_', get(hObject, 'String'), '.txt']);

% If Vie_LSM and Vie_GLOB are run both
if get(handles.checkbox_run_outputDirectories_runVieGlob,'Value') && ...
(get(handles.checkbox_run_outputDirectories_runVieLsmScanwiseUpdate,'Value') || get(handles.checkbox_run_outputDirectories_runVieLsm,'Value'))
	% Use different output Sub-directories
	if get(handles.checkbox_run_outDirs_diffSubs, 'Value')
		str_lsm_subDir=get(handles.edit_run_outDirs_level2, 'String');
	% Use one output Sub-directory	
	else
		str_lsm_subDir=get(handles.edit_run_outDirs_oneSub, 'String');
	end
	set(handles.edit_run_outDirs_glob_level2Sub, 'String', str_lsm_subDir)
end
	
% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_run_outDirs_oneSub_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_run_outDirs_oneSub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_run_outDirs_level0_Callback(hObject, eventdata, handles)
% hObject    handle to edit_run_outDirs_level0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_run_outDirs_level0 as text
%        str2double(get(hObject,'String')) returns contents of edit_run_outDirs_level0 as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_run_outDirs_level0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_run_outDirs_level0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_run_outDirs_level1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_run_outDirs_level1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_run_outDirs_level1 as text
%        str2double(get(hObject,'String')) returns contents of edit_run_outDirs_level1 as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_run_outDirs_level1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_run_outDirs_level1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_run_outDirs_level3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_run_outDirs_level3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_run_outDirs_level3 as text
%        str2double(get(hObject,'String')) returns contents of edit_run_outDirs_level3 as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_run_outDirs_level3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_run_outDirs_level3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_run_outDirs_level2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_run_outDirs_level2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_run_outDirs_level2 as text
%        str2double(get(hObject,'String')) returns contents of edit_run_outDirs_level2 as a double

% If Vie_LSM and Vie_GLOB are run both
if get(handles.checkbox_run_outputDirectories_runVieGlob,'Value') && ...
(get(handles.checkbox_run_outputDirectories_runVieLsmScanwiseUpdate,'Value') || get(handles.checkbox_run_outputDirectories_runVieLsm,'Value'))
	% Use different output Sub-directories
	if get(handles.checkbox_run_outDirs_diffSubs, 'Value')
		str_lsm_subDir=get(handles.edit_run_outDirs_level2, 'String');
	% Use one output Sub-directory	
	else
		str_lsm_subDir=get(handles.edit_run_outDirs_oneSub, 'String');
	end
	set(handles.edit_run_outDirs_glob_level2Sub, 'String', str_lsm_subDir)
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_run_outDirs_level2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_run_outDirs_level2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function newestAutoSavedParameterFile=getNewestAutoSavedParameterFile(hObject, handles)
% This function returns the newest automatically saved parameter file in
% PARAMETERS/, where all auto saved files should be stored to. This
% function is used when the GUI-parameters should be saved to a textfile.

allAutoSavedParameterFiles=dir('PARAMETERS/auto_save_fromGUI_*.mat');
[sortedNames, sortInd]=sort({allAutoSavedParameterFiles.name});
allAutoSavedParameterFiles=allAutoSavedParameterFiles(sortInd);
newestAutoSavedParameterFile=allAutoSavedParameterFiles(end);

% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% save gui parameters (again)
auto_save_parameterfile(hObject, handles)
 
% prepare everything for calling vie_batch
save_runp(hObject, handles)

% call vie_batch
vie_batch




function edit106_Callback(hObject, eventdata, handles)
% hObject    handle to edit106 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit106 as text
%        str2double(get(hObject,'String')) returns contents of edit106 as a double


% --- Executes during object creation, after setting all properties.
function edit106_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit106 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function popupmenu_parameters_refFrames_otherCRF_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_refFrames_otherCRF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of popupmenu_parameters_refFrames_otherCRF as text
%        str2double(get(hObject,'String')) returns contents of popupmenu_parameters_refFrames_otherCRF as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu_parameters_refFrames_otherCRF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_refFrames_otherCRF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu_parameters_refFrames_otherTRF_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_refFrames_otherTRF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of popupmenu_parameters_refFrames_otherTRF as text
%        str2double(get(hObject,'String')) returns contents of popupmenu_parameters_refFrames_otherTRF as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu_parameters_refFrames_otherTRF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_refFrames_otherTRF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popupmenu_parameters_iono_ext_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_iono_ext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of popupmenu_parameters_iono_ext as text
%        str2double(get(hObject,'String')) returns contents of popupmenu_parameters_iono_ext as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu_parameters_iono_ext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_iono_ext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function popupmenu_parameters_eop_aPriori_other_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_eop_aPriori_other (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of popupmenu_parameters_eop_aPriori_other as text
%        str2double(get(hObject,'String')) returns contents of popupmenu_parameters_eop_aPriori_other as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu_parameters_eop_aPriori_other_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_eop_aPriori_other (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_setInput_browseForProcLists.
function pushbutton_setInput_browseForProcLists_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_setInput_browseForProcLists (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName, PathName] = uigetfile('*.mat','Select process lists', './PROCESSLIST/', 'multiselect', 'on');

fileWasChosen=1;

if isempty(FileName)
    fileWasChosen=0;
elseif ~iscell(FileName)
    if FileName == 0
        fileWasChosen=0;
    end
end
       
if fileWasChosen
    % get number of files
    if iscell(FileName)
        nFiles=size(FileName,2);
    else
        nFiles=1;
    end
    
    % for all files
    for iFile=1:nFiles
        % load process list
        
        % if we have a cell == if we have more than 1 process_list selected
        if iscell(FileName)
            load([PathName, FileName{iFile}])
        else
            load([PathName, FileName])
        end
        
        % if we have now the process_list variable
        if exist('process_list', 'var')
            % get current listbox entries
            curContent=get(handles.listbox_setInput_processList, 'String');
            
            % delete last entry of process_list when it is ''
            if strcmp(process_list(size(process_list,1)), ' ')
                process_list(size(process_list,1),:)=[];
            end
            
            % make cellstr out of process_list
            process_list_cellstr=cellstr(process_list);
            
            % if listbox is empty: just take chosen process_list
            if isempty(curContent)
                newContent=process_list_cellstr;
            else % else: make unique
                newContent=unique([curContent; process_list_cellstr]);
            end

            % update listbox
            set(handles.listbox_setInput_processList, 'String', newContent);
        end
        
    end

    % save all selected sessions to handles struct
    handles.allSelectedFiles=newContent;
    
    % save changes to handles struct
    guidata(hObject, handles);
    
    % save parameter file automatically 
    auto_save_parameterfile(hObject, handles)
end
    


% --- Executes on button press in pushbutton_setInput_clearProcList.
function pushbutton_setInput_clearProcList_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_setInput_clearProcList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% get current content of listbox
curContent=get(handles.listbox_setInput_processList, 'String');

if ~isempty(curContent)
    % delete all selected entries entry
    
    curContent(get(handles.listbox_setInput_processList, 'Value'))=[];
    
    % get current value
    curValue=get(handles.listbox_setInput_processList, 'Value');
    
    % set the current value either to the current value or to the latest
    % (if the former latest was deleted). If no session is there anymore,
    % take 1.
    set(handles.listbox_setInput_processList, 'Value', ...
        min([max(curValue), max(length(curContent),1)]));

    % update listbox
    set(handles.listbox_setInput_processList, 'String', curContent);

    % and also save it to handles struct
    handles.allSelectedFiles=curContent;

    % save changes to handles struct
    guidata(hObject, handles);

    % save parameter file automatically 
    auto_save_parameterfile(hObject, handles)
end



% --- Executes on button press in checkbox_parameters_statCorr_oceanPoleTides.
function checkbox_parameters_statCorr_oceanPoleTides_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameters_statCorr_oceanPoleTides (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameters_statCorr_oceanPoleTides

if get(hObject, 'Value')
    set(handles.radiobutton_parameters_statCorr_poleModel_lin, 'Enable', 'on')
    set(handles.radiobutton_parameters_statCorr_poleModel_cub, 'Enable', 'on')
    set(handles.radiobutton_parameters_statCorr_poleModel_iers2015, 'Enable', 'on')
else
    % only set to disable when also other checkbox is 0
    if get(handles.checkbox_parameters_statCorr_poleTides, 'Value')==0
        set(handles.radiobutton_parameters_statCorr_poleModel_lin, 'Enable', 'off')
        set(handles.radiobutton_parameters_statCorr_poleModel_cub, 'Enable', 'off')
        set(handles.radiobutton_parameters_statCorr_poleModel_iers2015, 'Enable', 'off')
    end
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


function popupmenu_parameters_statCorr_tidalOceanLoad_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_statCorr_tidalOceanLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of popupmenu_parameters_statCorr_tidalOceanLoad as text
%        str2double(get(hObject,'String')) returns contents of popupmenu_parameters_statCorr_tidalOceanLoad as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu_parameters_statCorr_tidalOceanLoad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_statCorr_tidalOceanLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function popupmenu_parameters_statCorr_tidalAtmoOceanLoad_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_statCorr_tidalAtmoOceanLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of popupmenu_parameters_statCorr_tidalAtmoOceanLoad as text
%        str2double(get(hObject,'String')) returns contents of popupmenu_parameters_statCorr_tidalAtmoOceanLoad as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu_parameters_statCorr_tidalAtmoOceanLoad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_statCorr_tidalAtmoOceanLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function popupmenu_parameters_statCorr_nonTidalAtmoOceanLoad_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_statCorr_nonTidalAtmoOceanLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of popupmenu_parameters_statCorr_nonTidalAtmoOceanLoad as text
%        str2double(get(hObject,'String')) returns contents of popupmenu_parameters_statCorr_nonTidalAtmoOceanLoad as a double

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu_parameters_statCorr_nonTidalAtmoOceanLoad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_statCorr_nonTidalAtmoOceanLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes when selected object is changed in uipanel88.
function uipanel88_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel88 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(hObject, 'Tag')
    case 'radiobutton_parameters_refFrames_otherTRF'
        if get(hObject, 'Value')
            set(handles.popupmenu_parameters_refFrames_otherTRF, 'Enable', 'on')
            set(handles.pushbutton_parameters_refFrames_superstationTRF_chose, 'Enable', 'off')
            set(handles.popupmenu_parameters_refFrames_superstationTRF, 'Enable', 'off')
        end
    case 'radiobutton_parameters_refFrames_superstationTRF'
        if get(hObject, 'Value')
            set(handles.popupmenu_parameters_refFrames_otherTRF, 'Enable', 'off')
            set(handles.pushbutton_parameters_refFrames_superstationTRF_chose, 'Enable', 'on')
            set(handles.popupmenu_parameters_refFrames_superstationTRF, 'Enable', 'on')
        end
    otherwise
        set(handles.popupmenu_parameters_refFrames_otherTRF, 'Enable', 'off')
        set(handles.pushbutton_parameters_refFrames_superstationTRF_chose, 'Enable', 'off')
        set(handles.popupmenu_parameters_refFrames_superstationTRF, 'Enable', 'off')
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)

% --- Executes when selected object is changed in uipanel89.
function uipanel89_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel89 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(hObject, 'Tag')
    case 'radiobutton_parameters_refFrames_otherCRF'
        if get(hObject, 'Value')
            set(handles.popupmenu_parameters_refFrames_otherCRF, 'Enable', 'on')
            set(handles.popupmenu_parameters_refFrames_supersourceCRF, 'Enable', 'Off')
            set(handles.pushbutton_parameters_refFrames_superstationCRF_chose, 'Enable', 'Off')
        end
    otherwise
        set(handles.popupmenu_parameters_refFrames_otherCRF, 'Enable', 'off')
        set(handles.pushbutton_parameters_refFrames_superstationCRF_chose, 'Enable', 'On')
        set(handles.popupmenu_parameters_refFrames_supersourceCRF, 'Enable', 'On')
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes when selected object is changed in uipanel99.
function uipanel99_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel99 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% defines if the popupmenu (containing all folders of TRP) should be
% updated -> is needed when external tropospheric files are created in the
% meantime
updatePopupmenu=0;

switch get(hObject, 'Tag')
    case 'radiobutton_parameters_iono_ext'
        set(handles.popupmenu_parameters_iono_ext, 'Enable', 'on')
        set(handles.pushbutton_parameters_iono_create, 'Enable', 'on')
        updatePopupmenu=1;
    otherwise
        set(handles.popupmenu_parameters_iono_ext, 'Enable', 'off')
        set(handles.pushbutton_parameters_iono_create, 'Enable', 'off')
end

% update popupmenu
if updatePopupmenu==1
    curContent=get(handles.popupmenu_parameters_iono_ext, 'String');
    if ~iscell(curContent)
        curContent = {curContent};
    end
    curSelected=curContent{get(handles.popupmenu_parameters_iono_ext, 'Value')};
    
    % get directory content
    dirsInIonFolder=dir('../ION/FILES/');
    % delete '.', '..', and all files
    dirsInIonFolder(strcmp({dirsInIonFolder.name}, '.')|strcmp({dirsInIonFolder.name}, '..')|~[dirsInIonFolder.isdir])=[];
    
    % write folders to popupmenu
    if isempty(dirsInIonFolder)
        set(handles.popupmenu_parameters_iono_ext, 'String', ' ')
    else
        set(handles.popupmenu_parameters_iono_ext, 'String', {dirsInIonFolder.name})
    end
    
    % try to find previous selected to select this again
    valueOfPrevSelectedFolder=find(strcmp({dirsInIonFolder.name}, curSelected));
    if isempty(valueOfPrevSelectedFolder)
        set(handles.popupmenu_parameters_iono_ext, 'Value', 1)
    else
        set(handles.popupmenu_parameters_iono_ext, 'Value', valueOfPrevSelectedFolder)
    end
    
end %update popupmenu
% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes when selected object is changed in panel_models_troposphere_zhd.
function panel_models_troposphere_zhd_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in panel_models_troposphere_zhd 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)

% % defines if the popupmenu (containing all folders of TRP) should be  
% % updated -> is needed when external tropospheric files are created in the
% % meantime
% updatePopupmenu=0;
% 
% switch get(hObject, 'Tag')
%     case 'radiobutton_parameters_troposphere_externalFile'
%         set(handles.popupmenu_parameters_tropo_externalFile, 'Enable', 'on')
%         set(handles.pushbutton_parameters_troposphere_externalsCreate, 'Enable', 'on')
%         % enable also (again?) mapping functions and gradients
%         set(handles.radiobutton_parameters_troposphere_mfh_VMF3, 'Enable', 'Off')
%         set(handles.radiobutton_parameters_troposphere_mfh_VMF1, 'Enable', 'Off')
%         set(handles.radiobutton_parameters_troposphere_mfh_GPT3, 'Enable', 'Off')
%         set(handles.radiobutton_parameters_troposphere_mfw_VMF3, 'Enable', 'Off')
%         set(handles.radiobutton_parameters_troposphere_mfw_VMF1, 'Enable', 'Off')
%         set(handles.radiobutton_parameters_troposphere_mfw_GPT3, 'Enable', 'Off')
%         set(handles.radiobutton_parameters_troposphere_gradients_h_no, 'Enable', 'Off')
%         set(handles.radiobutton_parameters_troposphere_gradients_h_GRAD, 'Enable', 'Off')
%         set(handles.radiobutton_parameters_troposphere_gradients_h_GPT3, 'Enable', 'Off')
%         set(handles.radiobutton_parameters_troposphere_gradients_h_DAO, 'Enable', 'Off')
%         set(handles.radiobutton_parameters_troposphere_gradients_w_no, 'Enable', 'Off')
%         set(handles.radiobutton_parameters_troposphere_gradients_w_GRAD, 'Enable', 'Off')
%         set(handles.radiobutton_parameters_troposphere_gradients_w_GPT3, 'Enable', 'Off')
%         updatePopupmenu=1;
%     otherwise
%         set(handles.popupmenu_parameters_tropo_externalFile, 'Enable', 'off')
%         set(handles.pushbutton_parameters_troposphere_externalsCreate, 'Enable', 'off')
%         % enable also (again?) mapping functions and gradients
%         set(handles.radiobutton_parameters_troposphere_mfh_VMF3, 'Enable', 'On')
%         set(handles.radiobutton_parameters_troposphere_mfh_VMF1, 'Enable', 'On')
%         set(handles.radiobutton_parameters_troposphere_mfh_GPT3, 'Enable', 'On')
%         set(handles.radiobutton_parameters_troposphere_mfw_VMF3, 'Enable', 'On')
%         set(handles.radiobutton_parameters_troposphere_mfw_VMF1, 'Enable', 'On')
%         set(handles.radiobutton_parameters_troposphere_mfw_GPT3, 'Enable', 'On')
%         set(handles.radiobutton_parameters_troposphere_gradients_h_no, 'Enable', 'On')
%         set(handles.radiobutton_parameters_troposphere_gradients_h_GRAD, 'Enable', 'On')
%         set(handles.radiobutton_parameters_troposphere_gradients_h_GPT3, 'Enable', 'On')
%         set(handles.radiobutton_parameters_troposphere_gradients_h_DAO, 'Enable', 'On')
%         set(handles.radiobutton_parameters_troposphere_gradients_w_no, 'Enable', 'On')
%         set(handles.radiobutton_parameters_troposphere_gradients_w_GRAD, 'Enable', 'On')
%         set(handles.radiobutton_parameters_troposphere_gradients_w_GPT3, 'Enable', 'On')
% end
% 
% % update popupmenu ?
% if updatePopupmenu==1
%     curContent=get(handles.popupmenu_parameters_tropo_externalFile, 'String');
%     curSelected=curContent{get(handles.popupmenu_parameters_tropo_externalFile, 'Value')};
%     
%     % get directory content
%     dirsInTrpFolder=dir('../TRP/OUTPUT_DATA/');
%     % delete '.', '..', and all files
%     dirsInTrpFolder(strcmp({dirsInTrpFolder.name}, '.')|strcmp({dirsInTrpFolder.name}, '..')|~[dirsInTrpFolder.isdir])=[];
%     
%     % write folders to popupmenu
%     set(handles.popupmenu_parameters_tropo_externalFile, 'String', {dirsInTrpFolder.name})
%     
%     % try to find previous selected to select this again
%     valueOfPrevSelectedFolder=find(strcmp({dirsInTrpFolder.name}, curSelected));
%     if isempty(valueOfPrevSelectedFolder)
%         set(handles.popupmenu_parameters_tropo_externalFile, 'Value', 1)
%     else
%         set(handles.popupmenu_parameters_tropo_externalFile, 'Value', valueOfPrevSelectedFolder)
%     end
%     
% end % update popupmenu
    

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)

% --- Executes when selected object is changed in panel_models_troposphere_zwd.
function panel_models_troposphere_zwd_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in panel_models_troposphere_zwd 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes when selected object is changed in uipanel116.
function uipanel116_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel116 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(hObject, 'Tag')
    case 'radiobutton_parameters_eop_aPriori_other'
        set(handles.popupmenu_parameters_eop_aPriori_other, 'Enable', 'on')
    otherwise
        set(handles.popupmenu_parameters_eop_aPriori_other, 'Enable', 'off')
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox_parameters_eop_inclHf_oceanTidesOther.
function checkbox_parameters_eop_inclHf_oceanTidesOther_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameters_eop_inclHf_oceanTidesOther (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameters_eop_inclHf_oceanTidesOther
if get(hObject, 'Value')
    set(handles.popupmenu_parameters_eop_oceanTideModel, 'enable', 'on')
else
    set(handles.popupmenu_parameters_eop_oceanTideModel, 'enable', 'off')
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes when selected object is changed in btnGrp_parameters_statCorr_meanPoleModel.
function btnGrp_parameters_statCorr_meanPoleModel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in btnGrp_parameters_statCorr_meanPoleModel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes when selected object is changed in btnGrp_parameters_statCorr_thermalDef_temp.
function btnGrp_parameters_statCorr_thermalDef_temp_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in btnGrp_parameters_statCorr_thermalDef_temp 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)



% --- Executes when selected object is changed in uipanel91.
function uipanel91_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel91 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes when selected object is changed in uipanel119.
function uipanel119_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel119 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes when selected object is changed in uipanel123.
function uipanel123_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel123 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes when selected object is changed in panel_models_troposphere_mfh.
function panel_models_troposphere_mfh_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in panel_models_troposphere_mfh 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)

% --- Executes when selected object is changed in panel_models_troposphere_mfw.
function panel_models_troposphere_mfw_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in panel_models_troposphere_mfw 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes when selected object is changed in panel_models_troposphere_gradients_h.
function panel_models_troposphere_gradients_h_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in panel_models_troposphere_gradients_h 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% set en/disable
set(handles.radiobutton_parameters_troposphere_gradients_w_no, 'Enable', 'On')
set(handles.radiobutton_parameters_troposphere_gradients_w_GRAD, 'Enable', 'On')
set(handles.radiobutton_parameters_troposphere_gradients_w_GPT3, 'Enable', 'On')

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)

% --- Executes when selected object is changed in panel_models_troposphere_gradients_w.
function panel_models_troposphere_gradients_w_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in panel_models_troposphere_gradients_w 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes when selected object is changed in uipanel_estimation_leastSquares_clock.
function uipanel_estimation_leastSquares_clock_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_estimation_leastSquares_clock 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes when selected object is changed in uipanel43.
function uipanel43_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel43 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes when selected object is changed in uipanel46.
function uipanel46_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel46 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes when selected object is changed in uipanel45.
function uipanel45_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel45 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes when selected object is changed in uipanel44.
function uipanel44_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel44 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes when selected object is changed in uipanel47.
function uipanel47_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel47 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes when selected object is changed in uipanel48.
function uipanel48_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel48 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes when selected object is changed in uipanel131.
function uipanel131_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel131 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes when user attempts to close figure_vievs2.
function figure_vievs2_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure_vievs2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

auto_save_parameterfile(hObject, handles)

% Hint: delete(hObject) closes the figure
delete(hObject);



function popupmenu_plot_folder2_subfolder_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder2_subfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of popupmenu_plot_folder2_subfolder as text
%        str2double(get(hObject,'String')) returns contents of popupmenu_plot_folder2_subfolder as a double


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_folder2_subfolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder2_subfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_plot_folder2_sessions.
function popupmenu_plot_folder2_sessions_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder2_sessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_folder2_sessions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_folder2_sessions

handles=update_antennaParameterPopupmenus_plotParameters(handles,2);

handles=plotToAxes(hObject, handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_plot_folder2_sessions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder2_sessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_plot_folder2_param.
function popupmenu_plot_folder2_param_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder2_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_folder2_param contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_folder2_param

handles=plotToAxes(hObject, handles);

chosenPanel=2;

% write unit to text (as info for the user)
unitPerParam={'pwclk', 'cm'; 'rclk', 'cm/day'; 'qclk', 'cm/(day^2)'; ...
    'rqclk', 'cm/(day^2) warning: rate+quadr'; 'zwd', 'cm'; 'ngr', 'cm'; 'egr', 'cm';...
    'xpol', 'mas'; 'ypol', 'mas'; 'dut1', 'ms'; 'nutdx', 'mas'; 'nutdy',...
    'mas'; 'coorx', 'cm'; 'coory', 'cm'; 'coorz', 'cm'; 'soude', 'mas'; 'soura', 'mas';...
    'sat_pos1', 'cm'; 'sat_pos2', 'cm'; 'sat_pos3', 'cm'};
% get selected parameter
allParamInPopupmenu=get(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_param']), 'String');
curSelParam=allParamInPopupmenu{get(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_param']), 'Value')};
strForText=unitPerParam{strcmp(unitPerParam(:,1), curSelParam),2};
if isempty(strForText)
    set(eval(['handles.text_plot_folder', num2str(chosenPanel), '_unit']), 'String', '')
else
    set(eval(['handles.text_plot_folder', num2str(chosenPanel), '_unit']), 'String', ['[', strForText, ']'])
end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_folder2_param_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder2_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_plot_folder2_stat.
function popupmenu_plot_folder2_stat_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder2_stat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_folder2_stat contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_folder2_stat

handles=plotToAxes(hObject, handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_folder2_stat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder2_stat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function popupmenu_plot_folder1_subfolder_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder1_subfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of popupmenu_plot_folder1_subfolder as text
%        str2double(get(hObject,'String')) returns contents of popupmenu_plot_folder1_subfolder as a double


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_folder1_subfolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder1_subfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_plot_folder1_sessions.
function popupmenu_plot_folder1_sessions_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder1_sessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_folder1_sessions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_folder1_sessions

handles=update_antennaParameterPopupmenus_plotParameters(handles,1);

handles=plotToAxes(hObject, handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_plot_folder1_sessions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder1_sessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_plot_folder1_param.
function popupmenu_plot_folder1_param_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder1_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_folder1_param contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_folder1_param

handles=plotToAxes(hObject, handles);

chosenPanel=1;

% write unit to text (as info for the user)
unitPerParam={'pwclk', 'cm'; 'rclk', 'cm/day'; 'qclk', 'cm/(day^2)'; ...
    'rqclk', 'cm/(day^2) warning: rate+quadr'; 'zwd', 'cm'; 'ngr', 'cm'; 'egr', 'cm';...
    'xpol', 'mas'; 'ypol', 'mas'; 'dut1', 'ms'; 'nutdx', 'mas'; 'nutdy',...
    'mas'; 'coorx', 'cm'; 'coory', 'cm'; 'coorz', 'cm'; 'soude', 'mas'; 'soura', 'mas';...
    'sat_pos1', 'cm'; 'sat_pos2', 'cm'; 'sat_pos3', 'cm'; 'scale', 'ppb'; 'bdclko','cm'};
% get selected parameter
allParamInPopupmenu=get(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_param']), 'String');
curSelParam=allParamInPopupmenu{get(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_param']), 'Value')};
strForText=unitPerParam{strcmp(unitPerParam(:,1), curSelParam),2};
if isempty(strForText)
    set(eval(['handles.text_plot_folder', num2str(chosenPanel), '_unit']), 'String', '')
else
    set(eval(['handles.text_plot_folder', num2str(chosenPanel), '_unit']), 'String', ['[', strForText, ']'])
end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_folder1_param_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder1_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_plot_folder1_stat.
function popupmenu_plot_folder1_stat_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder1_stat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_folder1_stat contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_folder1_stat

handles=plotToAxes(hObject, handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_folder1_stat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder1_stat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_plot_folder1_sources.
function popupmenu_plot_folder1_sources_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder1_stat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=plotToAxes(hObject, handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_plot_folder2_sources.
function popupmenu_plot_folder2_sources_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder2_sources (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=plotToAxes(hObject, handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in popupmenu_plot_folder3_sources.
function popupmenu_plot_folder3_sources_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder3_sources (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=plotToAxes(hObject, handles);

% Update handles structure
guidata(hObject, handles);



function popupmenu_plot_folder3_subfolder_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder3_subfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of popupmenu_plot_folder3_subfolder as text
%        str2double(get(hObject,'String')) returns contents of popupmenu_plot_folder3_subfolder as a double


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_folder3_subfolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder3_subfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_plot_folder3_sessions.
function popupmenu_plot_folder3_sessions_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder3_sessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_folder3_sessions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_folder3_sessions

handles=update_antennaParameterPopupmenus_plotParameters(handles,3);

handles=plotToAxes(hObject, handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_folder3_sessions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder3_sessions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_plot_folder3_param.
function popupmenu_plot_folder3_param_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder3_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_folder3_param contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_folder3_param

handles=plotToAxes(hObject, handles);

chosenPanel=3;

% write unit to text (as info for the user)
unitPerParam={'pwclk', 'cm'; 'rclk', 'cm/day'; 'qclk', 'cm/(day^2)'; ...
    'rqclk', 'cm/(day^2) warning: rate+quadr'; 'zwd', 'cm'; 'ngr', 'cm'; 'egr', 'cm';...
    'xpol', 'mas'; 'ypol', 'mas'; 'dut1', 'ms'; 'nutdx', 'mas'; 'nutdy',...
    'mas'; 'coorx', 'cm'; 'coory', 'cm'; 'coorz', 'cm'; 'soude', 'mas'; 'soura', 'mas';...
    'sat_pos1', 'cm'; 'sat_pos2', 'cm'; 'sat_pos3', 'cm'};
% get selected parameter
allParamInPopupmenu=get(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_param']), 'String');
curSelParam=allParamInPopupmenu{get(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_param']), 'Value')};
strForText=unitPerParam{strcmp(unitPerParam(:,1), curSelParam),2};
if isempty(strForText)
    set(eval(['handles.text_plot_folder', num2str(chosenPanel), '_unit']), 'String', '')
else
    set(eval(['handles.text_plot_folder', num2str(chosenPanel), '_unit']), 'String', ['[', strForText, ']'])
end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_folder3_param_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder3_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_plot_folder3_stat.
function popupmenu_plot_folder3_stat_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder3_stat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_folder3_stat contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_folder3_stat

handles=plotToAxes(hObject, handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_folder3_stat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_folder3_stat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_reset_figure in Plotting parameters
function pushbutton_reset_figure_1_Callback(hObject, eventdata, handles)

if length(handles.plotting_parameter_data) >= 1
    if isvalid(handles.plotting_parameter_data{1})
        delete(handles.plotting_parameter_data{1});
        delete(handles.plotting_parameter_current_figure_save{1});
        delete(handles.plotting_parameter_xlabel{1});
        delete(handles.plotting_parameter_current_figure_save_xlabel{1});
        handles.popupmenu_plot_folder1_param.Value=1;
        handles.popupmenu_plot_folder1_param.String=char.empty(1,0);

        handles=plotToAxes(hObject, handles);
        % Update handles structure
        guidata(hObject, handles)
    end
end

function pushbutton_reset_figure_2_Callback(hObject, eventdata, handles)

if length(handles.plotting_parameter_data) >= 2 
    if isvalid(handles.plotting_parameter_data{2})
        delete(handles.plotting_parameter_data{2});
        delete(handles.plotting_parameter_current_figure_save{2});
        delete(handles.plotting_parameter_xlabel{2});
        delete(handles.plotting_parameter_current_figure_save_xlabel{2});
        handles.popupmenu_plot_folder2_param.Value=1;
        handles.popupmenu_plot_folder2_param.String=char.empty(1,0);

        handles=plotToAxes(hObject, handles);
        % Update handles structure
        guidata(hObject, handles)
    end
end


function pushbutton_reset_figure_3_Callback(hObject, eventdata, handles)

if length(handles.plotting_parameter_data) == 3 
    if isvalid(handles.plotting_parameter_data{3})
        delete(handles.plotting_parameter_data{3});
        delete(handles.plotting_parameter_current_figure_save{3});
        delete(handles.plotting_parameter_xlabel{3});
        delete(handles.plotting_parameter_current_figure_save_xlabel{3});
        handles.popupmenu_plot_folder3_param.Value=1;
        handles.popupmenu_plot_folder3_param.String=char.empty(1,0);

        handles=plotToAxes(hObject, handles);
        % Update handles structure
        guidata(hObject, handles)
    end
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over listbox_setInput_processList.
function listbox_setInput_processList_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to listbox_setInput_processList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if there is an entry in the listbox at all
if ~isempty(get(handles.listbox_setInput_processList, 'String'))

    % create the menu at proper position
    curMousePosition=get(handles.figure_vievs2, 'CurrentPoint');
    cmenu = uicontextmenu('Parent',handles.figure_vievs2,'Position',curMousePosition);
    
    % create function handles
    %fHandles_openOPTfile=@{openOPTfile, handles};
    allSessionsInList=get(handles.listbox_setInput_processList, 'String');
    session=allSessionsInList{get(handles.listbox_setInput_processList, 'Value')};
    
   
    % create entries with the pointer to the callbacks
    item1 = uimenu(cmenu, 'Label', 'Open/Create OPT file', 'Callback', {@openOPTfile,hObject,handles});
    item2 = uimenu(cmenu, 'Label', 'Open outlier file', 'Callback', {@openOutlierFile,hObject,handles});
    % is the selected dataset a netCDF file?
%     if isempty(strfind(session,'/')) && isempty(strfind(session,'\'))
    if strfind(session, ' [vgosDB]')
        item3 = uimenu(cmenu, 'Label', 'Analyse netCDF file', 'Callback', {@analyseNetcdfFile,hObject,handles});
    end    
        
    % set visible
    set(cmenu, 'Visible', 'on');
end

function analyseNetcdfFile(src,eventdata,hObject,handles)
% This function opens the GUI to visualise the content of netCDF datasets

% get selected file
allSessionsInList=get(handles.listbox_setInput_processList, 'String');
session=allSessionsInList{get(handles.listbox_setInput_processList, 'Value')};
vgosdb_path_str = ['../DATA/vgosDB/', session(1 : (strfind(session, ' [vgosDB]')-1)), '/'];
% Open GUI
analyseNetcdf(vgosdb_path_str)



    
function openOPTfile(src,eventdata,hObject,handles)
% This function opens the OPT file (empty=create or open existing)
% of the currently selected session

% get selected file
allSessionsInList=get(handles.listbox_setInput_processList, 'String');
session=allSessionsInList{get(handles.listbox_setInput_processList, 'Value')};

% get currently selected OPT directory
OPTdirs=get(handles.popupmenu_setInput_optDir, 'String');
selectedOPTdir=OPTdirs{get(handles.popupmenu_setInput_optDir, 'value')};

% Get file format of input data file:
datatype_str = 'ngs';
if strfind(session, ' [vgosDB]')
    datatype_str = 'vgosdb'; 
elseif strfind(session, ' [VSO]')
    datatype_str = 'vso'; 
end

% Get the OPT file name and path:
switch(datatype_str)
    case 'ngs'
        optFileName = [session(1 : end-5), '.OPT'];
    case 'vso'
        optFileName = [session(1 : (strfind(session, ' [VSO]')-1)), '.OPT'];
    case 'vgosdb'
        optFileName = [session(1 : (strfind(session, ' [vgosDB]')-1)), '.OPT'];
end % switch(datatype_str)

wantedOPTfile = ['../../VLBI_OPT/', selectedOPTdir, '/', optFileName];
yrFolder = ['../../VLBI_OPT/', selectedOPTdir, '/', session(1:4)];

% if the year folder does not exist - create it
if ~exist(yrFolder, 'dir')
    mkdir(yrFolder);
end

% if the OPT file does not exist - create it
if ~exist(wantedOPTfile, 'file')
    writeNewOptFile(wantedOPTfile);
end
    
% open OPT file
if ispc
    matlabVersion = ver('MATLAB');
    if str2double(matlabVersion.Version)<=8
        dos(['start wordpad ',wantedOPTfile]);
    else
        % does not work for 7.11 (R2010b):
        winopen(wantedOPTfile) % Open with default editor which is set for *.opt files
    end
elseif isunix
    % system(['xterm -e ''vi ',wantedOPTfile '''']);
    system(['gedit ',wantedOPTfile]); % Open gedit
else % 
    warning('Operating system unknown: not able to open OPT file in text editor.');
end
    
    
    

function openOutlierFile(src,eventdata,hObject,handles)
% This function opens the outlier file of the currently selected session
    
% get selected file
allSessionsInList=get(handles.listbox_setInput_processList, 'String');
session=allSessionsInList{get(handles.listbox_setInput_processList, 'Value')};

% get currently selected Outlier directory
allOutDirInList = get(handles.popupmenu_setInput_outDir, 'String');
selOutDir = allOutDirInList{get(handles.popupmenu_setInput_outDir, 'Value')};

% Get file format of input data file:
datatype_str = 'ngs';
if strfind(session, ' [vgosDB]')
    datatype_str = 'vgosdb'; 
elseif strfind(session, ' [VSO]')
    datatype_str = 'vso'; 
end

% Get the Outlier file name and path:
switch(datatype_str)
    case 'ngs'
        sess_name_str = session(1 : end);
    case 'vso'
        sess_name_str = session(1 : (strfind(session, ' [VSO]')-1));
    case 'vgosdb'
        sess_name_str = session(1 : (strfind(session, ' [vgosDB]')-1));
end % switch(datatype_str)


if isempty(selOutDir)
    wantedOutlierFile=['../DATA/OUTLIER/', sess_name_str, '.OUT'];
else
    wantedOutlierFile=['../DATA/OUTLIER/', selOutDir,'/', sess_name_str, '.OUT'];
end

if exist(wantedOutlierFile, 'file')
    % open 
    % open(wantedOutlierFile) % Open with MATLAB Editor
    winopen(wantedOutlierFile) % Open with default editor which is set for *.OUT files
else
    msgbox('No outlier file exists yet', 'No outlier file found', 'help')
end

if exist(wantedOutlierFile, 'file')
    % open Outlier file
    if ispc
        matlabVersion = ver('MATLAB');
        if str2double(matlabVersion.Version)<=8
            dos(['start wordpad ',wantedOPTfile]);
        else
            % does not work for 7.11 (R2010b):
            winopen(wantedOutlierFile) % Open with default editor which is set for *.OUT files
        end
    elseif isunix
        % system(['xterm -e ''vi ',wantedOutlierFile '''']);
        system(['gedit ',wantedOutlierFile]); % Open gedit
    else % 
        warning('Operating system unknown: not able to open Outlier file in text editor.');
    end
else
    msgbox('Outlier file does not exists yet', 'No outlier file found', 'help')
end



% --- Executes on button press in pushbutton_plot_folder3_load.
function pushbutton_plot_folder3_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_folder3_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles=loadX_folder(hObject, handles);


% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_plot_folder2_load.
function pushbutton_plot_folder2_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_folder2_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles=loadX_folder(hObject, handles);


% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_plot_folder1_load.
function pushbutton_plot_folder1_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_folder1_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles=loadX_folder(hObject, handles);


% Update handles structure
guidata(hObject, handles);


% --- Executes when selected object is changed in uipanel_plot_folder1.
function uipanel_plot_folder1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_plot_folder1 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% update popupmenus
handles=update_antennaParameterPopupmenus_plotParameters(handles,1);

if isequal(hObject, handles.radiobutton_plot_folder1_oneSess)
    set(handles.popupmenu_plot_folder1_sessions, 'Enable', 'on')
else
    set(handles.popupmenu_plot_folder1_sessions, 'Enable', 'off')
end

handles=plotToAxes(hObject, handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes when selected object is changed in uipanel_plot_folder2.
function uipanel_plot_folder2_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_plot_folder2 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% update popupmenus
handles=update_antennaParameterPopupmenus_plotParameters(handles,2);

if isequal(hObject, handles.radiobutton_plot_folder2_oneSess)
    set(handles.popupmenu_plot_folder2_sessions, 'Enable', 'on')
else
    set(handles.popupmenu_plot_folder2_sessions, 'Enable', 'off')
end

handles=plotToAxes(hObject, handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes when selected object is changed in uipanel_plot_folder3.
function uipanel_plot_folder3_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_plot_folder3 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% update popupmenus
handles=update_antennaParameterPopupmenus_plotParameters(handles,3);

if isequal(hObject, handles.radiobutton_plot_folder3_oneSess)
    set(handles.popupmenu_plot_folder3_sessions, 'Enable', 'on')
else
    set(handles.popupmenu_plot_folder3_sessions, 'Enable', 'off')
end

handles=plotToAxes(hObject, handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_plot_data2ascii.
function pushbutton_plot_data2ascii_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_data2ascii (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% This function writes all plotted values to a txt file to a folder
% specified by ythe user (uiputfile)

% variable to define if data should be output
doWriteAscii=0;

% see if there are data
if isfield(handles, 'plot')
    if isfield(handles.plot, 'data')
        doWriteAscii=1;
    end
end

% make the "final" decision wether to output ascii file
if doWriteAscii==0;
    msgbox('No data available to be written to file', 'Error', 'warn');
else
    % get file to be specified by user
    [filename, pathname, ~] = uiputfile( ...
        {'*.txt', 'Textfile (*.txt)'},...
         'Save as',...
         '../OUT/');
     
     % see if something was chosen
    if ~isempty(filename)
        if ~strcmp('0', num2str(filename))
            % open file for writing
            fid=fopen([pathname, filename], 'w');
            
            % write header info
            c=clock;
            fprintf(fid, '# This file contains plotted data values from the Vienna VLBI Software (VieVS).\n# File was created on %02.0f.%02.0f.%04.0f (%02.0f:%02.0f:%02.1f)\n#\n', c(3), c(2), c(1), c(4), c(5), c(6));
            
            % get units for all parameters:
            unitPerParam={'pwclk', 'cm'; 'rclk', 'cm/day'; 'qclk', 'cm/(day^2)'; ...
                'rqclk', 'cm/(day^2) warning: rate+quadr'; 'zwd', 'cm'; 'ngr', 'cm'; 'egr', 'cm';...
                'xpol', 'mas'; 'ypol', 'mas'; 'dut1', 'ms'; 'nutdx', 'mas'; 'nutdy',...
                'mas'; 'coorx', 'cm'; 'coory', 'cm'; 'coorz', 'cm'; 'soude', 'mas'; 'soura', 'mas'};
            
            
            % go through all folders
            for iFolder=1:size(handles.plot.data.vals,2)
                % if there was something plotted in current folder
                if ~isempty(handles.plot.data.vals{iFolder})
                    
                    % get chosen parameter and antenna
                    allParameters=eval(['get(handles.popupmenu_plot_folder', num2str(iFolder), '_param, ''String'')']);
                    if ~isempty(allParameters) % if data were removed from plot window by clicking Reset
                        allAntennas=eval(['get(handles.popupmenu_plot_folder', num2str(iFolder), '_stat, ''String'')']);
                        allSubfolders=eval(['get(handles.popupmenu_plot_folder', num2str(iFolder),'_subfolder, ''String'')']);
                        chosenParameter=allParameters{eval(['get(handles.popupmenu_plot_folder', num2str(iFolder), '_param, ''Value'')'])};
                    
                        chosenAntenna=allAntennas{eval(['get(handles.popupmenu_plot_folder', num2str(iFolder), '_stat, ''Value'')'])};
                        chosenSubfolder=allSubfolders{eval(['get(handles.popupmenu_plot_folder', num2str(iFolder), '_subfolder, ''Value'')'])};
                    
                        % get unit of param
                        indChosenParamInUnitCell=find(strcmp(unitPerParam(:,1),chosenParameter));
                        if isempty(indChosenParamInUnitCell)
                            chosenParameterUnit='unknown';
                        else
                            chosenParameterUnit=unitPerParam{indChosenParamInUnitCell,2};
                        end
                    
                    
                        fprintf(fid, '# ========================================\n# Folder %1.0f (subfolder: %s)\n# parameter: %s (unit: %s), antenna: %s\n#  mjd (%%16.7e)   value (%%16.7e)    std (%%16.7e)\n', ...
                            iFolder, chosenSubfolder, ...
                            chosenParameter, chosenParameterUnit, chosenAntenna);
                    
                        % number of values
                        nVals=length(handles.plot.data.vals{iFolder});
                        for iVal=1:nVals
                            fprintf(fid, '%16.7e%16.7e%16.7e\n', handles.plot.data.mjds{iFolder}(iVal), ...
                                handles.plot.data.vals{iFolder}(iVal), ...
                                handles.plot.data.mxs{iFolder}(iVal));
                        end
                    end 
                end
            end
        end
    end
    
     
end


% --------------------------------------------------------------------
function menu_file_parameterFiles_saveProcessList_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_parameterFiles_saveProcessList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get a file
[filename, pathname, filterindex] = uiputfile( ...
{'process_list.mat', 'Matlab Binary Format (*.mat)'},...
 'Save as',...
 '../WORK/PROCESSLIST/');

% see if something was chosen
if ~isempty(filename)
    if ~strcmp('0', num2str(filename))
        % call function for saving a parameter file
        saveProcessList(hObject, handles, [pathname, filename])
        
        % write message box
        msgbox('Process list successfully saved!', 'Save process list', 'help');
    end
end


% VIE_SIM start -----------------------------------------------------------
% --------------------------------------------------------------------
function menu_simulation_Callback(hObject, eventdata, handles)
% hObject    handle to menu_simulation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton1.
function pushbutton_vie_sim_saveParams_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

saveSimParameters(hObject, handles)


function edit_vie_sim_nDaysSim_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_nDaysSim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sim_nDaysSim as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sim_nDaysSim as a double

% --- Executes during object creation, after setting all properties.
function edit_vie_sim_nDaysSim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_nDaysSim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vie_sim_startInd_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_startInd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sim_startInd as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sim_startInd as a double

% --- Executes during object creation, after setting all properties.
function edit_vie_sim_startInd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_startInd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox_vie_sim_setRefClk0.
function checkbox_vie_sim_setRefClk0_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sim_setRefClk0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sim_setRefClk0

function edit_vie_sim_tropo_dhseg_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_tropo_dhseg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sim_tropo_dhseg as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sim_tropo_dhseg as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sim_tropo_dhseg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_tropo_dhseg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vie_sim_tropo_dh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_tropo_dh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sim_tropo_dh as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sim_tropo_dh as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sim_tropo_dh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_tropo_dh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listbox_vie_sim_paramFile.
function listbox_vie_sim_paramFile_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_vie_sim_paramFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_vie_sim_paramFile contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_vie_sim_paramFile

% --- Executes during object creation, after setting all properties.
function listbox_vie_sim_paramFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_vie_sim_paramFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vie_sim_noise_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sim_noise as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sim_noise as a double

% --- Executes during object creation, after setting all properties.
function edit_vie_sim_noise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vie_sim_clock_asd_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_clock_asd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sim_clock_asd as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sim_clock_asd as a double

% --- Executes during object creation, after setting all properties.
function edit_vie_sim_clock_asd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_clock_asd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vie_sim_clock_asdInt_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_clock_asdInt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sim_clock_asdInt as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sim_clock_asdInt as a double

% --- Executes during object creation, after setting all properties.
function edit_vie_sim_clock_asdInt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_clock_asdInt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vie_sim_tropo_vn_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_tropo_vn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sim_tropo_vn as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sim_tropo_vn as a double

% --- Executes during object creation, after setting all properties.
function edit_vie_sim_tropo_vn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_tropo_vn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vie_sim_tropo_ve_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_tropo_ve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sim_tropo_ve as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sim_tropo_ve as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sim_tropo_ve_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_tropo_ve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vie_sim_tropo_cn_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_tropo_cn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sim_tropo_cn as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sim_tropo_cn as a double

% --- Executes during object creation, after setting all properties.
function edit_vie_sim_tropo_cn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_tropo_cn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vie_sim_tropo_H_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_tropo_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sim_tropo_H as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sim_tropo_H as a double

% --- Executes during object creation, after setting all properties.
function edit_vie_sim_tropo_H_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_tropo_H (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vie_sim_tropo_wzd0_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_tropo_wzd0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sim_tropo_wzd0 as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sim_tropo_wzd0 as a double

% --- Executes during object creation, after setting all properties.

function edit_vie_sim_tropo_wzd0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_tropo_wzd0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox_vie_sim_writeNgsFile.
function checkbox_vie_sim_writeNgsFile_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sim_writeNgsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sim_writeNgsFile

% VIE_SIM end -------------------------------------------------------------




function edit_vie_sched_sundistmin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_sundistmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sched_sundistmin as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sched_sundistmin as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sched_sundistmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_sundistmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vie_sched_cutel_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_cutel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sched_cutel as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sched_cutel as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sched_cutel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_cutel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vie_sched_flux_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_flux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sched_flux as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sched_flux as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sched_flux_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_flux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_vie_sched_edsnrxband_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_edsnrxband (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sched_edsnrxband as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sched_edsnrxband as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sched_edsnrxband_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_edsnrxband (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vie_sched_edsnrsband_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_edsnrsband (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sched_edsnrsband as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sched_edsnrsband as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sched_edsnrsband_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_edsnrsband (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_vie_sched_days_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_days (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sched_days as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sched_days as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sched_days_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_days (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vie_sched_mons_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_mons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sched_mons as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sched_mons as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sched_mons_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_mons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vie_sched_years_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_years (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sched_years as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sched_years as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sched_years_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_years (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vie_sched_hours_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_hours (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sched_hours as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sched_hours as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sched_hours_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_hours (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vie_sched_mins_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_mins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sched_mins as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sched_mins as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sched_mins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_mins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vie_sched_secs_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_secs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sched_secs as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sched_secs as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sched_secs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_secs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vie_sched_duration_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sched_duration as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sched_duration as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sched_duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listbox_vie_sched_selectsta.
function listbox_vie_sched_selectsta_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_vie_sched_selectsta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_vie_sched_selectsta contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_vie_sched_selectsta

% common part
%global staname stanetstr;
allStatsInList=cellstr(get(handles.listbox_vie_sched_selectsta, 'String'));
selectedStat=allStatsInList(get(handles.listbox_vie_sched_selectsta, 'Value'));

% make the selected stations 8 digits long
if length(selectedStat{1})~=8
    selectedStat{1}(length(selectedStat{1})+1:8)=' ';
end

% set value to one if needed (otherwise: matlab error!)
if isempty(get(handles.listbox_vie_sched_stanet, 'Value'))
    set(handles.listbox_vie_sched_stanet, 'value', 1)
end
% update listbox
set(handles.listbox_vie_sched_stanet, 'String', ...
    unique([get(handles.listbox_vie_sched_stanet, 'String'); selectedStat]))

% Choose default command line output for vie_setup
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% staid = get(handles.listbox_vie_sched_selectsta,'Value');
% station(1:8) = staname(staid,1:8);
% stanetstr = strcat(stanetstr,station,'|');
% set(handles.listbox_vie_sched_stanet,'String',stanetstr);

% --- Executes during object creation, after setting all properties.
function listbox_vie_sched_selectsta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_vie_sched_selectsta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% -listbox_vie_sched_selectsta-- Executes on selection change in listbox_vie_sched_stanet.
function listbox_vie_sched_stanet_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_vie_sched_stanet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_vie_sched_stanet contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_vie_sched_stanet

% common part
% global stanetstr;
% get stations in listbox
% stanetstr=get(handles.listbox_vie_sched_stanet, 'String');
% 
% if ~isempty(stanetstr)
% 
%     % get chosen index
%     delstaid = get(handles.listbox_vie_sched_stanet,'Value');
%     
%     % delete entry
%     stanetstr(delstaid)=[];
%     
%     % set value to one
%     set(handles.listbox_vie_sched_stanet,'Value',1);
%     
%     % update listbox
%     set(handles.listbox_vie_sched_stanet,'String',stanetstr);
%     
%     % Choose default command line output for vie_setup
%     handles.output = hObject;
% 
%     % Update handles structure
%     guidata(hObject, handles);
% endv

% --- Executes during object creation, after setting all properties.
function listbox_vie_sched_stanet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_vie_sched_stanet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listbox_vie_sched_prenet.
function listbox_vie_sched_prenet_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_vie_sched_prenet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_vie_sched_prenet contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_vie_sched_prenet

allLists=cellstr(get(handles.listbox_vie_sched_prenet, 'String'));
curList=allLists{get(handles.listbox_vie_sched_prenet, 'Value')};

curListContent=load(['../WORK/STATIONLIST/', curList]);
fieldn=fieldnames(curListContent);
curListContent=curListContent.(fieldn{1});

set(handles.listbox_vie_sched_stanet, 'String', curListContent);

% % common part
% %global stanetstr;
% % define stations of default lists
% netvlbi2010={'HART2010';'KOKE2010';'NYAL2010';'GOLD2010';'URUM2010';'WARK2010';'TSUK2010';'AREC2010';'FORT2010';'GGAO2010';'KATH2010';'MTPL2010';'YARR2010';'KORE2010';'AZOR2010';'CNAR2010';'WETZ2010';'YEBE2010'};
% netivsr1   ={'WESTFORD';'TIGOCONC';'TSUKUB32';'SESHAN25';'SVETLOE ';'ZELENCHK';'HARTRAO ';'WETTZELL';'NYALES20'};
% netivsr4   ={'TIGOCONC';'KOKEE   ';'ZELENCHK';'MEDICINA';'WETTZELL';'BADARY  ';'SVETLOE ';'NYALES20'};
% netint1    ={'KOKEE   ';'WETTZELL'};
%     
% % get value of selected entry
% netid = get(handles.listbox_vie_sched_prenet,'Value');
% allStatListsInListbox=get(handles.listbox_vie_sched_prenet, 'String');
% 
% set(handles.listbox_vie_sched_stanet, 'Value', 1)
% 
% % if VLBI2010 was selected
% if strcmp(allStatListsInListbox{netid}, 'VLBI2010')
%     set(handles.listbox_vie_sched_stanet,'String',...
%         unique([get(handles.listbox_vie_sched_stanet, 'String'); netvlbi2010]));
% end
% % if IVSR1 was selected
% if strcmp(allStatListsInListbox{netid}, 'IVSR1')
%     set(handles.listbox_vie_sched_stanet,'String',...
%         unique([get(handles.listbox_vie_sched_stanet, 'String'); netivsr1]));
% end
% % if IVSR4 was selected
% if strcmp(allStatListsInListbox{netid}, 'IVSR4')
%     set(handles.listbox_vie_sched_stanet,'String',...
%         unique([get(handles.listbox_vie_sched_stanet, 'String'); netivsr4]));
% end
% % if INT1 was selected
% if strcmp(allStatListsInListbox{netid}, 'INT1')
%     set(handles.listbox_vie_sched_stanet,'String',...
%         unique([get(handles.listbox_vie_sched_stanet, 'String'); netint1]));
% end

% --- Executes during object creation, after setting all properties.
function listbox_vie_sched_prenet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_vie_sched_prenet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_vie_sched_srcnum_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_srcnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sched_srcnum as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sched_srcnum as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sched_srcnum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_srcnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_vie_sched_saveParams.
function pushbutton_vie_sched_saveParams_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_vie_sched_saveParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% run the function to save the GUI sched parameters
[vie_sched_error_code] = saveSchedParameters(hObject, handles);
if vie_sched_error_code ~= 0
    fprintf(1, 'VIE_SCHED input parameter error! Error code: %g\n', vie_sched_error_code);
end

% --- Executes when selected object is changed in uipanel208.
function uipanel208_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel208 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% if "specify now" was ticked
if strcmp(get(hObject, 'Tag'), 'radiobutton_vie_sim_specifyNow')
    set(handles.listbox_vie_sim_paramFile, 'Enable', 'off')
    % set enabling
    newState='on';
else
    newState='off';
    set(handles.listbox_vie_sim_paramFile, 'Enable', 'on')
end
set(handles.uipanel_vie_sim_text1, 'Enable', newState);
set(handles.edit_vie_sim_tropo_cn, 'Enable', newState);
set(handles.uipanel_vie_sim_text3, 'Enable', newState);
set(handles.edit_vie_sim_tropo_H, 'Enable', newState);
set(handles.uipanel_vie_sim_text5, 'Enable', newState);
set(handles.edit_vie_sim_tropo_ve, 'Enable', newState);
set(handles.uipanel_vie_sim_text14, 'Enable', newState);
set(handles.edit_vie_sim_tropo_dh, 'Enable', newState);
set(handles.uipanel_vie_sim_text12, 'Enable', newState);
set(handles.edit_vie_sim_tropo_wzd0, 'Enable', newState);
set(handles.uipanel_vie_sim_text4, 'Enable', newState);
set(handles.edit_vie_sim_tropo_vn, 'Enable', newState);
set(handles.uipanel_vie_sim_text13, 'Enable', newState);
set(handles.edit_vie_sim_tropo_dhseg, 'Enable', newState);
set(handles.uipanel_vie_sim_text9, 'Enable', newState);
set(handles.edit_vie_sim_clock_asd, 'Enable', newState);
set(handles.uipanel_vie_sim_text10, 'Enable', newState);
set(handles.edit_vie_sim_clock_asdInt, 'Enable', newState);
set(handles.uipanel_vie_sim_text16, 'Enable', newState);
set(handles.uipanel_vie_sim_text11, 'Enable', newState);
set(handles.edit_vie_sim_noise, 'Enable', newState);
set(handles.text_vie_sim_wn_sat, 'Enable', newState);
set(handles.edit_vie_sim_noise_sat, 'Enable', newState);


% --- Executes when selected object is changed in uipanel_vie_sched_rbgstrategy.
function uipanel_vie_sched_rbgstrategy_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_vie_sched_rbgstrategy 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% if "source based strategy" is ticked
if strcmp(get(hObject, 'Tag'), get(handles.radiobutton_vie_sched_rbsource, 'Tag'))
    newState='On';
    set(handles.checkbox_vie_sched_distribute, 'Enable', 'Off');
    % enable sched analyser
    set(handles.checkbox_vie_sched_openSchedAnalyser,'Enable','on');
	
% if "station based strategy" is ticked
elseif strcmp(get(hObject, 'Tag'), get(handles.radiobutton_vie_sched_rbstation, 'Tag'))
    newState='Off';
    set(handles.checkbox_vie_sched_distribute, 'Enable', 'On');
    % enable sched analyser
    set(handles.checkbox_vie_sched_openSchedAnalyser,'Enable','on');
else 
    set(handles.checkbox_vie_sched_distribute, 'Enable', 'Off');
	newState='Off';
end

% set enabling
set(handles.text146, 'Enable', newState)
set(handles.edit_vie_sched_srcnum, 'Enable', newState)
set(handles.text147, 'Enable', newState)

% if "Satellite obbservation" is ticked
if strcmp(get(hObject, 'Tag'), get(handles.radiobutton_vie_sched_satellites, 'Tag'))
    % Set satellite scheduling panel visible
	set(handles.uipanel_vie_sched_sat_sched,'Visible','On')
	
	% --- get available TLE file names ---
    path_tle_dir = '../ORBIT/TLE/';    % TLE data folder
	[tle_str] = getTLEFileNames(path_tle_dir);
	
	% Write TLE file names to popup menu:
	set(handles.popupmenu_vie_sched_select_tle,'String',tle_str);
    
    % disable sched analyser
    set(handles.checkbox_vie_sched_openSchedAnalyser,'Value',0)
    set(handles.checkbox_vie_sched_openSchedAnalyser,'Enable','off')
else
	% Set satellite scheduling panel invisible
    set(handles.uipanel_vie_sched_sat_sched,'Visible','Off')
end

% if "manually" is ticked
if strcmp(get(hObject, 'Tag'), get(handles.radiobutton_vie_sched_manually, 'Tag'))
    % enable sched analyser
    set(handles.checkbox_vie_sched_openSchedAnalyser,'Enable','on')
end
% --------------------------------------------------------------------
function menu_estimation_globalSolution_Callback(hObject, eventdata, handles)
% hObject    handle to menu_estimation_globalSolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_estimation_globalSolution, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in checkbox_run_estParameters.
function checkbox_run_estParameters_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_estParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_estParameters



function edit_run_outDirs_glob_outDirVieGlob_Callback(hObject, eventdata, handles)
% hObject    handle to edit_run_outDirs_glob_outDirVieGlob (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_run_outDirs_glob_outDirVieGlob as text
%        str2double(get(hObject,'String')) returns contents of edit_run_outDirs_glob_outDirVieGlob as a double


% --- Executes during object creation, after setting all properties.
function edit_run_outDirs_glob_outDirVieGlob_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_run_outDirs_glob_outDirVieGlob (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_glob_param_other_maxRMS_Callback(hObject, eventdata, handles)
% hObject    handle to edit_glob_param_other_maxRMS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_glob_param_other_maxRMS as text
%        str2double(get(hObject,'String')) returns contents of edit_glob_param_other_maxRMS as a double


% --- Executes during object creation, after setting all properties.
function edit_glob_param_other_maxRMS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_glob_param_other_maxRMS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_global_param_other_back.
function checkbox_global_param_other_back_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_global_param_other_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_global_param_other_back



function edit_run_outDirs_glob_level2Sub_Callback(hObject, eventdata, handles)
% hObject    handle to edit_run_outDirs_glob_level2Sub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_run_outDirs_glob_level2Sub as text
%        str2double(get(hObject,'String')) returns contents of edit_run_outDirs_glob_level2Sub as a double


% --- Executes during object creation, after setting all properties.
function edit_run_outDirs_glob_level2Sub_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_run_outDirs_glob_level2Sub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_run_outDirs_glob_pathToLevel2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_run_outDirs_glob_pathToLevel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_run_outDirs_glob_pathToLevel2 as text
%        str2double(get(hObject,'String')) returns contents of edit_run_outDirs_glob_pathToLevel2 as a double


% --- Executes during object creation, after setting all properties.
function edit_run_outDirs_glob_pathToLevel2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_run_outDirs_glob_pathToLevel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_run_set_default_path.
function pushbutton_run_set_default_path_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_run_set_default_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get default to edit textfield
set(handles.edit_run_outDirs_glob_pathToLevel2, 'String', '../DATA/LEVEL2/')


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_vie_glob, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_vie_glob_TRFCRFparam, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in checkbox_vie_sim_simSwd.
function checkbox_vie_sim_simSwd_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sim_simSwd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sim_simSwd


% --- Executes on button press in checkbox_vie_sim_simClock.
function checkbox_vie_sim_simClock_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sim_simClock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sim_simClock


% --- Executes on button press in checkbox_vie_sim_simNoise.
function checkbox_vie_sim_simNoise_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sim_simNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sim_simNoise


% --- Executes on button press in pushbutton_plot_save_saveAs.
function pushbutton_plot_save_saveAs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_save_saveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% only save if invisible figure exists
if isfield(handles, 'figure_save')
    
    allSaveAsEndingsOrig={'*.fig', 'Matlab Figure (*.fig)';...
    '*.jpg', 'JPEG 24-bit (*.jpg)'; ...
    '*.bmp', 'Bitmap 24-bit (*.bmp)';...
    '*.png', 'PNG 24-bit (*.png';...
    '*.eps', 'EPS color (Vector) (*.eps)';...
    '*.eps', 'EPS black and white (Vector) (*.eps)';...
    '*.pdf', 'PDF color (Vector) (*.pdf)'};

    % let user chose file
    [filename, pathname, filterindex] = uiputfile( ...
    allSaveAsEndingsOrig,...
     'Save as');
 
     % see if something was chosen
    if ~isempty(filename)
        if ~strcmp('0', num2str(filename))
            allTypes={' ', '-djpeg', '-dbmp', '-dpng', '-depsc', '-deps', '-dpdf'};
            
            % set size of figure
            curWidth=8; %[cm]
            curHeight=6;
            curResolution=150; %[dpi]
            set(handles.figure_save, 'PaperUnits', 'centimeters', ...
                'PaperPosition', [0 0 curWidth curHeight],...
                'PaperSize',     [curWidth, curHeight]);
            
            % set limits (important when figure is zoomed)
            curXlim=get(get(handles.figure_vievs2, 'currentaxes'), 'xlim');
            curYlim=get(get(handles.figure_vievs2, 'currentaxes'), 'ylim');
            set(get(handles.figure_save, 'currentaxes'), 'xlim', curXlim);
            set(get(handles.figure_save, 'currentaxes'), 'ylim', curYlim);
            
            % print
            if strcmp(allSaveAsEndingsOrig{filterindex,1}, '*.fig') ==1
                set(handles.figure_save, 'Visible', 'On');
                saveas(handles.figure_save, [pathname, '\', filename]);
                set(handles.figure_save, 'Visible', 'Off');
            else
                print(handles.figure_save, allTypes{filterindex}, [pathname, '\', filename], ['-r', num2str(curResolution)]);
            end
            
        end
    end
    
    
else
    msgbox('Nothing has been printed yet', 'Plot something', 'warn');
end

% --- Executes on button press in pushbutton_plot_save_runOwnCode.
function pushbutton_plot_save_runOwnCode_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_save_runOwnCode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if we have the invisible saving-figure
if isfield(handles, 'figure_save')
    % make saving-figure current figure (since user uses gcf)
    set(0,'CurrentFigure',handles.figure_save);
    
    % set limits (important when figure is zoomed)
    curXlim=get(get(handles.figure_vievs2, 'currentaxes'), 'xlim');
    curYlim=get(get(handles.figure_vievs2, 'currentaxes'), 'ylim');
    set(get(handles.figure_save, 'currentaxes'), 'xlim', curXlim);
    set(get(handles.figure_save, 'currentaxes'), 'ylim', curYlim);

    %  try to run that code
    try 
        eval(get(handles.edit_plot_save_ownCode, 'String'))
    catch exception
        msgbox(exception.message, 'Following matlab error occured', 'warn');
    end

    % make main figure the current figure again
    set(0,'CurrentFigure',handles.figure_vievs2);
else
    msgbox('Plot something first', 'Nothing plotted yet', 'warn')
end



function edit_plot_save_ownCode_Callback(hObject, eventdata, handles)
% hObject    handle to edit_plot_save_ownCode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_plot_save_ownCode as text
%        str2double(get(hObject,'String')) returns contents of edit_plot_save_ownCode as a double


% --- Executes during object creation, after setting all properties.
function edit_plot_save_ownCode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_plot_save_ownCode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_run_outputDirectories_runVieGlob.
function checkbox_run_outputDirectories_runVieGlob_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_outputDirectories_runVieGlob (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_outputDirectories_runVieGlob

% If Vie_GLOB and Vie_LSM (scanwise update) are run both:
if get(handles.checkbox_run_outputDirectories_runVieGlob,'Value') && ...
(get(handles.checkbox_run_outputDirectories_runVieLsmScanwiseUpdate,'Value') || get(handles.checkbox_run_outputDirectories_runVieLsm,'Value'))
    
		% Use different output Sub-directories
		if get(handles.checkbox_run_outDirs_diffSubs, 'Value')
			str_lsm_subDir=get(handles.edit_run_outDirs_level2, 'String');
		% Use one output Sub-directory	
		else
			str_lsm_subDir=get(handles.edit_run_outDirs_oneSub, 'String');
		end
		
		set(handles.edit_run_outDirs_glob_pathToLevel2, 'Enable', 'off');
		set(handles.edit_run_outDirs_glob_pathToLevel2, 'String', '../DATA/LEVEL2/');
		set(handles.edit_run_outDirs_glob_level2Sub, 'Enable', 'off');
		set(handles.edit_run_outDirs_glob_level2Sub, 'String', str_lsm_subDir);
		
		msgbox('If Vie_LSM and Vie_GLOB are both run, the Vie_LSM output directory for N-matrices have to be equal to the Vie_GLOB input directory!', 'Warning (subdirectories not equal)');
else
	set(handles.edit_run_outDirs_glob_pathToLevel2, 'Enable', 'on');
	set(handles.edit_run_outDirs_glob_level2Sub, 'Enable', 'on');
end



% --- Executes on button press in checkbox_run_outputDirectories_runVieSched.
function checkbox_run_outputDirectories_runVieSched_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_outputDirectories_runVieSched (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_outputDirectories_runVieSched


% --- Executes on button press in checkbox_run_outputDirectories_runVieSim.
function checkbox_run_outputDirectories_runVieSim_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_outputDirectories_runVieSim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_outputDirectories_runVieSim


% --- Executes on button press in checkbox_run_outputDirectories_runVieInit.
function checkbox_run_outputDirectories_runVieInit_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_outputDirectories_runVieInit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_outputDirectories_runVieInit


% --- Executes on button press in pushbutton_run_runOptions_saveAllParamsToAscii.
function pushbutton_run_runOptions_saveAllParamsToAscii_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_run_runOptions_saveAllParamsToAscii (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get a file
[filename, pathname, filterindex] = uiputfile( ...
{'*.txt', 'Textfile (*.txt)'},...
 'Save as',...
 '../WORK/PARAMETERS/');

% see if something was chosen
if ~isempty(filename)
    if ~strcmp('0', num2str(filename))
        % first auto save parameters (again)
        auto_save_parameterfile(hObject, handles);
        
        % get newest auto saved parameter file
        newestAutoSavedParameterFile=getNewestAutoSavedParameterFile(hObject, handles);
        
        % load newest parameter file
        load(['PARAMETERS/', newestAutoSavedParameterFile.name]);
        
        % write parameter file to text output
        create_input_protocol(parameter, [pathname, filename]);
        
        msgbox(sprintf('GUI options have been saved to\n%s', [pathname, filename]), 'Saving completed', 'help');
    end
end

% --- Executes on button press in checkbox_run_outputDirectories_runVieMod.
function checkbox_run_outputDirectories_runVieMod_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_outputDirectories_runVieMod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_outputDirectories_runVieMod


% --- Executes on button press in checkbox_run_outputDirectories_runVieLsm.
function checkbox_run_outputDirectories_runVieLsm_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_outputDirectories_runVieLsm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_outputDirectories_runVieLsm

if get(hObject, 'Value')
	% if checkbox is ticked -> untick scanwise update
    set(handles.checkbox_run_outputDirectories_runVieLsmScanwiseUpdate, 'Value', 0);
	
    % if also Vie_GLOB is ticked
	if get(handles.checkbox_run_outputDirectories_runVieGlob, 'value')
    
		% Use different output Sub-directories
		if get(handles.checkbox_run_outDirs_diffSubs, 'Value')
			str_lsm_subDir=get(handles.edit_run_outDirs_level2, 'String');
		% Use one output Sub-directory	
		else
			str_lsm_subDir=get(handles.edit_run_outDirs_oneSub, 'String');
		end
		
		%str_glob_subDir=get(handles.edit_run_outDirs_glob_level2Sub, 'String');
		
		set(handles.edit_run_outDirs_glob_pathToLevel2, 'Enable', 'off');
		set(handles.edit_run_outDirs_glob_pathToLevel2, 'String', '../DATA/LEVEL2/');
		set(handles.edit_run_outDirs_glob_level2Sub, 'Enable', 'off');
		set(handles.edit_run_outDirs_glob_level2Sub, 'String', str_lsm_subDir);
		
		msgbox('If Vie_LSM and Vie_GLOB are both run, the Vie_LSM output directory for N-matrices have to be equal to the Vie_GLOB input directory!', 'Warning (subdirectories not equal)');
		
	end
	
else
	set(handles.edit_run_outDirs_glob_pathToLevel2, 'Enable', 'on');
	set(handles.edit_run_outDirs_glob_level2Sub, 'Enable', 'on');
end

% --- Executes on selection change in popupmenu_vie_glob_trfCrf_trf_stations4Nnt.
function popupmenu_vie_glob_trfCrf_trf_stations4Nnt_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_vie_glob_trfCrf_trf_stations4Nnt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_vie_glob_trfCrf_trf_stations4Nnt contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_vie_glob_trfCrf_trf_stations4Nnt


% --- Executes during object creation, after setting all properties.
function popupmenu_vie_glob_trfCrf_trf_stations4Nnt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_vie_glob_trfCrf_trf_stations4Nnt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_vie_glob_trfCrf_trf_positDisc.
function popupmenu_vie_glob_trfCrf_trf_positDisc_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_vie_glob_trfCrf_trf_positDisc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_vie_glob_trfCrf_trf_positDisc contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_vie_glob_trfCrf_trf_positDisc


% --- Executes during object creation, after setting all properties.
function popupmenu_vie_glob_trfCrf_trf_positDisc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_vie_glob_trfCrf_trf_positDisc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_vie_glob_trfCrf_trf_reduceStats.
function checkbox_vie_glob_trfCrf_trf_reduceStats_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_glob_trfCrf_trf_reduceStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_glob_trfCrf_trf_reduceStats

if get(hObject, 'Value')
    newState='on';
else
    newState='off';
end

set(handles.popupmenu_vie_glob_trfCrf_trf_sessionWiseRedStat, 'Enable', newState);
set(handles.text278, 'Enable', newState);

% --- Executes on button press in checkbox_vie_glob_trfCrf_trf_keepConstVel.
function checkbox_vie_glob_trfCrf_trf_keepConstVel_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_glob_trfCrf_trf_keepConstVel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_glob_trfCrf_trf_keepConstVel

if get(hObject, 'Value')
    newState='on';
else
    newState='off';
end

set(handles.popupmenu_vie_glob_trfCrf_trf_keepConstVel, 'Enable', newState);
set(handles.text281, 'Enable', newState);


% --- Executes on button press in checkbox_vie_glob_trfCrf_trf_useVelTies.
function checkbox_vie_glob_trfCrf_trf_useVelTies_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_glob_trfCrf_trf_useVelTies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_glob_trfCrf_trf_useVelTies

if get(hObject, 'Value')
    newState='on';
else
    newState='off';
end

set(handles.popupmenu_vie_glob_trfCrf_trf_velTies, 'Enable', newState);
set(handles.text280, 'Enable', newState);


% --- Executes on selection change in popupmenu_vie_glob_trfCrf_trf_sessionWiseRedStat.
function popupmenu_vie_glob_trfCrf_trf_sessionWiseRedStat_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_vie_glob_trfCrf_trf_sessionWiseRedStat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_vie_glob_trfCrf_trf_sessionWiseRedStat contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_vie_glob_trfCrf_trf_sessionWiseRedStat


% --- Executes during object creation, after setting all properties.
function popupmenu_vie_glob_trfCrf_trf_sessionWiseRedStat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_vie_glob_trfCrf_trf_sessionWiseRedStat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_vie_glob_trfCrf_trf_keepConstVel.
function popupmenu_vie_glob_trfCrf_trf_keepConstVel_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_vie_glob_trfCrf_trf_keepConstVel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_vie_glob_trfCrf_trf_keepConstVel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_vie_glob_trfCrf_trf_keepConstVel


% --- Executes during object creation, after setting all properties.
function popupmenu_vie_glob_trfCrf_trf_keepConstVel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_vie_glob_trfCrf_trf_keepConstVel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_vie_glob_trfCrf_trf_velTies.
function popupmenu_vie_glob_trfCrf_trf_velTies_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_vie_glob_trfCrf_trf_velTies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_vie_glob_trfCrf_trf_velTies contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_vie_glob_trfCrf_trf_velTies


% --- Executes during object creation, after setting all properties.
function popupmenu_vie_glob_trfCrf_trf_velTies_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_vie_glob_trfCrf_trf_velTies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_vie_glob_trfCrf_crf_nnrTransInDecl.
function checkbox_vie_glob_trfCrf_crf_nnrTransInDecl_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_glob_trfCrf_crf_nnrTransInDecl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_glob_trfCrf_crf_nnrTransInDecl
if get(hObject, 'value')
    set(handles.checkbox_vie_glob_trfCrf_crf_onlyNnrOnSources, 'Value', 0)
    set(handles.popupmenu_vie_glob_trfCrf_crf_sourcesForNnr, 'Enable', 'on');
    set(handles.text283, 'Enable', 'on');
else
    % if also 2nd checkbox is unticked -> make popupmenu disabled
    if get(handles.checkbox_vie_glob_trfCrf_crf_onlyNnrOnSources, 'Value')==0
        set(handles.popupmenu_vie_glob_trfCrf_crf_sourcesForNnr, 'Enable', 'off');
        set(handles.text283, 'Enable', 'off');
    end
end


% --- Executes on button press in checkbox_vie_glob_trfCrf_crf_onlyNnrOnSources.
function checkbox_vie_glob_trfCrf_crf_onlyNnrOnSources_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_glob_trfCrf_crf_onlyNnrOnSources (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_glob_trfCrf_crf_onlyNnrOnSources

if get(hObject, 'value')
    set(handles.checkbox_vie_glob_trfCrf_crf_nnrTransInDecl, 'Value', 0)
    set(handles.popupmenu_vie_glob_trfCrf_crf_sourcesForNnr, 'Enable', 'on');
    set(handles.text283, 'Enable', 'on');
else
    % if also 2nd checkbox is unticked -> make popupmenu disabled
    if get(handles.checkbox_vie_glob_trfCrf_crf_nnrTransInDecl, 'Value')==0
        set(handles.popupmenu_vie_glob_trfCrf_crf_sourcesForNnr, 'Enable', 'off');
        set(handles.text283, 'Enable', 'off');
    end
end


% --- Executes on button press in checkbox_vie_glob_trfCrf_crf_fixSomeSources.
function checkbox_vie_glob_trfCrf_crf_fixSomeSources_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_glob_trfCrf_crf_fixSomeSources (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_glob_trfCrf_crf_fixSomeSources
if get(hObject, 'Value')
    newState='on';
else
    newState='off';
end

set(handles.popupmenu_vie_glob_trfCrf_crf_fixedSources, 'Enable', newState);
set(handles.text284, 'Enable', newState);


% --- Executes on button press in checkbox_vie_glob_trfCrf_crf_sessionWiseReduceSources.
function checkbox_vie_glob_trfCrf_crf_sessionWiseReduceSources_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_glob_trfCrf_crf_sessionWiseReduceSources (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_glob_trfCrf_crf_sessionWiseReduceSources

if get(hObject, 'Value')
    newState='on';
else
    newState='off';
end

set(handles.popupmenu_vie_glob_trfCrf_crf_sessionWiseRedSources, 'Enable', newState);
set(handles.text285, 'Enable', newState);


% --- Executes on selection change in popupmenu_vie_glob_trfCrf_crf_sourcesForNnr.
function popupmenu_vie_glob_trfCrf_crf_sourcesForNnr_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_vie_glob_trfCrf_crf_sourcesForNnr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_vie_glob_trfCrf_crf_sourcesForNnr contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_vie_glob_trfCrf_crf_sourcesForNnr


% --- Executes during object creation, after setting all properties.
function popupmenu_vie_glob_trfCrf_crf_sourcesForNnr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_vie_glob_trfCrf_crf_sourcesForNnr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_vie_glob_trfCrf_crf_fixedSources.
function popupmenu_vie_glob_trfCrf_crf_fixedSources_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_vie_glob_trfCrf_crf_fixedSources (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_vie_glob_trfCrf_crf_fixedSources contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_vie_glob_trfCrf_crf_fixedSources


% --- Executes during object creation, after setting all properties.
function popupmenu_vie_glob_trfCrf_crf_fixedSources_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_vie_glob_trfCrf_crf_fixedSources (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_vie_glob_trfCrf_crf_sessionWiseRedSources.
function popupmenu_vie_glob_trfCrf_crf_sessionWiseRedSources_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_vie_glob_trfCrf_crf_sessionWiseRedSources (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_vie_glob_trfCrf_crf_sessionWiseRedSources contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_vie_glob_trfCrf_crf_sessionWiseRedSources


% --- Executes during object creation, after setting all properties.
function popupmenu_vie_glob_trfCrf_crf_sessionWiseRedSources_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_vie_glob_trfCrf_crf_sessionWiseRedSources (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_vie_glob_trfCrf_saveGlobParams.
function pushbutton_vie_glob_trfCrf_saveGlobParams_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_vie_glob_trfCrf_saveGlobParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

saveGlobParameters(hObject, handles)


% --- Executes when selected object is changed in uipanel29.
function uipanel29_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel29 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% if source coordinates are estimated -> enable parts in trf/crf globe gui
if strcmp(get(hObject, 'Tag'), 'radiobuttonedit_glob_param_estSrcCoords')
    set(handles.checkbox_vie_glob_trfCrf_crf_nnrTransInDecl, 'Enable', 'on');
    set(handles.text282, 'Enable', 'on');
    set(handles.checkbox_vie_glob_trfCrf_crf_onlyNnrOnSources, 'Enable', 'on');
    set(handles.checkbox_vie_glob_trfCrf_crf_fixSomeSources, 'Enable', 'on');
    set(handles.checkbox_vie_glob_trfCrf_crf_sessionWiseReduceSources, 'Enable', 'on');
    set(handles.checkbox_vie_glob_trfCrf_crf_nnrTransInDecl, 'Enable', 'on');
    set(handles.checkbox_vie_glob_trfCrf_crf_nnrTransInDecl, 'Enable', 'on');
    set(handles.checkbox_vie_glob_trfCrf_crf_nnrTransInDecl, 'Enable', 'on');
    if get(handles.checkbox_vie_glob_trfCrf_crf_nnrTransInDecl, 'Value')+...
            get(handles.checkbox_vie_glob_trfCrf_crf_onlyNnrOnSources, 'Value')==0
        set(handles.popupmenu_vie_glob_trfCrf_crf_sourcesForNnr, 'Enable', 'off')
        set(handles.text283, 'Enable', 'off')
    else
        set(handles.popupmenu_vie_glob_trfCrf_crf_sourcesForNnr, 'Enable', 'on')
        set(handles.text283, 'Enable', 'on')
    end
    
    if get(handles.checkbox_vie_glob_trfCrf_crf_fixSomeSources, 'value')
        set(handles.popupmenu_vie_glob_trfCrf_crf_fixedSources, 'Enable', 'on')
        set(handles.text284, 'Enable', 'on')
    else
        set(handles.popupmenu_vie_glob_trfCrf_crf_fixedSources, 'Enable', 'off')
        set(handles.text284, 'Enable', 'off')
    end
    
    if get(handles.checkbox_vie_glob_trfCrf_crf_sessionWiseReduceSources, 'Value')
        set(handles.popupmenu_vie_glob_trfCrf_crf_sessionWiseRedSources, 'Enable', 'on')
        set(handles.text285, 'Enable', 'on')
    else
        set(handles.popupmenu_vie_glob_trfCrf_crf_sessionWiseRedSources, 'Enable', 'off')
        set(handles.text285, 'Enable', 'off')
    end
else
    set(handles.checkbox_vie_glob_trfCrf_crf_nnrTransInDecl, 'Enable', 'off')
    set(handles.text282, 'Enable', 'off')
    set(handles.text282, 'Enable', 'off')
    set(handles.checkbox_vie_glob_trfCrf_crf_onlyNnrOnSources, 'Enable', 'off')
    set(handles.checkbox_vie_glob_trfCrf_crf_fixSomeSources, 'Enable', 'off')
    set(handles.checkbox_vie_glob_trfCrf_crf_sessionWiseReduceSources, 'Enable', 'off')
    set(handles.popupmenu_vie_glob_trfCrf_crf_sourcesForNnr, 'Enable', 'off')
    set(handles.text283, 'Enable', 'off')
    set(handles.popupmenu_vie_glob_trfCrf_crf_fixedSources, 'Enable', 'off')
    set(handles.text284, 'Enable', 'off')
    set(handles.popupmenu_vie_glob_trfCrf_crf_sessionWiseRedSources, 'Enable', 'off')
    set(handles.text285, 'Enable', 'off')
end 



% --- Executes on button press in pushbutton_parameters_iono_create.
function pushbutton_parameters_iono_create_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_parameters_iono_create (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set "from NGS" enabled again, because when we create the
% external iono files, we should update the popupmenu -> this should be
% done when clicking on "external files"
set(handles.radiobutton_parameters_iono_fromNGS, 'Value', 1)
set(handles.popupmenu_parameters_iono_ext, 'Enable', 'off')
set(handles.pushbutton_parameters_iono_create, 'Enable', 'off')

% save automatically 
auto_save_parameterfile(hObject, handles)

% add folder to path
addpath('../ION/PROGRAM/GUI/');

% open GUI to create external ionospheric files
start_createExtIonoFilesGUI


% --- Executes on selection change in popupmenu_parameters_refFrames_superstationTRF.
function popupmenu_parameters_refFrames_superstationTRF_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_refFrames_superstationTRF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_parameters_refFrames_superstationTRF contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_parameters_refFrames_superstationTRF

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu_parameters_refFrames_superstationTRF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_refFrames_superstationTRF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
   

% --- Executes on button press in pushbutton_parameters_refFrames_superstationTRF_chose.
function pushbutton_parameters_refFrames_superstationTRF_chose_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_parameters_refFrames_superstationTRF_chose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get file from explorer
[FileName, PathName] = uigetfile({'*.mat', 'Matlab Binary Format (*.mat)'}, 'Select a superstation file', '../TRF/');

if ~isempty(FileName)
    if ~strcmp('0', num2str(FileName))
        handles=loadSuperstationFile(hObject, handles, [PathName, FileName]);
        handles.data.superstationFile=[PathName, FileName];
        set(handles.text_parameters_refFrames_selected_superstation_file, 'String', handles.data.superstationFile);
        % write message box for information for user
        msgbox('Superstations file has been successfully loaded.', 'Load superstations file', 'help')
    end
end

% save automatically 
auto_save_parameterfile(hObject, handles)

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% save gui parameters (again)
auto_save_parameterfile(hObject, handles)
 
% prepare everything for calling vie_batch
save_runp(hObject, handles)

msgbox('Variable runp saved! Please run vie_batch', 'Saving finished', 'help');



function save_runp(hObject, handles)
% This function saves the runp variable which is needed for running
% vie_batch.

% create textfile containing all infos
% first - get latest auto saved file
newestAutoSavedParameterFile=getNewestAutoSavedParameterFile(hObject, handles);

% load parameter struct
load(['PARAMETERS/', newestAutoSavedParameterFile.name]);

create_input_protocol(parameter);

% SAVE PATHS for init, mod and lsm in any case
if get(handles.checkbox_run_outDirs_diffSubs, 'Value')
	runp.init_path=get(handles.edit_run_outDirs_level0, 'String');
    runp.mod_path=get(handles.edit_run_outDirs_level1, 'String');
    runp.lsm_path=get(handles.edit_run_outDirs_level3, 'String');
    runp.glob_path=get(handles.edit_run_outDirs_level2, 'String');
else
	runp.init_path=get(handles.edit_run_outDirs_oneSub, 'String');
    runp.mod_path=runp.init_path;
    runp.lsm_path=runp.init_path;
    runp.glob_path=runp.init_path;
end

% VIE_SCHED - saving parameters
% ==================================
runp.sched=get(handles.checkbox_run_outputDirectories_runVieSched, 'Value');
if runp.sched
    [vie_sched_error_code] = saveSchedParameters(hObject, handles);
    %If an error occurs => do not start vie_sched!
    if vie_sched_error_code ~= 0
        fprintf(1, 'VIE_SCHED input parameter error! Error code: %g\n', vie_sched_error_code);
        runp.sched=0;
    end
end


% VIEVS
% ================================
curProcessList=get(handles.listbox_setInput_processList, 'String');

% if session was selected -> save (1|0) to run init, mod, lsm.
if isempty(curProcessList)
    % if at least one of init,mod,lsm was selected
    if get(handles.checkbox_run_outputDirectories_runVieInit, 'Value') + ...
            get(handles.checkbox_run_outputDirectories_runVieMod, 'Value')+ ...
            get(handles.checkbox_run_outputDirectories_runVieLsm, 'Value')>0 ...
            && get(handles.checkbox_run_outputDirectories_runVieSched, 'Value') == 0
        msgbox('Select session(s) first! Cannot run vie_init, vie_mod or vie_lsm!', 'No session selected', 'warn');
    end  
else
    if get(handles.checkbox_run_outputDirectories_runVieSched, 'Value') == 1
        msgbox('vie_sched will override your process list!', 'No session selected', 'warn')
    end
    % save parameter file for each session
    % save first parameter file
    findSlashInPl = strfind(curProcessList{1}, '/');
    if isnan(str2double(curProcessList{1}(1))) % if absolute path is taken (first element is character) - use SIM directory!
        runp.ngsdir='SIM';
    end
    saveParamFile(hObject, handles, 'guiparameter.mat');
end

% save parameter and process_list to WORK directory to be loaded from
% vie_batch
% saveParamFile(hObject, handles, 'currentf_pa.mat');
saveProcessList(hObject, handles, 'process_list')

% create runp struct
runp.init = get(handles.checkbox_run_outputDirectories_runVieInit, 'Value');
runp.mod = get(handles.checkbox_run_outputDirectories_runVieMod, 'Value');
runp.lsm = get(handles.checkbox_run_outputDirectories_runVieLsm, 'Value');
runp.lsm_scanwise = get(handles.checkbox_run_outputDirectories_runVieLsmScanwiseUpdate, 'Value');
runp.parallel = get(handles.checkbox_run_runOptions_parallel, 'Value');
runp.error = get(handles.checkbox_run_error_routine, 'Value');
allNCoreEntries = get(handles.popupmenu_run_runOptions_nCores, 'String');
runp.nCores = allNCoreEntries(get(handles.popupmenu_run_runOptions_nCores,...
    'Value'));

   
    
% VIE_SIM - saving parameters (if SIM should be run)
% =================================================
if get(handles.checkbox_run_outputDirectories_runVieSim, 'Value')
    saveSimParameters(hObject, handles)
    runp.sim=1;
    runp.ngsdir='NGS';
end

% does not work better yet


   

% VIE_GLOB - saving parameters (if SIM should be run)
% =================================================

runp.glob=get(handles.checkbox_run_outputDirectories_runVieGlob, 'Value');
if runp.glob
    saveGlobParameters(hObject, handles)
end

save('runp.mat', 'runp'); 


function saveGlobParameters(hObject, handles)
% saves global parameters as given in interface

% save global parameter-paths
items_datumsta=get(handles.popupmenu_vie_glob_trfCrf_trf_stations4Nnt, 'String');
items_posDiscont=get(handles.popupmenu_vie_glob_trfCrf_trf_positDisc, 'String');
items_axisOffset=get(handles.listbox_vie_glob_axisOffsets, 'String');
items_stSeasPos=get(handles.listbox_vie_glob_stseaspos, 'String');
items_APLrg=get(handles.listbox_vie_glob_APLrg, 'String');

pathGS.path_in=get(handles.edit_run_outDirs_glob_pathToLevel2, 'String');
pathGS.L2=get(handles.edit_run_outDirs_glob_level2Sub, 'String');
pathGS.out=get(handles.edit_run_outDirs_glob_outDirVieGlob, 'String');
pathGS.bckwrdsol=get(handles.checkbox_global_param_other_back, 'Value');
if isempty(items_axisOffset) % listbox is empty
    pathGS.ao='';
else
    pathGS.ao=items_axisOffset{get(handles.listbox_vie_glob_axisOffsets, 'Value')};
end

if isempty(items_stSeasPos) % listbox is empty
    pathGS.stseason='';
else
    pathGS.stseason=items_stSeasPos{get(handles.listbox_vie_glob_stseaspos, 'Value')};
end

if isempty(items_APLrg) % listbox is empty
    pathGS.aplrg='';
else
    pathGS.aplrg=items_APLrg{get(handles.listbox_vie_glob_APLrg, 'Value')};
end

pathGS.datumsta=items_datumsta{get(handles.popupmenu_vie_glob_trfCrf_trf_stations4Nnt, 'Value')};
if get(handles.checkbox_vie_glob_trfCrf_trf_reduceStats, 'Value')==0
    pathGS.antred='';
else
    items_antred=get(handles.popupmenu_vie_glob_trfCrf_trf_sessionWiseRedStat, 'String');
    pathGS.antred=items_antred{get(handles.popupmenu_vie_glob_trfCrf_trf_sessionWiseRedStat, 'Value')};
end

pathGS.discont=items_posDiscont{get(handles.popupmenu_vie_glob_trfCrf_trf_positDisc, 'Value')};

if get(handles.checkbox_vie_glob_trfCrf_trf_keepConstVel, 'Value')==0
    pathGS.velconst='';
else
    items_velconst=get(handles.popupmenu_vie_glob_trfCrf_trf_keepConstVel, 'String');
    pathGS.velconst=items_velconst{get(handles.popupmenu_vie_glob_trfCrf_trf_keepConstVel, 'Value')};
end

if get(handles.checkbox_vie_glob_trfCrf_trf_useVelTies, 'Value')==0
    pathGS.velties='';
else
    items_velties=get(handles.popupmenu_vie_glob_trfCrf_trf_velTies, 'String');
    pathGS.velties=items_velties{get(handles.popupmenu_vie_glob_trfCrf_trf_velTies, 'Value')};
end

if (get(handles.checkbox_vie_glob_trfCrf_crf_nnrTransInDecl, 'Value')==0) && ...
       (get(handles.checkbox_vie_glob_trfCrf_crf_onlyNnrOnSources , 'Value')==0) || ...
       (strcmpi(get(handles.checkbox_vie_glob_trfCrf_crf_nnrTransInDecl, 'Enable'), 'Off'))
    pathGS.datumsou='';
else
    items_datumsou=get(handles.popupmenu_vie_glob_trfCrf_crf_sourcesForNnr, 'String');
    pathGS.datumsou=items_datumsou{get(handles.popupmenu_vie_glob_trfCrf_crf_sourcesForNnr, 'Value')};
end

if get(handles.checkbox_vie_glob_trfCrf_crf_fixSomeSources, 'Value')==0
    pathGS.fixedsou='';
else
    items_fixedsou=get(handles.popupmenu_vie_glob_trfCrf_crf_fixedSources, 'String');
    pathGS.fixedsou=items_fixedsou{get(handles.popupmenu_vie_glob_trfCrf_crf_fixedSources, 'Value')};
end

if get(handles.checkbox_vie_glob_trfCrf_crf_sessionWiseReduceSources, 'Value')==0
    pathGS.soured='';
else
    items_soured=get(handles.popupmenu_vie_glob_trfCrf_crf_sessionWiseRedSources, 'String');
    pathGS.soured=items_soured{get(handles.popupmenu_vie_glob_trfCrf_crf_sessionWiseRedSources, 'Value')};
end

pathGS.datumst_scale=num2str(get(handles.radiobutton2_vie_glob_trfCrf_trf_7HelmertParams, 'Value'));
pathGS.datumsou_dz=num2str(get(handles.checkbox_vie_glob_trfCrf_crf_nnrTransInDecl, 'Value'));
pathGS.maxRMS=str2double(get(handles.edit_glob_param_other_maxRMS, 'String'));


% save global parameters
baseName='radiobuttonedit_glob_param_';
parNames={'Clk', 'Zwd', 'TropoGrad', 'AntCoords', 'AntVel', 'SrcCoords', ...
    'Xpol', 'Ypol', 'Dut1', 'DX', 'DY', 'AxisOffsets','StSeasPos','hpole','lpole','APLrg'};
structNames={'pwck' 'rqck' 'zwd' 'ngr' 'egr'   'ant_x' 'ant_y'    'ant_z'  'vx' 'vy'    'vz'   'sra' 'sde' ...
    'xpol'    'ypol'    'dut1' 'dX'    'dY' 'ao' ...
    'stseaspos_Acr' 'stseaspos_Ace' 'stseaspos_Acn' 'stseaspos_Asr' 'stseaspos_Ase' 'stseaspos_Asn' ...
    'hpole' 'lpole' 'aplrg'};

curStructInd=1;

for k=1:length(parNames)
    
    % get 0,1,2
    if eval(['get(handles.', baseName, 'est', parNames{k}, ', ''Value'')'])==1
        par=1;
    elseif eval(['get(handles.', baseName, 'fix', parNames{k}, ', ''Value'')'])==1
        par=0;
    else
        par=2;
    end
    
    parGS(curStructInd).name=structNames{curStructInd};
    parGS(curStructInd).id=par;
    
    if strcmp(parNames{k},'Clk' )
        parGS(curStructInd+1).name=structNames{curStructInd+1};
        parGS(curStructInd+1).id=par;
        curStructInd=curStructInd+1;
    elseif strcmp(parNames{k},'TropoGrad') 
        parGS(curStructInd+1).name=structNames{curStructInd+1};
        parGS(curStructInd+1).id=par; 
        curStructInd=curStructInd+1;
    elseif strcmp(parNames{k},'AntCoords') 
        parGS(curStructInd+1).name=structNames{curStructInd+1};
        parGS(curStructInd+2).name=structNames{curStructInd+2};
        parGS(curStructInd+1).id=par; 
        parGS(curStructInd+2).id=par;
        curStructInd=curStructInd+2;
    elseif strcmp(parNames{k},'AntVel') 
        parGS(curStructInd+1).name=structNames{curStructInd+1};
        parGS(curStructInd+2).name=structNames{curStructInd+2};
        parGS(curStructInd+1).id=par;
        parGS(curStructInd+2).id=par;
        curStructInd=curStructInd+2;
    elseif strcmp(parNames{k},'SrcCoords')
         parGS(curStructInd+1).name=structNames{curStructInd+1};
        parGS(curStructInd+1).id=par; 
         curStructInd=curStructInd+1; 
    elseif strcmp(parNames{k},'StSeasPos' )
        parGS(curStructInd+1).name=structNames{curStructInd+1};
        parGS(curStructInd+2).name=structNames{curStructInd+2};
        parGS(curStructInd+3).name=structNames{curStructInd+3};
        parGS(curStructInd+4).name=structNames{curStructInd+4};
        parGS(curStructInd+5).name=structNames{curStructInd+5};
        
        parGS(curStructInd+1).id=par;
        parGS(curStructInd+2).id=par;
        parGS(curStructInd+3).id=par;
        parGS(curStructInd+4).id=par;
        parGS(curStructInd+5).id=par;
        
        if par==1
            if ~get(handles.checkbox_global_param_stseaspos_R, 'Value')
                parGS(curStructInd).id=0;
                parGS(curStructInd+3).id=0;
            end
            if ~get(handles.checkbox_global_param_stseaspos_E, 'Value')
                parGS(curStructInd+1).id=0;
                parGS(curStructInd+4).id=0;
            end
            if ~get(handles.checkbox_global_param_stseaspos_N, 'Value')
                parGS(curStructInd+2).id=0;
                parGS(curStructInd+5).id=0;
            end
        end
        curStructInd=curStructInd+5;
    end
    curStructInd=curStructInd+1;
end

% tidal ERP variations
parGS(curStructInd).name = 'tidpm';
parGS(curStructInd).id = 0; 
parGS(curStructInd+1).name = 'tidut'; 
parGS(curStructInd+1).id = 0; 
if get(handles.checkbox_global_param_tidERPvar, 'Value')
    parGS(curStructInd).id = 1; 
    parGS(curStructInd+1).id = 1; 
end

% save parameter and paths for global solution and 
save('../DATA/GLOB/parGS.mat', 'parGS')
save('../DATA/GLOB/pathGS.mat', 'pathGS')

% --- Executes on button press in pushbutton_setInput_addPreviousProcList.
function pushbutton_setInput_addPreviousProcList_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_setInput_addPreviousProcList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


proc_list_file='../WORK/process_list.mat';

% if exist -> load it
if exist(proc_list_file, 'file')
    % load
    load(proc_list_file);
    
    % put already selected sessions and last process_list together (but
    % unique)
    newSessions=unique([get(handles.listbox_setInput_processList, 'String');...
        cellstr(process_list)]);
    
    % clear empty cells
    newSessions(cellfun(@isempty, newSessions))=[];
    
    % write sessions to popupmenu
    set(handles.listbox_setInput_processList, 'String', newSessions)
    
    % save to handles struct
    handles.allSelectedFiles=newSessions;
else
    msgbox('Process list does not exist in WORK directory. Probably VieVS was never run before.', 'process_list not found', 'warn');
    
    
end

% Choose default command line output for vie_setup
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in checkbox_parameters_statCorr_hydroLoading.
function checkbox_parameters_statCorr_hydroLoading_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameters_statCorr_hydroLoading (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameters_statCorr_hydroLoading
if get(hObject, 'Value')
    set(handles.popupmenu_parameters_statCorr_hydroLoading, 'Enable', 'on')
else
    set(handles.popupmenu_parameters_statCorr_hydroLoading, 'Enable', 'off')
end


% --- Executes on selection change in popupmenu_parameters_statCorr_hydroLoading.
function popupmenu_parameters_statCorr_hydroLoading_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_statCorr_hydroLoading (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_parameters_statCorr_hydroLoading contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_parameters_statCorr_hydroLoading


% --- Executes during object creation, after setting all properties.
function popupmenu_parameters_statCorr_hydroLoading_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_statCorr_hydroLoading (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menu_plotting_Callback(hObject, eventdata, handles)
% hObject    handle to menu_plotting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_plotting_residuals_Callback(hObject, eventdata, handles)
% hObject    handle to menu_plotting_residuals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% shift the whole window a little bit down because the menu bar appears
curPos=get(handles.figure_vievs2, 'Position');
curPos(2)=curPos(2)-27;
set(handles.figure_vievs2, 'Position', curPos);

% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_plot_residuals, 'Visible', 'On');

% set bar for plotting options to visible
set(handles.uitoolbar_plotOptions, 'Visible', 'On');

% update the folders (since there could have been created antoher subfolder)
dirsInDataFolder=dir('../DATA/LEVEL3/');
dirsInDataFolder(strcmp({dirsInDataFolder.name}, '.')|strcmp({dirsInDataFolder.name}, '..')|~[dirsInDataFolder.isdir])=[];
set(handles.popupmenu_plot_residuals_folder, 'String', ['/', {dirsInDataFolder.name}])

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton53.
function pushbutton53_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton54.
function pushbutton54_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton55.
function pushbutton55_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton56.
function pushbutton56_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit230_Callback(hObject, eventdata, handles)
% hObject    handle to edit230 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit230 as text
%        str2double(get(hObject,'String')) returns contents of edit230 as a double


% --- Executes during object creation, after setting all properties.
function edit230_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit230 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu47.
function popupmenu47_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu47 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu47


% --- Executes during object creation, after setting all properties.
function popupmenu47_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu47 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu48.
function popupmenu48_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu48 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu48


% --- Executes during object creation, after setting all properties.
function popupmenu48_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu49.
function popupmenu49_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu49 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu49


% --- Executes during object creation, after setting all properties.
function popupmenu49_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu50.
function popupmenu50_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu50 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu50


% --- Executes during object creation, after setting all properties.
function popupmenu50_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton57.
function pushbutton57_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu51.
function popupmenu51_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu51 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu51


% --- Executes during object creation, after setting all properties.
function popupmenu51_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu51 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu52.
function popupmenu52_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu52 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu52


% --- Executes during object creation, after setting all properties.
function popupmenu52_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu53.
function popupmenu53_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu53 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu53


% --- Executes during object creation, after setting all properties.
function popupmenu53_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu54.
function popupmenu54_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu54 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu54


% --- Executes during object creation, after setting all properties.
function popupmenu54_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton58.
function pushbutton58_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu55.
function popupmenu55_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu55 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu55


% --- Executes during object creation, after setting all properties.
function popupmenu55_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu56.
function popupmenu56_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu56 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu56


% --- Executes during object creation, after setting all properties.
function popupmenu56_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu56 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu57.
function popupmenu57_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu57 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu57


% --- Executes during object creation, after setting all properties.
function popupmenu57_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu57 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu58.
function popupmenu58_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu58 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu58


% --- Executes during object creation, after setting all properties.
function popupmenu58_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu58 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton59.
function pushbutton59_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton59 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in pushbutton_plot_residuals_clockBreak.
function pushbutton_plot_residuals_clockBreak_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_residuals_clockBreak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if there is a session selected
if ~isempty(get(handles.popupmenu_plot_residuals_session, 'String'))

    % only for station-wise residuals that should be working - so check if
    % that's chosen
    if get(handles.radiobutton_plot_residuals_perStat, 'Value')
        % let user select the mjd from the plot
        [x,y] = ginput(1);
        
        writeClockbreakToOptFile(handles,x);
    else
        msgbox('Selecting clock breaks is only working for station-wise residuals', 'Warning', 'warn');
    end
else
    msgbox('Select session first (click load above)', 'No session selected', 'warn');
end

% --- Executes on selection change in popupmenu_plot_residuals_refClock.
function popupmenu_plot_residuals_refClock_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_residuals_refClock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_residuals_refClock contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_residuals_refClock


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_residuals_refClock_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_residuals_refClock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in popupmenu_plot_residuals_folder.
function popupmenu_plot_residuals_folder_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_residuals_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_residuals_folder contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_residuals_folder


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_residuals_folder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_residuals_folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_plot_residuals_session.
function popupmenu_plot_residuals_session_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_residuals_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_residuals_session contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_residuals_session

% update popupmenus in panel
handles=updatePopupmenusInResidualPlot(hObject, handles);

% plot right away
handles=plotResidualsToAxes(handles);


% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_residuals_session_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_residuals_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function togglebutton_plot_residuals_selectData_Callback(hObject, eventdata, handles)

if get(hObject, 'Value')
    set(handles.axes_plot_residuals,        'ButtonDownFcn', {@startSelectingOutliers,handles})
	set(handles.figure_vievs2,              'WindowButtonUpFcn', {@endSelectingOutliers, handles})

%     set(handles.togglebutton_plot_residuals_SelectedData_output, 'Enable', 'On')
%     set(handles.radiobutton_unit_plot, 'Enable', 'On')
%     set(handles.radiobutton_unit_UTC, 'Enable', 'On')
%     set(handles.radiobutton_unit_MJD, 'Enable', 'On')
%     set(handles.edit_plot_residuals_interval_input_1_unit, 'Enable', 'On')
%     set(handles.edit_plot_residuals_interval_input_2_unit, 'Enable', 'On')
    
    set(handles.pushbutton_plot_residuals_removeOutliers, 'Enable', 'Off')
    if  get(handles.radiobutton_plot_residuals_perAll, 'Value') || get(handles.radiobutton_plot_residuals_perBasel, 'Value') 
        set(handles.pushbutton_selectedData_writeOPT,'Enable','Off');
    else
        set(handles.pushbutton_selectedData_writeOPT,'Enable','On');
    end
else
	set(handles.axes_plot_residuals,                'ButtonDownFcn', '')
    set(handles.figure_vievs2,              'WindowButtonUpFcn', '')
%     set(handles.togglebutton_plot_residuals_SelectedData_output, 'Enable', 'Off')
    set(handles.pushbutton_selectedData_writeOPT,'Enable','Off');
% 	set(handles.radiobutton_unit_plot, 'Enable', 'Off')
%     set(handles.radiobutton_unit_UTC, 'Enable', 'Off')
%     set(handles.radiobutton_unit_MJD, 'Enable', 'Off')
%     set(handles.edit_plot_residuals_interval_input_1_unit, 'Enable', 'Off')
%     set(handles.edit_plot_residuals_interval_input_2_unit, 'Enable', 'Off')
%     set(handles.edit_plot_residuals_interval_show,'String','xx.xx - xx.xx')

% 	set(handles.edit_plot_residuals_interval_input_1_unit, 'String', 'hh')
%     set(handles.edit_plot_residuals_interval_input_2_unit, 'String', 'hh')

    % remove all black object (box and crosses)
    allLineHandles=findobj(handles.axes_plot_residuals,'Type','line', 'color', 'k');
    delete(allLineHandles);
    allLineHandles=findobj(handles.axes_plot_residuals,'Type','line', ...
        'color', [0 0 0.04]);
    delete(allLineHandles);
end


function togglebutton_plot_residuals_SelectedData_output_Callback(hObject, eventdata, handles)

if get(hObject, 'Value')
    if get(handles.radiobutton_unit_MJD,'Value')
        set(handles.radiobutton_unit_MJD,'Value',0);
        set(handles.radiobutton_unit_plot,'Value',1);
    end
    if strcmp(get(handles.edit_plot_residuals_interval_input_1_unit,'String'),'hh') || ...
        strcmp(get(handles.edit_plot_residuals_interval_input_2_unit,'String'),'hh')
            h = msgbox('Please, choose any real number instead of ''hh'' in both windows and switch between Plot Units and UTC to specify dimension for ''hh'' values', 'Warning','warn');
        set(hObject, 'Value',0);
    else
        outlierSelection(hObject, eventdata, handles);
        %set(handles.radiobutton_unit_plot,'Enable','Off');
        %set(handles.radiobutton_unit_UTC,'Enable','Off');
        set(handles.radiobutton_unit_MJD,'Enable','Off');
        set(hObject, 'Value',0);
    end
end
   

function radiobutton_unit_plot_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    set(handles.radiobutton_unit_UTC,'Value',0);
    set(handles.radiobutton_unit_MJD,'Value',0);
    if ~handles.data.plot.foget_unit && handles.data.plot.intervalSelectedatPlot
        set(handles.edit_plot_residuals_interval_show, 'String', ...
            sprintf('%2.3f - %2.3f',handles.data.plot.unit_plot(1),handles.data.plot.unit_plot(2)))
    elseif ~handles.data.plot.foget_unit
        set(handles.edit_plot_residuals_interval_show, 'String', ...
            sprintf(num2str(handles.data.plot.unit_plot)))
    end
else
    set(handles.radiobutton_unit_plot,'Value',1)
end


function radiobutton_unit_UTC_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    set(handles.radiobutton_unit_plot,'Value',0);
    set(handles.radiobutton_unit_MJD,'Value',0);
    if ~handles.data.plot.foget_unit && handles.data.plot.intervalSelectedatPlot
        set(handles.edit_plot_residuals_interval_show, 'String', ...
            sprintf('%2.3f - %2.3f',handles.data.plot.unit_UTC(1),handles.data.plot.unit_UTC(2)))
    elseif ~handles.data.plot.foget_unit
        set(handles.edit_plot_residuals_interval_show, 'String', ...
            sprintf(num2str(handles.data.plot.unit_UTC)))
    end
else
    set(handles.radiobutton_unit_UTC,'Value',1)
end

function radiobutton_unit_MJD_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    set(handles.radiobutton_unit_plot,'Value',0)
    set(handles.radiobutton_unit_UTC,'Value',0)
    if ~handles.data.plot.foget_unit && handles.data.plot.intervalSelectedatPlot
        set(handles.edit_plot_residuals_interval_show, 'String', ...
            sprintf('%2.3f - %2.3f',handles.data.plot.unit_MJD(1),handles.data.plot.unit_MJD(2)))
    elseif ~handles.data.plot.foget_unit
        set(handles.edit_plot_residuals_interval_show, 'String', ...
            sprintf(num2str(handles.data.plot.unit_MJD)))
    end
else
    set(handles.radiobutton_unit_MJD,'Value',1)
end


function pushbutton_selectedData_writeOPT_Callback(hObject, eventdata, handles)

allLineHandles=findobj(handles.axes_plot_residuals,'Type','line', 'color', 'k');
if length(allLineHandles == 2)
    curSession=get(handles.popupmenu_plot_residuals_session, 'Value');
    SessionStartTimeMJD =  handles.data.plot.res(curSession).mjd(1);

    timeSpan = SessionStartTimeMJD+sort([allLineHandles(1).XData(1) allLineHandles(2).XData(1)])/24;
    
    if get(handles.radiobutton_plot_residuals_perStat, 'Value')
        WriteSelectDatatobeExcluded(handles, handles.data.plot.currentStation, timeSpan)
    else
        WriteSelectDatatobeExcluded(handles, handles.data.plot.currentSources, timeSpan)
    end
else
    
end
 

function togglebutton_plot_residuals_selectOutliers_Callback(hObject, eventdata, handles)

% button was selected
if get(hObject, 'Value')
% 	set(handles.radiobutton_unit_plot, 'Enable', 'Off')
%     set(handles.radiobutton_unit_UTC, 'Enable', 'Off')
%     set(handles.radiobutton_unit_MJD, 'Enable', 'Off')
    
    set(handles.axes_plot_residuals,        'ButtonDownFcn', {@startSelectingOutliers,handles})
    set(handles.figure_vievs2,              'WindowButtonUpFcn', {@endSelectingOutliers, handles})
    set(handles.pushbutton_plot_residuals_removeOutliers, 'Enable', 'On')
    

else
    set(handles.axes_plot_residuals,                'ButtonDownFcn', '')
    set(handles.figure_vievs2,              'WindowButtonUpFcn', '')
    set(handles.pushbutton_plot_residuals_removeOutliers, 'Enable', 'Off')
    set(handles.pushbutton_plot_residuals_removeOutliers, 'String', 'Remove Outliers')

    % remove all black object (box and crosses)
    allLineHandles=findobj(handles.axes_plot_residuals,'Type','line', 'color', 'k');
    delete(allLineHandles);
    allLineHandles=findobj(handles.axes_plot_residuals,'Type','line', ...
        'color', [0 0 0.04]);
    delete(allLineHandles);
end


% --- Executes when selected object is changed in uipanel27.
function uipanel27_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel27 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

if strcmp(get(hObject, 'Tag'), 'radiobuttonedit_glob_param_estAntCoords')
    % enable options
    set(handles.radiobutton2_vie_glob_trfCrf_trf_6HelmertParams, 'Enable', 'On')
    set(handles.radiobutton2_vie_glob_trfCrf_trf_7HelmertParams, 'Enable', 'On')
    set(handles.text275, 'Enable', 'On')
    set(handles.popupmenu_vie_glob_trfCrf_trf_stations4Nnt, 'Enable', 'On')
    set(handles.checkbox_vie_glob_trfCrf_trf_reduceStats, 'Enable', 'On')
    if get(handles.checkbox_vie_glob_trfCrf_trf_reduceStats, 'Value')
        set(handles.popupmenu_vie_glob_trfCrf_trf_sessionWiseRedStat, 'Enable', 'On')
        set(handles.text278, 'Enable', 'On')
    end
    set(handles.text277, 'Enable', 'On')
    set(handles.checkbox_vie_glob_trfCrf_trf_keepConstVel, 'Enable', 'On')
    if get(handles.checkbox_vie_glob_trfCrf_trf_keepConstVel, 'Value')
        set(handles.popupmenu_vie_glob_trfCrf_trf_keepConstVel, 'Enable', 'on')
        set(handles.text281, 'Enable', 'on')
    end
    set(handles.checkbox_vie_glob_trfCrf_trf_useVelTies, 'Enable', 'On')
    if get(handles.checkbox_vie_glob_trfCrf_trf_useVelTies, 'Value')
        set(handles.popupmenu_vie_glob_trfCrf_trf_velTies, 'Enable', 'On')
        set(handles.text280, 'Enable', 'On')
    end
    set(handles.text299, 'Enable', 'On')
    
    
else % disable the options in the second interface
    set(handles.radiobutton2_vie_glob_trfCrf_trf_6HelmertParams, 'Enable', 'Off')
    set(handles.radiobutton2_vie_glob_trfCrf_trf_7HelmertParams, 'Enable', 'Off')
    set(handles.text275, 'Enable', 'Off')
    set(handles.popupmenu_vie_glob_trfCrf_trf_stations4Nnt, 'Enable', 'Off')
    set(handles.checkbox_vie_glob_trfCrf_trf_reduceStats, 'Enable', 'Off')
    set(handles.popupmenu_vie_glob_trfCrf_trf_sessionWiseRedStat, 'Enable', 'Off')
    set(handles.text278, 'Enable', 'Off')
    set(handles.text277, 'Enable', 'Off')
    set(handles.checkbox_vie_glob_trfCrf_trf_keepConstVel, 'Enable', 'Off')
    set(handles.popupmenu_vie_glob_trfCrf_trf_keepConstVel, 'Enable', 'Off')
    set(handles.text281, 'Enable', 'Off')
    set(handles.checkbox_vie_glob_trfCrf_trf_useVelTies, 'Enable', 'Off')
    set(handles.popupmenu_vie_glob_trfCrf_trf_velTies, 'Enable', 'Off')
    set(handles.text280, 'Enable', 'Off')
end
    


% --- Executes on button press in pushbutton_plot_residuals_load.
function pushbutton_plot_residuals_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_residuals_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=loadResFiles(hObject, handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes when selected object is changed in uipanel_plot_residuals_plotPer.
function uipanel_plot_residuals_plotPer_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_plot_residuals_plotPer 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% if "station is ticked
if strcmp(get(hObject, 'Tag'), 'radiobutton_plot_residuals_perStat')
    set(handles.popupmenu_plot_residuals_station, 'Enable', 'On')
    set(handles.popupmenu_plot_residuals_baseline, 'Enable', 'Off')
    set(handles.popupmenu_plot_residuals_source, 'Enable', 'Off')
    
elseif strcmp(get(hObject, 'Tag'), 'radiobutton_plot_residuals_perBasel')
    set(handles.popupmenu_plot_residuals_baseline, 'Enable', 'On')
    set(handles.popupmenu_plot_residuals_station, 'Enable', 'Off')
    set(handles.popupmenu_plot_residuals_source, 'Enable', 'Off')
elseif strcmp(get(hObject, 'Tag'), 'radiobutton_plot_residuals_perSource')
    set(handles.popupmenu_plot_residuals_source, 'Enable', 'On')
    set(handles.popupmenu_plot_residuals_baseline, 'Enable', 'Off')
    set(handles.popupmenu_plot_residuals_station, 'Enable', 'Off')
else
    set(handles.popupmenu_plot_residuals_baseline, 'Enable', 'Off')
    set(handles.popupmenu_plot_residuals_station, 'Enable', 'Off')
    set(handles.popupmenu_plot_residuals_source, 'Enable', 'Off')
end

handles=plotResidualsToAxes(handles);

% Update handles structure
guidata(hObject, handles);


function removeBlackObjectsOnAxes(axes)
% This function removes black Objects (lines, markers) from the axes given
% as input
allBlackObjects=findobj(handles.axes_plot_residuals,'Type','line', 'color', 'k');
    delete(allBlackObjects);
    

% --- Executes on selection change in popupmenu_plot_residuals_station.
function popupmenu_plot_residuals_station_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_residuals_station (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_residuals_station contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_residuals_station

% plot residuals to axes
handles=plotResidualsToAxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_residuals_station_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_residuals_station (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_plot_residuals_baseline.
function popupmenu_plot_residuals_baseline_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_residuals_baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_residuals_baseline contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_residuals_baseline

% plot residuals to axes
handles=plotResidualsToAxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_residuals_baseline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_residuals_baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel_plot_residuals_firstOrMainSolution.
function uipanel_plot_residuals_firstOrMainSolution_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_plot_residuals_firstOrMainSolution 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% "clear" outlier selection
% set(handles.axes_plot_residuals,                'ButtonDownFcn', '')
% set(handles.figure_vievs2,              'WindowButtonUpFcn', '')
set(handles.pushbutton_plot_residuals_removeOutliers, 'Enable', 'Off')
set(handles.pushbutton_plot_residuals_removeOutliers, 'String', 'Remove Outliers')

% remove all black object (box and crosses)
allLineHandles=findobj(handles.axes_plot_residuals,'Type','line', 'color', 'k');
delete(allLineHandles);


handles=plotResidualsToAxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in checkbox_run_sinex_sources.
function checkbox_run_sinex_sources_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_sinex_sources (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_sinex_sources

if get(hObject, 'Value')
    set(handles.text_run_sinex_sources, 'Enable', 'On')
    set(handles.radiobutton_run_sinex_sources_incl, 'Enable', 'On')
else
    set(handles.text_run_sinex_sources, 'Enable', 'Off')
    set(handles.radiobutton_run_sinex_sources_incl, 'Enable', 'Off')
end


% --- Executes on button press in checkbox_vie_sched_ngs.
function checkbox_vie_sched_ngs_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sched_ngs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sched_ngs


% --- Executes on button press in checkbox_vie_sched_skd.
function checkbox_vie_sched_skd_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sched_skd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sched_skd


% --- Executes when uipanel171 is resized.
function uipanel171_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel171 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_vie_sched_clearStations.
function pushbutton_vie_sched_clearStations_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_vie_sched_clearStations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get current content of listbox
curContent=get(handles.listbox_vie_sched_stanet, 'String');

if ~isempty(curContent)
    % delete all selected entries entry
    
    curContent(get(handles.listbox_vie_sched_stanet, 'Value'))=[];
    
    % get current value
    curValue=get(handles.listbox_vie_sched_stanet, 'Value');
    
    set(handles.listbox_vie_sched_stanet, 'Value', ...
        min([max(curValue), max(length(curContent),1)]));

    % update listbox
    set(handles.listbox_vie_sched_stanet, 'String', curContent);

    % and also save it to handles struct
    % handles.allSelectedFiles=curContent;

    % save changes to handles struct
    guidata(hObject, handles);

    % save parameter file automatically 
    auto_save_parameterfile(hObject, handles)
end


% --- Executes on selection change in popupmenu_plot_residuals_source.
function popupmenu_plot_residuals_source_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_residuals_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_residuals_source contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_residuals_source

% plot residuals to axes
handles=plotResidualsToAxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_residuals_source_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_residuals_source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togglebutton5.
function togglebutton5_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton5



% --------------------------------------------------------------------
function menu_scheduling_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to menu_scheduling_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_vie_sched, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);



% --------------------------------------------------------------------
function menu_scheduling_Callback(hObject, eventdata, handles)
% hObject    handle to menu_scheduling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_run_Callback(hObject, eventdata, handles)
% hObject    handle to menu_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_simulation_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to menu_simulation_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_vie_sim, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);




% --- Executes on selection change in popupmenu_plot_sessionAnalysis_subfolder.
function popupmenu_plot_sessionAnalysis_subfolder_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_sessionAnalysis_subfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_sessionAnalysis_subfolder contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_sessionAnalysis_subfolder


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_sessionAnalysis_subfolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_sessionAnalysis_subfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_plot_sessionAnalysis_session.
function popupmenu_plot_sessionAnalysis_session_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_sessionAnalysis_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_sessionAnalysis_session contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_sessionAnalysis_session

% update popupmenus in panel
handles=updatePopupmenusInSessionAnalysisPlot(hObject, handles);

% plot right away
handles=plotSessionAnalysisToAxes(handles);
    
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_sessionAnalysis_session_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_sessionAnalysis_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_plot_sessionAnalysis_load.
function pushbutton_plot_sessionAnalysis_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_sessionAnalysis_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load
handles=loadSessionAnalysisData(hObject, handles);

% update popupmenus in panel
handles=updatePopupmenusInSessionAnalysisPlot(hObject, handles);

% plot right away
handles=plotSessionAnalysisToAxes(handles);
    
guidata(hObject, handles);
    
% i need following files: 
% x_ (level3)      ... baseline lnegth rep
% antenna (level3) ... map network
% atpa (level3)    ... corrmat


% --- Executes on selection change in popupmenu_plot_sessionAnalysis_corrMat_firstPar.
function popupmenu_plot_sessionAnalysis_corrMat_firstPar_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_sessionAnalysis_corrMat_firstPar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_sessionAnalysis_corrMat_firstPar contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_sessionAnalysis_corrMat_firstPar


handles=plotSessionAnalysisToAxes(handles);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_sessionAnalysis_corrMat_firstPar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_sessionAnalysis_corrMat_firstPar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_plot_sessionAnalysis_corrMat_secondPar.
function popupmenu_plot_sessionAnalysis_corrMat_secondPar_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_sessionAnalysis_corrMat_secondPar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_sessionAnalysis_corrMat_secondPar contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_sessionAnalysis_corrMat_secondPar


handles=plotSessionAnalysisToAxes(handles);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_sessionAnalysis_corrMat_secondPar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_sessionAnalysis_corrMat_secondPar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menu_plotting_sessionAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to menu_plotting_sessionAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% shift the whole window a little bit down because the menu bar appears
curPos=get(handles.figure_vievs2, 'Position');
curPos(2)=curPos(2)-27;
set(handles.figure_vievs2, 'Position', curPos);

% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_plot_sessionAnalysis, 'Visible', 'On');

% set bar for plotting options to visible
set(handles.uitoolbar_plotOptions, 'Visible', 'On');

% update the folders (since there could have been created antoher subfolder)
dirsInDataFolder=dir('../DATA/LEVEL3/');
dirsInDataFolder(strcmp({dirsInDataFolder.name}, '.')|strcmp({dirsInDataFolder.name}, '..')|~[dirsInDataFolder.isdir])=[];
set(handles.popupmenu_plot_sessionAnalysis_subfolder, 'String', ['/', {dirsInDataFolder.name}])

% Update handles structure
guidata(hObject, handles);
    
            
% --- Executes when selected object is changed in uipanel_plot_sessionAnalysis_options.
function uipanel_plot_sessionAnalysis_options_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_plot_sessionAnalysis_options 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% session network
if get(handles.radiobutton_plot_sessionAnalysis_network, 'Value')
    set(handles.popupmenu_plot_sessionAnalysis_corrMat_firstPar, 'Enable', 'Off')
    set(handles.popupmenu_plot_sessionAnalysis_corrMat_secondPar, 'Enable', 'Off') 
    set(handles.checkbox_plot_sessionAnalysis_baselineNames, 'Enable', 'Off')

% baseline length repeatability  
elseif get(handles.radiobutton_plot_sessionAnalysis_corMatrix, 'Value')
    set(handles.popupmenu_plot_sessionAnalysis_corrMat_firstPar, 'Enable', 'On')
    set(handles.popupmenu_plot_sessionAnalysis_corrMat_secondPar, 'Enable', 'On')
    
    set(handles.checkbox_plot_sessionAnalysis_baselineNames, 'Enable', 'Off')

% correlation matrix
elseif get(handles.radiobutton_plot_sessionAnalysis_baselLeRep, 'Value')
    set(handles.checkbox_plot_sessionAnalysis_baselineNames, 'Enable', 'On')
    
    set(handles.popupmenu_plot_sessionAnalysis_corrMat_firstPar, 'Enable', 'Off')
    set(handles.popupmenu_plot_sessionAnalysis_corrMat_secondPar, 'Enable', 'Off')
    

end


handles=plotSessionAnalysisToAxes(handles);

guidata(hObject, handles);



% --- Executes on button press in checkbox_plot_sessionAnalysis_baselineNames.
function checkbox_plot_sessionAnalysis_baselineNames_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_sessionAnalysis_baselineNames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_sessionAnalysis_baselineNames


handles=plotSessionAnalysisToAxes(handles);

guidata(hObject, handles);


% --- Executes on selection change in popupmenu73.
function popupmenu73_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu73 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu73 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu73


% --- Executes during object creation, after setting all properties.
function popupmenu73_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu73 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu74.
function popupmenu74_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu74 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu74 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu74


% --- Executes during object creation, after setting all properties.
function popupmenu74_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu74 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox263.
function checkbox263_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox263 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox263


% --- Executes on selection change in popupmenu75.
function popupmenu75_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu75 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu75 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu75


% --- Executes during object creation, after setting all properties.
function popupmenu75_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu75 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu76.
function popupmenu76_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu76 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu76 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu76


% --- Executes during object creation, after setting all properties.
function popupmenu76_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu76 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton73.
function pushbutton73_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton73 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_plotting_scheduling_Callback(hObject, eventdata, handles)
% hObject    handle to menu_plotting_scheduling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sched_analyser();



% --- Executes on button press in pushbutton_parameters_refFrames_createSupersourceFile.
function pushbutton_parameters_refFrames_createSupersourceFile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_parameters_refFrames_createSupersourceFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

createSupersourcesFile


% --- Executes on button press in pushbutton_parameters_refFrames_superstationTRF_create.
function pushbutton_parameters_refFrames_superstationTRF_create_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_parameters_refFrames_superstationTRF_create (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

createSuperstationsFile



% --- Executes on selection change in popupmenu_parameters_refFrames_supersourceCRF.
function popupmenu_parameters_refFrames_supersourceCRF_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_refFrames_supersourceCRF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_parameters_refFrames_supersourceCRF contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_parameters_refFrames_supersourceCRF


% --- Executes during object creation, after setting all properties.
function popupmenu_parameters_refFrames_supersourceCRF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_refFrames_supersourceCRF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_parameters_refFrames_superstationCRF_chose.
function pushbutton_parameters_refFrames_superstationCRF_chose_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_parameters_refFrames_superstationCRF_chose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% get file from explorer
[FileName, PathName] = uigetfile({'*.mat', 'Matlab Binary Format (*.mat)'}, 'Select a supersource file', '../CRF/');

if ~isempty(FileName)
    if ~strcmp('0', num2str(FileName))
        handles=loadSupersourceFile(hObject, handles, [PathName, FileName]);
        handles.data.supersourceFile=[PathName, FileName];
        set(handles.text_parameters_refFrames_selected_supersource_file, 'String', handles.data.supersourceFile);
        % write message box for information for user
        msgbox('Supersource file has been successfully loaded.', 'Load supersource file', 'help')
    end
end

% save automatically 
auto_save_parameterfile(hObject, handles)

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in checkbox_plot_residuals_showStatNumbers.
function checkbox_plot_residuals_showStatNumbers_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_residuals_showStatNumbers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_residuals_showStatNumbers

handles=plotResidualsToAxes(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in checkbox_plot_residuals_showSourceNames.
function checkbox_plot_residuals_showSourceNames_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_residuals_showSourceNames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_residuals_showSourceNames

handles=plotResidualsToAxes(handles);

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in checkbox_run_outputDirectories_runVieLsmScanwiseUpdate.
function checkbox_run_outputDirectories_runVieLsmScanwiseUpdate_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_outputDirectories_runVieLsmScanwiseUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_outputDirectories_runVieLsmScanwiseUpdate

if get(hObject, 'Value')
	% if checkbox is ticked -> untick normal LSM
    set(handles.checkbox_run_outputDirectories_runVieLsm, 'Value', 0);
	
    % if also Vie_GLOB is ticked
	if get(handles.checkbox_run_outputDirectories_runVieGlob, 'value')
    
		% Use different output Sub-directories
		if get(handles.checkbox_run_outDirs_diffSubs, 'Value')
			str_lsm_subDir=get(handles.edit_run_outDirs_level2, 'String');
		% Use one output Sub-directory	
		else
			str_lsm_subDir=get(handles.edit_run_outDirs_oneSub, 'String');
		end
		
		%str_glob_subDir=get(handles.edit_run_outDirs_glob_level2Sub, 'String');
		
		set(handles.edit_run_outDirs_glob_pathToLevel2, 'Enable', 'off');
		set(handles.edit_run_outDirs_glob_pathToLevel2, 'String', '../DATA/LEVEL2/');
		set(handles.edit_run_outDirs_glob_level2Sub, 'Enable', 'off');
		set(handles.edit_run_outDirs_glob_level2Sub, 'String', str_lsm_subDir);
		
		msgbox('If Vie_LSM and Vie_GLOB are both run, the Vie_LSM output directory for N-matrices have to be equal to the Vie_GLOB input directory!', 'Warning (subdirectories not equal)');
		
	end
	
else
	set(handles.edit_run_outDirs_glob_pathToLevel2, 'Enable', 'on');
	set(handles.edit_run_outDirs_glob_level2Sub, 'Enable', 'on');
end





% --- Executes on button press in checkbox_vie_sched_sum.
function checkbox_vie_sched_sum_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sched_sum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sched_sum



% --- Executes on button press in radiobutton_parameters_troposphere_indModeling.
function radiobutton_parameters_troposphere_indModeling_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_parameters_troposphere_indModeling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_parameters_troposphere_indModeling

% set en/disable
set(handles.radiobutton_parameters_troposphere_zhd_no, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_zhd_fromInSitu, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_zhd_VMF3, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_zhd_VMF1, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_zhd_GPT3, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_zwd_no, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_zwd_fromInSitu, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_zwd_VMF3, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_zwd_VMF1, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_zwd_GPT3, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_mfh_VMF3, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_mfh_VMF1, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_mfh_GPT3, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_mfw_VMF3, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_mfw_VMF1, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_mfw_GPT3, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_gradients_h_no, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_gradients_h_GRAD, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_gradients_h_GPT3, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_gradients_h_DAO, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_gradients_w_no, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_gradients_w_GRAD, 'Enable', 'on')
set(handles.radiobutton_parameters_troposphere_gradients_w_GPT3, 'Enable', 'on')

auto_save_parameterfile(hObject, handles)


% --- Executes on button press in radiobutton_parameters_troposphere_raytr.
function radiobutton_parameters_troposphere_raytr_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_parameters_troposphere_raytr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_parameters_troposphere_raytr

% set en/disable
set(handles.radiobutton_parameters_troposphere_zhd_no, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_zhd_fromInSitu, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_zhd_VMF3, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_zhd_VMF1, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_zhd_GPT3, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_zwd_no, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_zwd_fromInSitu, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_zwd_VMF3, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_zwd_VMF1, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_zwd_GPT3, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_mfh_VMF3, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_mfh_VMF1, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_mfh_GPT3, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_mfw_VMF3, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_mfw_VMF1, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_mfw_GPT3, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_gradients_h_no, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_gradients_h_GRAD, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_gradients_h_GPT3, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_gradients_h_DAO, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_gradients_w_no, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_gradients_w_GRAD, 'Enable', 'off')
set(handles.radiobutton_parameters_troposphere_gradients_w_GPT3, 'Enable', 'off')

auto_save_parameterfile(hObject, handles)


% --- Executes on button press in radiobutton_parameters_troposphere_gradients_h_DAO.
function radiobutton_parameters_troposphere_gradients_h_DAO_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_parameters_troposphere_gradients_h_DAO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_parameters_troposphere_gradients_h_DAO

% set en/disable
set(handles.radiobutton_parameters_troposphere_gradients_w_no, 'Enable', 'Off')
set(handles.radiobutton_parameters_troposphere_gradients_w_GRAD, 'Enable', 'Off')
set(handles.radiobutton_parameters_troposphere_gradients_w_GPT3, 'Enable', 'Off')

% --- Executes on button press in radiobutton_parameters_eop_interp_lin.
function radiobutton_parameters_eop_interp_lin_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_parameters_eop_interp_lin 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_parameters_eop_interp_lin

% set en/disable
set(handles.checkbox_parameters_eop_interp_lin48h, 'Enable', 'On')


% --- Executes on selection change in popupmenu_plot_sessionAnalysis_subfolder2.
function popupmenu_plot_sessionAnalysis_subfolder2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_sessionAnalysis_subfolder2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_sessionAnalysis_subfolder2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_sessionAnalysis_subfolder2


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_sessionAnalysis_subfolder2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_sessionAnalysis_subfolder2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_plot_sessionAnalysis_load2.
function pushbutton_plot_sessionAnalysis_load2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_sessionAnalysis_load2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load
handles=loadSessionAnalysisData(hObject, handles);

% plot right away
handles=plotSessionAnalysisToAxes(handles);
    
guidata(hObject, handles);


% --- Executes on selection change in popupmenu_plot_sessionAnalysis_session2.
function popupmenu_plot_sessionAnalysis_session2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_sessionAnalysis_session2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_sessionAnalysis_session2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_sessionAnalysis_session2

% plot right away
handles=plotSessionAnalysisToAxes(handles);
    
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_sessionAnalysis_session2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_sessionAnalysis_session2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_plot_sessionAnalysis_add2.
function checkbox_plot_sessionAnalysis_add2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_sessionAnalysis_add2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_sessionAnalysis_add2

% plot right away
handles=plotSessionAnalysisToAxes(handles);
    
guidata(hObject, handles);


% --- Executes on selection change in popupmenu_plot_sessionAnalysis_subfolder4.
function popupmenu_plot_sessionAnalysis_subfolder4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_sessionAnalysis_subfolder4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_sessionAnalysis_subfolder4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_sessionAnalysis_subfolder4


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_sessionAnalysis_subfolder4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_sessionAnalysis_subfolder4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_plot_sessionAnalysis_load4.
function pushbutton_plot_sessionAnalysis_load4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_sessionAnalysis_load4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load
handles=loadSessionAnalysisData(hObject, handles);

% plot right away
handles=plotSessionAnalysisToAxes(handles);
    
guidata(hObject, handles);


% --- Executes on selection change in popupmenu_plot_sessionAnalysis_session4.
function popupmenu_plot_sessionAnalysis_session4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_sessionAnalysis_session4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_sessionAnalysis_session4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_sessionAnalysis_session4

% plot right away
handles=plotSessionAnalysisToAxes(handles);
    
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_sessionAnalysis_session4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_sessionAnalysis_session4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_plot_sessionAnalysis_add4.
function checkbox_plot_sessionAnalysis_add4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_sessionAnalysis_add4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_sessionAnalysis_add4

% plot right away
handles=plotSessionAnalysisToAxes(handles);
    
guidata(hObject, handles);


% --- Executes on selection change in popupmenu_plot_sessionAnalysis_subfolder3.
function popupmenu_plot_sessionAnalysis_subfolder3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_sessionAnalysis_subfolder3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_sessionAnalysis_subfolder3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_sessionAnalysis_subfolder3


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_sessionAnalysis_subfolder3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_sessionAnalysis_subfolder3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_plot_sessionAnalysis_load3.
function pushbutton_plot_sessionAnalysis_load3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_sessionAnalysis_load3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load
handles=loadSessionAnalysisData(hObject, handles);

% plot right away
handles=plotSessionAnalysisToAxes(handles);
    
guidata(hObject, handles);


% --- Executes on selection change in popupmenu_plot_sessionAnalysis_session3.
function popupmenu_plot_sessionAnalysis_session3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_sessionAnalysis_session3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_sessionAnalysis_session3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_sessionAnalysis_session3

% plot right away
handles=plotSessionAnalysisToAxes(handles);
    
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_plot_sessionAnalysis_session3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_sessionAnalysis_session3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_plot_sessionAnalysis_add3.
function checkbox_plot_sessionAnalysis_add3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_sessionAnalysis_add3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_sessionAnalysis_add3

% plot right away
handles=plotSessionAnalysisToAxes(handles);
    
guidata(hObject, handles);


% --- Executes on button press in checkbox_run_runOptions_parallel.
function checkbox_run_runOptions_parallel_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_runOptions_parallel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_runOptions_parallel

if get(hObject, 'Value')
    newState='On';
else
    newState='Off';
end
set(handles.text328, 'Enable', newState);
set(handles.popupmenu_run_runOptions_nCores, 'Enable', newState);
set(handles.text329, 'Enable', newState);

% --- Executes on selection change in popupmenu_run_runOptions_nCores.
function popupmenu_run_runOptions_nCores_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_run_runOptions_nCores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_run_runOptions_nCores contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_run_runOptions_nCores


% --- Executes during object creation, after setting all properties.
function popupmenu_run_runOptions_nCores_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_run_runOptions_nCores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton84.
function pushbutton84_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton84 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

object_handles = findall(handles.figure_vievs2);
nObj=length(object_handles);
for iObj=1:nObj
    % available unit?
    unitAv=1;
    try get(object_handles(iObj),'Units');
    catch
        fprintf('No unit for ''%s''\n',get(object_handles(iObj), 'Tag'));
        unitAv=0;
    end
    if unitAv==1
        % do we have to delete object (one seems to be not needed and
        % tag=''
        deletedObj=0;
        try
            oldUnit=get(handles.(get(object_handles(iObj),'Tag')), 'Units');
        catch
            delete(object_handles(iObj));
            deletedObj=1;
        end
        if deletedObj==0
            set(handles.(get(object_handles(iObj),'Tag')), 'Units', 'pixels');
            fprintf(' Object %3.0f/%3.0f done (%s --> %s) (%-s)\n', iObj,nObj,...
                oldUnit,get(handles.(get(object_handles(iObj),'Tag')), 'Units'),...
                get(object_handles(iObj), 'Tag'));
        end
    end
    
end

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
fprintf('DONE!\n');


% --- Executes on button press in checkbox_run_globalPram_axisOffset.
function checkbox_run_globalPram_axisOffset_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_globalPram_axisOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_globalPram_axisOffset


% --- Executes on selection change in listbox_vie_glob_axisOffsets.
function listbox_vie_glob_axisOffsets_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_vie_glob_axisOffsets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_vie_glob_axisOffsets contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_vie_glob_axisOffsets


% --- Executes during object creation, after setting all properties.
function listbox_vie_glob_axisOffsets_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_vie_glob_axisOffsets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel260.
function uipanel260_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel260 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% fr
% change enable of listboxom estimate to fix (or otherwise)
if get(handles.radiobuttonedit_glob_param_estAxisOffsets, 'Value')
    % estimate was selected
    state='On';
else
    state='Off';
end

set(handles.listbox_vie_glob_axisOffsets, 'Enable', state);



function edit_parameter_obsRestr_jetang_Callback(hObject, eventdata, handles)
% hObject    handle to edit_parameter_obsRestr_jetang (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_parameter_obsRestr_jetang as text
%        str2double(get(hObject,'String')) returns contents of edit_parameter_obsRestr_jetang as a double


% --- Executes during object creation, after setting all properties.
function edit_parameter_obsRestr_jetang_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_parameter_obsRestr_jetang (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_parameters_ss_catalog.
function popupmenu_parameters_ss_catalog_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_ss_catalog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_parameters_ss_catalog contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_parameters_ss_catalog


% --- Executes during object creation, after setting all properties.
function popupmenu_parameters_ss_catalog_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_ss_catalog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_parameters_ss_applyss.
function checkbox_parameters_ss_applyss_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameters_ss_applyss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameters_ss_applyss


% --- Executes on button press in checkbox_parameters_ss_writejet.
function checkbox_parameters_ss_writejet_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameters_ss_writejet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameters_ss_writejet


% --- Executes on selection change in listbox_vie_sim_ss.
function listbox_vie_sim_ss_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_vie_sim_ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_vie_sim_ss contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_vie_sim_ss


% --- Executes during object creation, after setting all properties.
function listbox_vie_sim_ss_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_vie_sim_ss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_vie_sim_simSS.
function checkbox_vie_sim_simSS_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sim_simSS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sim_simSS
if get(handles.checkbox_vie_sim_simSS,'Value')==1
    status='on';
else
    status='off';
end 
    set(handles.text332,'Enable',status);
    set(handles.listbox_vie_sim_ss,'Visible',status);


% --------------------------------------------------------------------
function menu_parameters_sourcestructure_Callback(hObject, eventdata, handles)
% hObject    handle to menu_parameters_sourcestructure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setAllPanelsToInvisible(hObject,handles);
set(handles.uipanel_parameters_sourcestructure,'Visible','on');

% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in uipanel_plot_eopOut_write.
function uipanel_plot_eopOut_write_Callback(hObject, eventdata, handles)
% hObject    handle to uipanel_plot_eopOut_write (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out_basEop(hObject, handles);

function edit_plot_eopOut_outFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_plot_eopOut_outFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_plot_eopOut_outFile as text
%        str2double(get(hObject,'String')) returns contents of edit_plot_eopOut_outFile as a double


% --- Executes during object creation, after setting all properties.
function edit_plot_eopOut_outFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_plot_eopOut_outFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_plot_eopOut_outFile_getFile.
function pushbutton_plot_eopOut_outFile_getFile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_eopOut_outFile_getFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% get a file
[filename, pathname, filterindex] = uiputfile( ...
{'*.txt', 'Text file (*.txt)'},...
 'Save as',...
 '../OUT/');

% see if something was chosen
if ~isempty(filename)
    if ~strcmp('0', num2str(filename))
        % write new filename to edit textbox
        set(handles.edit_plot_eopOut_outFile, 'String', [pathname, filename]);
    end
end

% --- Executes on selection change in popupmenu_plot_eopOut_subfolder.
function popupmenu_plot_eopOut_subfolder_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_eopOut_subfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_eopOut_subfolder contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_eopOut_subfolder


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_eopOut_subfolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_eopOut_subfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_plot_eopOut_pl.
function popupmenu_plot_eopOut_pl_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_eopOut_pl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_eopOut_pl contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_eopOut_pl


% --- Executes during object creation, after setting all properties.
function popupmenu_plot_eopOut_pl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_eopOut_pl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel_plot_eopOut_pl.
function uipanel_plot_eopOut_pl_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_plot_eopOut_pl 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if strcmpi(get(hObject, 'Tag'), 'radiobutton_plot_eopOut_pl_current') || ...
        strcmpi(get(hObject, 'Tag'), 'radiobutton_use_subfolder_data')
	newState='Off';
else
	newState='On';
end
set(handles.popupmenu_plot_eopOut_pl, 'Enable', newState);


% --------------------------------------------------------------------
function menu_plotting_eopOut_Callback(hObject, eventdata, handles)
% hObject    handle to menu_plotting_eopOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_plot_eopOut, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);


% --- Executes when selected object is changed in uipanel_plot_eopOut_subfolder.
function uipanel_plot_eopOut_subfolder_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_plot_eopOut_subfolder 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if strcmpi(get(hObject, 'Tag'), 'radiobutton_plot_eopOut_subfolder_current')
	newState='Off';
else
	newState='On';
end
set(handles.popupmenu_plot_eopOut_subfolder, 'Enable', newState);


% --- Executes when selected object is changed in uipanel_plot_eopOut_outFile.
function uipanel_plot_eopOut_outFile_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_plot_eopOut_outFile 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
if strcmpi(get(hObject, 'Tag'), 'radiobutton_plot_eopOut_outFile_default')
	newState='Off';
else
	newState='On';
end
set(handles.edit_plot_eopOut_outFile, 'Enable', newState);
set(handles.pushbutton_plot_eopOut_outFile_getFile, 'Enable', newState);


% --------------------------------------------------------------------
function menu_plotting_basOut_Callback(hObject, eventdata, handles)
% hObject    handle to menu_plotting_basOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in uipanel_plot_eopOut_writeInt.
function uipanel_plot_eopOut_writeInt_Callback(hObject, eventdata, handles)
% hObject    handle to uipanel_plot_eopOut_writeInt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out_basEop(hObject, handles);


% --- Executes on button press in uipanel_plot_basOut_write.
function uipanel_plot_basOut_write_Callback(hObject, eventdata, handles)
% hObject    handle to uipanel_plot_basOut_write (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

out_basEop(hObject, handles);


% --- Executes on button press in pushbutton_plot_residuals_clockRef_get.
function pushbutton_plot_residuals_clockRef_get_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_residuals_clockRef_get (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

OptFileExist=0;

% read OPT file if exist
% get selected file
allSessionsInList=get(handles.popupmenu_plot_residuals_session, 'String');
session=allSessionsInList{get(handles.popupmenu_plot_residuals_session, 'Value')};
if str2double(session(1:2))>75
    yearStr=['19', session(1:2)];
else
    yearStr=['20', session(1:2)];
end

% get currently selected OPT directory
OPTdirs=get(handles.popupmenu_setInput_optDir, 'String');
selectedOPTdir=OPTdirs{get(handles.popupmenu_setInput_optDir, 'value')};

% check if file exist
wantedOPTfile=['../../VLBI_OPT/', selectedOPTdir, '/', yearStr, '/', session(1:9), '.OPT'];
% if the year folder does not exist - create it
if exist(wantedOPTfile, 'file')
    OptFileExist=1;
end

if OptFileExist
    [ini_opt, bas_excl]=readOPT(wantedOPTfile);
    if isempty(ini_opt.refclock)
        newValPopup=1;
    else
        allAntennas=get(handles.popupmenu_plot_residuals_refClock, 'String');
        newValPopup=find(strcmpi(allAntennas,ini_opt.refclock));
    end
else
    newValPopup=1;
end


set(handles.popupmenu_plot_residuals_refClock, 'Value', newValPopup);


% --- Executes on button press in pushbutton_plot_residuals_clockRef_set.
function pushbutton_plot_residuals_clockRef_set_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_residuals_clockRef_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

allAntennas=get(handles.popupmenu_plot_residuals_refClock, 'String');
newRefclockAntenna=allAntennas{get(handles.popupmenu_plot_residuals_refClock, 'Value')};

writeClockReferenceToOptfile(handles,newRefclockAntenna);


% --- Executes on button press in pushbutton_vie_sched_loadStations.
function pushbutton_vie_sched_loadStations_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_vie_sched_loadStations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_vie_sched_saveStations.
function pushbutton_vie_sched_saveStations_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_vie_sched_saveStations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

saveSchedStatList(hObject,eventdata,handles);


% --- Executes on button press in pushbutton_vie_sched_load_sched_param.
function pushbutton_vie_sched_load_sched_param_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_vie_sched_load_sched_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

loadSchedParameters(hObject,eventdata,handles);


% --- Executes on button press in pushbutton_vie_sched_load_catalog.
function pushbutton_vie_sched_load_catalog_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_vie_sched_load_catalog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
upd_catalogs();

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton_vie_sched_load_catalog.
function pushbutton_vie_sched_load_catalog_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_vie_sched_load_catalog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_run_bslDepWeights.
function checkbox_run_bslDepWeights_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_bslDepWeights (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_bslDepWeights
auto_save_parameterfile(hObject, handles)




% --- Executes on button press in checkbox_vie_sched_distribute.
function checkbox_vie_sched_distribute_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sched_distribute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sched_distribute


% --- Executes on button press in checkbox_setInput_useOptFiles.
function checkbox_setInput_useOptFiles_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_setInput_useOptFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_setInput_useOptFiles

% Enable/Disable popupmenu to select the OPT file
if get(handles.checkbox_setInput_useOptFiles, 'Value') == 0
    set(handles.popupmenu_setInput_optDir, 'Enable', 'off')
    set(handles.radiobutton_estimation_leastSquares_basdepClockoff_OPT, 'Enable', 'off')
else
    set(handles.popupmenu_setInput_optDir, 'Enable', 'on')
    if get(handles.checkbox_estimation_leastSquares_clocksBasDepOffset, 'Value')
        set(handles.radiobutton_estimation_leastSquares_basdepClockoff_OPT, 'Enable', 'on')
    end
end

auto_save_parameterfile(hObject, handles)


% --- Executes on key press with focus on listbox_setInput_processList and none of its controls.
function listbox_setInput_processList_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to listbox_setInput_processList (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_run_sinex_changeAnalystsName.
function checkbox_run_sinex_changeAnalystsName_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_sinex_changeAnalystsName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_sinex_changeAnalystsName
if get(hObject, 'Value')
    newState='On';
else
    newState='Off';
end

set(handles.edit_run_sinex_firstname, 'Enable', newState);
set(handles.edit_run_sinex_lastname, 'Enable', newState);
set(handles.edit_run_sinex_email, 'Enable', newState);


function edit_run_sinex_firstname_Callback(hObject, eventdata, handles)
% hObject    handle to edit_run_sinex_firstname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_run_sinex_firstname as text
%        str2double(get(hObject,'String')) returns contents of edit_run_sinex_firstname as a double


% --- Executes during object creation, after setting all properties.
function edit_run_sinex_firstname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_run_sinex_firstname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_run_sinex_lastname_Callback(hObject, eventdata, handles)
% hObject    handle to edit_run_sinex_lastname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_run_sinex_lastname as text
%        str2double(get(hObject,'String')) returns contents of edit_run_sinex_lastname as a double


% --- Executes during object creation, after setting all properties.
function edit_run_sinex_lastname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_run_sinex_lastname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_run_sinex_email_Callback(hObject, eventdata, handles)
% hObject    handle to edit_run_sinex_email (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_run_sinex_email as text
%        str2double(get(hObject,'String')) returns contents of edit_run_sinex_email as a double


% --- Executes during object creation, after setting all properties.
function edit_run_sinex_email_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_run_sinex_email (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_run_sinex_addSuffix.
function checkbox_run_sinex_addSuffix_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_sinex_addSuffix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_sinex_addSuffix
if get(hObject, 'Value')
    newState='On';
else
    newState='Off';
end

set(handles.edit_run_sinex_suffix, 'Enable', newState);


function edit_run_sinex_suffix_Callback(hObject, eventdata, handles)
% hObject    handle to edit_run_sinex_suffix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_run_sinex_suffix as text
%        str2double(get(hObject,'String')) returns contents of edit_run_sinex_suffix as a double


% --- Executes during object creation, after setting all properties.
function edit_run_sinex_suffix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_run_sinex_suffix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --------------------------------------------------------------------
function menu_minor_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to menu_minor_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_vie_sched_minor_parameters, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);




function edit_vie_sched_experiment_code_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_experiment_code (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sched_experiment_code as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sched_experiment_code as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sched_experiment_code_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_experiment_code (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_vie_sched_experiment_description_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_experiment_description (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sched_experiment_description as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sched_experiment_description as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sched_experiment_description_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_experiment_description (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_vie_sched_print_orbit_data_to_CW.
function checkbox_vie_sched_print_orbit_data_to_CW_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sched_print_orbit_data_to_CW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sched_print_orbit_data_to_CW


% --- Executes on button press in checkbox_vie_sched_create_sky_plot.
function checkbox_vie_sched_create_sky_plot_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sched_create_sky_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.checkbox_vie_sched_create_sky_plot, 'Value');
	set(handles.edit_vie_sched_sky_plot_int, 'Enable', 'on');
else 
	set(handles.edit_vie_sched_sky_plot_int, 'Enable', 'off');
end
% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sched_create_sky_plot




% --- Executes on button press in checkbox_vie_sched_create_elevation_plot.
function checkbox_vie_sched_create_elevation_plot_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sched_create_elevation_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sched_create_elevation_plot



function edit_vie_sched_sky_plot_int_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_sky_plot_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sched_sky_plot_int as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sched_sky_plot_int as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sched_sky_plot_int_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_sky_plot_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_vie_sched_init_prop_interval_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_init_prop_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sched_init_prop_interval as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sched_init_prop_interval as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sched_init_prop_interval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_init_prop_interval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_vie_sched_clear_sat_selection.
function pushbutton_vie_sched_clear_sat_selection_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_vie_sched_clear_sat_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get current content of listbox
curContent=get(handles.listbox_vie_sched_selected_satellites, 'String');


if ~isempty(curContent)
    % delete all selected entries
    curContent(get(handles.listbox_vie_sched_selected_satellites, 'Value'))=[];
    
    % get current value
    curValue=get(handles.listbox_vie_sched_selected_satellites, 'Value');
    
    set(handles.listbox_vie_sched_selected_satellites, 'Value', ...
        min([max(curValue), max(length(curContent),1)]));

    % update listbox
    set(handles.listbox_vie_sched_selected_satellites, 'String', curContent);

    % save changes to handles struct
    guidata(hObject, handles);
end


% --- Executes on selection change in listbox_vie_sched_available_satellites.
function listbox_vie_sched_available_satellites_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_vie_sched_available_satellites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

all_tle_datasets_in_list = cellstr(get(handles.listbox_vie_sched_available_satellites, 'String'));
selected_tle_dataset = all_tle_datasets_in_list(get(handles.listbox_vie_sched_available_satellites, 'Value'));

% make the selected tle datasets 24 digits long
if length(selected_tle_dataset{1})~=24
    selected_tle_dataset{1}(length(selected_tle_dataset{1})+1:24)=' ';
end

% set value to one if needed (otherwise: matlab error!)
if isempty(get(handles.listbox_vie_sched_selected_satellites, 'Value'))
    set(handles.listbox_vie_sched_selected_satellites, 'value', 1)
end
% update listbox
set(handles.listbox_vie_sched_selected_satellites, 'String', ...
    unique([get(handles.listbox_vie_sched_selected_satellites, 'String'); selected_tle_dataset]))

% Choose default command line output for vie_setup
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function listbox_vie_sched_available_satellites_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_vie_sched_available_satellites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_vie_sched_selected_satellites.
function listbox_vie_sched_selected_satellites_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_vie_sched_selected_satellites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_vie_sched_selected_satellites contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_vie_sched_selected_satellites


% --- Executes during object creation, after setting all properties.
function listbox_vie_sched_selected_satellites_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_vie_sched_selected_satellites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_vie_sched_open_scheduler_interface
function checkbox_vie_sched_open_scheduler_interface_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sched_open_scheduler_interface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sched_open_scheduler_interface


% --- Executes on selection change in popupmenu_vie_sched_select_tle.
function popupmenu_vie_sched_select_tle_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_vie_sched_select_tle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_vie_sched_select_tle contents as cell arra
%        contents{get(hObject,'Value')} returns selected item from popupmenu_vie_sched_select_tle

path_tle_folder = '../ORBIT/TLE/';    % TLE data folder

num_selected_element = get(handles.popupmenu_vie_sched_select_tle, 'Value');    % Number of selected element in list.
temp_str = get(handles.popupmenu_vie_sched_select_tle, 'String');

tle_filename = temp_str(num_selected_element,:);

[tle, error_code, error_msg] = loadTLEData([path_tle_folder, tle_filename], 1);

if (error_code == 0)    % No errors occured during tle data loading
    num_tle_datasets = length(tle);
    temp_str = '';

    for i = 1 : (num_tle_datasets - 1)
        temp_str = [temp_str, tle(i).sat_name, '|'];
    end
    temp_str = [temp_str, tle(num_tle_datasets).sat_name];

    set(handles.listbox_vie_sched_available_satellites, 'Value', 1);
	
    set(handles.listbox_vie_sched_available_satellites, 'String', temp_str);
    fprintf(1, '\n');
    fprintf(1, 'Loaded TLE data from file:     %s\n', tle_filename);
    fprintf(1, 'Total number of TLE datasets:  %d\n', num_tle_datasets);
    fprintf(1, '\n');
    
else    % in case of error
    set(handles.listbox_vie_sched_available_satellites, 'String', ' ');
    set(handles.listbox_vie_sched_selected_satellites, 'String', ' ');
    fprintf(1, 'ERROR (load_tle.m): %s', error_msg);
end

% Clear Selection Listbox
set(handles.listbox_vie_sched_selected_satellites, 'String', '');

% --- Executes during object creation, after setting all properties.
function popupmenu_vie_sched_select_tle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_vie_sched_select_tle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_vie_sched_update_tle.
function pushbutton_vie_sched_update_tle_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_vie_sched_update_tle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

downloadAndUpdateTLEData();

% --- Executes on selection change in listbox_vie_sched_sat_obs_network_available.
function listbox_vie_sched_sat_obs_network_available_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_vie_sched_sat_obs_network_available (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

allStatsInList=cellstr(get(handles.listbox_vie_sched_sat_obs_network_available, 'String'));
selectedStat=allStatsInList(get(handles.listbox_vie_sched_sat_obs_network_available, 'Value'));

% make the selected stations 8 digits long
if length(selectedStat{1})~=8
    selectedStat{1}(length(selectedStat{1})+1:8)=' ';
end

% set value to one if needed (otherwise: matlab error!)
if isempty(get(handles.listbox_vie_sched_sat_obs_network_selected, 'Value'))
    set(handles.listbox_vie_sched_sat_obs_network_selected, 'value', 1)
end
% update listbox
set(handles.listbox_vie_sched_sat_obs_network_selected, 'String', ...
    unique([get(handles.listbox_vie_sched_sat_obs_network_selected, 'String'); selectedStat]))

% Choose default command line output for vie_setup
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function listbox_vie_sched_sat_obs_network_available_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_vie_sched_sat_obs_network_available (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in listbox_vie_sched_sat_obs_network_selected.
function listbox_vie_sched_sat_obs_network_selected_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_vie_sched_sat_obs_network_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_vie_sched_sat_obs_network_selected contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_vie_sched_sat_obs_network_selected


% --- Executes during object creation, after setting all properties.
function listbox_vie_sched_sat_obs_network_selected_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_vie_sched_sat_obs_network_selected (see GCBO
% eventdata  reserved - to be defined in a future version of MATLA
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_sched_sat_obs_network_clear.
function pushbutton_sched_sat_obs_network_clear_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_sched_sat_obs_network_clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get current content of listbox
curContent=get(handles.listbox_vie_sched_sat_obs_network_selected, 'String');

if ~isempty(curContent)
    % delete all selected entry
    curContent(get(handles.listbox_vie_sched_sat_obs_network_selected, 'Value'))=[];
    
    % get current value
    curValue=get(handles.listbox_vie_sched_sat_obs_network_selected, 'Value');
    
    set(handles.listbox_vie_sched_sat_obs_network_selected, 'Value', ...
        min([max(curValue), max(length(curContent),1)]));

    % update listbox
    set(handles.listbox_vie_sched_sat_obs_network_selected, 'String', curContent);

    % save changes to handles struct
    guidata(hObject, handles);
end



% --- Executes on button press in checkbox_vie_sched_separate_stations_for_satellites.
function checkbox_vie_sched_separate_stations_for_satellites_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sched_separate_stations_for_satellites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLA
% handles    structure with handles and user data (see GUIDATA)

if get(handles.checkbox_vie_sched_separate_stations_for_satellites, 'Value');
	set(handles.text353, 'Enable', 'on');
	set(handles.text354, 'Enable', 'on');
	set(handles.listbox_vie_sched_sat_obs_network_available, 'Enable', 'on');
	set(handles.listbox_vie_sched_sat_obs_network_selected, 'Enable', 'on');
	set(handles.pushbutton_sched_sat_obs_network_clear, 'Enable', 'on');
else 
	set(handles.text353, 'Enable', 'off');
	set(handles.text354, 'Enable', 'off');
	set(handles.listbox_vie_sched_sat_obs_network_available, 'Enable', 'off');
	set(handles.listbox_vie_sched_sat_obs_network_selected, 'Enable', 'off');
	set(handles.pushbutton_sched_sat_obs_network_clear, 'Enable', 'off');
end
% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sched_separate_stations_for_satellites


% --- Executes on button press in checkbox_vie_sched_use_new_source_file.
function checkbox_vie_sched_use_new_source_file_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sched_use_new_source_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sched_use_new_source_file


% --- Executes when selected object is changed in uipanel_run_modules.
function uipanel_run_modules_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_run_modules 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_run_outDirs_diffSubs.
function checkbox_run_outDirs_diffSubs_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_outDirs_diffSubs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.checkbox_run_outDirs_diffSubs, 'Value')
% Different output sub-directories
	set(handles.edit_run_outDirs_oneSub, 'Enable', 'off')
	set(handles.text_run_outDirs_diff_level0, 'Enable', 'on')
	set(handles.edit_run_outDirs_level0, 'Enable', 'on')
	set(handles.text_run_outDirs_diff_level1, 'Enable', 'on')
	set(handles.edit_run_outDirs_level1, 'Enable', 'on')
	set(handles.text_run_outDirs_diff_level3, 'Enable', 'on')
	set(handles.edit_run_outDirs_level3, 'Enable', 'on')
	set(handles.text_run_outDirs_diff_level2, 'Enable', 'on')
	set(handles.edit_run_outDirs_level2, 'Enable', 'on')
	set(handles.text_run_one_subDir, 'Enable', 'off')
else
% One output sub-directory
	set(handles.edit_run_outDirs_oneSub, 'Enable', 'on')
	set(handles.text_run_outDirs_diff_level0, 'Enable', 'off')
	set(handles.edit_run_outDirs_level0, 'Enable', 'off')
	set(handles.text_run_outDirs_diff_level1, 'Enable', 'off')
	set(handles.edit_run_outDirs_level1, 'Enable', 'off')
	set(handles.text_run_outDirs_diff_level3, 'Enable', 'off')
	set(handles.edit_run_outDirs_level3, 'Enable', 'off')
	set(handles.text_run_outDirs_diff_level2, 'Enable', 'off')
	set(handles.edit_run_outDirs_level2, 'Enable', 'off')
	set(handles.text_run_one_subDir, 'Enable', 'on')
end

% save parameter file automatically 
auto_save_parameterfile(hObject, handles) 



function edit_vie_sched_out_sub_dir_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_out_sub_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sched_out_sub_dir as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sched_out_sub_dir as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sched_out_sub_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_out_sub_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_vie_sched_open_paramtxt.
function pushbutton_vie_sched_open_paramtxt_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_vie_sched_open_paramtxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parameter_file = '../CATALOGS/param.txt';
if exist(parameter_file, 'file')
    % open 
    if ispc
        current_work_dir = pwd;
		% Open file with appropriate Windows applications:
		winopen([current_work_dir(1:(end-4)), 'CATALOGS\param.txt'])
    else
        system(['xterm -e ''vi ', parameter_file, '''']);
    end
else
    msgbox(['Could not find ../CATALOGS/param.txt!'], 'File not found!', 'help')
end


% --- Executes on button press in checkbox_run_error_routine.
function checkbox_run_error_routine_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_error_routine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_error_routine


% --- Executes on button press in checkbox_plot_eopOut_write_detailed_eop_data.
function checkbox_plot_eopOut_write_detailed_eop_data_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_eopOut_write_detailed_eop_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_eopOut_write_detailed_eop_data


% --- Executes on button press in checkbox_plot_eopOut_write_sorted_eop_data.
function checkbox_plot_eopOut_write_sorted_eop_data_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_eopOut_write_sorted_eop_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_eopOut_write_sorted_eop_data


% --- Executes on button press in checkbox_plot_eopOut_write_vievs_eop_data.
function checkbox_plot_eopOut_write_vievs_eop_data_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_eopOut_write_vievs_eop_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_eopOut_write_vievs_eop_data


% --- Executes on button press in checkbox_run_globalPram_tidERPvar.
function checkbox_run_globalPram_tidERPvar_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_globalPram_tidERPvar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_globalPram_tidERPvar



% --- Executes on button press in checkbox_global_param_tidERPvar.
function checkbox_global_param_tidERPvar_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_global_param_tidERPvar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_global_param_tidERPvar


% --- Executes on button press in checkbox_plot_eopOut_basRepOptions_writeBasOut.
function checkbox_plot_eopOut_basRepOptions_writeBasOut_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_eopOut_basRepOptions_writeBasOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_eopOut_basRepOptions_writeBasOut
if get(hObject,'Value')
    newState='on';
else
    newState='off';
end
set(handles.edit_plot_eopOut_basRepOptions_writeBasOutFname, 'Enable', newState)


function edit_plot_eopOut_basRepOptions_minBasObs_Callback(hObject, eventdata, handles)
% hObject    handle to edit_plot_eopOut_basRepOptions_minBasObs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_plot_eopOut_basRepOptions_minBasObs as text
%        str2double(get(hObject,'String')) returns contents of edit_plot_eopOut_basRepOptions_minBasObs as a double
defaultVal=10;

if isnan(str2double(get(hObject,'String')))
    set(hObject, 'String', num2str(defaultVal));
end

% --- Executes during object creation, after setting all properties.
function edit_plot_eopOut_basRepOptions_minBasObs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_plot_eopOut_basRepOptions_minBasObs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_plot_eopOut_basRepOptions_simplePlot.
function checkbox_plot_eopOut_basRepOptions_simplePlot_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_eopOut_basRepOptions_simplePlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_eopOut_basRepOptions_simplePlot


% --- Executes on button press in checkbox_plot_eopOut_basRepOptions_commandWindowOutput.
function checkbox_plot_eopOut_basRepOptions_commandWindowOutput_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_eopOut_basRepOptions_commandWindowOutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_eopOut_basRepOptions_commandWindowOutput



function edit_plot_eopOut_basRepOptions_writeBasOutFname_Callback(hObject, eventdata, handles)
% hObject    handle to edit_plot_eopOut_basRepOptions_writeBasOutFname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_plot_eopOut_basRepOptions_writeBasOutFname as text
%        str2double(get(hObject,'String')) returns contents of edit_plot_eopOut_basRepOptions_writeBasOutFname as a double


% --- Executes during object creation, after setting all properties.
function edit_plot_eopOut_basRepOptions_writeBasOutFname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_plot_eopOut_basRepOptions_writeBasOutFname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_setInput_browseVgos.
function pushbutton_setInput_browseVgos_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_setInput_browseVgos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

vgosDir='../DATA/vgosDB/';
if ~exist(vgosDir,'dir')
    mkdir(vgosDir);
end

out = uipickfiles('FilterSpec', vgosDir);


if iscell(out)
    out=out';
    
    % Format: yyyy/<session_name> [vgosDB]
    
    for iF = 1 : length(out)
        curSlash = sort([strfind(out{iF},'/'), strfind(out{iF},'\')]);
        tgzDot = sort(strfind(out{iF},'.'));
        if length(tgzDot) > 2
            out{iF} = [out{iF}(curSlash(end-1)+1 : tgzDot(3)-1), ' [vgosDB]'];
        else
            out{iF} = [out{iF}(curSlash(end-1)+1 : end), ' [vgosDB]'];
        end
    end

    updateInputFilesBox(hObject, eventdata,handles,out)   
end





% --- Executes on button press in checkbox_estimation_leastSquares_sources_ICRF2_def.
function checkbox_estimation_leastSquares_sources_ICRF2_def_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_estimation_leastSquares_sources_ICRF2_def (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_estimation_leastSquares_sources_ICRF2_def
auto_save_parameterfile(hObject, handles)



function edit_estimation_leastSquares_sources_abs_constr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_sources_abs_constr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_estimation_leastSquares_sources_abs_constr as text
%        str2double(get(hObject,'String')) returns contents of edit_estimation_leastSquares_sources_abs_constr as a double


% --- Executes during object creation, after setting all properties.
function edit_estimation_leastSquares_sources_abs_constr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_estimation_leastSquares_sources_abs_constr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menu_models_space_crafts_Callback(hObject, eventdata, handles)
% hObject    handle to uipanel_models_space_crafts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set all uipanels to invisibel
setAllPanelsToInvisible(hObject, handles)

% set the one panel to visible
set(handles.uipanel_models_space_crafts, 'Visible', 'On');

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbutton_models_sc_browse_for_sp3.
function pushbutton_models_sc_browse_for_sp3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_models_sc_browse_for_sp3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Select Sp3 file:
[FileName, PathName] = uigetfile('*.*','Select SP3 file', '../ORBIT', 'multiselect', 'off');

if ischar(FileName) && ischar(PathName)
    set(handles.edit_models_sc_sp3_file, 'String', [PathName, FileName])
    
    % save parameter file automatically 
    auto_save_parameterfile(hObject, handles)
end
    
    

function edit_models_sc_sp3_file_Callback(hObject, eventdata, handles)
% hObject    handle to edit_models_sc_sp3_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_models_sc_sp3_file as text
%        str2double(get(hObject,'String')) returns contents of edit_models_sc_sp3_file as a double


% --- Executes during object creation, after setting all properties.
function edit_models_sc_sp3_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_models_sc_sp3_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_vie_sim_rng_use_seed.
function checkbox_vie_sim_rng_use_seed_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sim_rng_use_seed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sim_rng_use_seed
% Enable/disable edit field to enter the seed number for the random number generator:

if(get(handles.checkbox_vie_sim_rng_use_seed,'Value'))
    set(handles.text_vie_sim_seed, 'Enable', 'on') 
    set(handles.edit_vie_sim_rng_seed, 'Enable', 'on')
else
     set(handles.text_vie_sim_seed, 'Enable', 'off')
     set(handles.edit_vie_sim_rng_seed, 'Enable', 'off')
end


function edit_vie_sched_pathCatalogs_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_pathCatalogs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sched_pathCatalogs as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sched_pathCatalogs as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sched_pathCatalogs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_pathCatalogs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_vie_sim_rng_seed_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_rng_seed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vie_sim_rng_seed as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sim_rng_seed as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sim_rng_seed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_rng_seed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_vie_sched_browsCatalogs.
function pushbutton_vie_sched_browsCatalogs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_vie_sched_browsCatalogs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
folder_name = uigetdir('../');
if folder_name ~= 0
    handles.edit_vie_sched_pathCatalogs.String = [folder_name '\'];
end


function edit_vie_sim_noise_sat_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_noise_sat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit_vie_sim_noise_sat as text
%        str2double(get(hObject,'String')) returns contents of edit_vie_sim_noise_sat as a double


% --- Executes during object creation, after setting all properties.
function edit_vie_sched_bandnum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_bandnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit_vie_sim_noise_sat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sim_noise_sat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uitable_vie_sched_bandInfo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable_vie_sched_bandInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
d = cell(2,5);
d{1,1} = true;
d{1,2} = 'X';
d{1,3} = 0.034896;
d{1,4} = 10;
d{1,5} = 20;
d{2,1} = true;
d{2,2} = 'S';
d{2,3} = 3.800;
d{2,4} = 6;
d{2,5} = 15;
hObject.Data = d;


% --- Executes on button press in checkbox_vie_sched_disp_sky_coverage.
function checkbox_vie_sched_disp_sky_coverage_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sched_disp_sky_coverage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sched_disp_sky_coverage


% --- Executes on button press in checkbox_vie_sim_write_vso_file.
function checkbox_vie_sim_write_vso_file_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sim_write_vso_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sim_write_vso_file


% --- Executes on button press in checkbox_vie_sched_multiSched.
function checkbox_vie_sched_multiSched_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sched_multiSched (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value')
    handles.checkbox_vie_sched_saveOutput.Value = 1;
    handles.checkbox_vie_sched_saveOutput.Enable = 'off';
    handles.checkbox_vie_sched_dispSkyCoverage = 0;
    handles.checkbox_vie_sched_saveOutput.Enable = 'off';
else
    handles.checkbox_vie_sched_saveOutput.Value = 0;
    handles.checkbox_vie_sched_saveOutput.Enable = 'on';
    handles.checkbox_vie_sched_dispSkyCoverage = 1;
    handles.checkbox_vie_sched_saveOutput.Enable = 'on';
end
% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sched_multiSched



% --- Executes on button press in pushbutton_multiSched.
function pushbutton_multiSched_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_multiSched (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parameter_file = '../CATALOGS/multiSched.txt';
if exist(parameter_file, 'file')
    % open 
    if ispc
        current_work_dir = pwd;
		% Open file with appropriate Windows applications:
		winopen([current_work_dir(1:(end-4)), 'CATALOGS\multiSched.txt'])
    else
        system(['xterm -e ''vi ', parameter_file, '''']);
    end
else
    msgbox(['Could not find ../CATALOGS/multiSched.txt!'], 'File not found!', 'help')
end


% --- Executes on key press with focus on edit_vie_sched_init_prop_interval and none of its controls.
function edit_vie_sched_init_prop_interval_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to edit_vie_sched_init_prop_interval (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_run_ElevDepNoise.
function checkbox_run_ElevDepNoise_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_ElevDepNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of checkbox_run_ElevDepNoise
if get(handles.checkbox_run_ElevDepNoise, 'Value')
    set(handles.edit_run_add_noise, 'Enable', 'off'); 
    set(handles.text_run_add_noise, 'Enable', 'off');
else
    set(handles.edit_run_add_noise, 'Enable', 'on');
    set(handles.text_run_add_noise, 'Enable', 'on');
end
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox_run_globalPram_stseaspos.
function checkbox_run_globalPram_stseaspos_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_globalPram_stseaspos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_globalPram_stseaspos
auto_save_parameterfile(hObject, handles)


% --- Executes on button press in checkbox_run_globalPram_hlpole.
function checkbox_run_globalPram_hlpole_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_globalPram_hlpole (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_globalPram_hlpole
auto_save_parameterfile(hObject, handles)


% --- Executes on selection change in listbox_vie_glob_stseaspos.
function listbox_vie_glob_stseaspos_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_vie_glob_stseaspos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_vie_glob_stseaspos contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_vie_glob_stseaspos


% --- Executes during object creation, after setting all properties.
function listbox_vie_glob_stseaspos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_vie_glob_stseaspos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_global_param_stseaspos_N.
function checkbox_global_param_stseaspos_N_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_global_param_stseaspos_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_global_param_stseaspos_N


% --- Executes on button press in checkbox_global_param_stseaspos_E.
function checkbox_global_param_stseaspos_E_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_global_param_stseaspos_E (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_global_param_stseaspos_E


% --- Executes on button press in checkbox_global_param_stseaspos_R.
function checkbox_global_param_stseaspos_R_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_global_param_stseaspos_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_global_param_stseaspos_R




% --- Executes when selected object is changed in uibuttongroup10.
function uibuttongroup31_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup10 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.radiobuttonedit_glob_param_estStSeasPos, 'Value')
    % estimate was selected
    state='On';
else
    state='Off';
end

set(handles.listbox_vie_glob_stseaspos, 'Enable', state);
set(handles.checkbox_global_param_stseaspos_R, 'Enable', state);
set(handles.checkbox_global_param_stseaspos_E, 'Enable', state);
set(handles.checkbox_global_param_stseaspos_N, 'Enable', state);







% --- Executes on button press in checkbox_parameters_statCorr_temp_fromInSitu.
function checkbox_parameters_statCorr_temp_fromInSitu_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameters_statCorr_temp_fromInSitu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameters_statCorr_temp_fromInSitu


% --- Executes on button press in checkbox_parameters_statCorr_temp_GPT3.
function checkbox_parameters_statCorr_temp_GPT3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameters_statCorr_temp_GPT3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameters_statCorr_temp_GPT3


% --- Executes on button press in checkbox_parameters_statCorr_APLrg.
function checkbox_parameters_statCorr_APLrg_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameters_statCorr_APLrg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameters_statCorr_APLrg
if get(hObject, 'Value')
    newState1='on';
    newState2='off';
else
    newState1='off';
    newState2='on';
end

set(handles.popupmenu_parameters_statCorr_APLrg, 'Enable', newState1)

set(handles.checkbox_parameters_statCorr_tidalAtmoLoad, 'Enable', newState2)
set(handles.popupmenu_parameters_statCorr_tidalAtmoOceanLoad, 'Enable', newState2)
set(handles.checkbox_parameters_statCorr_nonTidalAtmoLoad, 'Enable', newState2)
set(handles.popupmenu_parameters_statCorr_nonTidalAtmoOceanLoad, 'Enable', newState2)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)


% --- Executes on selection change in popupmenu_parameters_statCorr_APLrg.
function popupmenu_parameters_statCorr_APLrg_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_statCorr_APLrg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_parameters_statCorr_APLrg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_parameters_statCorr_APLrg


% --- Executes during object creation, after setting all properties.
function popupmenu_parameters_statCorr_APLrg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_statCorr_APLrg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_parameters_statCorr_GIA.
function checkbox_parameters_statCorr_GIA_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_parameters_statCorr_GIA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_parameters_statCorr_GIA
if get(hObject, 'Value')
    newState='on';
else
    newState='off';
end

set(handles.popupmenu_parameters_statCorr_GIA, 'Enable', newState)

% save parameter file automatically 
auto_save_parameterfile(hObject, handles)




% --- Executes on selection change in popupmenu_parameters_statCorr_GIA.
function popupmenu_parameters_statCorr_GIA_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_statCorr_GIA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_parameters_statCorr_GIA contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_parameters_statCorr_GIA


% --- Executes during object creation, after setting all properties.
function popupmenu_parameters_statCorr_GIA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_parameters_statCorr_GIA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in listbox_vie_glob_APLrg.
function listbox_vie_glob_APLrg_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_vie_glob_APLrg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_vie_glob_APLrg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_vie_glob_APLrg


% --- Executes during object creation, after setting all properties.
function listbox_vie_glob_APLrg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_vie_glob_APLrg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uibuttongroup32.
function uibuttongroup32_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup32 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.radiobuttonedit_glob_param_estAPLrg, 'Value')
    % estimate was selected
    state='On';
else
    state='Off';
end

set(handles.listbox_vie_glob_APLrg, 'Enable', state);



% --- Executes on button press in checkbox_run_globalPram_APLrg.
function checkbox_run_globalPram_APLrg_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_run_globalPram_APLrg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_run_globalPram_APLrg
auto_save_parameterfile(hObject, handles)



function edit_run_add_noise_Callback(hObject, eventdata, handles)
% hObject    handle to edit_run_add_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_run_add_noise as text
%        str2double(get(hObject,'String')) returns contents of edit_run_add_noise as a double
auto_save_parameterfile(hObject, handles)


% --- Executes during object creation, after setting all properties.
function edit_run_add_noise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_run_add_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_vie_sched_conditions.
function checkbox_vie_sched_conditions_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_vie_sched_conditions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if hObject.Value
    handles.popupmenu_vie_sched_conditions.Enable = 'on';
    handles.pushbutton_vie_sched_conditions.Enable = 'on';
    handles.text_vie_sched_conditions.Enable = 'on';
else
    handles.popupmenu_vie_sched_conditions.Enable = 'off';
    handles.pushbutton_vie_sched_conditions.Enable = 'off';
    handles.text_vie_sched_conditions.Enable = 'off';
    handles.text_vie_sched_conditions.String = '';
end

% Hint: get(hObject,'Value') returns toggle state of checkbox_vie_sched_conditions


% --- Executes on selection change in popupmenu_vie_sched_conditions.
function popupmenu_vie_sched_conditions_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_vie_sched_conditions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_vie_sched_conditions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_vie_sched_conditions


% --- Executes during object creation, after setting all properties.
function popupmenu_vie_sched_conditions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_vie_sched_conditions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_vie_sched_conditions.
function pushbutton_vie_sched_conditions_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_vie_sched_conditions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
con = sched_optimization;
handles.text_vie_sched_conditions.String = con;



% --- Executes on button press in checkbox_plot_show_errorbars.
function checkbox_plot_show_errorbars_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_plot_show_errorbars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_plot_show_errorbars
handles=plotToAxes(hObject, handles);
% Update handles structure
guidata(hObject, handles)


% --- Executes on selection change in popupmenu_plot_select_time_ref_format.
function popupmenu_plot_select_time_ref_format_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_select_time_ref_format (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_plot_select_time_ref_format contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_plot_select_time_ref_format
handles=plotToAxes(hObject, handles);
% Update handles structure
guidata(hObject, handles)



% --- Executes during object creation, after setting all properties.
function popupmenu_plot_select_time_ref_format_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_plot_select_time_ref_format (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function rerun_vie_lsm(hObject, eventdata, handles)
try
    idx = handles.popupmenu_plot_residuals_session.Value;

    if handles.radiobutton_plot_residuals_perAll.Value
        flag = 0;
    elseif handles.radiobutton_plot_residuals_perBasel.Value
        flag = 1;
        flag_idx = handles.popupmenu_plot_residuals_baseline.Value;
    elseif handles.radiobutton_plot_residuals_perStat.Value
        flag = 2;
        flag_idx = handles.popupmenu_plot_residuals_station.Value;
    else
        flag = 3;
        flag_idx = handles.popupmenu_plot_residuals_source.Value;
    end
    ax = gca;
    xlim = ax.XLim;
    ylim = ax.YLim;

    name = handles.popupmenu_plot_residuals_session.String{idx};
    eliminOutliers = handles.checkbox_setInput_eliminOutliers.Value;
    if eliminOutliers
        allOutDirInList = get(handles.popupmenu_setInput_outDir, 'String');
        outlierDirectory = allOutDirInList{get(handles.popupmenu_setInput_outDir, 'Value')};
    else
        outlierDirectory = 'none';
    end
    if handles.checkbox_run_simpleOutlierTest.Value == 0 && handles.checkbox_run_normalOutlierTest.Value == 0
        outlierTest = 'none';
    elseif handles.checkbox_run_simpleOutlierTest.Value
        outlierTest = 'simple';
    elseif handles.checkbox_run_normalOutlierTest.Value
        outlierTest = 'normal';
    end
    
    vie_batch(name, outlierTest, outlierDirectory);
    handles=loadResFiles(hObject, handles);
    handles.popupmenu_plot_residuals_session.Value = idx;
    handles=updatePopupmenusInResidualPlot(hObject, handles);

    if flag == 0
        handles.radiobutton_plot_residuals_perAll.Value = 1;
    elseif flag == 1
        handles.radiobutton_plot_residuals_perBasel.Value = 1;
        handles.popupmenu_plot_residuals_baseline.Value = flag_idx;
        handles.popupmenu_plot_residuals_baseline.Enable = 'on';
    elseif flag == 2
        handles.radiobutton_plot_residuals_perStat.Value = 1;
        handles.popupmenu_plot_residuals_station.Value = flag_idx;
        handles.popupmenu_plot_residuals_station.Enable = 'on';
    else
        handles.radiobutton_plot_residuals_perSource.Value = 1;
        handles.popupmenu_plot_residuals_source.Value = flag_idx;
        handles.popupmenu_plot_residuals_source.Enable = 'on';
    end
    handles=plotResidualsToAxes(handles);
    ax.XLim = xlim;
    ax.YLim = ylim;
    guidata(hObject, handles)
catch ex
    warning('ERROR while rerunning vie_lsm!\nMessage: %s',ex.message);
end


function reference(hObject, eventdata, handles)
try
    web('http://vievswiki.geo.tuwien.ac.at/doku.php?id=public:vievs_publications:publications','-browser')
catch
    warning('An error occured when opening: http://vievswiki.geo.tuwien.ac.at/doku.php?id=public:vievs_publications:publications')
end

function wiki(hObject, eventdata, handles)
try
    web('http://vievswiki.geo.tuwien.ac.at/doku.php','-browser')
catch
    warning('An error occured when opening: http://vievswiki.geo.tuwien.ac.at/doku.php')
end

function github(hObject, eventdata, handles)
try
    web('https://github.com/TUW-VieVS','-browser')
catch
    warning('An error occured when opening: https://github.com/TUW-VieVS')
end

function pushbutton_vie_sim_browseStatistics_Callback(hObject,eventdata,handles)
    [file,path] = uigetfile('*.csv','Browse for VieSched++ statistics.csv file');
    handles.edit_vie_sim_statistics_csv.String = [path file];
