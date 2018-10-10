% #########################################################################
% #     saveSchedParameters
% #########################################################################
%
% DESCRITPION
% This function saves the scheduling parameters from the GUI to the LEVEL5
% folder. Those files are then used within the sched program.
%
% AUTHOR 
%   Matthias Madzak (?)
%
% INPUT
%   handles      structure from the GUI (also containing data, e.g. residuals)
%   hObject      ...unused so far.
%
% OUTPUT
%   error_code   = 0, if no error occurs.
%
% CHANGES
%   2015-20-01, A. Hellerschmied: Exctact code from "vievs2_3.m" file and
%     create this file.
%   2015-02-02, A. Hellerschmied: added fields related to GUI update.
%   2015-02-11, A. Hellerschmied: PARA.EXPER, PARA.DESCRIPTION now enterd via GUI.
%   2016-09-06  M. Schartner: changes for sched_manually and sched_analyser
%   2016-10-05, M. Schartner: browse catalog added
%   2016-10-06, M. Schartner: Bandinfo added
%   2016-10-19, M. Schartner: changes for disp Sky Coverage
%   2016-10-19, M. Schartner: no longer saves subdirectory name
%   2017-02-22, M. Schartner: new parameters for conditions


function [error_code] = saveSchedParameters(hObject, handles)

% #### Init ####
error_code = 0;  % No error


% (1) ##### station list #####
% check if listbox of stations is empty
stanetstr=char(get(handles.listbox_vie_sched_stanet, 'String'));
if ~isempty(stanetstr)
    save ('../DATA/LEVEL5/stanetstr.mat','stanetstr');
else 
    msgbox('Select station network for scheduling in Scheduling menu!', 'No station for scheduling selected!', 'warn');
    error_code = 1;
    return;
end


% (2) ##### scheduling parameters #####

% (2.1) ### Get logical values from GUI ###
PARA.USE_NEW_SOURCE_File    = get(handles.checkbox_vie_sched_use_new_source_file, 'Value');

% (2.2) ### Get session timing parameters from GUI ###
PARA.year_s   = str2double(get(handles.edit_vie_sched_years,'string'));
PARA.month_s  = str2double(get(handles.edit_vie_sched_mons,'string'));
PARA.day_s    = str2double(get(handles.edit_vie_sched_days,'string'));
PARA.hour_s   = str2double(get(handles.edit_vie_sched_hours,'string'));
PARA.minute_s = str2double(get(handles.edit_vie_sched_mins,'string'));
PARA.second_s = str2double(get(handles.edit_vie_sched_secs,'string'));
PARA.duration = str2double(get(handles.edit_vie_sched_duration,'string'));

if (isnan(PARA.year_s)||isnan(PARA.month_s)||isnan(PARA.day_s)||isnan(PARA.duration)||isnan(PARA.hour_s) ||isnan(PARA.minute_s)||isnan(PARA.second_s))
        msgbox('Please check start and duration of the session!', 'Invalid session time options!', 'warn');
        error_code = 2;
        return;
end

% (2.3) ### Get Band values from GUI ###
data = get(handles.uitable_vie_sched_bandInfo,'DATA');

c = 1;
for i = 1:2
    if data{i,1}
        if i==1
            PARA.BAND(c)     = 'X';   % 1-X
        else
            PARA.BAND(c)     = 'S';   % 2-S
        end
        PARA.MIN_SNR(c) = data{i,5};
        PARA.WAVEL(c) = data{i,3};
        PARA.CHANUM(c) = data{i,4};
        c = c+1;
    end
end

PARA.MAX_BANDNUM = c-1;

% (2.4) ### Get further parameters from GUI ###
PARA.MIN_SUNDIST = str2double(get(handles.edit_vie_sched_sundistmin,'String'))* pi / 180.0;
if (isnan(PARA.MIN_SUNDIST))
    msgbox('Invalid sun distance value (NaN)!', 'Invalid VIE_SCHED parameter input!', 'warn');
    error_code = 8;
    return;
end
PARA.MIN_CUTEL = str2double(get(handles.edit_vie_sched_cutel,'String'))* pi / 180.0;
if (isnan(PARA.MIN_CUTEL) || (PARA.MIN_CUTEL < 0) || (PARA.MIN_CUTEL > (90 * pi / 180.0)))
    msgbox('Invalid cut-off elevation value (valid between 0 and 90 deg)!', 'Invalid VIE_SCHED parameter input!', 'warn');
    error_code = 9;
    return;
end
PARA.MIN_FLUX = str2double(get(handles.edit_vie_sched_flux,'String'));
if (isnan(PARA.MIN_FLUX))
    msgbox('Invalid minimum flux density value (NaN)!', 'Invalid VIE_SCHED parameter input!', 'warn');
    error_code = 10;
    return;
end
PARA.EXPER = get(handles.edit_vie_sched_experiment_code,'String');
PARA.DESCRIPTION = get(handles.edit_vie_sched_experiment_description,'String');


% (2.5) ### Observing mode for twin telescopes ###
if (get(handles.radiobutton_vie_sched_twinm1,'Value'))
    PARA.TWINMODE = 1;
end
if (get(handles.radiobutton_vie_sched_twinm2,'Value'))
    PARA.TWINMODE = 2;
end
if (get(handles.radiobutton_vie_sched_twinm3,'Value'))
    PARA.TWINMODE = 3;
end

% (2.6) ### Get satellite scheduling parameters from GUI ###

% Get delta_t for satelite orbit propagation
PARA.INIT_PROP_INTERVAL = str2double(get(handles.edit_vie_sched_init_prop_interval, 'String'));
if (isnan(PARA.INIT_PROP_INTERVAL))
    msgbox('Invalid initial orbit propagation interval!', 'Invalid VIE_SCHED parameter input!', 'warn');
    error_code = 7;
    return;
end

% Get auxiliary output parameters
PARA.PRINT_ORBIT_TIMESERIES = get(handles.checkbox_vie_sched_print_orbit_data_to_CW, 'Value');
PARA.CREATE_SKYPLOTS = get(handles.checkbox_vie_sched_create_sky_plot, 'Value');
PARA.SKYPLOT_MARKER_INT = str2double(get(handles.edit_vie_sched_sky_plot_int, 'String'));
if (isnan(PARA.SKYPLOT_MARKER_INT))
    msgbox('Invalid skyplot marker interval!', 'Invalid VIE_SCHED parameter input!', 'warn');
    error_code = 6;
    return;
end
PARA.CREATE_ELEV_PLOT = get(handles.checkbox_vie_sched_create_elevation_plot, 'Value');

% Open Scheduler interface and create VEX files
PARA.SCHEDULE_WRITE_VEX = get(handles.checkbox_vie_sched_open_scheduler_interface, 'Value');

% Get separate statio network for satellite observations
PARA.USE_SEPARATE_STAT_NWW_FOR_SATS = get(handles.checkbox_vie_sched_separate_stations_for_satellites, 'Value');
if get(handles.checkbox_vie_sched_separate_stations_for_satellites, 'Value');
    % Get station names and save them to a *.mat file
    stanet_sat_str = char(get(handles.listbox_vie_sched_sat_obs_network_selected, 'String'));
    if ~isempty(stanet_sat_str)
        save ('../DATA/LEVEL5/stanet_sat_str.mat','stanet_sat_str');
    else
        msgbox('Please select a station network for satellite observbations!', 'Error!', 'warn');
        error_code = 6;
        return;
    end
end


% (2.6) ### Get output options from GUI ###
PARA.ngs                    = get(handles.checkbox_vie_sched_ngs, 'Value');
PARA.skdsum                 = get(handles.checkbox_vie_sched_sum, 'Value');
PARA.skd                    = get(handles.checkbox_vie_sched_skd, 'Value');

% Putput Subdirectory in /DATA/SCHED/ 
PARA.OUTPUT_SUBDIR = get(handles.edit_run_outDirs_oneSub, 'String');
if ~isempty(PARA.OUTPUT_SUBDIR)
    PARA.USE_OUTPUT_SUBDIR = 1;
else
    PARA.USE_OUTPUT_SUBDIR = 0;
end


% (3) ##### Selection of scheduling mode/strategy ##### 

% ### (3.1) Source based strategy ###
if (get(handles.radiobutton_vie_sched_rbsource,'Value'))
    PARA.OPTIMIZATION = 1;
    if (length(get(handles.edit_vie_sched_srcnum,'string')) == 0)
        msgbox('Please select the number of sources observed simultaneously!', 'Error!', 'warn')
        error_code = 4;
        return;
    else
        PARA.SRCNUM = str2double(get(handles.edit_vie_sched_srcnum,'string'));
        if isnan(PARA.SRCNUM) || ~(PARA.SRCNUM == 1 || PARA.SRCNUM == 2 || PARA.SRCNUM == 4)
            msgbox('Number of sources selected simultaneously: invalid input!', 'Error!', 'warn')
            error_code = 4;
            return;
        end
    end
end

% ### (3.2) Station based strategy ###
if (get(handles.radiobutton_vie_sched_rbstation,'Value'))
    PARA.OPTIMIZATION = 2;
end
PARA.DISTRIBUTE=get(handles.checkbox_vie_sched_distribute, 'Value');

% ### (3.3) satellite scheduling ###
if (get(handles.radiobutton_vie_sched_satellites,'Value'))
    PARA.OPTIMIZATION = 3;
 
    % satellite parameter:
    satellite_str = char(get(handles.listbox_vie_sched_selected_satellites, 'String'));
    if isempty(satellite_str)
        msgbox('Please select satellites!', 'No satellites selected', 'warn');
        error_code = 5;
        return;
    end

    % Save list of selected satellites to *.mat file 
    save ('../DATA/LEVEL5/satellite_str.mat','satellite_str');

    % Get TLE Filename
    temp_str = get(handles.popupmenu_vie_sched_select_tle, 'String');
    PARA.TLE_FILENAME = temp_str(get(handles.popupmenu_vie_sched_select_tle, 'Value'), :);

end
 
if (get(handles.radiobutton_vie_sched_manually,'Value'))
    PARA.OPTIMIZATION = 4;
end


if (get(handles.checkbox_vie_sched_openSchedAnalyser,'Value'))
	PARA.openAnalyser = 1;
else
	PARA.openAnalyser = 0;
end
 
if (get(handles.checkbox_vie_sched_disp_sky_coverage,'Value'))
	PARA.disp_sky_coverage = 1;
else
	PARA.disp_sky_coverage = 0;
end

% ### (3.4) Error case ###
if(PARA.OPTIMIZATION == 0)
    msgbox('Please select the optimization criteria!', 'Error!', 'warn');
    error_code = 7;
    return; 
end

if (get(handles.checkbox_vie_sched_multiSched,'Value'))
	PARA.MULTISCHED = 1;
else
	PARA.MULTISCHED = 0;
end

if (get(handles.checkbox_vie_sched_saveOutput,'Value'))
	PARA.SAVEOUTPUT = 1;
else
	PARA.SAVEOUTPUT = 0;
end


PARA.optimization_type = handles.popupmenu_vie_sched_conditions.String(handles.popupmenu_vie_sched_conditions.Value);
PARA.optimization_condition = handles.text_vie_sched_conditions.String;

PARA.pinfile = handles.edit_vie_sched_pathCatalogs.String;

PARA.parallel = handles.checkbox_run_runOptions_parallel.Value;
PARA.nCores = handles.popupmenu_run_runOptions_nCores.String{handles.popupmenu_run_runOptions_nCores.Value};

% ##### (4) Save selected parameters to *.mat file ##### 
save ('../DATA/LEVEL5/schedparam.mat','PARA');
