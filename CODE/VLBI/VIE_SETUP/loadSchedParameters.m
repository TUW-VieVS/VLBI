% #########################################################################
% #     loadSchedParameters
% #########################################################################
%
% DESCRITPION
% This function loads the scheduling parameters from the *.mat file defined 
% below.
%
% AUTHOR 
%   (?) 9.4.2014
%
% INPUT
%   handles      structure from the GUI (also containing data, e.g. residuals)
%   hObject      ...unused so far.
%   eventdata    ...unused so far.
%
% OUTPUT
%
% CHANGES
%   2015-02-02, A. Hellerschmied: added fields related to GUI update.
%   2015-02-11, A. Hellerschmied: PARA.EXPER, PARA.DESCRIPTION now enterd via GUI.
%   2016-09-12, M. Schartner: changes for sched manually and sched analyser
%   2016-10-05, M. Schartner: browse catalog added
%   2016-10-06, M. Schartner: snr information deleted
%   2016-10-19, M. Schartner: changes for disp Sky Coverage
%   2016-10-19, M. Schartner: no longer loads subdirectory name
%   2017-02-22, M. Schartner: new parameters for conditions




function loadSchedParameters(hObject,eventdata,handles)

    % ##### (1.) Load Parameter File #####
    schedPFile='../DATA/LEVEL5/schedparam.mat';

    if ~exist(schedPFile, 'file')
        msgbox('schedparam.mat does not exist', 'warning', 'warn');
        return;
    end

    schedP=load(schedPFile);
    fields=fieldnames(schedP);
    schedP=schedP.(fields{1});


    % ##### (2.) Set scheduling parameters in GUI #####

    % (2.1) ### Set logical values###
    set(handles.checkbox_vie_sched_use_new_source_file, 'Value', schedP.USE_NEW_SOURCE_File);

    % (2.2) ### Set session timing parameters###
    set(handles.edit_vie_sched_years,'string', num2str(schedP.year_s));
    set(handles.edit_vie_sched_mons,'string', num2str(schedP.month_s));
    set(handles.edit_vie_sched_days,'string', num2str(schedP.day_s));
    set(handles.edit_vie_sched_hours,'string', num2str(schedP.hour_s));
    set(handles.edit_vie_sched_mins,'string', num2str(schedP.minute_s));
    set(handles.edit_vie_sched_secs,'string', num2str(schedP.second_s));
    set(handles.edit_vie_sched_duration,'string', num2str(schedP.duration));

    % (2.3) ### Set SNR values ###
%     set(handles.edit_vie_sched_edsnrxband,'String', num2str(schedP.MIN_SNR(1)));
%     set(handles.edit_vie_sched_edsnrsband,'String', num2str(schedP.MIN_SNR(2)));

    % (2.4) ### Set further parameters ###
    set(handles.edit_vie_sched_sundistmin,'String', num2str(schedP.MIN_SUNDIST*180/pi));
    set(handles.edit_vie_sched_cutel,'String', num2str(schedP.MIN_CUTEL*180/pi));
    set(handles.edit_vie_sched_flux,'String', num2str(schedP.MIN_FLUX));
    set(handles.edit_vie_sched_experiment_code,'String', schedP.EXPER);
    set(handles.edit_vie_sched_experiment_description,'String', schedP.DESCRIPTION);

    % (2.5) ### Set observing mode for twin telescopes ###
    switch schedP.TWINMODE
        case 1
            set(handles.radiobutton_vie_sched_twinm1,'Value', 1);
        case 2
            set(handles.radiobutton_vie_sched_twinm2, 'Value', 1);
        case 3
            set(handles.radiobutton_vie_sched_twinm3, 'Value', 1);
    end

    % (2.6) ### Set satellite scheduling parameters ###
    set(handles.edit_vie_sched_init_prop_interval,'String', num2str(schedP.INIT_PROP_INTERVAL));
    set(handles.checkbox_vie_sched_print_orbit_data_to_CW, 'Value', schedP.PRINT_ORBIT_TIMESERIES);
    set(handles.checkbox_vie_sched_create_sky_plot, 'Value', schedP.CREATE_SKYPLOTS);
    set(handles.edit_vie_sched_sky_plot_int,'String', num2str(schedP.SKYPLOT_MARKER_INT));
    set(handles.checkbox_vie_sched_create_elevation_plot, 'Value', schedP.CREATE_ELEV_PLOT);
    set(handles.checkbox_vie_sched_open_scheduler_interface, 'Value', schedP.SCHEDULE_WRITE_VEX);
    set(handles.checkbox_vie_sched_separate_stations_for_satellites, 'Value', schedP.USE_SEPARATE_STAT_NWW_FOR_SATS);
    
    if schedP.USE_SEPARATE_STAT_NWW_FOR_SATS;
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
    
    % (2.6) ### Set output options ###
    set(handles.checkbox_vie_sched_ngs, 'Value', schedP.ngs);
    set(handles.checkbox_vie_sched_sum, 'Value', schedP.skdsum);
    set(handles.checkbox_vie_sched_skd, 'Value', schedP.skd);


    % (3) ##### Selection of scheduling mode/strategy ##### 
    switch schedP.OPTIMIZATION

        % ### (3.1) Source based strategy ###
        case 1
            set(handles.radiobutton_vie_sched_rbsource, 'Value', 1);
            set(handles.radiobutton_vie_sched_rbstation, 'Value', 0);
            set(handles.radiobutton_vie_sched_satellites, 'Value', 0);
            set(handles.radiobutton_vie_sched_manually, 'Value', 0);
            set(handles.uipanel_vie_sched_sat_sched, 'Visible', 'off');
            set(handles.edit_vie_sched_srcnum,'string', num2str(schedP.SRCNUM));
            set(handles.text146, 'Enable', 'on');
            set(handles.edit_vie_sched_srcnum, 'Enable', 'on');
            set(handles.text147, 'Enable', 'on');
            set(handles.checkbox_vie_sched_distribute, 'Enable', 'Off');

        % ### (3.2) Station based strategy ###
        case 2
            set(handles.radiobutton_vie_sched_rbsource, 'Value', 0);
            set(handles.radiobutton_vie_sched_rbstation, 'Value', 1);
            set(handles.radiobutton_vie_sched_satellites, 'Value', 0);
            set(handles.radiobutton_vie_sched_manually, 'Value', 0);
            set(handles.uipanel_vie_sched_sat_sched, 'Visible', 'off');
            set(handles.text146, 'Enable', 'off');
            set(handles.edit_vie_sched_srcnum, 'Enable', 'off');
            set(handles.text147, 'Enable', 'off');
            set(handles.checkbox_vie_sched_distribute, 'Enable', 'on');
            set(handles.checkbox_vie_sched_distribute, 'Value', schedP.DISTRIBUTE);


        % ### (3.3) satellite scheduling ###
        case 3
            set(handles.radiobutton_vie_sched_rbsource, 'Value', 0);
            set(handles.radiobutton_vie_sched_rbstation, 'Value', 0);
            set(handles.radiobutton_vie_sched_satellites, 'Value', 1);
            set(handles.radiobutton_vie_sched_manually, 'Value', 0);
            set(handles.uipanel_vie_sched_sat_sched, 'Visible', 'on');
            set(handles.text146, 'Enable', 'off');
            set(handles.edit_vie_sched_srcnum, 'Enable', 'off');
            set(handles.text147, 'Enable', 'off');
            set(handles.checkbox_vie_sched_distribute, 'Enable', 'off');

            % ?? Add function to load a list of selected satellites here ??
            % Load from "../DATA/LEVEL5/satellite_str.mat"
        case 4
            set(handles.radiobutton_vie_sched_rbsource, 'Value', 0);
            set(handles.radiobutton_vie_sched_rbstation, 'Value', 0);
            set(handles.radiobutton_vie_sched_satellites, 'Value', 0);
            set(handles.radiobutton_vie_sched_manually, 'Value', 1);
            set(handles.uipanel_vie_sched_sat_sched, 'Visible', 'off');
            set(handles.text146, 'Enable', 'off');
            set(handles.edit_vie_sched_srcnum, 'Enable', 'off');
            set(handles.text147, 'Enable', 'off');
            set(handles.checkbox_vie_sched_distribute, 'Enable', 'off');
            set(handles.checkbox_vie_sched_distribute, 'Value', 0);
    end % switch schedP.OPTIMIZATION
        
    if schedP.openAnalyser == 1
        set(handles.checkbox_vie_sched_openSchedAnalyser, 'Value', 1);
    else
        set(handles.checkbox_vie_sched_openSchedAnalyser, 'Value', 0);
    end
	
    if schedP.disp_sky_coverage == 1
        set(handles.checkbox_vie_sched_disp_sky_coverage, 'Value', 1);
    else
        set(handles.checkbox_vie_sched_disp_sky_coverage, 'Value', 0);
    end
    
    
    if schedP.MULTISCHED == 1
        set(handles.checkbox_vie_sched_multiSched, 'Value', 1);
    else
        set(handles.checkbox_vie_sched_multiSched, 'Value', 0);
    end
	
    if schedP.SAVEOUTPUT == 1
        set(handles.checkbox_vie_sched_saveOutput, 'Value', 1);
    else
        set(handles.checkbox_vie_sched_saveOutput, 'Value', 0);
    end

    if ~isempty(schedP.optimization_condition)
        handles.checkbox_vie_sched_conditions.Value = 1;
        if strcmp(schedP.optimization_type,'quick fillin')
        	handles.popupmenu_vie_sched_conditions.Value = 1;
        else
            handles.popupmenu_vie_sched_conditions.Value = 2;
        end
        handles.text_vie_sched_conditions.String = schedP.optimization_condition;
        handles.popupmenu_vie_sched_conditions.Enable = 'on';
        handles.pushbutton_vie_sched_conditions.Enable = 'on';
        handles.text_vie_sched_conditions.Enable = 'on';
    else
        handles.checkbox_vie_sched_conditions.Value = 0;
        handles.popupmenu_vie_sched_conditions.Enable = 'off';
        handles.pushbutton_vie_sched_conditions.Enable = 'off';
        handles.tex_vie_sched_conditions.Enable = 'off';
        handles.text_vie_sched_conditions.String = '';
    end
	set(handles.edit_vie_sched_pathCatalogs,'String',schedP.pinfile)

end % function

  