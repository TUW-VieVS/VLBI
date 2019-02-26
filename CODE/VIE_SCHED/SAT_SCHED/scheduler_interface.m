% -------------------------------------------------------------------------
%
%                              scheduler_interface.m
%
%   This procedure handles the command window communication (in- and output) for satellite scheduling.
%
%   Author: 
%       Andreas Hellerschmied, 19.01.2014
%   
%   changes       : 
%   - 2015-06-29, A. Hellerschmied: Cable wrap sector can now be written to the VEX file
%   - 2015-07-08, A. Hellerschmied: Changed name from "manual_scheduling.m" to "scheduler_interface.m"
%   - 2015-08-20, A. Hellerschmied: Bug fix: input check of t_start for station based scheduling (input format <yyyy-mm-dd HH:MM:SS>) corrected! Crashed if there was no scan before.
%   - 2015-08-21, A. Hellerschmied: Functions added to save and load scheduling data via the interface, which gives the possibility to generate VEX files from saved schedules at a later time
%   - 2015-08-24, A. Hellerschmied: Added "source" and "stat_data" structures to output arguments. This is required, when existing scheduling data was loaded for VEX file creation.
%   - 2015-08-28, A. Hellerschmied: Added possibility to write a schedule summary to the Matlab CW or to a *.txt file
%   - 2015-10-25, A. Hellerschmied: Submenu for automtic scheduling of satellites and qusasars added
%   - 2015-11-26, A. Hellerschmied: Option added to create sky-plots, which visialize the creted schedule.
%   - 2016-04-??, A. Hellerschmied: Added possibility to change the station network (select sub-net) for sat. observations.  sched_data.stat(1).sat(31).obs_times
%   - 2016-04-28, A. Hellerschmied: Available obs. times per sat are now saved to "sched_data.stat(i_stat).sat(i_sat).obs_times" from "obs_data"
%   - 2016-04-28, A. Hellerschmied: Added to "sched_data": horizontal mask data, TLE data for each satellite, Ell. stat coordinates (long, lat, altitude)
%   - 2016-05-02, A. Hellerschmied: Added to "sched_data": quasar coordinates and names (fields: ra_rad, de_rad, name, info)
%   - 2016-05-03, A. Hellerschmied: Added functions to create sky plots for the data in "sched_data" (create_skyplot_sum.m)
%   - 2016-05-09, A. Hellerschmied: - Save sky_plot_sum figures as .fig and/or .pdf files
%                                   - Save sched_sum file and vievs2tie schedule file to the subdir. defined in PARA.OUTPUT_SUBDIR
%   - 2016-05-11, A. Hellerschmied: - Update sky/elev. plots at end of auto_sched
%                                   - Added subroutine to check time epoch input (check_t_input_19char)
%                                   - Possibility to start auto_sched for defined time windows only
%   - 2016-06-16, A. Hellerschmied: Various changes applied to main menu:
%                                   - Sub-menu to load/save schedule data
%                                   - Sub-men to generate output
%   - 2016-06-20, A. Hellerschmied: When laoding a session with the "PARA" structure, the current sub-directory set in the GUI is kept. 
%   - 2016-06-21, A. Hellerschmied: More options to change the station network for satellite observations (recover orig. station list, etc.)
%   - 2016-06-29: A. Hellerschmied: - Errors at inputting obs. end times via time tags fixed.
%                                   - Option for debug output added (Set PARA.DEBUG_FLAG = true) 
%   - 2016-07-13: A. Hellerschmied: Output earliest possible scan start time to CW, when the entered scan start time in not valid.
%   - 2016-08-02: A. Hellerschmied: Added the "commandwindow" command in front of all "input" commands to set the CW active (was not done automatically any more in MATLAB R2016a)
%   - 2016-10-28: A. Hellerschmied: Changes required for new "station_based_quasar_scheduling.m" function
%   - 2016-11-02: A. Hellerschmied: Option added to select alternative cable wrap section for defined antennas manually
%   - 2016-11-14: A. Hellerschmied: - Option added to disable scan start time checks for single scans (Only for the next scan!!!)
%                                   - station_based_quasar_scheduling.m => The user is now asked, if the scans should be added to the schedule
%   - 2016-11-25: A. Hellerschmied: - Option added to set Time constants added between slewing and obs. start to zero
%   - 2016-11-25: A. Hellerschmied: - New option to change scan start time for quasars manually
%   - 2016-11-30: A. Hellerschmied: - New option to save sched data to LEVEL5 dir. using version number
%   - 2017-01-23: A. Hellerschmied: - New function added to write VSO files (write_vso_file.m)
%   - 2017-02-07: A. Hellerschmied: - Option added to write VSO files with geocentric baselines 
%   - 2018-12-20: A. Corbin       : - highlight selected satellite
%									- empty string in input of scan start time does no longer lead to crash
%                                   - new option for start time of new satellite scan (peak)
%									- default values for some user inputs
%                                                                
%            
%
%   inputs        :
%   - obs_data          : observation data structure used during scheduling process (temp data is stored there)
%   - stat_data         : station data structure (for VieVS satellite scheduling)
%   - source            : Source (quasar) data structure
%   - sched_data        : Scheduling data structure (for VieVS satellite scheduling)
%   - sched_handles     : Figure handles structure for VieVS satellite scheduling
%   - station           : Station structure (stations defined via GUI "General setup")
%   - PARA              : General scheduling parameter structure
%   - station           : station data structure from VIE_SCHED
%
%   outputs       :
%   - error_code        : Error Code (0 = no erros occured)
%   - error_msg         : Error Message (empty, if no errors occured)
%   - sched_data        : Scheduling data structure (for VieVS satellite scheduling) - updated!
%   - stat_data         : station data structure (for VieVS satellite scheduling)   - Updated, if existing scheduling data was loaded from /DATA/LEVEL5/<experiment_name>/
%   - source            : Source (quasar) data structure                            - Updated, if existing scheduling data was loaded from /DATA/LEVEL5/<experiment_name>/
% 
%
% 
%
%   coupling      :
%   - update_sky_plots_pointer
%   - update_elevation_plot_pointer
%   - update_sky_plots
%   - update_sky_plots_highlight_quasars
%   - check_t_start
%   - calc_t_start
%   - invjday
%   - check_scan_duration_and_calc_repos_epochs
%   - check_t_end
%   - save_scan
%   - preselect_observable_quasars
%   - find_closest_quasars_and_calc_scan
%   - select_quasar_in_skyplot
%   - calc_quasar_scan
%   - station_based_quasar_scheduling
%   - write_string_to_file
%   - write_sched_sum_str
%   - create_skyplot_sum
%   - vie_sched_sat_o.m
%   - write_vso_file.m
%
%   sub-routines    :
%   - check_t_input_19char
%   
%
%   references    :
%
%-------------------------------------------------------------------------



function [sched_data, stat_data, source, error_code, error_msg] = scheduler_interface(obs_data, stat_data, source, sched_data, sched_handles, PARA, station)

% ##### Init #####
error_code = 0;
error_msg = '';
flag_sat_network_id_list_changed = 0;


% ##### Constants #####
mjd2jd = 2.400000500000000e+006;


% ##### Set options #####
filepath_sched_sum_file_out = ['../DATA/SCHED/', PARA.OUTPUT_SUBDIR, '/'];
filepath_vso_out            = ['../DATA/SCHED/', PARA.OUTPUT_SUBDIR, '/'];

% Debug output:
PARA.DEBUG_FLAG = false;

PARA.CHECK_T_START      = true;
PARA.DISABLE_SLEW_CONST = false;


% ##### Prepare station lists #####
% List of station ID, referring to stat_data.stat fields
station_id_list = 1 : stat_data.number_of_stations;
i_quasar_stat = 0;
i_sat_stat = 0;
for i_stat = 1 : stat_data.number_of_stations 

    switch(stat_data.stat(i_stat).obs_type)
       
        % Quasar network:
        case{'quasar_only', 'twin_quasar'}
            i_quasar_stat = i_quasar_stat + 1;
            quasar_network_id_list(i_quasar_stat) = i_stat;
            
        % Satellite network:
        case{'sat_only', 'twin_sat'}
            i_sat_stat = i_sat_stat + 1;
            sat_network_id_list(i_sat_stat) = i_stat;
            quasar_network_id_list = [];
            
        % Satellite & Quasar network:
        case{'sat_and_quasar'}
            i_sat_stat = i_sat_stat + 1;
            sat_network_id_list(i_sat_stat) = i_stat;
            i_quasar_stat = i_quasar_stat + 1;
            quasar_network_id_list(i_quasar_stat) = i_stat;
            
        % Error case
        otherwise
            error_code = 1;
            error_msg = ['Observation type of station ', stat_data.stat(i_stat).name, ' is unknown.'];
            return;
    end
    
    % ##### Save station properties to "sched_data" structure #####
    % This data is needed for the generation of VEX files later on (in function "write_vex_file.m")
    sched_data.stat(i_stat).name                    = stat_data.stat(i_stat).name;                      % [8 char]
    sched_data.stat(i_stat).trf_x                   = stat_data.stat(i_stat).location.TRF.x;            % [m]        
    sched_data.stat(i_stat).trf_y                   = stat_data.stat(i_stat).location.TRF.y;            % [m] 
    sched_data.stat(i_stat).trf_z                   = stat_data.stat(i_stat).location.TRF.z;            % [m] 
    sched_data.stat(i_stat).trf_vx                  = stat_data.stat(i_stat).location.TRF.vx;           % [m/year]
    sched_data.stat(i_stat).trf_vy                  = stat_data.stat(i_stat).location.TRF.vy;           % [m/year]
    sched_data.stat(i_stat).trf_vz                  = stat_data.stat(i_stat).location.TRF.vz;           % [m/year]
    sched_data.stat(i_stat).trf_epoch               = stat_data.stat(i_stat).location.TRF.epoch;        % [MJD]
    sched_data.stat(i_stat).trf_label               = stat_data.stat(i_stat).location.TRF.TRF_label;    % [string]
    sched_data.stat(i_stat).long                    = stat_data.stat(i_stat).location.ellipsoid.long;   % [deg]
    sched_data.stat(i_stat).lat                     = stat_data.stat(i_stat).location.ellipsoid.lat;    % [deg]
    sched_data.stat(i_stat).altitude                = stat_data.stat(i_stat).location.ellipsoid.altitude; % [m]
    sched_data.stat(i_stat).label                   = stat_data.stat(i_stat).label;                     % [2 char]
    sched_data.stat(i_stat).axis_type               = stat_data.stat(i_stat).axis_type;                 % [4 char]
    sched_data.stat(i_stat).obs_type                = stat_data.stat(i_stat).obs_type;                  % [string]
    sched_data.stat(i_stat).axis_1_slew_rate        = stat_data.stat(i_stat).max_axis1_rate * 60;       % [deg/min]  
    sched_data.stat(i_stat).axis_2_slew_rate        = stat_data.stat(i_stat).max_axis2_rate * 60;       % [deg/min]  
    sched_data.stat(i_stat).axis_1_settling_time    = stat_data.stat(i_stat).c1;                        % [sec] 
    sched_data.stat(i_stat).axis_2_settling_time    = stat_data.stat(i_stat).c2;                        % [sec] 
    sched_data.stat(i_stat).axis_1_acc              = stat_data.stat(i_stat).max_axis1_acc;             % [deg/s^2]
    sched_data.stat(i_stat).axis_2_acc              = stat_data.stat(i_stat).max_axis2_acc;             % [deg/s^2]
    sched_data.stat(i_stat).axis_offset             = stat_data.stat(i_stat).axis_offset;               % [m] 
    sched_data.stat(i_stat).az_n1                   = stat_data.stat(i_stat).az_n1;                     % Cable wrap neutral pointing sector limit 1 [rad]
    sched_data.stat(i_stat).az_n2                   = stat_data.stat(i_stat).az_n2;                     % Cable wrap neutral pointing sector limit 2 [rad]
    sched_data.stat(i_stat).az_cw1                  = stat_data.stat(i_stat).az_cw1;                    % Cable wrap clockwise pointing sector limit 1 [rad]
    sched_data.stat(i_stat).az_cw2                  = stat_data.stat(i_stat).az_cw2;                    % Cable wrap clockwise pointing sector limit 2 [rad]
    sched_data.stat(i_stat).az_ccw1                 = stat_data.stat(i_stat).az_ccw1;                   % Cable wrap counter-clockwise pointing sector limit 1 [rad]
    sched_data.stat(i_stat).az_ccw2                 = stat_data.stat(i_stat).az_ccw2;                   % Cable wrap counter-clockwise pointing sector limit 2 [rad]
    sched_data.stat(i_stat).az_ccw1                 = stat_data.stat(i_stat).az_ccw1;                   % Cable wrap counter-clockwise pointing sector limit 1 [rad]
    sched_data.stat(i_stat).az_ccw2                 = stat_data.stat(i_stat).az_ccw2;                   % Cable wrap counter-clockwise pointing sector limit 2 [rad]
    sched_data.stat(i_stat).lim11                   = stat_data.stat(i_stat).lim11;                     % Axis 1, limit 1 (Limit from antenna.cat incl. marge from param.txt)[rad]
    sched_data.stat(i_stat).lim12                   = stat_data.stat(i_stat).lim12;                     % Axis 1, limit 2 (Limit from antenna.cat incl. marge from param.txt)[rad]
    sched_data.stat(i_stat).lim21                   = stat_data.stat(i_stat).lim21;                     % Axis 2, limit 1 (Limit from antenna.cat incl. marge from param.txt)[rad]
    sched_data.stat(i_stat).lim22                   = stat_data.stat(i_stat).lim22;                     % Axis 2, limit 2 (Limit from antenna.cat incl. marge from param.txt)[rad]
    sched_data.stat(i_stat).horizontal_mask         = stat_data.stat(i_stat).horizontal_mask;           % horizontal mask
    sched_data.stat(i_stat).horizontal_mask_num     = stat_data.stat(i_stat).horizontal_mask_num;       % horizontal mask: Number of elements
    
    
    % ##### Save available "obs_times" of satellites for particular stations to "sched_data" struct #####
    for i_sat = 1 : length(obs_data.sat)
        if(strcmp(sched_data.stat(i_stat).obs_type, 'sat_only') || strcmp(sched_data.stat(i_stat).obs_type, 'twin_sat') || strcmp(sched_data.stat(i_stat).obs_type, 'sat_and_quasar'))
            sched_data.stat(i_stat).sat(i_sat).obs_times = obs_data.sat(i_sat).stat(i_stat).obs_times;
        else
            sched_data.stat(i_stat).sat(i_sat).obs_times = [];
        end
    end

end % for i_stat = 1 : stat_data.number_of_stations 

sat_network_id_list_orig = sat_network_id_list;


% ##### Save satellite data to "sched_data" struct #####
% stored at "stat_data" in calc_pointing_data.m => Same TLE info for every station!
for i_sat = 1 : length(stat_data.stat(1).sat)
    sched_data.sat(i_sat).TLE_data.TLE_header_line  = stat_data.stat(1).sat(i_sat).TLE_data.TLE_header_line;
    sched_data.sat(i_sat).TLE_data.sat_number       = stat_data.stat(1).sat(i_sat).TLE_data.sat_number;
end

% ##### Save quasar data to "sched_data" struct #####
for i_quasar = 1 : length(source)
    sched_data.quasar(i_quasar).name    = source(i_quasar).name;
    sched_data.quasar(i_quasar).ra_rad  = source(i_quasar).ra;
    sched_data.quasar(i_quasar).de_rad  = source(i_quasar).de;
    sched_data.quasar(i_quasar).info    = source(i_quasar).info;
end


% ##### Preallocation ##### 

% "obs_data.source"
if (length(obs_data.quasars) == 1) && isempty(obs_data.quasars(1).number_of_scans)
    for i = 1 : length(source)
       obs_data.quasars(i).last_obs_jd      = 0;
       obs_data.quasars(i).number_of_scans  = 0;
    end
end
% "obs_data.sat" 
if (length(obs_data.sat) == 1) && isempty(obs_data.sat(1).number_of_scans)
    for i = 1 : length(obs_data.sat)
       obs_data.sat(i).last_obs_jd          = 0;
       obs_data.sat(i).number_of_scans      = 0;
    end
end

fprintf(1, '\n');
fprintf(1, '###### Scheduler interface for VieVS satellite scheduling ######\n');
fprintf(1, '\n');

        

% ##### Main loop #####

% loop init:
sched_state = 0;
flag_exit_user_input = 0;

add_scan_state = 0;
flag_exit_add_scan = 0;

while(~flag_exit_user_input)

    switch(sched_state)
        
        
        % ##################################
        % ##### Choose Experiment Name #####
        % ##################################
        case 0 
            fprintf(1, '#### Type in an Experiment Name ####\n');
            fprintf(1,'  => Input Length: Between 1 and 4 characters\n');
            fprintf(1,'  => Legal characters: "A-Z", "a-z", "0-9", "_" and "-"\n');
            
            commandwindow
            input_str = input(' Experiment name: ', 's');
            len_input_str = length(input_str);
            
            if (    (   (len_input_str >= 1) && (len_input_str <= 4) )                                          && ...
                    (   (sum((input_str(1:len_input_str) <= '9') & (input_str(1:len_input_str) >= '0')) + ...
                         sum((input_str(1:len_input_str) <= 'z') & (input_str(1:len_input_str) >= 'a')) + ...
                         sum((input_str(1:len_input_str) <= 'Z') & (input_str(1:len_input_str) >= 'A')) + ...
                         sum(input_str(1:len_input_str) == '_') + ...
                         sum(input_str(1:len_input_str) == '-') ) == len_input_str)                             && ...
                    (~strcmp(input_str, 'exit'))                                                                && ...
                    (~strcmp(input_str, 'help'))                                                                    )
                         
                
                sched_data.exper_name = input_str;
                sched_state = 1;
                fprintf(1,'\n');
        

            elseif (strcmp(input_str, 'exit'))
                sched_state = 6;
            elseif (strcmp(input_str, 'help'))
                disp(' ');
                disp('  => Type in an Experiment Name:');
                disp('     Input Length: Between 1 and 4 characters');
                disp('     Legal characters: "A-Z", "a-z", "0-9", "_" and "-"');
                disp(' ');
                disp('  => exit - Cancel operator input and exit satallite schedule procedure.');
                disp(' ');
            else
                disp('  => ERROR: Invalid input! Type "help" to see input options.')
            end
            
            
           
        % #######################################
        % ##### Main menu: Choose an action #####
        % #######################################
        case 1 % Main Menu 
            fprintf(1, '#### Main menu ####\n');
            fprintf(1, 'use "back" and "exit" to escape from submenus \n');
            fprintf(1, 'Choose an action\n\n');
            fprintf(1, ' 1   - Add a scan to the the current schedule (append)\n');
            fprintf(1, ' 2   - Create output\n');
            fprintf(1, ' 3   - Load/save data\n');
            fprintf(1, ' 4   - Edit current schedule\n');
            fprintf(1, ' 6   - Create schedule summary and statistics\n');
            fprintf(1, ' 7   - Automatic scheduling (satellites and quasars)\n');
            fprintf(1, ' 8   - Change settings\n');
            fprintf(1, ' 9   - Exit\n');
            fprintf(1,'\n');
%             fprintf(1, ' 2   - Finish user input and create VEX file\n');
%             fprintf(1, ' 3   - Create schedule file in the VieVS2tie format\n');
%             fprintf(1, ' 4   - Save scheduling data to /DATA/LEVEL5/%s/\n', sched_data.exper_name);
%             fprintf(1, ' 5   - Load scheduling data from /DATA/LEVEL5/, create VEX file and exit\n');
%             fprintf(1, ' 6   - Create schedule summary and statistics\n');
%             fprintf(1, ' 7   - Automatic scheduling (satellites and quasars)\n');
%             fprintf(1, ' 8   - Change settings\n');
%             fprintf(1, ' 9   - Exit\n');
%             fprintf(1,'\n');
            sched_state = 1.1;
            
        case 1.1 % Main Menu input
            commandwindow
            input_str = input(' Please select: ', 's');
            switch(input_str)
                case '1'                % Add a scan to the the current schedule
                    sched_state = 2;
                case '2'                % Create output
                    sched_state = 5;
                case '3'                % Load/save data
                    sched_state = 8;
                case '4'                % Edit current schedule
                    sched_state = 1;
                    fprintf(1, '\nSorry, this function is not available by now!\n')
%                 case '5'
%                     sched_state = 9;
                case '6'                % Create schedule summary and statistics
                    sched_state = 10;
                case '7'                % Automatic scheduling
                    sched_state = 11;
                case '8'                % Change settings
                    sched_state = 20;
                case {'9', 'exit'}      % Exit
                    sched_state = 6;
                otherwise
                    fprintf(1, ' ERROR: Invalid input, try again!\n');
            end
            
            
            
        % #################################################  
        % ##### Add scan to current Schedule (append) #####
        % #################################################
        case 2 
            
            % Init.:
            flag_exit_add_scan = 0;
            add_scan_state = 0;
            
            while(~flag_exit_add_scan)
                
                switch(add_scan_state)
                    
                    case 0
                        fprintf(1, '#### Add scan to current Schedule (append) ####\n');
                        fprintf(1, ' 1   - Add a satellite scan\n');
                        fprintf(1, ' 2   - Add a quasar scan\n');
                        % fprintf(1, ' 3   - Add satellite scan with subnetting\n');
                        
                        [sched_handles, error_code, error_msg] = update_sky_plots_pointer(sched_handles, station_id_list, obs_data, 4, PARA);
                        if error_code > 0
                            error_msg = ['update_sky_plots_pointer: ', error_msg];
                            fprintf(1, [' ERROR: ', error_msg, '\n']);
                            add_scan_state = 0;
                        else 
                            [sched_handles, error_code, error_msg] = update_elevation_plot_pointer(sched_handles, obs_data, 4, [], PARA, station_id_list);
                            if error_code > 0
                                error_msg = ['update_elevation_plot_pointer: ', error_msg];
                                fprintf(1, [' ERROR: ', error_msg, '\n']);
                                add_scan_state = 0;
                            else
                                if isempty(obs_data.end_of_last_scan_jd)
                                    t_epoch_jd_temp = PARA.startmjd + 2.400000500000000e+006; % first scan of session
                                else
                                    t_epoch_jd_temp = obs_data.end_of_last_scan_jd;
                                end
                                flag_plot_sources = 0;
                                [sched_handles, error_code, error_msg] = update_sky_plots(sched_handles, stat_data, station_id_list, source, PARA, obs_data, flag_plot_sources, t_epoch_jd_temp);
                                if error_code > 0
                                    error_msg = ['update_sky_plots: ', error_msg];
                                    fprintf(1, [' ERROR: ', error_msg, '\n']);
                                    add_scan_state = 0;
                                else
                                    % Plot quasar markers (highlight):
                                    [sched_handles, error_code, error_msg] = update_sky_plots_highlight_quasars(sched_handles, stat_data, station_id_list, source, t_epoch_jd_temp, [], 2);
                                    if error_code > 0
                                         error_msg = ['update_sky_plots_highlight_quasars: ', error_msg];
                                         fprintf(1, [' ERROR: ', error_msg, '\n']);
                                         add_scan_state = 0;
                                    else
                                        add_scan_state = 0.1;
                                    end
                                end
                            end
                        end
                        
                        
                    case 0.1 % Choose satellite/quasar as target
                        commandwindow
                        input_str = input(' Please select: ', 's');
                        switch(input_str)
                            case '1'
                                add_scan_state = 1; % satellite
                            case '2'
                                add_scan_state = 2; % quasar
%                             case '3'
%                                 add_scan_state = 4; % add satellite scan with subnetting 
                            case 'exit'
                                add_scan_state = 9;
                            case 'back'
                                add_scan_state = 9;
                            otherwise
                                fprintf(1, ' ERROR: Invalid input, please try again!\n');
                        end
                        
                        
                        
                    % ##########################
                    % ##### Satellite scan #####
                    % ##########################
                    
                    % ##### Choose a satellite scan #####
                    case 1 
                        fprintf(1, ' \n++++ Choose a satellite: ++++ \n');
                        add_scan_state = 1.1;
                        
                    case 1.1 % Choose a satellite
                        commandwindow
                        input_str = input(' Satellite number : ', 's');
                        [satellite_id, status] = str2num(input_str);
                        if (status == 1) % Conversion successfull? yes => numeric value
                            if (satellite_id <= obs_data.number_of_sat)
                                    sat_name_str = obs_data.sat(satellite_id).name(1 : end-2);
                                    fprintf(1, ' Choosen satellite: %s\n\n', sat_name_str);
                                    add_scan_state = 1.3; % Got satellite name!
                                    
                                    % Highligh selected stellite
                                    highlight_satellite( sched_handles, station_id_list, satellite_id, true );
                                    
                            else
                                fprintf(1, ' ERROR: Entered satellite number is not available. Try another one!\n');
                            end
                        else % not a numeric value
                            switch(input_str)
                                case 'back'
                                    add_scan_state = 0;
                                case 'exit'
                                    add_scan_state = 9;
                                otherwise
                                    fprintf(1, ' ERROR: Invalid input, try again!\n');
                            end
                        end
                        
                    case 1.3 % Scan start time: how to enter it?
                        fprintf(1, ' \n++++ Scan start time: ++++ \n');
                        fprintf(1, ' 1   - Enter it manually\n');
                        fprintf(1, ' 2   - Next possible start\n');
                        fprintf(1, ' 3   - Next peak\n');
                        add_scan_state = 1.31;
                        % Init.:
                        t_start_jd_temp = 0;
                        t_start_jd = 0;
                        t_start_jd_min = 0;
                        
                    case 1.31 % Scan start time: how to enter it?
                        commandwindow
                        input_str = input(' Please select: ', 's');
                        switch(input_str)
                            case '1'
                                add_scan_state = 1.4; % Enter it manually
                            case '2'
                                add_scan_state = 1.5; % Next possible start
                            case '3'
                                add_scan_state = 1.55; % Next peak
                            case 'exit'
                               highlight_satellite( sched_handles, station_id_list, satellite_id, false );
                                add_scan_state = 9;
                            case 'back'
                               highlight_satellite( sched_handles, station_id_list, satellite_id, false );
                                add_scan_state = 1;
                            otherwise
                                fprintf(1, ' ERROR: Invalid input, please try again!\n');
                        end

                    case 1.4 % #### Scan start time: Enter it manually ####
                        fprintf(1, ' \n++++ Enter the start date/time of the scan - t_start ++++\n');
                        fprintf(1,'  => Input Format: <Time Tag Label>  or <yyyy-mm-dd HH:MM:SS>\n');
                        add_scan_state = 1.41;
                        
                    case 1.41 % Scan start time: Enter it manually
                        commandwindow
                        input_str = input(' t_start = ', 's');
                        len_input_str = length(input_str);
                                    
                        if isempty(len_input_str) || len_input_str == 0
                             fprintf(1, ' ERROR: Ivalid input. Please try again!\n');
                        % <Time Tag label>
                        elseif ((input_str(1) == 't') && (len_input_str >= 2) && ...  % Check input Format!
                            (sum((input_str(2:len_input_str) <= '9') & (input_str(2:len_input_str) >= '0')) == (len_input_str - 1)) )
                            
                            % Search <Time Tag label> in obs_data.sat.obs_times get the according time [JD]
                            t_num = str2double(input_str(2:len_input_str));
                            t_start_jd_temp = obs_data.sat(satellite_id).obs_times(logical(obs_data.sat(satellite_id).obs_times(:, 3) == t_num), 1);
                            if isempty(t_start_jd_temp)
                                fprintf(1, ' ERROR: Time tag label cannot be found!\n');
                            else
                                add_scan_state = 1.42;
                            end
                                    
                        % <yyyy-mm-dd HH:MM:SS>  
                        elseif (len_input_str == 19)
                            if (    (sum((input_str(1:4) <= '9') & (input_str(1:4) >= '0')) == 4)       && ... % Check input Format!
                                    (sum((input_str(6:7) <= '9') & (input_str(6:7) >= '0')) == 2)       && ...
                                    (sum((input_str(9:10) <= '9') & (input_str(9:10) >= '0')) == 2)     && ...
                                    (sum((input_str(12:13) <= '9') & (input_str(12:13) >= '0')) == 2)   && ...
                                    (sum((input_str(15:16) <= '9') & (input_str(15:16) >= '0')) == 2)   && ...
                                    (sum((input_str(18:19) <= '9') & (input_str(18:19) >= '0')) == 2)   && ...
                                    (input_str(5) == '-')                                               && ...
                                    (input_str(8) == '-')                                               && ...
                                    (input_str(11) == ' ')                                              && ...
                                    (input_str(14) == ':')                                              && ...
                                    (input_str(17) == ':')                                                      )   
                                
                                % Convert input to JD:
                                yr = str2num(input_str(1:4)); %#ok<*ST2NM>
                                mon = str2num(input_str(6:7));
                                day = str2num(input_str(9:10));
                                hr = str2num(input_str(12:13));
                                min = str2num(input_str(15:16));
                                sec = str2num(input_str(18:19));
                                t_start_jd_temp = jday(yr, mon, day, hr, min, sec);
                                add_scan_state = 1.42;
                            else
                                fprintf(1, ' ERROR: Ivalid input. Please try again!\n');
                            end
                            
                        % other commands    
                        else
                            switch(input_str)
                                case 'exit'
                                   highlight_satellite( sched_handles, station_id_list, satellite_id, false );
                                    add_scan_state = 9;
                                case 'back'
                                    add_scan_state = 1.3;
                                case 'help'
                                    disp(' ');
                                    disp('  => Choose t_start: Two input formats are possible:');
                                    disp('      (a) <Time Tag Label> , e.g. "t2"');
                                    disp('          The input time tag label has to correspond with an available label of the choosen satellite.');
                                    disp('      (b) <yyyy-mm-dd HH:MM:SS> , e.g. "2000-03-01 15:45:17"');
                                    disp('          The input date and time has to be within an available observation periode!');
                                    disp(' ');

                                otherwise
                                    fprintf(1, ' ERROR: Invalid input. Type "help" to see input options!\n');
                            end
                        end

                        
                    % Check, if t_start is OK:
                    case 1.42 
                        source_quasar = [];
                        if PARA.CHECK_T_START
                        	[flag_t_start_ok, t_start_jd_min, obs_data, error_code, error_msg] = check_t_start(stat_data, sat_network_id_list, obs_data, PARA, source_quasar, satellite_id, t_start_jd_temp, flag_sat_network_id_list_changed);
                            if error_code > 0
                                error_msg = ['check_t_start: ', error_msg];
                                fprintf(1, [' ERROR: ', error_msg, '\n']);
                                add_scan_state = 0;
                            end
                        else % Do not check t_start_jd_temp, just calc. un_az at begin of scan
                            fprintf('==> WARNING: Checks for scan start time is disabled manually for this scan only! <==\n')
                            PARA.CHECK_T_START = true;
                            fprintf(' - Scan start checks are re-activated now!\n')
                            [obs_data, error_code, error_msg] = calc_un_az_begin_of_scan(stat_data, sat_network_id_list, obs_data, PARA, source_quasar, satellite_id, t_start_jd_temp);
                            if error_code > 0
                                error_msg = ['calc_un_az_begin_of_scan: ', error_msg];
                                fprintf(1, [' ERROR: ', error_msg, '\n']);
                                add_scan_state = 0;
                            end
                            flag_t_start_ok = 1;
                        end
                        if error_code == 0
                            if flag_t_start_ok % NEXT STATE! t_start found!
                                % Update graphics:
                               highlight_satellite( sched_handles, station_id_list, satellite_id, false );
                                [sched_handles, error_code, error_msg] = update_sky_plots_pointer(sched_handles, sat_network_id_list, obs_data, 2, PARA);
                                if error_code > 0
                                    error_msg = ['update_sky_plots_pointer: ', error_msg];
                                    fprintf(1, [' ERROR: ', error_msg, '\n']);
                                    add_scan_state = 0;
                                else
                                    flag_plot_sources = 0;
                                    [sched_handles, error_code, error_msg] = update_sky_plots(sched_handles, stat_data, sat_network_id_list, source, PARA, obs_data, flag_plot_sources, t_start_jd_temp);
                                    if error_code > 0
                                        error_msg = ['update_sky_plots: ', error_msg];
                                        fprintf(1, [' ERROR: ', error_msg, '\n']);
                                        add_scan_state = 0;
                                    else
                                        [sched_handles, error_code, error_msg] = update_elevation_plot_pointer(sched_handles, obs_data, 2, [], PARA, sat_network_id_list);
                                        if error_code > 0
                                            error_msg = ['update_elevation_plot_pointer: ', error_msg];
                                            fprintf(1, [' ERROR: ', error_msg, '\n']);
                                            add_scan_state = 0;
                                        else
                                            % Assign final t_start:
                                            t_start_jd = t_start_jd_temp;
                                            add_scan_state = 1.6; 
                                        end
                                    end
                                end
                            else
                                fprintf(1, ' ERROR: Scan start time is not valid. Try again!\n');
                                if t_start_jd_min ~= 0
                                    fprintf(1, '  => Earliest possible scan start: %s\n', jd2datestr(t_start_jd_min));
                                end
                                add_scan_state = 1.4;
                            end
                        end

                        
                        
                    % #### Scan start time: Next possible start ####
                    case  1.5 
                        fprintf(1, ' \n++++ Calculating the earliest possible scan start time ++++\n');
                        source_quasar = [];
                        [t_start_jd_temp, flag_t_start_found, obs_data, error_code, error_msg] = calc_t_start(stat_data, sat_network_id_list, obs_data, PARA, source_quasar, satellite_id);                        
                        if error_code > 0
                            error_msg = ['calc_t_start: ', error_msg];
                            fprintf(1, [' ERROR: ', error_msg, '\n']);
                            add_scan_state = 0;
                        else
                            if flag_t_start_found
                                % Print scan start time:
                                [year, mon, day, hr, min, sec] = invjday (t_start_jd_temp); % Conversion
                                fprintf(1, ' Calculated t_start_min: %4.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\n', year, mon, day, hr, min, sec);    % <yyyy-mm-dd HH:MM:SS> 
                                % Update graphics:
                                flag_plot_sources = 0;
                                [sched_handles, error_code, error_msg] = update_sky_plots(sched_handles, stat_data, sat_network_id_list, source, PARA, obs_data, flag_plot_sources, t_start_jd_temp);
                                if error_code > 0
                                    error_msg = ['update_sky_plots: ', error_msg];
                                    fprintf(1, [' ERROR: ', error_msg, '\n']);
                                    add_scan_state = 0;
                                else
                                    [sched_handles, error_code, error_msg] = update_elevation_plot_pointer(sched_handles, obs_data, 9, t_start_jd_temp, PARA, sat_network_id_list);
                                    if error_code > 0
                                        error_msg = ['update_elevation_plot: ', error_msg];
                                        fprintf(1, [' ERROR: ', error_msg, '\n']);
                                        add_scan_state = 0;
                                    else
                                        add_scan_state = 1.51;
                                    end
                                end
                            else
                                fprintf(1, ' ERROR: It is not possible to calculate a valid scan start time for the current setup!\n');
                                add_scan_state = 1.3;
                            end
                        end
                        
                    case  1.51 % Scan start time: Calculate it; input
                        commandwindow
                        input_str = input(' Take it? (y=yes, n=no, 9=back to main) <y>: ', 's');
                        switch(input_str)
                            case {'n' 'no'}
                                add_scan_state = 1.3; % Scan start time: how to enter it?
                            case {'9'}
                                highlight_satellite( sched_handles, station_id_list, satellite_id, false );
                                add_scan_state = 9;
                            case {'y' 'yes', ''}
                                [sched_handles, error_code, error_msg] = update_sky_plots_pointer(sched_handles, sat_network_id_list, obs_data, 2, PARA);
                                if error_code > 0
                                    error_msg = ['update_sky_plots_pointer: ', error_msg];
                                    fprintf(1, [' ERROR: ', error_msg, '\n']);
                                    add_scan_state = 0;
                                else
                                    [sched_handles, error_code, error_msg] = update_elevation_plot_pointer(sched_handles, obs_data, 2, t_start_jd_temp, PARA, sat_network_id_list);
                                    if error_code > 0
                                        error_msg = ['update_elevation_plot_pointer: ', error_msg];
                                        fprintf(1, [' ERROR: ', error_msg, '\n']);
                                        add_scan_state = 0;
                                    else
                                        flag_plot_sources = 0;
                                        [sched_handles, error_code, error_msg] = update_sky_plots(sched_handles, stat_data, sat_network_id_list, source, PARA, obs_data, flag_plot_sources, t_start_jd_temp);
                                        if error_code > 0
                                            error_msg = ['update_sky_plots: ', error_msg];
                                            fprintf(1, [' ERROR: ', error_msg, '\n']);
                                            add_scan_state = 0;
                                        else
                                            t_start_jd = t_start_jd_temp;
                                            add_scan_state = 1.6; % #### Scan stop time ####
                                        end
                                    end
                                end
                            case 'exit'
                                highlight_satellite( sched_handles, station_id_list, satellite_id, false );
                                add_scan_state = 9;
                            case 'back'
                                add_scan_state = 1.3;
                            otherwise
                                fprintf(1, ' ERROR: Invalid input, please try again!\n');
                        end
                    % #### Scan start time: Next peak ####
                    case 1.55
                        commandwindow
                        offset = str2num( input(' offset to peak in sec <0>: ', 's') );
                        if( isempty(offset) || ~isnumeric(offset) )
                            offset = 0;
                        end
                        offset = offset/86400;
                        
                        peaksta = str2num( input(' station <1>: ', 's') );
                        if( isempty(peaksta) || ~isnumeric(peaksta) )
                            peaksta = 1;
                        end
                        
                        if peaksta < 1 || peaksta > stat_data.number_of_stations
                            fprintf(1, ' ERROR: invalid station id \n');
                             add_scan_state = 1.55;
                        else
                            fprintf(1, '   %1.0f - %s \n', i_stat, stat_data.stat(peaksta).name);
                            source_quasar = [];
                            [t_earliest, obs_data, error_code, error_msg] = calc_start_of_scan(stat_data, station_id_list, obs_data, PARA, source_quasar, satellite_id);
                            t_start_jd_temp = nan;
                            fprintf(1, ' \n++++ Calculating peak ++++\n');
                            for i_overp = 1 : stat_data.stat(1).sat(satellite_id).number_of_overpasses

                                for i_peak = 1 : stat_data.stat(peaksta).sat(satellite_id).overpass(i_overp).number_of_peaks
                                    t_peak = stat_data.stat(peaksta).sat(satellite_id).overpass(i_overp).peak(i_peak).jd;
                                    if t_peak - offset > t_earliest 
                                        t_start_jd_temp = t_peak - offset;
                                        break;
                                    end
                                end
                                if ~isnan(t_start_jd_temp)
                                    break
                                end
                            end

                            if (~isnan(t_start_jd_temp) )

                                 [obs_data, error_code, error_msg] = calc_un_az_begin_of_scan(stat_data, station_id_list, obs_data, PARA, source_quasar, satellite_id, t_start_jd_temp);
                                 [year, mon, day, hr, min, sec] = invjday (t_start_jd_temp); % Conversion
                                 fprintf(1, ' Calculated t_start_min: %4.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\n', year, mon, day, hr, min, sec);
                                 % Update graphics:
                                    flag_plot_sources = 0;
                                    [sched_handles, error_code, error_msg] = update_sky_plots(sched_handles, stat_data, sat_network_id_list, source, PARA, obs_data, flag_plot_sources, t_start_jd_temp);
                                    if error_code > 0
                                        error_msg = ['update_sky_plots: ', error_msg];
                                        fprintf(1, [' ERROR: ', error_msg, '\n']);
                                        add_scan_state = 0;
                                    else
                                        [sched_handles, error_code, error_msg] = update_elevation_plot_pointer(sched_handles, obs_data, 9, t_start_jd_temp, PARA, sat_network_id_list);
                                        if error_code > 0
                                            error_msg = ['update_elevation_plot: ', error_msg];
                                            fprintf(1, [' ERROR: ', error_msg, '\n']);
                                            add_scan_state = 0;
                                        else
                                            add_scan_state = 1.51; % #### Scan stop time ####
                                        end
                                    end
                            else
                                fprintf(1, ' ERROR: No peak found\n');
                                add_scan_state = 1.3;
                            end
                           
                        end
                        
                   
                        
                    % ##### Scan end time/duration: how to enter it? #####
                    case  1.6 
                        fprintf(1, ' \n++++ Scan duration & scan end time: ++++ \n');
%                         fprintf(1, ' 1   - Enter scan duration or scan end time manually\n');
%                         fprintf(1, ' 2   - Calculate scan duration & scan end time (NOT AVAILABLE SO FAR!)\n');
                        % Init.: 
                        t_end_jd = 0;
                        t_end_jd_temp = 0;
                        scan_duration_sec_temp = 0;
                        % Update graphics (plot_opt = 5)
                        [sched_handles, error_code, error_msg] = update_sky_plots_pointer(sched_handles, sat_network_id_list, obs_data, 5, PARA);
                        if error_code > 0
                            error_msg = ['update_sky_plots_pointer: ', error_msg];
                            fprintf(1, [' ERROR: ', error_msg, '\n']);
                            add_scan_state = 0;
                        else 
                            [sched_handles, error_code, error_msg] = update_elevation_plot_pointer(sched_handles, obs_data, 5, [], PARA, sat_network_id_list);
                            if error_code > 0
                                error_msg = ['update_elevation_plot_pointer: ', error_msg];
                                fprintf(1, [' ERROR: ', error_msg, '\n']);
                                add_scan_state = 0;
                            else
                                [sched_handles, error_code, error_msg] = update_sky_plots(sched_handles, stat_data, sat_network_id_list, source, PARA, obs_data, flag_plot_sources, t_start_jd);
                                if error_code > 0
                                    error_msg = ['update_sky_plots: ', error_msg];
                                    fprintf(1, [' ERROR: ', error_msg, '\n']);
                                    add_scan_state = 0;
                                else
                                    add_scan_state = 1.7;
                                end
                            end
                        end
                       
                                      
                        
                    % #### Enter scan duration or scan end time manually ####
                    case 1.7 
                        fprintf(1, ' \n++++ Enter the end date/time of the scan or the scan duration ++++\n');
                        fprintf(1,'  => Input Format: <duration [sec]> or <Time Tag Label>  or <yyyy-mm-dd HH:MM:SS>\n');
                        fprintf(1,'  => The scan duration has to be a multiple of the antenna reposition time (%1.2f sec)\n', PARA.REPOS_INT);
                        add_scan_state = 1.71;
                        
                    case 1.71 % Scan end time / duration: Enter it manually
                        commandwindow
                        input_str = input(' t_end / duration = ', 's');
                        len_input_str = length(input_str);
                        
                        % <Time Tag label>
                        if ((input_str(1) == 't') && (len_input_str >= 2) && ...  % Check input Format!
                            (sum((input_str(2:len_input_str) <= '9') & (input_str(2:len_input_str) >= '0')) == (len_input_str - 1)) )
                            
                            % Search <Time Tag label> in obs_data.sat.obs_times get the according time [JD]
                            t_num = str2double(input_str(2:len_input_str));
                            t_end_jd_temp = obs_data.sat(satellite_id).obs_times(logical(obs_data.sat(satellite_id).obs_times(:, 4) == t_num), 2);
                            if isempty(t_end_jd_temp)
                                fprintf(1, ' Error: Time tag label is not available. Please try again!\n');
                            else
                                add_scan_state = 1.72;
                            end
                                    
                        % <yyyy-mm-dd HH:MM:SS>  
                        elseif (len_input_str == 19)
                            if (    (sum((input_str(1:4) <= '9') & (input_str(1:4) >= '0')) == 4)       && ... % Check input Format!
                                    (sum((input_str(6:7) <= '9') & (input_str(6:7) >= '0')) == 2)       && ...
                                    (sum((input_str(9:10) <= '9') & (input_str(9:10) >= '0')) == 2)     && ...
                                    (sum((input_str(12:13) <= '9') & (input_str(12:13) >= '0')) == 2)   && ...
                                    (sum((input_str(15:16) <= '9') & (input_str(15:16) >= '0')) == 2)   && ...
                                    (sum((input_str(18:19) <= '9') & (input_str(18:19) >= '0')) == 2)   && ...
                                    (input_str(5) == '-')                                               && ...
                                    (input_str(8) == '-')                                               && ...
                                    (input_str(11) == ' ')                                              && ...
                                    (input_str(14) == ':')                                              && ...
                                    (input_str(17) == ':')                                                      )   
                                
                                % Convert input to JD:
                                yr = str2num(input_str(1:4));
                                mon = str2num(input_str(6:7));
                                day = str2num(input_str(9:10));
                                hr = str2num(input_str(12:13));
                                min = str2num(input_str(15:16));
                                sec = str2num(input_str(18:19));
                                t_end_jd_temp = jday(yr, mon, day, hr, min, sec);
                                add_scan_state = 1.72;
                            else
                                fprintf(1, ' ERROR: Ivalid input. Please try again!\n');
                            end
                            
                        % <duration [sec]>  
                        elseif (    (len_input_str < 19) && ...                                                             % Check length...
                                    (sum((input_str(1:end) <= '9') & (input_str(1:end) >= '0')) == length(input_str)) ...   % Only numbers...
                                )
                            scan_duration_sec_temp = str2double(input_str);
                            % Calc scan end time:
                            t_end_jd_temp = t_start_jd + scan_duration_sec_temp / 86400;
                            [year, mon, day, hr, min, sec] = invjday (t_end_jd_temp); % Conversion
                            fprintf(1, ' Calculated t_end with a scan duration of %1.1f seconds: %4.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\n', scan_duration_sec_temp, year, mon, day, hr, min, sec);    % <yyyy-mm-dd HH:MM:SS> 
                            add_scan_state = 1.72;
                        % other commands    
                        else
                            switch(input_str)
                                case 'exit'
                                    add_scan_state = 9;
                                case 'back'
                                    add_scan_state = 1.6;
                                case 'help'
                                    disp(' ');
                                    disp('  => Choose t_end: Three input formats are possible:');
                                    disp('      (a) <Time Tag Label> , e.g. "t2"');
                                    disp('          The input time tag label has to correspond with an available label of the choosen satellite.');
                                    disp('      (b) <yyyy-mm-dd HH:MM:SS> , e.g. "2000-03-01 15:45:17"');
                                    disp('          The input date and time has to be within an available observation periode.');
                                    disp('      (c) <duration [sec]> , e.g. "15"');
                                    disp('          The resulting scan end time has to be within an available observation periode.');
                                    disp(' ');

                                otherwise
                                    fprintf(1, ' ERROR: Invalid input. Type "help" to see input options!\n');
                            end
                        end
                        
                    % Check the scan duration:
                    case 1.72 
                        [flag_scan_duration_ok, alt_t_end_jd, repos_epochs_jd, error_code, error_msg] = check_scan_duration_and_calc_repos_epochs(PARA, t_start_jd, t_end_jd_temp);
                        if error_code > 0
                            error_msg = ['check_scan_duration_and_calc_repos_epochs: ', error_msg];
                            fprintf(1, [' ERROR: ', error_msg, '\n']);
                            add_scan_state = 0;
                        else
                            if ~flag_scan_duration_ok % the scan duration is not a multiple of the antenna repos. interval...
                                add_scan_state = 1.721;
                            else % Scan duration is OK...
                                % Update graphics:
                                flag_plot_sources = 0;
                                [sched_handles, error_code, error_msg] = update_sky_plots(sched_handles, stat_data, sat_network_id_list, source, PARA, obs_data, flag_plot_sources, t_end_jd_temp);
                                if error_code > 0
                                    error_msg = ['update_sky_plots: ', error_msg];
                                    fprintf(1, [' ERROR: ', error_msg, '\n']);
                                    add_scan_state = 0;
                                else
                                    [sched_handles, error_code, error_msg] = update_elevation_plot_pointer(sched_handles, obs_data, 9, t_end_jd_temp, PARA, sat_network_id_list);
                                    if error_code > 0
                                        error_msg = ['update_elevation_plot_pointer: ', error_msg];
                                        fprintf(1, [' ERROR: ', error_msg, '\n']);
                                        add_scan_state = 0;
                                    else
                                        add_scan_state = 1.73;
                                    end
                                end 
                            end
                        end
                        
                    case 1.721 % Check the scan duration: NOT OK!
                        scan_duration_sec_temp = (t_end_jd_temp - t_start_jd) * 86400;
                        fprintf(1, ' ERROR: Scan duration (%1.1f sec) is not a multiple of the antenna repos. interval (%1.1f sec)!\n', scan_duration_sec_temp, PARA.REPOS_INT);
                        [year, mon, day, hr, min, sec] = invjday (alt_t_end_jd); % Conversion
                        scan_duration_sec_temp = (alt_t_end_jd - t_start_jd) * 86400;
                        % Update graphics:
                        flag_plot_sources = 0;
                        [sched_handles, error_code, error_msg] = update_sky_plots(sched_handles, stat_data, sat_network_id_list, source, PARA, obs_data, flag_plot_sources, alt_t_end_jd);
                        if error_code > 0
                            error_msg = ['update_sky_plots: ', error_msg];
                            fprintf(1, [' ERROR: ', error_msg, '\n']);
                            add_scan_state = 0;
                        else
                            [sched_handles, error_code, error_msg] = update_elevation_plot_pointer(sched_handles, obs_data, 9, alt_t_end_jd, PARA, sat_network_id_list);
                            if error_code > 0
                                error_msg = ['update_elevation_plot_pointer: ', error_msg];
                                fprintf(1, [' ERROR: ', error_msg, '\n']);
                                add_scan_state = 0;
                            else
                                fprintf(1, '  => Use a duration of %1.1f sec instead? (Resulting scan end: %4.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f)\n', scan_duration_sec_temp, year, mon, day, hr, min, sec);
                                add_scan_state = 1.722;
                            end
                        end

                    case  1.722 % Check the scan duration: Use alternaticve duration instead?
                        commandwindow
                        input_str = input(' Input(y=yes, n=no): ', 's');
                        switch(input_str)
                            case {'n' 'no'}
                                add_scan_state = 1.6; % Enter scan duration manually
                            case {'y' 'yes'}
                                t_end_jd_temp = alt_t_end_jd;
                                add_scan_state = 1.73; % Check scan end time
                            case 'exit'
                                add_scan_state = 9;
                            case 'back'
                                add_scan_state = 1.6; % Enter scan duration manually
                            otherwise
                                fprintf(1, ' ERROR: Invalid input, please try again!\n');
                        end
                        
                    % Check scan end time
                    case 1.73
                        source_quasar = [];
                        [flag_t_end_ok, obs_data, error_code, error_msg] = check_t_end(stat_data, sat_network_id_list, obs_data, PARA, source_quasar, satellite_id, t_end_jd_temp, t_start_jd, flag_sat_network_id_list_changed);
                        if error_code > 0
                            error_msg = ['check_t_end: ', error_msg];
                            fprintf(1, [' ERROR: ', error_msg, '\n']);
                            add_scan_state = 0;
                        else
                            if flag_t_end_ok % NEXT STATE! t_end found!
                                % Update graphics:
                                flag_plot_sources = 0;
                                [sched_handles, error_code, error_msg] = update_sky_plots(sched_handles, stat_data, sat_network_id_list, source, PARA, obs_data, flag_plot_sources, t_end_jd_temp);
                                if error_code > 0
                                    error_msg = ['update_sky_plots: ', error_msg];
                                    fprintf(1, [' ERROR: ', error_msg, '\n']);
                                    add_scan_state = 0;
                                else
                                    [sched_handles, error_code, error_msg] = update_elevation_plot_pointer(sched_handles, obs_data, 9, t_end_jd_temp, PARA, sat_network_id_list);
                                    if error_code > 0
                                        error_msg = ['update_elevation_plot_pointer: ', error_msg];
                                        fprintf(1, [' ERROR: ', error_msg, '\n']);
                                        add_scan_state = 0;
                                    else
                                        add_scan_state = 1.731; % Confirm input
                                    end
                                end
                            else
                                fprintf(1, ' ERROR: Scan end time is not valid. Try again!\n');
                                add_scan_state = 1.6;
                            end
                        end

                        
                    case 1.731 % Finally, confirm the entered/calculated scan end time (<duration> input, or calculated)
                        commandwindow
                        input_str = input(' Take it? (y=yes, n=no, 9=back) <y> : ', 's');
                        remove_highlight = true;
                        switch(input_str)
                            case {'n' 'no'}
                                add_scan_state = 1.6; % Scan end time: how to enter it?
                            case {'9'}
                                add_scan_state = 9;
                            case {'y' 'yes' ''}
                                remove_highlight = false;
                                [sched_handles, error_code, error_msg] = update_sky_plots_pointer(sched_handles, sat_network_id_list, obs_data, 3, PARA);
                                if error_code > 0
                                    error_msg = ['update_sky_plots_pointer: ', error_msg];
                                    fprintf(1, [' ERROR: ', error_msg, '\n']);
                                    add_scan_state = 0;
                                else
                                    [sched_handles, error_code, error_msg] = update_elevation_plot_pointer(sched_handles, obs_data, 3, [], PARA, sat_network_id_list);
                                    if error_code > 0
                                        error_msg = ['update_elevation_plot_pointer: ', error_msg];
                                        fprintf(1, [' ERROR: ', error_msg, '\n']);
                                        add_scan_state = 0;
                                    else
                                        t_end_jd = t_end_jd_temp;
                                        add_scan_state = 1.9; % #### Confirm scan ####
                                    end
                                end
                            case 'exit'
                                add_scan_state = 9;
                            case 'back'
                                add_scan_state = 1.6;
                            otherwise
                                remove_highlight = false;
                                fprintf(1, ' ERROR: Invalid input, please try again!\n');
                            
                        end
                        
                        if remove_highlight
                            highlight_satellite( sched_handles, station_id_list, satellite_id, false );
                        end
                        
                    % #### Confirm scan and save it ####
                    case  1.9 
                        [year_start, mon_start, day_start, hr_start, min_start, sec_start] = invjday(t_start_jd); % Conversion
                        [year_end, mon_end, day_end, hr_end, min_end, sec_end] = invjday(t_end_jd); % Conversion
                        scan_duration_sec_temp = (t_end_jd - t_start_jd) * 86400;
                        fprintf(1, ' \n++++ Confirm scan and save it: ++++ \n');
                        fprintf(1, '   => Satellite:      %s\n', sat_name_str);
                        fprintf(1, '   => Scan start:     %4.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\n',  year_start, mon_start, day_start, hr_start, min_start, sec_start);
                        fprintf(1, '   => Scan end:       %4.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\n',  year_end, mon_end, day_end, hr_end, min_end, sec_end);
                        fprintf(1, '   => Scan duration:  %1.1f seconds\n',  scan_duration_sec_temp);
                        fprintf(1, ' 1   - Save scan and append another one\n');
                        fprintf(1, ' 2   - Save scan and go back to main menu\n');
                        fprintf(1, ' 3   - Discard scan and go back to main menu\n');
                        add_scan_state = 1.91;
                        
                    case 1.91 % Confirm scan and save it
                        commandwindow
                        input_str = input(' Please select <1>: ', 's');
                        quasar_id = [];
                        remove_highlight = true;
                        switch(input_str)
                            
                            case {'1', ''}
                                % Save scan:
                                [sched_data, obs_data, error_code, error_msg] = save_scan(stat_data, sat_network_id_list, sched_data, obs_data, t_start_jd, t_end_jd, source, quasar_id, satellite_id, PARA);
                                if error_code > 0
                                    error_msg = ['save_scan: ', error_msg];
                                    fprintf(1, [' ERROR: ', error_msg, '\n']);
                                    add_scan_state = 0;
                                else
                                    fprintf(1, ' ...scan number %1.0f saved! \n\n', obs_data.number_of_scans);
                                    add_scan_state = 0; % Add scan to current Schedule (append)
                                end
                                 
                            case '2'
                                % Save scan here:
                                [sched_data, obs_data, error_code, error_msg] = save_scan(stat_data, sat_network_id_list, sched_data, obs_data, t_start_jd, t_end_jd, source, quasar_id, satellite_id, PARA);
                                if error_code > 0
                                     error_msg = ['save_scan: ', error_msg];
                                     fprintf(1, [' ERROR: ', error_msg, '\n']);
                                     add_scan_state = 0;
                                else
                                    fprintf(1, ' ...scan number %1.0f saved! \n\n', obs_data.number_of_scans);
                                    add_scan_state = 9; % Exit the "add scan" routine and go back to main menu
                                end
                                
                            case '3'
                                add_scan_state = 9; % Exit the "add scan" routine and go back to main menu
                            case 'exit'
                                add_scan_state = 9; % Exit the "add scan" routine and go back to main menu
                            case 'back'
                                add_scan_state = 1.6; % Enter the end date/time of the scan or the scan duration
                            otherwise
                                fprintf(1, ' ERROR: Invalid input, please try again!\n');
                                remove_highlight = false;
                        end
                        
                        if remove_highlight
                            highlight_satellite( sched_handles, station_id_list, satellite_id, false );
                        end

                        
                    % ##########################
                    % ##### Quasar scan    #####
                    % ##########################

                    % ##### How to choose a quasar? #####
                    case 2.0
                        fprintf(1, ' \n++++ How to choose a quasar: ++++ \n');
                        fprintf(1, ' 1   - Finde the closest one at the end of the last scan\n');
                        fprintf(1, ' 2   - Finde the n closest quasars at the end of the last scan\n');
                        fprintf(1, ' 3   - Enter quasar name\n');
                        fprintf(1, ' 4   - Select quasar in a sky-plot\n'); 
                        fprintf(1, ' 5   - Automatic station based scheduling for a defined time period\n');
                        
                        % Init.:
                        quasar_id = []; % id of selected quasar (refers to field in "source" structure)
                        t_start_jd = 0;
                        t_start_jd_temp = 0;
                        t_end_jd = 0;
                        t_end_jd_temp = 0;
                        
                        % Get the epoch for further calculations:
                        if obs_data.number_of_scans == 0 % first scan...
                            t_epoch_jd = PARA.startmjd + 2.400000500000000e+006;
                            [year, mon, day, hr, min, sec] = invjday (t_epoch_jd); % Conversion
                            fprintf(1, ' Session start: %4.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\n', year, mon, day, hr, min, sec);    % <yyyy-mm-dd HH:MM:SS> 
                        else
                            t_epoch_jd = obs_data.end_of_last_scan_jd;
                            [year, mon, day, hr, min, sec] = invjday (t_epoch_jd); % Conversion
                            fprintf(1, ' \nEnd of last scan: %4.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\n', year, mon, day, hr, min, sec);    % <yyyy-mm-dd HH:MM:SS> 
                        end
                        
                        % Pre-seelction of observable quasars at the end of the last scan:
                        [source_observable, flag_observable_quasars_list, error_code, error_msg] = preselect_observable_quasars(stat_data, quasar_network_id_list, source, t_epoch_jd);
                        if error_code > 0
                             error_msg = ['preselect_observable_quasars: ', error_msg];
                             fprintf(1, [' ERROR: ', error_msg, '\n']);
                             add_scan_state = 0;
                        else
                            
                            % Update graphics:
                            [sched_handles, error_code, error_msg] = update_elevation_plot_pointer(sched_handles, obs_data, 4, [], PARA, quasar_network_id_list);
                            if error_code > 0
                                 error_msg = ['update_elevation_plot_pointer: ', error_msg];
                                 fprintf(1, [' ERROR: ', error_msg, '\n']);
                                 add_scan_state = 0;
                            else
                                flag_plot_sources = 1;
                                [sched_handles, error_code, error_msg] = update_sky_plots(sched_handles, stat_data, quasar_network_id_list, source, PARA, obs_data, flag_plot_sources, t_epoch_jd, flag_observable_quasars_list);
                                if error_code > 0
                                     error_msg = ['update_sky_plots: ', error_msg];
                                     fprintf(1, [' ERROR: ', error_msg, '\n']);
                                     add_scan_state = 0;
                                else
                                    [sched_handles, error_code, error_msg] = update_sky_plots_pointer(sched_handles, quasar_network_id_list, obs_data, 4, PARA);
                                    if error_code > 0
                                         error_msg = ['update_sky_plots_pointer: ', error_msg];
                                         fprintf(1, [' ERROR: ', error_msg, '\n']);
                                         add_scan_state = 0;
                                    else
                                    [sched_handles, error_code, error_msg] = update_sky_plots_highlight_quasars(sched_handles, stat_data, quasar_network_id_list, source, t_epoch_jd, quasar_id, 2); % Delete highlighted sources (markers)
                                        if error_code > 0
                                             error_msg = ['update_sky_plots_pointer: ', error_msg];
                                             fprintf(1, [' ERROR: ', error_msg, '\n']);
                                             add_scan_state = 0;
                                        else
                                            add_scan_state = 2.01;
                                        end
                                    end
                                end
                            end
                        end
                        
                    case 2.01
                        commandwindow
                        input_str = input(' Please select: ', 's');
                        switch(input_str)
                            case '1'
                                number_quasars_output = 1;
                                add_scan_state = 2.1; % Finde the closest one
                            case '2'
                                fprintf(1, ' \n++++ Enter the number of quasars (n):  ++++ \n');
                                add_scan_state = 2.2; % Finde the closest n sources
                            case '3'
                                add_scan_state = 2.3; % Enter quasar name
                            case '4'
                                add_scan_state = 2.4; % Select quasar in a sky-plot
                            case '5'
                                add_scan_state = 3; % Automatic station based scheduling for a defined time period
                            case 'exit'
                                add_scan_state = 9;
                            case 'back'
                                add_scan_state = 0;
                            otherwise
                                fprintf(1, ' ERROR: Invalid input, please try again!\n');
                        end
                    

                    % ##### 1. Find the closest quasar #####
                    case 2.1
                        number_quasars_output = 1;
                        flag_print_quasar_infos = 1;
                        [quasars_info_list, obs_data, error_code, error_msg] = find_closest_quasars_and_calc_scan(stat_data, quasar_network_id_list, sat_network_id_list, source, obs_data, flag_observable_quasars_list, number_quasars_output, t_epoch_jd, flag_print_quasar_infos, PARA);
                        if error_code > 0
                             error_msg = ['find_closest_quasars_and_calc_scan: ', error_msg];
                             fprintf(1, [' ERROR: ', error_msg, '\n']);
                             add_scan_state = 0;
                        else
                            if size(quasars_info_list, 1) == 1;
                                % Plot quasar markers (highlight):
                                [sched_handles, error_code, error_msg] = update_sky_plots_highlight_quasars(sched_handles, stat_data, quasar_network_id_list, source, t_epoch_jd, quasars_info_list(:,1), 1);
                                if error_code > 0
                                     error_msg = ['update_sky_plots_highlight_quasars: ', error_msg];
                                     fprintf(1, [' ERROR: ', error_msg, '\n']);
                                     add_scan_state = 0;
                                else
                                    add_scan_state = 2.23;
                                end
                            else
                                fprintf(1, ' ERROR: No close quasar available!\n');
                                add_scan_state = 2.0;
                            end
                        end
                        

                    % ##### 2. Find the closest n quasars #####
                    case 2.20 % Choose number of sources "n"
                        commandwindow
                        input_str = input(' Enter n: ', 's');
                        if sum((input_str(1:end) <= '9') & (input_str(1:end) >= '0')) == length(input_str) % Only numbers...
                            number_quasars_output = str2double(input_str);
                            add_scan_state = 2.21;
                        else
                            switch(input_str)
                                case 'exit'
                                    add_scan_state = 9;
                                case 'back'
                                    add_scan_state = 2;
                                otherwise
                                    fprintf(1, ' ERROR: Invalid input, please try again!\n');
                            end
                        end
                    
                    
                    case 2.21
                        % Output options:
                        flag_print_quasar_infos = 1;
                        [quasars_info_list, obs_data, error_code, error_msg] = find_closest_quasars_and_calc_scan(stat_data, quasar_network_id_list, sat_network_id_list, source, obs_data, flag_observable_quasars_list, number_quasars_output, t_epoch_jd, flag_print_quasar_infos, PARA);
                        if error_code > 0
                             error_msg = ['find_closest_quasars_and_calc_scan: ', error_msg];
                             fprintf(1, [' ERROR: ', error_msg, '\n']);
                             add_scan_state = 0;
                        else
                            % Plot quasar markers (highlight):
                            [sched_handles, error_code, error_msg] = update_sky_plots_highlight_quasars(sched_handles, stat_data, quasar_network_id_list, source, t_epoch_jd, quasars_info_list(:,1), 1);
                            if error_code > 0
                                 error_msg = ['update_sky_plots_highlight_quasars: ', error_msg];
                                 fprintf(1, [' ERROR: ', error_msg, '\n']);
                                 add_scan_state = 0;
                            else
                                number_of_quasars_to_select = size(quasars_info_list, 1);
                                if number_of_quasars_to_select > 1
                                    add_scan_state = 2.22;
                                    fprintf(1, ' \n++++ Sellect one quasar:  ++++ \n');
                                elseif number_of_quasars_to_select == 1
                                    add_scan_state = 2.23;
                                    fprintf(1, ' \n++++ Observe this quasar?  ++++ \n');
                                else
                                    fprintf(1, ' ERROR: List of observable quasars is empty!\n');
                                    add_scan_state = 0;
                                end
                            end
                        end
                        
                        
                    case 2.22 % Select one quasar:
                        commandwindow
                        input_str = input(' Select one quasar by its number: ', 's');
                        if sum((input_str(1:end) <= '9') & (input_str(1:end) >= '0')) == length(input_str) % Only numbers...
                            quasar_selection_number = str2double(input_str);
                            if quasar_selection_number <= number_of_quasars_to_select
                                add_scan_state = 2.6;
                            else
                                fprintf(1, ' ERROR: Invalid input, please try again!\n');
                            end
                        else
                            switch(input_str)
                                case 'exit'
                                    add_scan_state = 9;
                                case 'back'
                                    add_scan_state = 2;
                                otherwise
                                    fprintf(1, ' ERROR: Invalid input, please try again!\n');
                            end
                        end
                        
                        
                    case 2.23 % Take the only available quasar?
                        commandwindow
                        input_str = input(' Observe this quasar? (y=yes, n=no) : ', 's');
                        switch(input_str)
                            case {'n' 'no'}
                                add_scan_state = 2;
                            case {'y' 'yes'}
                                quasar_selection_number = 1;
                                add_scan_state = 2.6;
                            case 'exit'
                                add_scan_state = 9;
                            case 'back'
                                add_scan_state = 2;
                            otherwise
                                fprintf(1, ' ERROR: Invalid input, please try again!\n');
                        end
                        

                        
                    % ##### 3. Enter quasar name #####
                    case 2.3 
                        commandwindow
                        input_str = input(' Enter quasar name: ', 's');
                        switch(input_str)
                            case 'exit'
                                add_scan_state = 9;
                            case 'back'
                                add_scan_state = 2;
                            otherwise
                                % Search for the enterd name in the list of observable quasars:
                                fag_found_quasar = 0;
                                for i_quasar = 1 : length(source_observable)
                                    if strcmp(source_observable(i_quasar).name, input_str)
                                        fag_found_quasar = 1;
                                        quasar_id_list = find(flag_observable_quasars_list);
                                        quasar_id = quasar_id_list(i_quasar);
                                        break;
                                    end
                                end
                                if fag_found_quasar
                                    add_scan_state = 2.5;
                                else
                                    fprintf(1, ' ERROR: The entered quasar name was not found in the list of observable quasars. Please try again!\n');
                                end
                        end
                        
  
                        
                    % ##### 4. Select quasar in a sky plot #####
                    case 2.4
                        % Print station list:
                        for i_stat = 1 : length(quasar_network_id_list)
                            station_id = quasar_network_id_list(i_stat);
                            fprintf(1, '   %1.0f - %s \n', i_stat, stat_data.stat(station_id).name);
                        end
                        % user input:
                        commandwindow
                        input_str = input(' Select skyplot by station number: ', 's');
                        if sum((input_str(1:end) <= '9') & (input_str(1:end) >= '0')) == length(input_str) % Only numbers...
                            i_stat = str2double(input_str);
                            station_id = quasar_network_id_list(i_stat);
                            if i_stat <= length(quasar_network_id_list)
                                % ##### Select quasar via mouse input in the selected skyplot: #####
                                [error_code, error_msg, quasar_id] = select_quasar_in_skyplot(sched_handles, station_id);
                                if error_code > 0
                                     error_msg = ['select_quasar_in_skyplot: ', error_msg];
                                     fprintf(1, [' ERROR: ', error_msg, '\n']);
                                     add_scan_state = 0;
                                else
                                    % ##### Highlight selected quasar in skyplot #####
                                    [sched_handles, error_code, error_msg] = update_sky_plots_highlight_quasars(sched_handles, stat_data, quasar_network_id_list, source, t_epoch_jd, quasar_id, 1);
                                    if error_code > 0
                                         error_msg = ['update_sky_plots_highlight_quasars: ', error_msg];
                                         fprintf(1, [' ERROR: ', error_msg, '\n']);
                                         add_scan_state = 0;
                                    else
                                        add_scan_state = 2.5;
                                    end
                                end
                            else
                                fprintf(1, ' ERROR: Invalid input, please try again!\n');
                            end
                        else
                            switch(input_str)
                                case 'exit'
                                    add_scan_state = 9;
                                case 'back'
                                    add_scan_state = 2;
                                otherwise
                                    fprintf(1, ' ERROR: Invalid input, please try again!\n');
                            end
                        end
                        
  
                    % ##### Calculate quasar scan #####
                    case 2.5
                        flag_print_quasar_infos = 1;
                        [quasars_info_list, obs_data, error_code, error_msg] = calc_quasar_scan(stat_data, quasar_network_id_list, source, obs_data, quasar_id, flag_print_quasar_infos, PARA);
                        if error_code > 0
                             error_msg = ['calc_quasar_scan: ', error_msg];
                             fprintf(1, [' ERROR: ', error_msg, '\n']);
                             add_scan_state = 0;
                        else
                            if isempty(quasars_info_list)
                                add_scan_state = 2.0;
                            else
                                quasar_selection_number = 1;
                                add_scan_state = 2.6; % Quasar selected successfully => Assign parameters
                            end
                        end
                        

                    % ##### Quasar selected successfully => Assign parameters: #####
                    % If one quasar is selected ("quasar_selection_number") => go to this case!!
                    case 2.6
                        quasar_id           = quasars_info_list(quasar_selection_number, 1);
                        t_start_jd_temp     = quasars_info_list(quasar_selection_number, 2);
                        t_end_jd_temp       = quasars_info_list(quasar_selection_number, 3);
                        number_of_stations = length(quasar_network_id_list);
                        for i_stat = 1 : number_of_stations
                            station_id = quasar_network_id_list(i_stat);
                            obs_data.stat(station_id).begin_of_new_obs_temp.jd      = quasars_info_list(quasar_selection_number, 3 + i_stat); % [JD]
                            obs_data.stat(station_id).begin_of_new_obs_temp.un_az   = quasars_info_list(quasar_selection_number, 3 + i_stat + number_of_stations); % [rad]
                            obs_data.stat(station_id).begin_of_new_obs_temp.el      = quasars_info_list(quasar_selection_number, 3 + i_stat + 2 * number_of_stations); % [rad]
                            obs_data.stat(station_id).end_of_new_obs_temp.jd        = quasars_info_list(quasar_selection_number, 3 + i_stat + 3 * number_of_stations);  % [JD]
                            obs_data.stat(station_id).end_of_new_obs_temp.un_az     = quasars_info_list(quasar_selection_number, 3 + i_stat + 4 * number_of_stations);  % [rad]
                            obs_data.stat(station_id).end_of_new_obs_temp.el        = quasars_info_list(quasar_selection_number, 3 + i_stat + 5 * number_of_stations);  % [rad]
                        end
                        fprintf(1, '  Selected quasar: %s\n', source(quasar_id).name);
                        
                        % Update graphics:
                        [sched_handles, error_code, error_msg] = update_elevation_plot_pointer(sched_handles, obs_data, 6, [], PARA, quasar_network_id_list);
                        if error_code > 0
                             error_msg = ['update_elevation_plot_pointer: ', error_msg];
                             fprintf(1, [' ERROR: ', error_msg, '\n']);
                             add_scan_state = 0;
                        else
                            flag_plot_sources = 1;
                            [sched_handles, error_code, error_msg] = update_sky_plots(sched_handles, stat_data, quasar_network_id_list, source, PARA, obs_data, flag_plot_sources, obs_data.stat(station_id).end_of_new_obs_temp.jd, flag_observable_quasars_list);
                            if error_code > 0
                                 error_msg = ['update_sky_plots: ', error_msg];
                                 fprintf(1, [' ERROR: ', error_msg, '\n']);
                                 add_scan_state = 0;
                            else
                                [sched_handles, error_code, error_msg] = update_sky_plots_pointer(sched_handles, quasar_network_id_list, obs_data, 6, PARA);
                                if error_code > 0
                                     error_msg = ['update_sky_plots_pointer: ', error_msg];
                                     fprintf(1, [' ERROR: ', error_msg, '\n']); 
                                     add_scan_state = 0;
                                else
                                    % Plot quasar markers (highlight):
                                    [sched_handles, error_code, error_msg] = update_sky_plots_highlight_quasars(sched_handles, stat_data, quasar_network_id_list, source, obs_data.stat(station_id).end_of_new_obs_temp.jd, quasars_info_list(:,1), 1);
                                    if error_code > 0
                                         error_msg = ['update_sky_plots_highlight_quasars: ', error_msg];
                                         fprintf(1, [' ERROR: ', error_msg, '\n']);
                                         add_scan_state = 0;
                                    else
                                        add_scan_state = 2.61;
                                    end
                                end
                            end
                        end

                    
                    % ##### Select options #####
                    case 2.61
                        fprintf(1, ' \n++++ Please select: ++++ \n');
                        fprintf(1, ' 1   - Take the quasar, scan start time and duration\n');
                        fprintf(1, ' 2   - Change the scan duration\n');
                        fprintf(1, ' 3   - Change scan start time and duration\n');
                        fprintf(1, ' 4   - Select another quasar\n');
                        add_scan_state = 2.611;

                    case 2.611 % Select options?
                        commandwindow
                        input_str = input(' Please select: ', 's');
                        switch(input_str)
                            case '1'
                                t_start_jd = t_start_jd_temp;
                                t_end_jd = t_end_jd_temp;
                                add_scan_state = 2.7; % Take the quasar, scan start time and duration => #### Confirm quasar scan and save it ####
                            case '2'
                                add_scan_state = 2.62; % Take the quasar and scan start time, but enter duration
                            case '3'
                                add_scan_state = 2.63; % Take the quasar and enter scan start- time and duration/scan end-time manually 
                            case '4'
                                add_scan_state = 2.0; %
                            case 'exit'
                                add_scan_state = 9;
                            case 'back'
                                add_scan_state = 2.0;
                            otherwise
                                fprintf(1, ' ERROR: Invalid input, please try again!\n');
                        end
                        

                    % ####  Change scan duration ####
                    case 2.62 
                        commandwindow
                        input_str = input(' Scan duration [sec] = ', 's');
                        % Check input:
                        if sum((input_str(1:end) <= '9') & (input_str(1:end) >= '0')) == length(input_str) % Only numbers...
                            scan_duration_sec_temp = str2double(input_str);
                            [quasars_info_list, obs_data, error_code, error_msg] = calc_quasar_scan(stat_data, quasar_network_id_list, source, obs_data, quasar_id, true, PARA, scan_duration_sec_temp);
                            if error_code > 0
                                 error_msg = ['calc_quasar_scan: ', error_msg];
                                 fprintf(1, [' ERROR: ', error_msg, '\n']);
                                 add_scan_state = 0;
                            else
                                if isempty(quasars_info_list) % Invalid scan
                                    add_scan_state = 2.61;
                                else
                                    quasar_selection_number = 1;
                                    add_scan_state = 2.6; % Quasar selected successfully => Assign parameters
                                end
                            end
                        else
                            switch(input_str)
                                case 'exit'
                                    add_scan_state = 9;
                                case 'back'
                                    add_scan_state = 2;
                                otherwise
                                    fprintf(1, ' ERROR: Invalid input, please try again!\n');
                            end
                        end
                            

                    % ####  Change scan duration and start time ####
                    case 2.63 
                        commandwindow
                        input_str = input(' Scan duration [sec] = ', 's');
                        % Check input:
                        if sum((input_str(1:end) <= '9') & (input_str(1:end) >= '0')) == length(input_str) % Only numbers...
                            scan_duration_sec_temp = str2double(input_str);
                            add_scan_state = 2.64; % select start time
                        else
                            switch(input_str)
                                case 'exit'
                                    add_scan_state = 9;
                                case 'back'
                                    add_scan_state = 2;
                                otherwise
                                    fprintf(1, ' ERROR: Invalid input, please try again!\n');
                            end
                        end
                        
                    case 2.64 % select start time
                        commandwindow
                        fprintf(1, '   Enter scan start time\n');
                        fprintf(1,'    => Input Format <yyyy-mm-dd HH:MM:SS>\n');
                        input_str = input(' t_start = ', 's');
                        if (length(input_str) == 19)
                            [flag_t_epoch_jd_ok, t_epoch_jd_ok] = check_t_input_19char(input_str, sched_data, PARA);
                            if flag_t_epoch_jd_ok
                                t_start_jd_temp = t_epoch_jd_ok;
                                add_scan_state = 2.65; % Check scan
                            end
                        else
                            switch(input_str)
                                case 'exit'
                                    add_scan_state = 1;
                                case 'back'
                                    add_scan_state = 4;
                                otherwise
                                    fprintf(1, ' ERROR: Invalid input. Please try again!\n');
                            end
                        end
                        
                    case 2.65 % Check scan
                        [quasars_info_list, obs_data, error_code, error_msg] = calc_quasar_scan(stat_data, quasar_network_id_list, source, obs_data, quasar_id, true, PARA, scan_duration_sec_temp, t_start_jd_temp);
                        if error_code > 0
                             error_msg = ['calc_quasar_scan: ', error_msg];
                             fprintf(1, [' ERROR: ', error_msg, '\n']);
                             add_scan_state = 0;
                        else
                            if isempty(quasars_info_list) % Invalid scan
                                add_scan_state = 2.61;
                            else
                                quasar_selection_number = 1;
                                add_scan_state = 2.6; % Quasar selected successfully => Assign parameters
                            end
                        end
                        

                    % #### Confirm quasar scan and save it ####
                    case  2.7 
                        [year_start, mon_start, day_start, hr_start, min_start, sec_start] = invjday(t_start_jd); % Conversion
                        [year_end, mon_end, day_end, hr_end, min_end, sec_end] = invjday(t_end_jd); % Conversion
                        scan_duration_sec_temp = (t_end_jd - t_start_jd) * 86400;
                        fprintf(1, ' \n++++ Confirm scan and save it: ++++ \n');
                        fprintf(1, '   => Quasar:         %s\n', source(quasar_id).name);
                        fprintf(1, '   => Scan start:     %4.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\n',  year_start, mon_start, day_start, hr_start, min_start, sec_start);
                        fprintf(1, '   => Scan end:       %4.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\n',  year_end, mon_end, day_end, hr_end, min_end, sec_end);
                        fprintf(1, '   => Scan duration:  %1.1f seconds\n',  scan_duration_sec_temp);
                        fprintf(1, ' 1   - Save scan and append another one\n');
                        fprintf(1, ' 2   - Save scan and go back to main menu\n');
                        fprintf(1, ' 3   - Discard scan and go back to main menu\n');
                        add_scan_state = 2.71;
                        

                    case 2.71 % Confirm quasar scan and save it
                        commandwindow
                        input_str = input(' Please select: ', 's');
                        satellite_id = [];
                        switch(input_str)
                            
                            case '1'
                                % Save scan:
                                [sched_data, obs_data, error_code, error_msg] = save_scan(stat_data, quasar_network_id_list, sched_data, obs_data, t_start_jd, t_end_jd, source, quasar_id, satellite_id, PARA);
                                if error_code > 0
                                    error_msg = ['save_scan: ', error_msg];
                                    fprintf(1, [' ERROR: ', error_msg, '\n']);
                                    add_scan_state = 0;
                                else
                                    fprintf(1, ' ...scan number %1.0f saved! \n\n', obs_data.number_of_scans);
                                    add_scan_state = 0; % Add scan to current Schedule (append)
                                end
                                 
                            case '2'
                                % Save scan here:
                                [sched_data, obs_data, error_code, error_msg] = save_scan(stat_data, quasar_network_id_list, sched_data, obs_data, t_start_jd, t_end_jd, source, quasar_id, satellite_id, PARA);
                                if error_code > 0
                                     error_msg = ['save_scan: ', error_msg];
                                     fprintf(1, [' ERROR: ', error_msg, '\n']);
                                     add_scan_state = 0;
                                else
                                    fprintf(1, ' ...scan number %1.0f saved! \n\n', obs_data.number_of_scans);
                                    add_scan_state = 9; % Exit the "add scan" routine and go back to main menu
                                end
                                
                            case '3'
                                add_scan_state = 9; % Exit the "add scan" routine and go back to main menu
                            case 'exit'
                                add_scan_state = 9; % Exit the "add scan" routine and go back to main menu
                            case 'back'
                                add_scan_state = 2.61; % ##### Select options #####
                            otherwise
                                fprintf(1, ' ERROR: Invalid input, please try again!\n');
                        end
                       
                    % ######################################################################
                    % #### Automatic station based scheduling for a defined time period ####
                    % ######################################################################
                    case 3
                        fprintf(1, ' \n++++ Station based scheduling for quasars ++++ \n');
                        
                        % Init.:
                        t_start_jd_temp = 0;
                        t_start_jd = 0;
                        fprintf(1, '   Enter start time\n');
                        fprintf(1,'    => Input Format <yyyy-mm-dd HH:MM:SS> or "y" to take the end of the last scan\n');
                        add_scan_state = 3.1;
                        
                    case 3.1 % Scan start time: Enter it
                        commandwindow
                        input_str = input(' t_start = ', 's');
                        len_input_str = length(input_str);

                        if (len_input_str == 19)
                            if (    (sum((input_str(1:4) <= '9') & (input_str(1:4) >= '0')) == 4)       && ... % Check input Format!
                                    (sum((input_str(6:7) <= '9') & (input_str(6:7) >= '0')) == 2)       && ...
                                    (sum((input_str(9:10) <= '9') & (input_str(9:10) >= '0')) == 2)     && ...
                                    (sum((input_str(12:13) <= '9') & (input_str(12:13) >= '0')) == 2)   && ...
                                    (sum((input_str(15:16) <= '9') & (input_str(15:16) >= '0')) == 2)   && ...
                                    (sum((input_str(18:19) <= '9') & (input_str(18:19) >= '0')) == 2)   && ...
                                    (input_str(5) == '-')                                               && ...
                                    (input_str(8) == '-')                                               && ...
                                    (input_str(11) == ' ')                                              && ...
                                    (input_str(14) == ':')                                              && ...
                                    (input_str(17) == ':')                                                      )   
                                % Convert input to JD:
                                yr = str2num(input_str(1:4)); %#ok<*ST2NM>
                                mon = str2num(input_str(6:7));
                                day = str2num(input_str(9:10));
                                hr = str2num(input_str(12:13));
                                min = str2num(input_str(15:16));
                                sec = str2num(input_str(18:19));
                                t_start_jd_temp = jday(yr, mon, day, hr, min, sec);
                                % Check start time:
                                if obs_data.number_of_scans > 0 % not the first scan
                                    if (t_start_jd_temp > obs_data.end_of_last_scan_jd) && (t_start_jd_temp < (PARA.endmjd + 2.400000500000000e+006))
                                        t_start_jd = t_start_jd_temp;
                                        add_scan_state = 3.2;
                                    else
                                        fprintf(1, ' ERROR: Start time is not valid. Please try again!\n');
                                    end
                                else % first scan => "obs_data.end_of_last_scan_jd" is empty!
                                    if (t_start_jd_temp >= (PARA.startmjd + 2.400000500000000e+006)) && (t_start_jd_temp < (PARA.endmjd + 2.400000500000000e+006))
                                        t_start_jd = t_start_jd_temp;
                                        add_scan_state = 3.2;
                                    else
                                        fprintf(1, ' ERROR: Start time is not valid. Please try again!\n');
                                    end
                                end
                                
                            else
                                fprintf(1, ' ERROR: Invalid input. Please try again!\n');
                            end 
                        else
                            switch(input_str)
                                case 'y'
                                    if obs_data.number_of_scans == 0 % First scan
                                        t_start_jd = PARA.startmjd + 2.400000500000000e+006;
                                    else % Not the first scan
                                        t_start_jd = obs_data.end_of_last_scan_jd;
                                    end
                                    [year_start, mon_start, day_start, hr_start, min_start, sec_start] = invjday(t_start_jd); % Conversion
                                    fprintf(1, '   => Start:     %4.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\n',  year_start, mon_start, day_start, hr_start, min_start, sec_start);
                                    add_scan_state = 3.2;
                                case 'exit'
                                    add_scan_state = 9;
                                case 'back'
                                    add_scan_state = 2;
                                otherwise
                                    fprintf(1, ' ERROR: Invalid input. Please try again!\n');
                            end
                        end
                        
                    case 3.2 % Enter duration of station based scheduling
                        commandwindow
                        input_str = input(' Duration [min] = ', 's');
                        % Check input:
                        if sum((input_str(1:end) <= '9') & (input_str(1:end) >= '0')) == length(input_str) % Only numbers...
                            stat_based_duration_min_temp = str2double(input_str);
                            t_end_jd_temp = t_start_jd + stat_based_duration_min_temp / (24*60);
                            if (t_end_jd_temp < t_start_jd) && (t_end_jd_temp > (PARA.endmjd + 2.400000500000000e+006))
                                fprintf(1, ' ERROR: Invalid input. Please try again!\n');
                            else
                                t_end_jd = t_end_jd_temp;
                                add_scan_state = 3.3;
                            end
                        else
                            switch(input_str)
                                case 'exit'
                                    add_scan_state = 9;
                                case 'back'
                                    add_scan_state = 2;
                                otherwise
                                    fprintf(1, ' ERROR: Invalid input, please try again!\n');
                            end
                        end
                        
                    case 3.3 % Confirm input (1)
                        [year_start, mon_start, day_start, hr_start, min_start, sec_start] = invjday(t_start_jd); % Conversion
                        [year_end, mon_end, day_end, hr_end, min_end, sec_end] = invjday(t_end_jd); % Conversion
                        duration_min_temp = (t_end_jd - t_start_jd) * (60*24);
                        fprintf(1, ' \n++++ Confirm input and start station based scheduling: ++++ \n');
                        fprintf(1, '   => Start:     %4.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\n',  year_start, mon_start, day_start, hr_start, min_start, sec_start);
                        fprintf(1, '   => End:       %4.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\n',  year_end, mon_end, day_end, hr_end, min_end, sec_end);
                        fprintf(1, '   => Duration:  %1.1f minutes\n',  duration_min_temp);
                        fprintf(1, ' 1   - Start station based scheduling\n');
                        fprintf(1, ' 2   - Change setting\n');
                        add_scan_state = 3.4;
                        
                    case 3.4 % Confirm input (2)
                        commandwindow
                        input_str = input(' Please select: ', 's');
                        satellite_id = [];
                        switch(input_str)
                            case '1' % start station based scheduling
                                add_scan_state = 3.5;
                            case '2' % Change setting
                                add_scan_state = 3.1;
                            case 'exit'
                                add_scan_state = 9; % Exit the "add scan" routine and go back to main menu
                            case 'back'
                                add_scan_state = 3.1; % ##### Select options #####
                            otherwise
                                fprintf(1, ' ERROR: Invalid input, please try again!\n');
                        end
                        
                    case 3.5 % start station based scheduling
                        [sched_data_tmp, obs_data_tmp, error_code, error_msg] = station_based_quasar_scheduling(source, station, PARA, sched_data, obs_data, sat_network_id_list, t_start_jd, t_end_jd);
%                         [sched_data, obs_data, error_code, error_msg] = station_based_quasar_scheduling(source, station, PARA, sched_data, obs_data, quasar_network_id_list, t_start_jd, t_end_jd);
                        if error_code > 0
                            error_msg = ['station_based_quasar_scheduling: ', error_msg];
                            fprintf(1, [' ERROR: ', error_msg, '\n']);
                            add_scan_state = 0;
                        else
                            fprintf(1, ' ...station based scheduling finished und data saved successfully (last scan number: %1.0f)! \n\n', obs_data.number_of_scans);
                            add_scan_state = 3.6; % Add scan to current Schedule (append)
                        end
                        
                    case 3.6 % Ask, if the scheduled scans should be added to the schedule:
                        fprintf(1, ' Do you want to add these scans to the observation schedule? \n');
                        commandwindow
                        input_str = input(' [y or yes] to add them, [n or no] to discard them: ', 's');
                        switch(input_str)
                            case {'y', 'yes'} % add scans
                                sched_data = sched_data_tmp;
                                obs_data = obs_data_tmp;
                                clear sched_data_tmp obs_data_tmp
                                fprintf(1, '  => Scans added! \n');
                                add_scan_state = 9;
                            case {'n', 'no'} % discard scans
                                add_scan_state = 9;
                                clear sched_data_tmp obs_data_tmp
                                fprintf(1, '  => Scans discarded! \n');
                            otherwise
                                fprintf(1, ' ERROR: Invalid input, please try again!\n');
                        end

                        
                    % #### EXIT the "add scan" routine and go back to main menu ####
                    case 9 
                        flag_exit_add_scan = 1; % Exit sub-loop
                        sched_state = 1; % Main menu

                        
                end % switch(case1_state)
                
            end % while(~flag_exit_case_2)
            
            
        % #############################################################
        % ##### Create output                                     #####
        % #############################################################
        case 5
            
            % init.: 
            output_state                = 0;
            flag_exit_output_state      = false;
            
            while(~flag_exit_output_state)
                
                switch(output_state)
                    
                    case 0
                        fprintf(1, '#### Create output for current schedule (exper. name: %s) ####\n', sched_data.exper_name);
                        fprintf(1, 'Nuber of scans: %d\n', sched_data.number_of_scans);
                        if sched_data.number_of_scans == 0
                            fprintf(1, 'WARNING: The current session does not contain any scans! Add scans, before creating output files!\n');
                        end
                        if isempty(sched_data.scan(1).stat(1).epoch(1).ra)
                            fprintf(1, 'WARNING: Tracking data has not been calculated yet!\n');
                        end
                        fprintf(1, ' 1   - Calculate/update satellite tarcking data (step-wise tracking)\n'); % [sched_data, error_code, error_msg] = calc_tracking_data(stat_data, source, sched_data, PARA)
                        fprintf(1, ' 2   - Station dependent VEX files\n');
                        fprintf(1, ' 3   - Combined VEX file\n');
                        fprintf(1, ' 4   - VSO file (station baseline, e.g. for simulations)\n');
                        fprintf(1, ' 5   - VSO file (geocentric baselines, e.g. for correlation)\n');
                        fprintf(1, ' 8   - Keyboard input\n');
                        fprintf(1, ' 9   - Back to main menu\n');
                        fprintf(1,'\n');
                        output_state = 0.1;
                        
                    case 0.1 % User input
                        commandwindow
                        input_str = input(' Please select: ', 's');
                        switch(input_str)
                            case '1'
                                output_state = 1; % Calculate/update satellite tarcking data
                            case '2'
                                output_state = 2; % station dependent VEX files
                            case '4'
                                output_state = 3; % VSO file (station baseline)
                                flag_geocentric_delay = false;
                                pre_sec = 0;
                                post_sec = 0;
                            case '5'
                                output_state = 3; % VSO file (geocentric baselines)
                                flag_geocentric_delay = true;
                            case '3'
                                output_state = 4; % comb. VEX file
                            case '8'
                                output_state = 8; % Keyboard input
                            case {'9', 'back', 'exit'}
                                sched_state = 1;
                                flag_exit_output_state = 1;
                            otherwise
                                fprintf(1, ' ERROR: Invalid input, try again!\n');
                        end
                        
                    % ### Calculate/update satellite tarcking data ###
                    case 1 
                        [sched_data, error_code, error_msg] = calc_tracking_data(stat_data, source, sched_data, PARA);
                        if error_code > 0
                             error_msg = ['calc_tracking_data: ', error_msg];
                             fprintf(1, [' ERROR: ', error_msg, '\n']);
                        else
                            fprintf(1,'   ...tracking data sucessfully calculated!\n\n ');
                        end
                        output_state = 0;

                        
                    % ### Station dependet VEX files ###  
                    case 2
                        % First, check if tracking data is available!
                        if isempty(sched_data.scan(1).stat(1).epoch(1).ra)
                            fprintf(1,'  ERROR: The current schdule file file does not contain tracking data! Calculated the tracking data first!\n\n ');
                        else
                            PARA.VEX_MODE = 'stat_dependent';
                            [error_code, error_msg] = vie_sched_sat_o(sched_data, PARA);
                            if error_code > 0
                                 error_msg = ['vie_sched_sat_o: ', error_msg];
                                 fprintf(1, [' ERROR: ', error_msg, '\n']);
                            else
                                fprintf(1,'   ...Station dependent VEX files successfully written!\n\n ');
                            end
                        end
                        output_state = 0;
                        
                    % ### VSO file ###  
                    case 3
                        fprintf(1, '#### Create VSO file ####\n');
                        if sched_data.number_of_scans > 0 % % Write file, if at least one scan was scheduled:
                            output_state = 3.1;
                            fprintf(1,' Enter the observation interval for satellites [sec]\n');
                            fprintf(1,'  - Each satellite track is fragmented into sub-scans in that interval\n'); 
                            fprintf(1,'  - Enter "0" to skip this option\n'); 
                        else
                            output_state = 0;
                            fprintf(1,'  ERROR: At least one scan has to be added to the observation schedule before an output file can be created!\n\n');
                        end
                        
                    case 3.1 % Enter obs. interval for satellite scans
                        commandwindow
                        input_str = input(' Sat. obs interval [integer sec]: ', 's');
                        [sat_int_sec, status] = str2num(input_str);
                        if (status == 1) % Conversion successfull? yes => numeric value
                            % Check integer value:
                            if mod(sat_int_sec, 1) == 0
                                if flag_geocentric_delay
                                    output_state = 3.2; % Enter pre/post times
                                else
                                    output_state = 3.5; % Write file 
                                end
                            else
                                fprintf(1, ' ERROR: Invalid input: Sat. obs interval has to be an integer value!\n');
                            end
                        else % not a numeric value
                            switch(input_str)
                                case {'back', 'exit'}
                                    output_state = 0; % Back to main menu
                                otherwise
                                    fprintf(1, ' ERROR: Invalid input, try again!\n');
                            end
                        end
                    
                        
                        
                    case 3.2 % Enter pre-observation time [sec]
                        fprintf(1,' Enter pre-observation time [sec] (subtracted from start times of satellite scans)\n');
                        commandwindow
                        input_str = input(' Pre-obs. [integer sec]: ', 's');
                        [pre_sec, status] = str2num(input_str);
                        if (status == 1) % Conversion successfull? yes => numeric value
                            % Check integer value:
                            if mod(pre_sec, 1) == 0
                                output_state = 3.3; % OK
                            else
                                fprintf(1, ' ERROR: Invalid input: Pre-obs. has to be an integer value!\n');
                            end
                        else % not a numeric value
                            switch(input_str)
                                case {'back', 'exit'}
                                    output_state = 0; % Back to main menu
                                otherwise
                                    fprintf(1, ' ERROR: Invalid input, try again!\n');
                            end
                        end
                    
                    case 3.3 % Enter post-observation time [sec]
                        fprintf(1,' Enter post-observation time [sec] (added to end times of satellite scans)\n');
                        commandwindow
                        input_str = input(' Post-obs. [integer sec]: ', 's');
                        [post_sec, status] = str2num(input_str);
                        if (status == 1) % Conversion successfull? yes => numeric value
                            % Check integer value:
                            if mod(post_sec, 1) == 0
                                output_state = 3.5; % OK
                            else
                                fprintf(1, ' ERROR: Invalid input: Pre-obs. has to be an integer value!\n');
                            end
                        else % not a numeric value
                            switch(input_str)
                                case {'back', 'exit'}
                                    output_state = 0; % Back to main menu
                                otherwise
                                    fprintf(1, ' ERROR: Invalid input, try again!\n');
                            end
                        end
                        
                    case 3.5 % Write file
                        filename_vso_out = [sched_data.exper_name, '.vso'];
                        [error_code, error_msg] = write_vso_file(sched_data, sat_int_sec, filepath_vso_out, filename_vso_out, flag_geocentric_delay, pre_sec, post_sec);
                        if error_code > 0
                             error_msg = ['write_vso_file: ', error_msg];
                             fprintf(1, [' ERROR: ', error_msg, '\n']);
                             output_state = 0;
                        else
                            output_state = 0;
                            fprintf(1,'   ...VSO file created: %s%s\n\n ', filepath_vso_out, filename_vso_out);
                        end
                        
                    % ### Comb. VEX file ###
                    case 4
                        % First, check if tracking data is available!
                        if isempty(sched_data.scan(1).stat(1).epoch(1).ra)
                            fprintf(1,'  ERROR: The current schdule file file does not contain tracking data! Calculated the tracking data first!\n\n ');
                        else
                            PARA.VEX_MODE = 'combined';
                            [error_code, error_msg] = vie_sched_sat_o(sched_data, PARA);
                            if error_code > 0
                                 error_msg = ['vie_sched_sat_o: ', error_msg];
                                 fprintf(1, [' ERROR: ', error_msg, '\n']);
                            else
                                fprintf(1,'   ...Combined VEX file successfully written!\n\n ');
                            end
                        end
                        output_state = 0;
                        
                    % ### Keyboard input ###
                    case 8
                        % ===================================================================================================================
                        % => Use this section e.g. to add functions for calculating tracking data accrding to your specific needs!
                        % ===================================================================================================================
                        
                        keyboard
                        
                        % #### Funcitons used to track apod with AzEl files (AuSCOPE) ####
                        
                        % Load/calc. AzEl values for APOD observations:
                        flag_get_apod_azel_values = 0;
                        if flag_get_apod_azel_values
                            keyboard
                            [sched_data, error_code, error_msg] = calc_tracking_data_azel_apod(stat_data, source, sched_data, PARA);  % function used to calculate AzEl tracking data based on TRF state vectors provided in ASCII files
                            if error_code > 0
                                 error_msg = ['calc_tracking_data_apod: ', error_msg];
                                 fprintf(1, [' ERROR: ', error_msg, '\n']);
                            else
                                fprintf(1,'   ...APOD tracking data (AzEl) sucessfully added to the sched_data structure!\n\n ');
                                check_azel(sched_data, PARA);                   % function used to check AzEl tracking data and to create various plots
%                           check_azel_jing(sched_data, PARA)                   % function used to compare AzEl values from sched_data structure with values provided by Jing for APOD
%                           sched_data = write_azel_to_sched_data(sched_data)   % function to directly load AzEl values from ASCII files and write them to the sched_data struct (with corrected cable wrap sections!)
                            end
                        end

                        % Write AzEl tracking file for AuSCOPE antennas:
                        flag_write_azel_tracking_file = 0;
                        if flag_write_azel_tracking_file
                            keyboard
                            [error_code, error_msg] = write_azel_file_auscope(sched_data, PARA);
                            if error_code > 0
                                 error_msg = ['write_azel_file_auscope: ', error_msg];
                                 fprintf(1, [' ERROR: ', error_msg, '\n']);
                            else
                                fprintf(1,'   ...AzEl tracking files written successfully!\n\n ');
                            end
                        end
 
                        
                        % #### Function used for tracking APOD in the 263x sessions (WzWnOn) ####
                        
                        % Load RaDec values from file for APOD observations:
                        flag_get_apod_radec_values = false;
                        if flag_get_apod_radec_values
                            [sched_data, error_code, error_msg] = calc_tracking_data_apod_263(stat_data, source, sched_data, PARA);
                            if error_code > 0
                                 error_msg = ['calc_tracking_data_apod: ', error_msg];
                                 fprintf(1, [' ERROR: ', error_msg, '\n']);
                            else
                                fprintf(1,'   ...APOD tracking data sucessfully added to the sched_data structure!\n\n ');
                            end
                        end
                        
                        % Write snp file:
                        flag_write_snp_file = false;
                        if flag_write_snp_file
                            [error_code, error_msg] = write_snp_file_On(sched_data, PARA);
                            if error_code > 0
                                 error_msg = ['write_snp_file: ', error_msg];
                                 fprintf(1, [' ERROR: ', error_msg, '\n']);
                            else
                                fprintf(1,'   ...APOD snp file sucessfully written!\n\n ');
                            end
                        end
                        
                        % Function to compare RaDecTLE with RaDecJing
%                        check_radec(sched_data, PARA)
                        output_state = 0;
                        
                end % switch(output_state)

            end % while(~flag_exit_output_state)
            

        % #############################################################
        % ##### Cancel current scheduling task                    #####
        % #############################################################
        case 6
            fprintf(1, '#### Cancel current scheduling task ####\n');
            flag_exit_user_input = 1;
        
           

            
        % #############################################################
        % ##### Load/save scheduling data from/to /DATA/LEVEL5/   #####
        % #############################################################
        case 8
            % init.: 
            load_save_state         = 1;
            flag_exit_load_save_state    = false;
            flag_load_old_PARA_file      = false;
            
            while(~flag_exit_load_save_state)
                
                switch(load_save_state)
                    
                    case 1
                        fprintf(1, '#### Load/save scheduling data ####\n');
                        fprintf(1, ' 1   - Load scheduling data from /DATA/LEVEL5/\n');
                        fprintf(1, ' 2   - Load scheduling data from /DATA/LEVEL5/, includign the old sched parameters (var.: PARA)\n');
                        fprintf(1, ' 3   - Save scheduling data to /DATA/LEVEL5/%s\n', sched_data.exper_name);
                        fprintf(1, ' 9   - Back to main menu\n');
                        fprintf(1,'\n');
                        load_save_state = 1.1;
                       
                    case 1.1 % User input
                        commandwindow
                        input_str = input(' Please select: ', 's');
                        switch(input_str)
                            case '1'
                                load_save_state = 2; % Load scheduling data
                            case '2'
                                load_save_state = 2; % Load scheduling data
                                flag_load_old_PARA_file = true;
                            case '3'
                                ver_num_str = '';
                                load_save_state = 3; % Save scheduling data
                            case {'9', 'back', 'exit'}
                                sched_state = 1;
                                flag_exit_load_save_state = 1;
                            otherwise
                                fprintf(1, ' ERROR: Invalid input, try again!\n');
                        end
                        
                        
                    % ### Load scheduling data from /DATA/LEVEL5/ ###
                    case 2 
                        % user input:
                        commandwindow
                        input_str = input(' Enter exper. name: ', 's');
                        switch(input_str)
                            case {'back', 'exit'}
                                load_save_state = 1;
                            otherwise
                                % Check, if a subfolder with the inputted session name exists in /DATA/LEVEL5/:
                                if exist([PARA.pfolder, input_str], 'file')
                                    % Check, if all files which should be loaded exist:
                                    flag_all_input_files_exist = 1;
                                    if ~exist([PARA.pfolder, input_str, '/sched_data.mat'], 'file')
                                        flag_all_input_files_exist = 0;
                                        fprintf(1, '  ERROR: %s does not exist!\n\n', [PARA.pfolder, input_str, '/sched_data.mat']);
                                    end
                                    if ~exist([PARA.pfolder, input_str, '/source.mat'], 'file')
                                        flag_all_input_files_exist = 0;
                                        fprintf(1, '  ERROR: %s does not exist!\n\n', [PARA.pfolder, input_str, '/source.mat']);
                                    end
                                    if ~exist([PARA.pfolder, input_str, '/obs_data.mat'], 'file')
                                        flag_all_input_files_exist = 0;
                                        fprintf(1, '  ERROR: %s does not exist!\n\n', [PARA.pfolder, input_str, '/obs_data.mat']);
                                    end
                                    if ~exist([PARA.pfolder, input_str, '/sched_handles.mat'], 'file')
                                        fprintf(1, '  WARNING: figure handle structure (%s) does not exist!\n', [PARA.pfolder, input_str, '/sched_handles.mat']);
                                    end
                                    if flag_load_old_PARA_file
                                        if ~exist([PARA.pfolder, input_str, '/PARA.mat'], 'file')
                                            flag_all_input_files_exist = 0;
                                            fprintf(1, '  ERROR: %s does not exist!\n\n', [PARA.pfolder, input_str, '/PARA.mat']);
                                        end
                                    end
                                    
                                    % load files required for calcualting the antenna pointing dataand for creating the file, if they exist:
                                    if flag_all_input_files_exist
                                        load([PARA.pfolder, input_str, '/sched_data.mat']);
                                        fprintf(1, '  ...%s loaded\n', [PARA.pfolder, input_str, '/sched_data.mat']);
                                        load([PARA.pfolder, input_str, '/source.mat']);
                                        fprintf(1, '  ...%s loaded\n', [PARA.pfolder, input_str, '/source.mat']);
                                        load([PARA.pfolder, input_str, '/stat_data.mat']);
                                        fprintf(1, '  ...%s loaded\n', [PARA.pfolder, input_str, '/stat_data.mat']);
                                        load([PARA.pfolder, input_str, '/obs_data.mat']);
                                        fprintf(1, '  ...%s loaded\n', [PARA.pfolder, input_str, '/obs_data.mat']);
                                        
                                        if flag_load_old_PARA_file
                                            old_sub_dir_str = PARA.OUTPUT_SUBDIR;
                                            load([PARA.pfolder, input_str, '/PARA.mat']);
                                            fprintf(1, '  ...%s loaded\n', [PARA.pfolder, input_str, '/PARA.mat']);
                                            PARA.OUTPUT_SUBDIR = old_sub_dir_str;
                                            fprintf(1, '     => The current output sub-directory is kept(%s)! \n', PARA.OUTPUT_SUBDIR);
                                        end
                                        
                                        % Close all open figures before loading the saved ones (if available):
                                        % el-plot:
                                        if ~isempty(sched_handles.el_plot)
                                            if ishghandle(sched_handles.el_plot.h_fig)
                                                close(sched_handles.el_plot.h_fig);
                                            end
                                        end
                                        % sky-plots:
                                        if ~isempty(sched_handles.sky_plots)
                                            for i_plot = 1 : length(sched_handles.sky_plots)
                                                if ishghandle(sched_handles.sky_plots(i_plot).h_fig)
                                                    close(sched_handles.sky_plots(i_plot).h_fig);
                                                end
                                            end
                                        end
                                        
                                        if exist([PARA.pfolder, input_str, '/sched_handles.mat'], 'file')
                                            load([PARA.pfolder, input_str, '/sched_handles.mat']);
                                            fprintf(1, '  ...%s loaded\n', [PARA.pfolder, input_str, '/sched_handles.mat']);
                                            if isempty(sched_handles.el_plot)
                                                fprintf(1, '    => No elevation plot available!\n');
                                            end
                                            if isempty(sched_handles.sky_plots)
                                                fprintf(1, '    => No sky plot(s) available!\n');
                                            end
                                        end
                                        
                                        fprintf(1, '  ...%d scans of session %s were loaded successfuly!\n\n', sched_data.number_of_scans, input_str);
                                        % ##### Finish scheduling and create schedule file #####
                                        load_save_state = 1;
                                    else % Error
                                        fprintf(1, '  ERROR: At least one of the input files does not exist in the directory: %s!\n\n', [PARA.pfolder, input_str, '/']);
                                        load_save_state = 1;
                                    end
                                else % Error
                                    fprintf(1, ' ERROR: Sub-folder with the inputted session name does not exist. Please try again!\n\n');
                                    load_save_state = 1;
                                end % exist([PARA.pfolder, input_str], 'file')

                        end % switch(input_str) 
                        
                        
                        
                    % ### Save scheduling data to /DATA/LEVEL5/ ###  
                    case 3
                        % Check, if the "session folder" already exists => If not, create it!
                        out_path_str = [PARA.pfolder, sched_data.exper_name, ver_num_str,'/'];
                        if ~exist(out_path_str, 'file')
                            [error_status] = mkdir(out_path_str);
                            if error_status == 0 % error
                                fprintf(1, 'ERROR: Creating the folder /DATA/LEVEL5/%s/ was not successful\n\n', sched_data.exper_name);
                            else 
                                load_save_state = 3.2;
                            end
                        else
                            load_save_state = 3.1;
                        end
                        
                    case 3.1
                        % Ask, if: 
                        % - The existing data should be overwritten
                        % - Add a version number to the directory name
                        fprintf(1, 'The directory %s already exist! Please select an option:\n', out_path_str);
                        fprintf(1, '  1 - Overwrite the existing data\n');
                        fprintf(1, '  2 - Append version number to new dir. and keep old data\n');
                        fprintf(1, '  9 - Back to main menu\n');
                        fprintf(1, '\n');
                        load_save_state = 3.11;
                        
                    case 3.11
                        input_str = input(' Please select: ', 's');
                        switch(input_str)
                            case '1'
                                load_save_state = 3.2; % Overwrite the existing data
                            case '2'
                                load_save_state = 3.12; % Append version number to new dir. and keep old data
                            case {'9', 'back', 'exit'}
                                load_save_state = 1; % Back to main menu
                            otherwise
                                fprintf(1, ' ERROR: Invalid input, try again!\n');
                        end
                        
                    case 3.12 % Input version number
                        input_str = input(' Please enter version number (integer): ', 's');
                        [ver_num, status] = str2num(input_str);
                        if (status == 1) % Conversion successfull? yes => numeric value
                            ver_num_str = sprintf('-%d', ver_num);
                            load_save_state = 3;
                        else % not a numeric value
                            switch(input_str)
                                case {'9', 'back', 'exit'}
                                    load_save_state = 1; % Back to main menu
                                otherwise
                                    fprintf(1, ' ERROR: Invalid input, try again!\n');
                            end
                        end
                        
                    case 3.2
%                             % Check, if the current schedule contains at least one scan:
%                             if sched_data.number_of_scans > 0
%                         out_path_str = [PARA.pfolder, sched_data.exper_name, ver_num_str,'/'];
                        save([out_path_str, 'sched_data.mat'], 'sched_data');
                        fprintf(1, '  ...%s saved\n', [out_path_str, 'sched_data.mat']);
                        save([out_path_str, 'source.mat'], 'source');
                        fprintf(1, '  ...%s saved\n', [out_path_str, 'source.mat']);
                        save([out_path_str, 'stat_data.mat'], 'stat_data');
                        fprintf(1, '  ...%s saved\n', [out_path_str, 'stat_data.mat']);
                        save([out_path_str, 'PARA.mat'], 'PARA');
                        fprintf(1, '  ...%s saved\n', [out_path_str, 'PARA.mat']);
                        save([out_path_str, 'obs_data.mat'], 'obs_data');
                        fprintf(1, '  ...%s saved\n', [out_path_str, 'obs_data.mat']);
                        save([out_path_str, 'sched_handles.mat'], 'sched_handles');
                        fprintf(1, '  ...%s saved\n', [out_path_str, 'sched_handles.mat']);
                        fprintf(1, '   ...files saved successfully!\n\n');
                        fprintf(1, '  ...%d scans of session %s were saved successfuly!\n\n', sched_data.number_of_scans, sched_data.exper_name);
%                             else
%                                 fprintf(1, 'ERROR: The current schedule has to contain at least one scan before it can be saved!\n\n');
%                             end
                        load_save_state = 1;
                        
                        
                end % switch(load_save_state)
                
            end % while(~flag_exit_load_save_state)
            
        
        % #################################################  
        % ##### Schedule summary and statistics       #####
        % #################################################
        case 10
            % Init.:
            sched_stats_state = 1;
            flag_exit_sched_stats = 0;
            
            while(~flag_exit_sched_stats)
                
                switch(sched_stats_state)
                    
                    case 1
                        fprintf(1, '#### Schedule summary and statistics ####\n');
                        fprintf(1, ' 1   - Write schedule summary to Matlab CW\n');
                        fprintf(1, ' 2   - Write schedule summary file\n');
                        fprintf(1, ' 3   - Print sky-plots\n');
                        fprintf(1, ' 9   - Back to main menu\n');
                        fprintf(1,'\n');
                        sched_stats_state = 1.1;
                        
                    case 1.1 % User input
                        commandwindow
                        input_str = input(' Please select: ', 's');
                        switch(input_str)
                            case '1'
                                sched_stats_state = 2; % Write schedule summary to Matlab CW
                                flag_exit_add_scan = 0;
                            case '2'
                                sched_stats_state = 3; % Write schedule summary file
                            case '3'
                                sched_stats_state = 4; % Print sky-plots
                            case '9'
                                sched_state = 1;
                                flag_exit_sched_stats = 1;
                            case 'back'
                                sched_state = 1;
                                flag_exit_sched_stats = 1;
                            case 'exit'
                                sched_state = 1;
                                flag_exit_sched_stats = 1;
                            otherwise
                                fprintf(1, ' ERROR: Invalid input, try again!\n');
                        end
                        
                    % ##### Write schedule summary to Matlab CW #####
                    case 2
                        % Create schedule summary string:
                        [sched_sum_str, error_code, error_msg] = write_sched_sum_str(sched_data);
                        if error_code > 0
                             error_msg = ['write_sched_sum_str: ', error_msg];
                             fprintf(1, [' ERROR: ', error_msg, '\n']);
                             sched_stats_state = 1.1;
                        else % no error while writing schedule summary string:
                            % Write string to Matlab CW:
                            fprintf(1, '\n');
                            fprintf(1, sched_sum_str);
                            fprintf(1, '\n');
                            sched_stats_state = 1;
                        end

                        
                    % ##### Write schedule summary file #####
                    case 3
                        % Create schedule summary string:
                        [sched_sum_str, error_code, error_msg] = write_sched_sum_str(sched_data);
                        if error_code > 0
                             error_msg = ['write_sched_sum_str: ', error_msg];
                             fprintf(1, [' ERROR: ', error_msg, '\n']);
                             sched_stats_state = 1.1;
                        else % no error while writing schedule summary string:
                            % Write string to file:
                            filename_sched_sum_file_out = [sched_data.exper_name, '_sum.txt'];
                            [error_code, error_msg] = write_string_to_file(sched_sum_str, filepath_sched_sum_file_out, filename_sched_sum_file_out);
                            if error_code > 0
                                 error_msg = ['write_string_to_file: ', error_msg];
                                 fprintf(1, [' ERROR: ', error_msg, '\n']);
                                 sched_stats_state = 1.1;
                            else % no error while writing schedule summary file:
                                fprintf(1,'   ...summary file created: %s%s\n\n ', filepath_sched_sum_file_out, filename_sched_sum_file_out);
                                sched_stats_state = 1;
                            end
                        end
                        
                        
                    % ##### Print sky-plots #####
                    case 4
                        fprintf(1, '#### Create sky plots ####\n');
                        fprintf(1, ' 1   - For the whole schedule\n');
                        fprintf(1, ' 2   - For a defined time window\n');
                        fprintf(1, ' 9   - Back\n');
                        fprintf(1,'\n');
                        sched_stats_state = 4.1;
                        
                    case 4.1 % User input
                        commandwindow
                        input_str = input(' Please select: ', 's');
                        switch(input_str)
                            case '1'
                                sched_stats_state = 4.5; % For the whole schedule
                            case '2'
                                sched_stats_state = 4.6;
                            case '9'
                                sched_stats_state = 1;
                            case 'back'
                                sched_stats_state = 1;
                            case 'exit'
                                sched_stats_state = 1;
                            otherwise
                                fprintf(1, ' ERROR: Invalid input, try again!\n');
                        end
                        
                    case 4.5 % For the whole schedule
                        t_min_jd          = sched_data.t_nominal_start_jd;
                        t_max_jd          = sched_data.t_nominal_end_jd;
                        sched_stats_state = 4.7;
                        
                    case 4.60 % For a defined time window
                        fprintf(1, '   Enter start time\n');
                        fprintf(1,'    => Input Format <yyyy-mm-dd HH:MM:SS>\n');
                        sched_stats_state = 4.61;
                        
                    case 4.61 % Scan start time: Enter it
                        commandwindow
                        input_str = input(' t_start = ', 's');
                        if (length(input_str) == 19)
                            [flag_t_epoch_jd_ok, t_epoch_jd_ok] = check_t_input_19char(input_str, sched_data, PARA);
                            if flag_t_epoch_jd_ok
                                t_min_jd = t_epoch_jd_ok;
                                sched_stats_state = 4.62;
                            end
                        else
                            switch(input_str)
                                case 'exit'
                                    sched_stats_state = 1;
                                case 'back'
                                    sched_stats_state = 4;
                                otherwise
                                    fprintf(1, ' ERROR: Invalid input. Please try again!\n');
                            end
                        end
                        
                    case 4.62 % For a defined time window
                        fprintf(1, '   Enter end time\n');
                        fprintf(1,'    => Input Format <yyyy-mm-dd HH:MM:SS>\n');
                        sched_stats_state = 4.63;
                        
                    case 4.63 % Scan end time: Enter it
                        commandwindow
                        input_str = input(' t_end = ', 's');
                        if (length(input_str) == 19)
                            [flag_t_epoch_jd_ok, t_epoch_jd_ok] = check_t_input_19char(input_str, sched_data, PARA);
                            if flag_t_epoch_jd_ok
                                t_max_jd = t_epoch_jd_ok;
                                sched_stats_state = 4.7;
                            end
                        else
                            switch(input_str)
                                case 'exit'
                                    sched_stats_state = 1;
                                case 'back'
                                    sched_stats_state = 4;
                                otherwise
                                    fprintf(1, ' ERROR: Invalid input. Please try again!\n');
                            end
                        end
                        
                    case 4.7 % Create Sky Plots
                        % Options
                        flag_print_sattracks    = true;
                        flag_print_skycov       = false;
                        flag_save_pdf           = true;
                        flag_save_fig           = true;
                        [error_code, error_msg] = create_skyplot_sum(sched_data, PARA, flag_print_sattracks, flag_print_skycov, flag_save_fig, flag_save_pdf, t_min_jd, t_max_jd);
                        if error_code > 0
                             error_msg = ['create_skyplot_sum: ', error_msg];
                             fprintf(1, [' ERROR: ', error_msg, '\n\n']);
                        end
                        sched_stats_state = 1;
  
                end % switch(sched_stats_state)
                
            end % while(~flag_exit_sched_stats)
            
            
            
        % #############################################################  
        % ##### Automatic Scheduling (satellites & quasars)       #####
        % #############################################################
        case 11  
            
            % ##### Init #####
            flag_exit_auto_sched  = 0;
            auto_sched_state      = 1;
            
            % ##### Main loop #####
            while(~flag_exit_auto_sched)
                
                switch(auto_sched_state)
                    
                    case 1            
                        fprintf(1, '#### Automatic scheduling of satellite & quasar scans ####\n');
                        fprintf(1, ' 1   - For the whole session time\n');
                        fprintf(1, ' 2   - For a defined time window\n');
                        fprintf(1, ' 9   - Back to main menu\n');
                        fprintf(1,'\n');
                        auto_sched_state = 1.1;
                        
                    case 1.1 % User input
                        commandwindow
                        input_str = input(' Please select: ', 's'); 
                        switch(input_str)
                            case '1'
                                auto_sched_state = 2; % For the whole session time
                            case '2'
                                auto_sched_state = 3; % For a defined time window
                            case '9'
                                sched_state = 1;
                                flag_exit_auto_sched = 1;
                            case 'back'
                                sched_state = 1;
                                flag_exit_auto_sched = 1;
                            case 'exit'
                                sched_state = 1;
                                flag_exit_auto_sched = 1;
                            otherwise
                                fprintf(1, ' ERROR: Invalid input, try again!\n');
                        end
                        
                    case 2 % For the whole session time
                        t_start_jd = PARA.startmjd + mjd2jd;
                        t_end_jd   = PARA.endmjd + mjd2jd;
                        auto_sched_state = 5;
                        
                    case 3 % For a defined time window (start time)
                        fprintf(1, '   Enter start time\n');
                        fprintf(1,'    => Input Format <yyyy-mm-dd HH:MM:SS>\n');
                        auto_sched_state = 3.1;
                        
                    case 3.1 % Schedule start time: Enter it
                        commandwindow
                        input_str = input(' t_start = ', 's');
                        if (length(input_str) == 19)
                            [flag_t_epoch_jd_ok, t_epoch_jd_ok] = check_t_input_19char(input_str, sched_data, PARA);
                            if flag_t_epoch_jd_ok
                                t_start_jd = t_epoch_jd_ok;
                                auto_sched_state = 3.2;
                            end
                        else
                            switch(input_str)
                                case 'exit'
                                    auto_sched_state = 1;
                                case 'back'
                                    auto_sched_state = 1;
                                otherwise
                                    fprintf(1, ' ERROR: Invalid input. Please try again!\n');
                            end
                        end
                        
                    case 3.2 % For a defined time window
                        fprintf(1, '   Enter end time\n');
                        fprintf(1,'    => Input Format <yyyy-mm-dd HH:MM:SS>\n');
                        auto_sched_state = 3.3;
                        
                    case 3.3 % Scan end time: Enter it
                        commandwindow
                        input_str = input(' t_end = ', 's');
                        if (length(input_str) == 19)
                            [flag_t_epoch_jd_ok, t_epoch_jd_ok] = check_t_input_19char(input_str, sched_data, PARA);
                            if flag_t_epoch_jd_ok
                                t_end_jd = t_epoch_jd_ok;
                                auto_sched_state = 5;
                            end
                        else
                            switch(input_str)
                                case 'exit'
                                    auto_sched_state = 1;
                                case 'back'
                                    auto_sched_state = 3.1;
                                otherwise
                                    fprintf(1, ' ERROR: Invalid input. Please try again!\n');
                            end
                        end
                        
                    % ##### Start auto. scheduling #####
                    case 5
                        [sched_data, obs_data, error_code, error_msg] = auto_sched(stat_data, source, PARA, sched_data, obs_data, t_start_jd, t_end_jd, quasar_network_id_list, sat_network_id_list);
                        if error_code > 0
                             error_msg = ['auto_sched: ', error_msg];
                             fprintf(1, [' ERROR: ', error_msg, '\n']);
                             auto_sched_state = 1.1;
                        else % no error => Write statud msg to CW
                            fprintf(1, '\n');
                            fprintf(1, '....auto_sched finished!');
                            fprintf(1, '\n');
                            auto_sched_state = 1;
                            
                            % ### Update Sky-plots and elevation plots ###
                            % Set all pointers to end of last obs.
                            [sched_handles, error_code, error_msg] = update_sky_plots_pointer(sched_handles, station_id_list, obs_data, 4, PARA);
                            if error_code > 0
                                error_msg = ['update_sky_plots_pointer: ', error_msg];
                                fprintf(1, [' ERROR: ', error_msg, '\n']);
                            else 
                                [sched_handles, error_code, error_msg] = update_elevation_plot_pointer(sched_handles, obs_data, 4, [], PARA, station_id_list);
                                if error_code > 0
                                    error_msg = ['update_elevation_plot_pointer: ', error_msg];
                                    fprintf(1, [' ERROR: ', error_msg, '\n']);
                                else
                                    if isempty(obs_data.end_of_last_scan_jd)
                                        t_epoch_jd_temp = PARA.startmjd + 2.400000500000000e+006; % first scan of session
                                    else
                                        t_epoch_jd_temp = obs_data.end_of_last_scan_jd;
                                    end
                                    flag_plot_sources = 0;
                                    [sched_handles, error_code, error_msg] = update_sky_plots(sched_handles, stat_data, station_id_list, source, PARA, obs_data, flag_plot_sources, t_epoch_jd_temp);
                                    if error_code > 0
                                        error_msg = ['update_sky_plots: ', error_msg];
                                        fprintf(1, [' ERROR: ', error_msg, '\n']);
                                    else
                                        % Plot quasar markers (highlight):
                                        [sched_handles, error_code, error_msg] = update_sky_plots_highlight_quasars(sched_handles, stat_data, station_id_list, source, t_epoch_jd_temp, [], 2);
                                        if error_code > 0
                                             error_msg = ['update_sky_plots_highlight_quasars: ', error_msg];
                                             fprintf(1, [' ERROR: ', error_msg, '\n']);
                                        end
                                    end
                                end
                            end
 
                        end % if error_code > 0
            
                end % switch(auto_sched_state)
                
            end % while(~flag_exit_sched_stats)
            
            
            
        % #############################################################  
        % ##### Change Settings                                   #####
        % #############################################################
        case 20
            
            % ##### Init #####
            flag_exit_settings       = 0;
            settings_state           = 0;
            
            % ##### Main loop #####
            while(~flag_exit_settings)
                
                switch(settings_state)
                    
                    case 0           
                        fprintf(1, '#### Change scheduling settings ####\n');
                        fprintf(1, ' 1   - Show/change current satellite observation station network\n');
                        fprintf(1, ' 2   - Recover original satellite observation station network\n');
                        fprintf(1, ' 3   - Deactivate to check the start time of the next scan\n');
                        fprintf(1, ' 4   - Use alternative cable wrap section for defined antennas\n');
                        fprintf(1, ' 5   - Set Time constants added between slewing and obs. start to zero\n');
%                         fprintf(1, ' 5   - Change antenna reposition interval (NOT AVAILABLE YET!)\n');
                        fprintf(1, ' 9   - Back to main menu\n');
                        fprintf(1,'\n');
                        settings_state = 0.1;
                        
                    case 0.1 % User input
                        commandwindow
                        input_str = input(' Please select: ', 's'); 
                        switch(input_str)
                            case '1'
                                settings_state = 1; 
                            case '2'
                                settings_state = 2;
                            case '3'
                                settings_state = 3;
                            case '4'
                                settings_state = 4;
                            case '5'
                                settings_state = 5;
                            case '9'
                                settings_state = 0;
                                flag_exit_settings = 1;
                                sched_state = 1;
                            case 'back'
                                settings_state = 0;
                                flag_exit_settings = 1;
                            case 'exit'
                                settings_state = 0;
                                flag_exit_settings = 1;
                                sched_state = 1;
                            otherwise
                                fprintf(1, ' ERROR: Invalid input, try again!\n');
                        end
                        
                    % ##### Choose new station-subnet for SATELLITE observations #####
                    case 1
                        fprintf(1, ' Current satellite obs. network [<ID>: <station name>]: \n');
                        for stat_id = sat_network_id_list
                            fprintf(1, '  - %3d: %s\n', stat_id, stat_data.stat(stat_id).name);
                        end                  

                        fprintf(1, ' Enter new station ID list (vector format: [<ID_2>, <ID_2>, ... ]) or enter "back": \n');
                        commandwindow
                        input_str = input('station ID vector: ', 's');
                        % Check input:
                        [sat_network_id_list_new, status] = str2num(input_str);
                        if (status == 1) % Conversion successfull? yes => numeric value
                            flag_error = false;
                            if (size(sat_network_id_list_new,2) <= size(sat_network_id_list_orig,2)) && (size(sat_network_id_list_new,1) == 1)
                                try
                                    fprintf(1, '\n');
                                    fprintf(1, '  => New satellite obs. network [<ID>: <station name>]: \n');
                                    for stat_id = sat_network_id_list_new
                                        fprintf(1, '    - %3d: %s\n', stat_id, stat_data.stat(stat_id).name);
                                    end
                                catch
                                    flag_error = true;
                                end
                            else
                                flag_error = true;
                            end
                            if ~flag_error
                                settings_state = 1.1;
                            else
                                fprintf(1, ' ERROR: Invalid input, please try again!\n');
                            end
                        else
                            switch(input_str)
                                case 'back'
                                    settings_state = 0;
                                otherwise
                                    fprintf(1, ' ERROR: Invalid input, please try again!\n');
                            end
                        end
                        
                    % Confirm new station ID list:
                    case 1.1
                        fprintf(1, '\n');
                        commandwindow
                        input_str = input(' Confirma new station ID list? (y=yes, n=no) : ', 's');
                        switch(input_str)
                            case {'n' 'no'}
                                settings_state = 0;
                            case {'y' 'yes'}
                                fprintf(1, ' New station ID list confirmed!\n');
                                sat_network_id_list = sat_network_id_list_new;
                                clear sat_network_id_list_new;
                                flag_sat_network_id_list_changed = 1;
                                settings_state = 0;
                            case 'exit'
                                settings_state = 0;
                            case 'back'
                                settings_state = 0;
                            otherwise
                                fprintf(1, ' ERROR: Invalid input, please try again!\n');
                        end
                        
                    % ##### Recover original satellite observation station network #####
                    case 2
                        sat_network_id_list = sat_network_id_list_orig;
                        fprintf(1, '\n');
                        fprintf(1, '  => Orig. satellite observation station network recovered: \n');
                        for stat_id = sat_network_id_list
                            fprintf(1, '    - %3d: %s\n', stat_id, stat_data.stat(stat_id).name);
                        end
                        fprintf(1, '\n');
                        settings_state = 0;
                        
                        
                    % ##### Deactivate to check the scan start time for the next scan #####
                    case 3
                        fprintf(1, ' Current status: \n');
                        if PARA.CHECK_T_START
                            fprintf(1, ' => Check start time of next scan: YES \n');
                        else
                            fprintf(1, ' => Check start time of next scan: NO \n');
                        end
                                        
                        fprintf(1, ' Disable check for next scan?\n');
                        commandwindow
                        input_str = input('Input y for YES, otherwise NO: ', 's');
                        % Check input:
                        switch(input_str)
                            case {'y' 'yes'}
                                PARA.CHECK_T_START = false;
                                fprintf(1, ' Scan start checks are disables for the next scan!\n');
                            otherwise
                                PARA.CHECK_T_START = true;
                                fprintf(1, ' Scan start checks are active!\n');
                        end
                        settings_state = 0;
                        
 
                    % ##### Select alternative cable wrap section for defined antennas #####
                    case 4
                        fprintf(1, ' Available antennas [<ID>: <station name>, => alt. cable wrap selected?]: \n');
                        for i_stat = 1:length(stat_data.stat)
                            if isfield(stat_data.stat(i_stat), 'flag_alt_cable_wrap')
                                if stat_data.stat(i_stat).flag_alt_cable_wrap
                                    tmp_str = '=> alternative cable wrap selected!';
                                else
                                    tmp_str = '';
                                end
                            else
                                tmp_str = '';
                            end
                            fprintf(1, '  - %3d: %8s, %s\n', i_stat, stat_data.stat(i_stat).name, tmp_str);
                        end                  

                        fprintf(1, ' Define stations where an alternative cable wrap section should be used by entering a station ID list\n');
                        fprintf(1, '    - vector format: [<ID_2>, <ID_2>, ... ])or enter "[]" to clear the current settings: \n');
                        commandwindow
                        input_str = input('station ID vector: ', 's');
                        % Check input:
                        [stat_id_list, status] = str2num(input_str);
                        if (status == 1) % Conversion successfull? yes => numeric value
                            flag_error = false;
                            if ((size(stat_id_list,2) <= length(stat_data.stat)) && (size(stat_id_list,1) == 1)) || isempty(stat_id_list)
                                try
                                    fprintf(1, '\n');
                                    fprintf(1, ' New status [<ID>: <station name>, => alt. cable wrap selected?]: \n');
                                    for i_stat = 1:length(stat_data.stat)
                                        if sum(i_stat == stat_id_list) == 1
                                            tmp_str = '=> alternative cable wrap selected!';
                                        else
                                            tmp_str = '';
                                        end
                                        fprintf(1, '  - %3d: %8s, %s\n', i_stat, stat_data.stat(i_stat).name, tmp_str);
                                    end 
                                catch
                                    flag_error = true;
                                end
                            else
                                flag_error = true;
                            end
                            if ~flag_error
                                settings_state = 4.1;
                            else
                                fprintf(1, ' ERROR: Invalid input, please try again!\n');
                            end
                        else
                            switch(input_str)
                                case 'back'
                                    settings_state = 0;
                                otherwise
                                    fprintf(1, ' ERROR: Invalid input, please try again!\n');
                            end
                        end
                        
                    case 4.1
                        fprintf(1, '\n');
                        fprintf(1, ' WARNING: When the cable wrap preference is changed manually, cable wrap information has to be written to the (VEX) schedule files!\n');
                        commandwindow
                        input_str = input(' Confirm new setting? (y=yes, n=no) : ', 's');
                        switch(input_str)
                            case {'n' 'no'}
                                settings_state = 0;
                            case {'y' 'yes'}
                                fprintf(1, ' New cable wrap setings confirmed!\n');
                                for i_stat = 1:length(stat_data.stat)
                                    if sum(i_stat == stat_id_list) == 1
                                        stat_data.stat(i_stat).flag_alt_cable_wrap = true;
                                    else
                                       stat_data.stat(i_stat).flag_alt_cable_wrap = false;
                                    end
                                    fprintf(1, '  - %3d: %8s, %s\n', i_stat, stat_data.stat(i_stat).name, tmp_str);
                                end 
                                settings_state = 0;
                            case 'exit'
                                settings_state = 0;
                            case 'back'
                                settings_state = 0;
                            otherwise
                                fprintf(1, ' ERROR: Invalid input, please try again!\n');
                        end
                        
                    % ##### Set Time constants added between slewing and obs. start to zero #####
                    % => PARA.SOURCE, PARA.TAPETM, PARA.IDLE, PARA.CALIBRATION
                    case 5
                        fprintf(1, ' Current status: \n');
                        if PARA.DISABLE_SLEW_CONST
                            fprintf(1, ' => PARA.SOURCE, PARA.TAPETM, PARA.IDLE, PARA.CALIBRATION set to zero: YES \n');
                        else
                            fprintf(1, ' => PARA.SOURCE, PARA.TAPETM, PARA.IDLE, PARA.CALIBRATION set to zero: NO \n');
                        end
                                        
                        fprintf(1, ' Set Time constants added between slewing and obs. start to zero?\n');
                        commandwindow
                        input_str = input('Input y for YES, otherwise NO: ', 's');
                        % Check input:
                        switch(input_str)
                            case {'y' 'yes'}
                                PARA.DISABLE_SLEW_CONST = true;
                                fprintf(1, ' Time constants added between slewing and obs. start set to zero!\n');
                            otherwise
                                PARA.DISABLE_SLEW_CONST = false;
                                fprintf(1, ' Time constants added between slewing and obs. start enabled (default)!\n');
                        end
                        settings_state = 0;
                    
                        
            
                end % switch(auto_sched_state)
                
            end % while(~flag_exit_settings)
            
     
    end % switch(sched_state)
    
end % while(~flag_scheduling_finished)

         
return;
end
         
    
%% ##### Sub routines #####

function [flag_t_epoch_jd_ok, t_epoch_jd_ok] = check_t_input_19char(input_str, sched_data, PARA)
% This function checks the time input with the format: <yyyy-mm-dd HH:MM:SS>
% Checks:
%  1.) Format (<yyyy-mm-dd HH:MM:SS>)
%  2.) Epoch intersect with existing schedule?
%  3.) Epoch within the defined session time?

    % ### Init.: ###
    mjd2jd = 2.400000500000000e+006;
    flag_t_epoch_jd_ok = false;
    flag_within_existing_schedule = false;
    t_epoch_jd_ok      = 0;

    % ##### 1.) Checks the input Format: <yyyy-mm-dd HH:MM:SS> #####
    if (    (sum((input_str(1:4) <= '9') & (input_str(1:4) >= '0')) == 4)       && ... % Check input Format!
            (sum((input_str(6:7) <= '9') & (input_str(6:7) >= '0')) == 2)       && ...
            (sum((input_str(9:10) <= '9') & (input_str(9:10) >= '0')) == 2)     && ...
            (sum((input_str(12:13) <= '9') & (input_str(12:13) >= '0')) == 2)   && ...
            (sum((input_str(15:16) <= '9') & (input_str(15:16) >= '0')) == 2)   && ...
            (sum((input_str(18:19) <= '9') & (input_str(18:19) >= '0')) == 2)   && ...
            (input_str(5) == '-')                                               && ...
            (input_str(8) == '-')                                               && ...
            (input_str(11) == ' ')                                              && ...
            (input_str(14) == ':')                                              && ...
            (input_str(17) == ':')                                                      )   
        % Convert input to JD:
        yr = str2num(input_str(1:4)); %#ok<*ST2NM>
        mon = str2num(input_str(6:7));
        day = str2num(input_str(9:10));
        hr = str2num(input_str(12:13));
        min = str2num(input_str(15:16));
        sec = str2num(input_str(18:19));
        t_epoch_jd_temp = jday(yr, mon, day, hr, min, sec);
        
        % ##### 2.) Check if t_epoch_jd_temp do not intersect with the existing schedule #####
        if ~isempty(sched_data.t_nominal_start_jd) && ~isempty(sched_data.t_nominal_end_jd) % If at least one scan has been scheduled already
            if (t_epoch_jd_temp >= sched_data.t_nominal_start_jd) && (t_epoch_jd_temp <= sched_data.t_nominal_end_jd)
                flag_within_existing_schedule = true;
                fprintf(1, ' Input epoch intersect with the existing schedule. Please try again!\n');
            end
        end
        
        % ##### 3.) Check, if t_epoch_jd_temp is within the defined session time #####
        if ~flag_within_existing_schedule
            if (t_epoch_jd_temp >= (PARA.startmjd + mjd2jd)) && (t_epoch_jd_temp <= (PARA.endmjd + mjd2jd))
                flag_t_epoch_jd_ok = true;
                t_epoch_jd_ok      = t_epoch_jd_temp;
            else
                fprintf(1, ' Input epoch is out of the scheduled session. Please try again!\n');
            end
        end 
    else
        fprintf(1, ' ERROR: Invalid input. Please try again!\n');
    end
    
end

function highlight_satellite( sched_handles, station_id_list, satellite_id, highlight )
     if highlight
         lw = 3.0;
     else
         lw = 0.5;
     end
     for sta_id = station_id_list
         sched_handles.el_plot.stat(sta_id).sat(satellite_id).h_sat_el.LineWidth = lw;
         sched_handles.el_plot.stat(sta_id).sat(satellite_id).h_sat_el_obs_times.LineWidth = lw;
         sched_handles.sky_plots(sta_id).sat(satellite_id).h_plot.LineWidth = lw;
         
     end
end
























