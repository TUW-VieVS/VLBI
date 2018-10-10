% #########################################################################
% #     vie_sched_sat
% #########################################################################
% DESCRIPTION
%   Calculate the schedule for VLBI satelite observations.
%
% CREATED  
%   2013 Okt     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
%   - load_tle_file
%   - istation
%   - get_trf_data
%   - tle_propagation
%   - get_station_data
%   - calc_pointing_data
%   - find_orbit_events
%   - find_observation_times
%   - print_available_sat_obs_periods
%   - overpass_softcopy
%   - create_elevation_plot
%   - create_skyplot
%   - scheduler_interface
%   - save_pdf_a4_landscape
%
%
% INPUT
% - station             - station data structure (for quasar observations)
% - source              - source data structure
% - PARA                - Global Parameter structure.
% - INFILE              - File path structure
% - obsmode             - Observation mode structure
%
%
% OUTPUT
% - sched_data          - scheduling data structure
%
% CHANGES:
%   - 2014-01-19: A. Hellerschmied: Error treatment added.
%   - 2014-01-19: A. Hellerschmied: Two procedures added:
%                   - manual_scheduling_cw_input
%                   - calc_stepwise_antenna_pointing
%   - 2014-01-30: A. Hellerschmied: case 777 and case 666 added
%   - 2014-01-30: A. Hellerschmied: Preallocation of sched_data within this
%                 procedure.
%   - 2015-02-11: A. Hellerschmied: Renamed it "vie_sched_sat" and revised
%       the function.
%   - 2015-06-19, A. Hellerschmied: Updated for VieVS 2.3
%   - 2015-07-20: A. Hellerschmied: function "load_tle_file.m" added. 
%   - 2015-08-18: A. Hellerschmied: Possibility added to plot violations of axis limits in the elevation plots (if flag PARA.PLOT_AXIS_LIMITS is set here)
%   - 2015-08-24, A. Hellerschmied: Function-call of "scheduler_interface" updated ("source" and "stat_data" added as output arguments)
%   - 2015-11-19, A. Hellerschmied: "sched_data" structure preallocation modified
%   - 2016-06-14, A. Hellerschmied: Calculation of tracking data ("calc_tracking_data.m") has been moved to the scheduler interface!  
%   - 2016-07-13, A. Hellerschmied: Print elev. plot to pdf file changed for MATLAB2016a ('-fillpage' added)
%   - 2016-11-25, A. Hellerschmied: - Additional input arguments for function "print_available_sat_obs_periods.m" added. 
%                                   - "sched_data.number_of_vie_sched_calls" added
%   - 2016-11-30, A. Hellerschmied: - New output argument: process _list
%                                   - Use "save_pdf_a4_landscape.m" to save elevation plot sched sub-dir.
%   - 2016-12-22, A. Hellerschmied: - PARA.INIT_PROP_INTERVAL in [sec] instead of [min]
%                                   - Changed TLE directory from "/CATALOGS/TLE/" to "/ORBIT/TLE/" in the VieVS root dir.


function [sched_data, process_list] = vie_sched_sat(station, source, PARA, INFILE, obsmode)

% Init:
sched_state = 0.0;
error_code = 0;     % 0 = no errrors
error_msg = '';     % '' = no errors
process_list = [];  %

% preallocation:
sched_data = struct('number_of_scans', [], 'TLE_filename', [], 'TLE_path', [], 't_nominal_start_jd', [], 't_nominal_end_jd', [], 'exper_name', [], 'flag_write_output_file', [], 'scan', [], 'ra_dec_epoch', [], 'number_of_vie_sched_calls', 0);
sched_data.stat = struct('name', [], 'trf_x', [], 'trf_y', [], 'trf_z', [], 'trf_vx', [], 'trf_vy', [], 'trf_vz', [], 'trf_epoch', [], 'trf_label', [], 'stat_label', []);
sched_data.scan = struct('obs_type', [], 'number_of_stat', [], 'stat', [], 'repos_int_sec', [], 'quasar_name', [], 'quasar_id', [], 'sat_id', [], 'sat_name', [], 'sat_number', [], 'flag_record_scan', [], 'number_of_epochs', [], 't_start_jd', [], 't_end_jd', []);
sched_data.scan.stat = struct('stat_id', [], 'epoch', [], 'end', [], 'start', [], 'last_active_scan', [], 'duration_sec', [], 'slew_time_sec', [], 'sky_cov_num', []);
sched_data.scan.stat.start = struct('jd', [], 'az', [], 'el', [], 'ha', [], 'dc', [], 'un_az', []);
sched_data.scan.stat.end = struct('jd', [], 'az', [], 'el', [], 'ha', [], 'dc', [], 'un_az', []);
sched_data.scan.stat.epoch = struct('jd', [], 'ra', [], 'dec', [], 'un_az', [], 'el', []);

sched_data.number_of_scans = 0;

sched_handles = struct('el_plot', [], 'sky_plots', [], 'color_markers_scan_start', [], 'color_markers_scan_end', [], 'color_markers_end_of_last_scan', []);

% Const.:
mjd2jd = 2.400000500000000e+006;

% Plot Parameters (Overpass Plot):
PARA.PLOT_SUN_DIST = 1;
PARA.PLOT_AZ_RATE = 1;
PARA.PLOT_EL_RATE = 1;
PARA.PLOT_AXIS_LIMITS = 1;

% observation mode:
PARA.obsmode = obsmode;

% Output options:
sched_data.flag_write_output_file = PARA.SCHEDULE_WRITE_VEX;

% Load data:
%   If a different station network for satellites was selected:
%   - station_quasar => for quasar observations
%   - station_sat => for stellite observations
load([PARA.pfolder 'satellite_str.mat'], 'satellite_str');
station_sat = station;
if PARA.USE_SEPARATE_STAT_NWW_FOR_SATS
    load([PARA.pfolder 'stanet_sat_str.mat'], 'stanet_sat_str');
    station_quasar = station_sat;
    clear station_sat;
    fprintf(1,'1. Read station data for satellite observation network\n');
    [station_sat] = istation(stanet_sat_str, INFILE, PARA); % Get station data from catalogs
else 
    station_quasar = station; 
end


% Color setting for Graphics:
sched_handles.color_markers_scan_start = [1 0 0]; % red
sched_handles.color_markers_scan_end = [0.5430 0 0]; % dark red 
sched_handles.color_markers_end_of_last_scan = [0 0.5 0]; 
sched_handles.color_current_epoch_marker = [0 0 0.5]; 
sched_handles.color_quasars = [1 0.5 0.2]; 
sched_handles.color_highlighted_quasars = [1 0.5 0.5]; 

% Output options:
sched_data.flag_write_output_file = PARA.SCHEDULE_WRITE_VEX;

% SGP4 prediction parameter setup:
PARA.SGP4_GRAV_CONST            = 72;  % WGS72
PARA.SGP4_WRITE_FILE            = 0;   % write output *.txt file (only ECI coordinates!)
PARA.SGP4_PATH_OUT              = '';
PARA.SGP4_FILENAME_OUT          = '';
PARA.SGP4_VERIFICATION_MODE     = 0;
%  grav_const = 72;    % WGS72
%  write_file = 0;     % write output *.txt file (only ECI coordinates!)
%  verification_mode = 0;
%  path_out = '';
%  filename_out = '';
flag_print_whole_overpass = PARA.PRINT_ORBIT_TIMESERIES;

% TLE setup:
PARA.TLE_FILEPATH = [PARA.pvievs, 'ORBIT/TLE/'];

sched_data.TLE_filename = PARA.TLE_FILENAME;
sched_data.TLE_path = PARA.TLE_FILEPATH;


% Time settings:
t_start = PARA.startmjd + mjd2jd;
t_stop = PARA.endmjd + mjd2jd;
delta_t = PARA.INIT_PROP_INTERVAL / 60; % [min]

             
fprintf(1,' ################# vie_sched_sat ##########################\n');

% ##### Start main loop #####
while (1)

    switch(sched_state)
        
        case 0 % ++++ Load TLE file ++++
            [PARA, error_code, error_msg] = load_tle_file(PARA);
            if (error_code == 0)
                sched_state = 1;
            else
                sched_state = 999;
                error_msg = ['load_tle_file: ', error_msg];
            end

        case 1 % ++++ Get TRF site data (position, velocity, ref.-epoch) (Same obs. network for satellites/quasars) ++++
            [station_sat, error_code, error_msg] = get_trf_data(station_sat, PARA);
            if (error_code == 0)
                % Separate obs. networks for satellites/quasars
                if PARA.USE_SEPARATE_STAT_NWW_FOR_SATS 
                    sched_state = 1.1;
                % Same obs. network for satellites/quasars
                else
                    sched_state = 2;
                end
            else
                sched_state = 999;
                error_msg = ['get_trf_data: ', error_msg];
            end
         
        case 1.1 % ++++ Get TRF site data (position, velocity, ref.-epoch) (Separate obs. networks for satellites/quasars) ++++
            [station_quasar, error_code, error_msg] = get_trf_data(station_quasar, PARA);
            if (error_code == 0)
                sched_state = 2;
            else
                sched_state = 999;
                error_msg = ['get_trf_data: ', error_msg];
            end

        case 2 % ++++ Initial Orbit Propagation for choosen satellites ++++
            [sat_data, error_code, error_msg] = tle_propagation(t_start, t_stop, delta_t, PARA, satellite_str);
            % [sat_data2, error_code, error_msg] = TLE_propagation_old(PARA.TLE_FILEPATH, PARA.TLE_FILENAME, t_start, t_stop, delta_t, grav_const, write_file, path_out, filename_out, verification_mode, satellite_str);
            if (error_code == 0)
                sched_state = 3;
            else
                sched_state = 999;
                error_msg = ['tle_propagation: ', error_msg];
            end
            
        case 3 % ++++ Get Station data and preallocate "stat_data" structure ++++
            [stat_data, error_code, error_msg] = get_station_data(station_sat, station_quasar, PARA, INFILE);
            if (error_code == 0)
                sched_state = 4;
            else
                sched_state = 999;
                error_msg = ['get_station_data: ', error_msg];
            end
            
        case 4 % ++++ Calculation of initial antenna pointing data ++++
            [stat_data, error_code, error_msg] = calc_pointing_data(sat_data, stat_data);
            if (error_code == 0)
                sched_state = 5;
            else
                sched_state = 999;
                error_msg = ['calc_pointing_data: ', error_msg];
            end
            
        case 5 % ++++ Find Satellite Orbit Events  ++++
            [stat_data, error_code, error_msg] = find_orbit_events(stat_data, PARA);
            if (error_code == 0)
                sched_state = 6;
            else
                sched_state = 999;
                error_msg = ['find_orbit_events: ', error_msg];
            end
            
        case 6 % ++++ Find possible Observations Times  ++++
            [obs_data, error_code, error_msg] = find_observation_times(stat_data, PARA);
            if (error_code == 0)
                sched_state = 6.1;
            else
                sched_state = 999;
                error_msg = ['observation_times: ', error_msg];
            end

        case 6.1 % ++++ Print available observbation periods for all satellites to the CW  ++++
            flag_from_network = true;
            flag_from_stations = true;
            [error_code, error_msg] = print_available_sat_obs_periods(obs_data, flag_from_network, flag_from_stations);
            if (error_code == 0)
                sched_state = 7;
            else
                sched_state = 999;
                error_msg = ['print_available_sat_obs_periods: ', error_msg];
            end
            
        case 7 % ++++ Print Overpass Data to Command Window  ++++
            [error_code, error_msg] = overpass_softcopy(stat_data, flag_print_whole_overpass);
            if (error_code == 0)
                if (PARA.CREATE_ELEV_PLOT == 1)   % Create elevation plot?
                    sched_state = 8;
                elseif (PARA.CREATE_SKYPLOTS == 1)   % Create Skyplots?
                    sched_state = 9;
                elseif (sched_data.flag_write_output_file == 1)   % Create VEX Output?
                        sched_state = 10;   
                else % sched_data.flag_write_output_file = 0
                        sched_state = 666;
                end
            else
                sched_state = 999;
                error_msg = ['overpass_softcopy: ', error_msg];
            end
            
        case 8 % ++++ Create Overpass Plot  ++++
            [sched_handles, error_code, error_msg] = create_elevation_plot(stat_data, obs_data, PARA, sched_handles);
            if (error_code == 0)
                
                % Save elevation plot to selected sub-folder in /DATA/SCHED/
                pdf_filename = mjd2datestr(PARA.startmjd);
                pdf_filename(strfind(pdf_filename, ':')) = '';
                pdf_filename(strfind(pdf_filename, ' ')) = '_';
                pdf_filename = [pdf_filename, '_elev.pdf'];
                save_pdf_a4_landscape(sched_handles.el_plot.h_fig, ['../DATA/SCHED/', PARA.OUTPUT_SUBDIR, '/'], pdf_filename)
                    
                if (PARA.CREATE_SKYPLOTS == 1)   % Create Skyplots?
                    sched_state = 9;
                elseif (sched_data.flag_write_output_file == 1)   % Create VEX Output?
                    sched_state = 10;
                else % sched_data.flag_write_output_file = 0
                    sched_state = 666;
                end
            else
                sched_state = 999;
                error_msg = ['create_elevation_plot: ', error_msg];
            end
            
        case 9 % ++++ Create Skyplots  ++++
            [error_code, error_msg, sched_handles] = create_skyplot(stat_data, obs_data, PARA, sched_handles);
            if (error_code == 0)
                if (sched_data.flag_write_output_file == 1)   % Create VEX Output?
                    sched_state = 10;
                else % sched_data.flag_write_output_file = 0
                    sched_state = 666;
                end
            else
                sched_state = 999;
                error_msg = ['create_skyplot: ', error_msg];
            end
            
        case 10 % ++++ Command Window In- and Output ++++
            [sched_data, stat_data, source, error_code, error_msg] = scheduler_interface(obs_data, stat_data, source, sched_data, sched_handles, PARA, station);
            if (error_code == 0)
                % Save mat-files:
                save([PARA.pfolder 'sched_data.mat'], 'sched_data');
                sched_state = 777;
            else
                sched_state = 999;
                error_msg = ['scheduler_interface: ', error_msg];
            end
           
            
            
        case 666 % ++++ Do not create an output file (VEX) +++
            fprintf('\nOption: No scheduling task.\n\n');
            break; % while loop
            
        case 777 % ++++ Exit Scheduling Process by CW input => End switch-case routine +++
            fprintf('\nScheduling Process => Cancelled! \n\n');
%             sched_data.flag_write_output_file = 0; % Do not create an output file
            break; % while loop

        case 888 % ++++ End switch-case routine +++
            fprintf('\nScheduling Process => Finished! \n\n');
            break; % while loop
  
        case 999 % ++++ ERROR CASE ++++
            fprintf('=> ERROR: %s.\n\n', error_msg);
            sched_data.flag_write_output_file = 0; % Do not create an output file
            break; % while loop

    end % switch(sched_state)
    
end % while (1)




