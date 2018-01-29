% #########################################################################
% #     calc_quasar_scan
% #########################################################################
%
% DESCRIPTION
%   This function calculates for a quasar scans and the defined station network:
%     - Earliest possible scan start time [JD]
%     - Scan duration [sec]
%     - Earliest possible scan end time [JD]
%     - un_az [rad] for calculated scan start/end
%     - el [rad] for calculated scan start/end
%
%   If the input parameter "duration_sec" is set, the scan duration is not calculated automatically and the inputted value is taken.
%
% CREATED  
%   2015-07-01     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
%   - calc_scan_duration_network
%   - calc_start_of_scan
%   - calc_un_az_begin_of_scan
%   - check_quasar_visibility_from_network
%   - calc_scan_duration_network
%   - check_t_end
%   - invjday
%
%
% INPUT
%   - stat_data                     : station data structure
%   - station_id_list               : List of IDs of the quasar observing stations (referring to the required entry of "stat_data.stat" and "obs_data.stat")
%   - source                        : Source structure
%   - obs_data                      : Observation data structure (preallocated in "find_observation_times.m")
%   - quasar_id                     : Flag list, to define observable quasars at "t_epoch_jd"
%   - flag_print_quasar_infos       : =1 => Print information to Command Window
%   - PARA                          : Global scheduling parameter strucutre
%   - duration_sec                  : Scan duration (optional), if [] => calc. scan duration
%   - t_start_jd_in                 : Scan start time [JD] (optional)
%
%
% OUTPUT
%   - quasars_info_list             : List of quasars (only one in case of this function...), including further sation dependent information
%   - obs_data                      : Observation data structure (preallocated in "find_observation_times.m")
%   - error_code                    : Error Code (0 = no erros occured)
%   - error_msg                     : Error Message (empty, if no errors occured)
%
% CHANGES:
% - A. Hellerschmied: Added possibility to defined the scan duration manually via input parameter "duration_sec" (optional input parameter)
% - 2016-11-29, A. Hellerschmied: Optional input argument to define scan start time ("t_start_jd") manually
%

function [quasars_info_list, obs_data, error_code, error_msg] = calc_quasar_scan(stat_data, station_id_list, source, obs_data, quasar_id, flag_print_quasar_infos, PARA, varargin)

    % ##### Init #####
    error_code = 0;
    error_msg = '';
    satellite_id = [];
    number_of_stations = length(station_id_list);
    quasars_info_list = zeros(1, number_of_stations * 6 + 3);
    flag_calc_scan_duration = 0;
    flag_calc_scan_start_time = true;
    
    
     % ##### Check input #####
    switch(nargin)
            
        case 7 % Calculate scan duration automatically
            flag_calc_scan_duration = 1;
            
        case 8 % Enter Scan duration manually (input parameter "duration_sec" is set)
            if isempty(varargin{1})
                flag_calc_scan_duration = 1;
            else
                flag_calc_scan_duration = 0;
                scan_duration_sec = varargin{1};
                % Check input:
                if (scan_duration_sec < 0)
                    error_code = 1;
                    error_msg = 'The scan duration has to be a positive value!';
                    return;
                end
            end
            
        case 9 % Enter Scan duration manually (input parameter "duration_sec" is set)
            if isempty(varargin{1})
                flag_calc_scan_duration = 1;
            else
                flag_calc_scan_duration = 0;
                scan_duration_sec = varargin{1};
                % Check input:
                if (scan_duration_sec < 0)
                    error_code = 1;
                    error_msg = 'The scan duration has to be a positive value!';
                    return;
                end
            end
            flag_calc_scan_start_time = false;
            t_start_jd_in = varargin{2};
                
        otherwise % Error
            error_code = 1;
            error_msg = 'Incorrect input!';
            return;
            
    end % switch(nargin)
    

    % ##### Calculations #####

    % Calc. earliest possible scan start time:
    [t_start_jd_min, obs_data, error_code, error_msg] = calc_start_of_scan(stat_data, station_id_list, obs_data, PARA, source(quasar_id), satellite_id);
    if error_code > 0
        error_msg = ['calc_start_of_scan: ', error_msg];
        return; 
    end
    if flag_calc_scan_start_time
        t_start_jd = t_start_jd_min;
    else 
        % ### Check, if entered t_start_jd_in is larger than the earliest possible scan start time ###
        if t_start_jd_in >= t_start_jd_min % OK
            t_start_jd = t_start_jd_in;
        else % error
            quasars_info_list = [];
            [year, mon, day, hr, min, sec] = invjday(t_start_jd_in); % Conversion
            fprintf(1,'Quasar %s is not observable at scan start time: %4.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\n', source(quasar_id).name , year, mon, day, hr, min, sec);
            [year, mon, day, hr, min, sec] = invjday(t_start_jd_min); % Conversion
            fprintf(1,' => Earlierst possible scan start time: %4.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\n', year, mon, day, hr, min, sec);
            return; 
        end
    end
    
    % Check scan start time:
    [flag_t_start_ok, error_code, error_msg] = check_quasar_visibility_from_network(stat_data, station_id_list, source(quasar_id), t_start_jd);
    if error_code > 0
        error_msg = ['check_quasar_visibility_from_network: ', error_msg];
        return; 
    end
    % If quasar is not observable => return empty "quasars_info_list"
    if ~flag_t_start_ok
        quasars_info_list = [];
        [year, mon, day, hr, min, sec] = invjday(t_start_jd); % Conversion
        fprintf(1,'Quasar %s is not observable at the calculated scan start time: %4.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\n', source(quasar_id).name , year, mon, day, hr, min, sec);
        return;
    end

    % Calc un_az for t_start:
    [obs_data, error_code, error_msg] = calc_un_az_begin_of_scan(stat_data, station_id_list, obs_data, PARA, source(quasar_id), satellite_id, t_start_jd);
    if error_code > 0
        error_msg = ['calc_un_az_begin_of_scan: ', error_msg];
        return; 
    end

    % Calc. scan duration [sec] (if it is not defined maually via input parameter)
    if flag_calc_scan_duration
        [scan_duration_sec, error_code, error_msg] = calc_scan_duration_network(stat_data, station_id_list, PARA, source(quasar_id), t_start_jd);
        if error_code > 0
            error_msg = ['calc_scan_duration_network: ', error_msg];
            return; 
        end
    end

    % Calc scan end time:     obs_data.stat(station_id).begin_of_new_obs_temp.un_az = [] !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    t_end_jd = t_start_jd + scan_duration_sec / 86400.0; % Scan end time = scan start + duration

    % Check scan end time (and calc. un_az and elevation):
    [flag_t_end_ok, obs_data, error_code, error_msg] = check_t_end(stat_data, station_id_list, obs_data, PARA, source(quasar_id), satellite_id, t_end_jd, t_start_jd);
    if error_code > 0
        error_msg = ['check_quasar_visibility_from_network: ', error_msg];
        return; 
    end
    % If quasar is not observable => return empty "quasars_info_list"
    if ~flag_t_end_ok
        quasars_info_list = [];
        [year, mon, day, hr, min, sec] = invjday(t_end_jd); % Conversion
        fprintf(1,'Quasar %s is not observable at the calculated scan end time: %4.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\n', source(quasar_id).name , year, mon, day, hr, min, sec);
        return;
    end


    % ##### Print informations to CW #####
    if flag_print_quasar_infos
        [year, mon, day, hr, min, sec] = invjday(t_start_jd); % Conversion
        fprintf(1, '\n');
        fprintf(1, '      Source name:                  %s\n', source(quasar_id).name);
        fprintf(1, '      Scan duration:                %1.1f sec\n', scan_duration_sec);  
        fprintf(1, '      Earliest possible scan start: %4.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\n', year, mon, day, hr, min, sec); 
        [year, mon, day, hr, min, sec] = invjday(t_end_jd); % Conversion
        fprintf(1, '      Resulting scan end:           %4.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\n', year, mon, day, hr, min, sec);
    end


    % ##### Assign output variable #####
    % (| quasar ID | t_start_jd_min | t_end_jd | jd_1 |...| jd_n | un_az_1 |...| un_az_n | el_1 |...| el_n | jd_1 |...| jd_n | un_az_1 |...| un_az_n | el_1 |...| el_n |) 
    quasars_info_list(1, 1) = quasar_id;
    quasars_info_list(1, 2) = t_start_jd;
    quasars_info_list(1, 3) = t_end_jd;
    for i_stat = 1 : number_of_stations
        station_id = station_id_list(i_stat);
        quasars_info_list(1, 3 + i_stat)                            = obs_data.stat(station_id).begin_of_new_obs_temp.jd;   % [JD]
        quasars_info_list(1, 3 + i_stat + number_of_stations)       = obs_data.stat(station_id).begin_of_new_obs_temp.un_az;% [rad]
        quasars_info_list(1, 3 + i_stat + 2 * number_of_stations)   = obs_data.stat(station_id).begin_of_new_obs_temp.el;   % [rad]
        quasars_info_list(1, 3 + i_stat + 3 * number_of_stations)   = obs_data.stat(station_id).end_of_new_obs_temp.jd;     % [JD]
        quasars_info_list(1, 3 + i_stat + 4 * number_of_stations)   = obs_data.stat(station_id).end_of_new_obs_temp.un_az;  % [rad]
        quasars_info_list(1, 3 + i_stat + 5 * number_of_stations)   = obs_data.stat(station_id).end_of_new_obs_temp.el;     % [rad]
    end


    % The following fields have to be assigned, after ONE quasar was chosen:
    % => Clear fields here!
    for i_stat = 1 : number_of_stations
        station_id = station_id_list(i_stat);
        obs_data.stat(station_id).begin_of_new_obs_temp.un_az   = [];  % [rad]
        obs_data.stat(station_id).begin_of_new_obs_temp.el      = [];  % [rad]
        obs_data.stat(station_id).begin_of_new_obs_temp.jd      = [];  % [JD]
        obs_data.stat(station_id).end_of_new_obs_temp.un_az     = [];  % [rad]
        obs_data.stat(station_id).end_of_new_obs_temp.el        = [];  % [rad]
        obs_data.stat(station_id).end_of_new_obs_temp.jd        = [];  % [JD]
    end

return;