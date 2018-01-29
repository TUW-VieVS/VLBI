% #########################################################################
% #     find_closest_quasars_and_calc_scan
% #########################################################################
%
% DESCRIPTION
%   This function finds the closest n quasars (in terms of separation angle)
%   to the antenna positions at "t_epoch_jd". 
%   Calculation of the separation angles are done for the station network defined in "station_id_list_sat" 
%   (find closest quasar to the last satellite which was observed, seen from the satellite station network).
%   The stations of the quasar network ("station_id_list_quasar") are used to check various observation conditions and
%   for the calculation of the scan duration, scan start/end, etc...
%
% CREATED  
%   2015-04-13     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
%   - zazel_s                       : Conversion Ra/Dec => Az/El (Simplified approach for scheduling)
%   - calc_separation_angle
%   - preselect_observable_quasars
%   - calc_scan_duration_network
%   - calc_start_of_scan
%   - check_quasar_visibility_from_network
%   - calc_scan_duration_network
%   - check_t_end
%   - calc_un_az_begin_of_scan
%   - invjday
%
%
% INPUT
%  
%   - stat_data                     : station data structure
%   - station_id_list_quasar        : List of IDs of the quasar observing stations (referring to the required entry of "stat_data.stat" and "obs_data.stat")
%   - station_id_list_sat           : List of IDs of the satellite observing stations (referring to the required entry of "stat_data.stat" and "obs_data.stat")
%   - source                        : Source structure
%   - obs_data                      : Observation data structure (preallocated in "find_observation_times.m")
%   - flag_observable_quasars_list  : Flag list, to define observable quasars at "t_epoch_jd"
%   - number_quasars_output         : Humber of the closest quasars, which should me returned
%   - t_epoch_jd                    : Calculation epoch
%   - flag_print_quasar_infos       : =1 => Print information to Command Window
%   - PARA                          : Global scheduling parameter strucutre
%
%
% OUTPUT
%   - quasars_info_list             : List of closest quasars, including further information
%   - obs_data                      : Observation data structure (preallocated in "find_observation_times.m")

%   - error_code                    : Error Code (0 = no erros occured)
%   - error_msg                     : Error Message (empty, if no errors occured)
%
% CHANGES:
% - 2015-07-01, A. Hellerschmied: satelliet network is used to calc. the separation angles; quasar network is used to calc. observation data (scan start/end, etc..)
%

function [quasars_info_list, obs_data, error_code, error_msg] = find_closest_quasars_and_calc_scan(stat_data, station_id_list_quasar, station_id_list_sat, source, obs_data, flag_observable_quasars_list, number_quasars_output, t_epoch_jd, flag_print_quasar_infos, PARA)

    % ##### Init #####
    error_code = 0;
    error_msg = '';

    count_quasars_output = 0;
    satellite_id = [];
    quasars_info_list = [];

    
    % ##### Pre-select quasars, which are observable from the current network #####
    % For the quasar station network
    if isempty(flag_observable_quasars_list)
        [obs_source, flag_observable_quasars_list, error_code, error_msg] = preselect_observable_quasars(stat_data, station_id_list_quasar, source, t_epoch_jd);
        if error_code > 0
            error_msg = ['preselect_observable_quasars: ', error_msg];
            return; 
        end
    else
        % Preselection => Exclude all sources, which are currently not observable:
        % obs_source = source(logical(flag_observable_quasars_list));
    end

    quasar_id_list = find(flag_observable_quasars_list);

    
    % ##### Init. #####
    number_of_observable_sources = sum(flag_observable_quasars_list);
    number_of_stations_sat = length(station_id_list_sat);
    
    
    % ##### Calculate the separation angle at t_epoch_jd and sort the quasars by this angle  #####
    % For the satellite station network
    % => Preselection of the nearest sources without iterative calculation of the exact slew times/separations angles ...
    
    % Preallocations:
    quasar_scan_info = zeros(number_of_observable_sources, (number_of_stations_sat + 2)); % | quasar_id | sum(sep_angle_i) | sep_angle_1 | sep_angle_2 | ... | sep_angle_n |; i=1:n ; n...number_of_stations_sat
    
    % loop over all observable sources:
    for i_src = 1 : number_of_observable_sources
        
        quasar_id = quasar_id_list(i_src);
        quasar_scan_info(i_src, 1) = quasar_id;

        % loop over all stations:
        for i_stat = 1 : number_of_stations_sat
            
            station_id = station_id_list_sat(i_stat);
            [az_src, el_src, ha_src, dc_src] = zazel_s(t_epoch_jd - 2400000.5, stat_data.stat(station_id).location.ellipsoid.long * pi/180, stat_data.stat(station_id).location.ellipsoid.lat * pi/180, source(quasar_id).ra, source(quasar_id).de);
            [quasar_scan_info(i_src, 2 + i_stat)] = calc_separation_angle(obs_data.stat(station_id).end_of_last_obs.un_az, obs_data.stat(station_id).end_of_last_obs.el, az_src, el_src);
            
        end
        quasar_scan_info(i_src, 2) = sum(quasar_scan_info(i_src, 3 : (number_of_stations_sat + 2)));
        
    end
    
    % Sort sources:
    quasar_scan_info = sortrows(quasar_scan_info,2);
    
    
    
    % ##### Further calculations #####
    
    quasars_info_list = zeros(1, number_of_stations_sat * 6 + 3);

    i_src = 0;
    while (count_quasars_output < number_quasars_output) && (i_src < number_of_observable_sources)
        i_src = i_src + 1;

        quasar_id = quasar_scan_info(i_src, 1);

        % Calc. earliest possible scan start time:
        [t_start_jd_min, obs_data, error_code, error_msg] = calc_start_of_scan(stat_data, station_id_list_quasar, obs_data, PARA, source(quasar_id), satellite_id);
        if error_code > 0
            error_msg = ['calc_start_of_scan: ', error_msg];
            return; 
        end
        [year, mon, day, hr, min, sec] = invjday(t_start_jd_min); % Conversion

        % Check scan start time:
        [flag_t_start_ok, error_code, error_msg] = check_quasar_visibility_from_network(stat_data, station_id_list_quasar, source(quasar_id), t_start_jd_min);
        if error_code > 0
            error_msg = ['check_quasar_visibility_from_network: ', error_msg];
            return; 
        end
        
        % Calc un_az for t_start:
        [obs_data, error_code, error_msg] = calc_un_az_begin_of_scan(stat_data, station_id_list_quasar, obs_data, PARA, source(quasar_id), satellite_id, t_start_jd_min);
        if error_code > 0
            error_msg = ['calc_un_az_begin_of_scan: ', error_msg];
            return; 
        end

        if flag_t_start_ok

            % Calc. scan duration [sec]
            [scan_duration_sec, error_code, error_msg] = calc_scan_duration_network(stat_data, station_id_list_quasar, PARA, source(quasar_id), t_start_jd_min);
            if error_code > 0
                error_msg = ['calc_scan_duration_network: ', error_msg];
                return; 
            end

            % Calc scan end time:     obs_data.stat(station_id).begin_of_new_obs_temp.un_az = [] !!!!!!!!!!!!!!!!!!!!!!!!!!!!
            t_end_jd_min = t_start_jd_min + scan_duration_sec / 86400.0; % Scan end time = scan start + duration


            % Check earliest possible scan end time (and calc. un_az and elevation):
            [flag_t_end_ok, obs_data, error_code, error_msg] = check_t_end(stat_data, station_id_list_quasar, obs_data, PARA, source(quasar_id), satellite_id, t_end_jd_min, t_start_jd_min);
            if error_code > 0
                error_msg = ['check_quasar_visibility_from_network: ', error_msg];
                return; 
            end

        end


        if (flag_t_start_ok) && (flag_t_end_ok)

            count_quasars_output = count_quasars_output + 1;

            % ##### Print informations to CW #####
            if flag_print_quasar_infos
                fprintf(1, '\n');
                fprintf(1, ' %1.0f - Source name:       %s\n', count_quasars_output, source(quasar_id).name);
                fprintf(1, '      Scan duration:     %1.1f sec\n', scan_duration_sec);  
                fprintf(1, '      Earliest possible scan start: %4.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\n', year, mon, day, hr, min, sec); 
                fprintf(1, '      Separation angles at end of last scan:\n');  
                for i_stat = 1 : number_of_stations_sat
                    station_id = station_id_list_sat(i_stat);
                    fprintf(1, '       - %8s:      %1.1f deg\n', stat_data.stat(station_id).name, quasar_scan_info(i_src, 2 + i_stat) * 180 / pi);           
                end
                fprintf(1, '      Separation angles at scan start:\n');  
                for i_stat = 1 : number_of_stations_sat
                    station_id = station_id_list_sat(i_stat);
                    [az_src, el_src, ha_src, dc_src] = zazel_s(t_start_jd_min - 2400000.5, stat_data.stat(station_id).location.ellipsoid.long * pi/180, stat_data.stat(station_id).location.ellipsoid.lat * pi/180, source(quasar_id).ra, source(quasar_id).de);
                    sep_angle_rad = calc_separation_angle(obs_data.stat(station_id).end_of_last_obs.un_az, obs_data.stat(station_id).end_of_last_obs.el, az_src, el_src);
                    fprintf(1, '       - %8s:      %1.1f deg\n', stat_data.stat(station_id).name, sep_angle_rad * 180 / pi);           
                end
            end

            % [obs_data, error_code, error_msg] = calc_un_az_begin_of_scan(stat_data, station_id_list_quasar, obs_data, PARA, source(quasar_id), satellite_id, t_start_jd_min);


            % Assign output variable (| quasar ID | t_start_jd_min | t_end_jd_min | jd_1 |...| jd_n | un_az_1 |...| un_az_n | el_1 |...| el_n | jd_1 |...| jd_n | un_az_1 |...| un_az_n | el_1 |...| el_n |);
            % For quasar station network
            quasars_info_list(count_quasars_output, 1) = quasar_id;
            quasars_info_list(count_quasars_output, 2) = t_start_jd_min;
            quasars_info_list(count_quasars_output, 3) = t_end_jd_min;
            for i_stat = 1 : length(station_id_list_quasar)
                station_id = station_id_list_quasar(i_stat);
                quasars_info_list(count_quasars_output, 3 + i_stat) = obs_data.stat(station_id).begin_of_new_obs_temp.jd;
                quasars_info_list(count_quasars_output, 3 + i_stat + number_of_stations_sat) = obs_data.stat(station_id).begin_of_new_obs_temp.un_az;
                quasars_info_list(count_quasars_output, 3 + i_stat + 2 * number_of_stations_sat) = obs_data.stat(station_id).begin_of_new_obs_temp.el;
                quasars_info_list(count_quasars_output, 3 + i_stat + 3 * number_of_stations_sat) = obs_data.stat(station_id).end_of_new_obs_temp.jd;
                quasars_info_list(count_quasars_output, 3 + i_stat + 4 * number_of_stations_sat) = obs_data.stat(station_id).end_of_new_obs_temp.un_az;
                quasars_info_list(count_quasars_output, 3 + i_stat + 5 * number_of_stations_sat) = obs_data.stat(station_id).end_of_new_obs_temp.el;
            end

        end

    end

    if (count_quasars_output < number_quasars_output)
        fprintf(1, '  \nWARNING: Only %1.0f observable quasars could be found!\n', count_quasars_output);
    end

    % The following fields have to be assigned, after ONE quasar was chosen:
    % For quasar station network
    % => Clear fields here!
    for i_stat = 1 : length(station_id_list_quasar)

        station_id = station_id_list_quasar(i_stat);

        obs_data.stat(station_id).begin_of_new_obs_temp.un_az = [];
        obs_data.stat(station_id).begin_of_new_obs_temp.el = [];
        obs_data.stat(station_id).begin_of_new_obs_temp.jd = [];
        obs_data.stat(station_id).end_of_new_obs_temp.un_az = [];  % [rad]
        obs_data.stat(station_id).end_of_new_obs_temp.el = [];  % [rad]
        obs_data.stat(station_id).end_of_new_obs_temp.jd = [];  % [JD]
    end

return;