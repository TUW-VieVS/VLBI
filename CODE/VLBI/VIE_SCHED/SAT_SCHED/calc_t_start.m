% #########################################################################
% #     calc_t_start
% #########################################################################
%
% DESCRIPTION
%   Calculation the earliest possible start time of a scan for a complete station network, 
%   taking into account the previous antenna postion and the antenna slew time.
%   Observations to quasars and to satellites are considered.
%
%   To specify an satellite observation:
%       - satellite_id...is not empty
%       - source_quasar...is empty
%
%   To specify an quasar observation:
%       - satellite_id...is empty
%       - source_quasar...is not empty
%
%
% CREATED  
%   2015-03-30     Andreas Hellerschmied
%
% REFERENCES
%p
%
% COUPLING
% - calc_start_of_scan
% - check_quasar_visibility_from_network
%
%
% INPUT
% - stat_data           - station data structure
% - station_id_list     - List of IDs of the observing stations (referring to the required entry of "stat_data.stat" and "obs_data.stat")
% - obs_data            - Observation data structure (preallocated in "find_observation_times.m")
% - PARA                - Global scheduling parameter strucutre
% - source_quasar       - source (quasars) data structure (ONE sub-set/field of the "source" structure, containing the required data)
% - satellite_id        - ID of the satellite to be observed (referring to an entry of "stat_data.stat(i_stat).sat(i_sat)" or "obs_data.sat(i_sat)")
%
%
% OUTPUT
% - error_code          - Error Code (0 = no erros occured)
% - error_msg           - Error Message (empty, if no errors occured)
% - t_start_jd          - Calculated scan start time for the current station network [JD]
% - flag_t_start_found  - If a valid t_start_jd was found => 1, otherwise 0.
% - obs_data            - Observation data structure (preallocated in "find_observation_times.m")
%
% CHANGES:
%

function [t_start_jd, flag_t_start_found, obs_data, error_code, error_msg] = calc_t_start(stat_data, station_id_list, obs_data, PARA, source_quasar, satellite_id)
 
    % ##### Init #####
    error_code = 0;
    error_msg = '';
    t_start_jd = 0;
    flag_t_start_ok = 0;
    flag_t_start_found = 0;
    

    % ##### Calculate the earliest possible scan start time (taking into account slew times, etc...) #####
    [t_start_jd_temp, obs_data, error_code, error_msg] = calc_start_of_scan(stat_data, station_id_list, obs_data, PARA, source_quasar, satellite_id);
    if error_code > 0
        error_msg = ['calc_start_of_scan: ', error_msg];
        return; 
    end
    
    
    
    % ##########################
    % ##### Satellite scan #####
    % ##########################
    if ~isempty(satellite_id) && isempty(source_quasar)
        
        
        % ##### 1.) Satellite is observable at t_start_jd_temp and at (t_start_jd_temp + PARA.MIN_SAT_SCAN_TIME) #####
        
        % Check, if the satellite is visible at t_start_jd_temp:
        for i_obs = 1 : size(obs_data.sat(satellite_id).obs_times, 1)
            if ((t_start_jd_temp >= obs_data.sat(satellite_id).obs_times(i_obs, 1)) && (t_start_jd_temp <= obs_data.sat(satellite_id).obs_times(i_obs, 2)))
                
                % Check, if the satellite is visible at (t_start_jd_temp + PARA.MIN_SAT_SCAN_TIME):
                for i_obs_2 = 1 : size(obs_data.sat(satellite_id).obs_times, 1)
                    if (((t_start_jd_temp + PARA.MIN_SAT_SCAN_TIME / 86400) >= obs_data.sat(satellite_id).obs_times(i_obs_2, 1)) && ((t_start_jd_temp + PARA.MIN_SAT_SCAN_TIME / 86400) <= obs_data.sat(satellite_id).obs_times(i_obs_2, 2)))

                        % Calc. un_az and el for t_start => Saved in "obs_data":
                        [obs_data, error_code, error_msg] = calc_un_az_begin_of_scan(stat_data, station_id_list, obs_data, PARA, source_quasar, satellite_id, t_start_jd_temp);
                        if error_code > 0
                            error_msg = ['calc_un_az_begin_of_scan: ', error_msg];
                            return;
                        else
                            % Assign output variables and exit the function:
                            flag_t_start_found = 1;
                            t_start_jd = t_start_jd_temp;
                            return;
                        end
                    end
                end
                
            end
        end
        
        
        % ##### 2.) Satellite is NOT observable at t_start_jd_temp and at (t_start_jd_temp + PARA.MIN_SAT_SCAN_TIME) #####
        
        if flag_t_start_found == 0
            
            % Get the beginning of the next available observation periods (store it in a list):
            scan_start_jd_list = [];
            i_scan = 0;
            for i_obs = 1 : size(obs_data.sat(satellite_id).obs_times, 1)
                if (    (t_start_jd_temp <= obs_data.sat(satellite_id).obs_times(i_obs, 1)) && ...
                        (((t_start_jd_temp + PARA.MIN_SAT_SCAN_TIME / 86400) <= obs_data.sat(satellite_id).obs_times(i_obs, 1)) && ((obs_data.sat(satellite_id).obs_times(i_obs, 1) + PARA.MIN_SAT_SCAN_TIME / 86400) <= obs_data.sat(satellite_id).obs_times(i_obs, 2))) ...  
                    ) 
                    i_scan = i_scan + 1;
                    scan_start_jd_list(i_scan) = obs_data.sat(satellite_id).obs_times(i_obs, 1);
                end
            end
            
            % Find the earliest start time in the list:
            if ~isempty(scan_start_jd_list)
                
                % Find t_start by sorting:
                scan_start_jd_list = sort(scan_start_jd_list);
                
                % Calc. un_az and el for t_start => Saved in "obs_data":
                [obs_data, error_code, error_msg] = calc_un_az_begin_of_scan(stat_data, station_id_list, obs_data, PARA, source_quasar, satellite_id, scan_start_jd_list(1));
                if error_code > 0
                    error_msg = ['calc_un_az_begin_of_scan: ', error_msg];
                    return; 
                end
                
                % Assign output variables and exit the function:
                t_start_jd = scan_start_jd_list(1);
                flag_t_start_found = 1;
                return;
                
            else % There is no subsequent observation period available:
                flag_t_start_found = 0;
                return;
            end

        end
        
        
        
    % ##########################
    % ##### Quasar scan    #####
    % ##########################
    elseif isempty(satellite_id) && ~isempty(source_quasar)
        
        % ##### 1.) Calculate the scan duration #####
        % Only a more or less rough estimation for the epoch defined in "t_start_jd_temp".
        
        % [scan_duration_sec, error_code, error_msg] = calc_scan_duration_network(stat_data, station_id_list, PARA, source_quasar, obsmode, t_start_jd_temp);
        [scan_duration_sec, error_code, error_msg] = calc_scan_duration_network(stat_data, station_id_list, PARA, source_quasar, t_start_jd_temp);
        if error_code > 0
            error_msg = ['calc_scan_duration_network: ', error_msg];
            return; 
        end
        
        
        % ###### 2.) Check, if (t_start_jd_temp + scan_duration_sec / 86400) is before session end #####
        if ((t_start_jd_temp + scan_duration_sec / 86400) >= (PARA.endmjd + 2400000.5)) 
            fprintf(1, '    => Session start time + min. scan duration exceeds the session end!');
            flag_t_start_found = 0;
            return;
        end
        
        
        % ##### 3.) Check, if the source is observable from all stations #####
        
        % Check is the source is observable at (t_start_jd_temp):
        [flag_t_start_ok, error_code, error_msg] = check_quasar_visibility_from_network(stat_data, station_id_list, source_quasar, t_start_jd_temp);
        if error_code > 0
            error_msg = ['check_quasar_visibility_from_network: ', error_msg];
            return; 
        end
        
        % Check if the source is observable at (t_start_jd_temp + scan_duration_sec / 86400):
        if flag_t_start_ok
            [flag_t_start_ok, error_code, error_msg] = check_quasar_visibility_from_network(stat_data, station_id_list, source_quasar, (t_start_jd_temp + scan_duration_sec / 86400));
            if error_code > 0
                error_msg = ['check_quasar_visibility_from_network: ', error_msg];
                return; 
            end
        end
        
        if flag_t_start_ok  % Source is observable at (t_start_jd_temp) and (t_start_jd_temp + scan_duration_sec / 86400) => Finished!
            
            % Calc. un_az and el for t_start => Saved in "obs_data":
            [obs_data, error_code, error_msg] = calc_un_az_begin_of_scan(stat_data, station_id_list, obs_data, PARA, source_quasar, satellite_id, t_start_jd_temp);
            if error_code > 0
                error_msg = ['calc_un_az_begin_of_scan: ', error_msg];
                return; 
            end
            
            % Assign output variables and exit the function:
            t_start_jd = t_start_jd_temp;
            flag_t_start_found = 1;
            
            
        else    % Source is NOT observable at (t_start_jd_temp) and (t_start_jd_temp + scan_duration_sec / 86400) => Start searching the next time epoch where the source is observable again:

            % Loop init.:
            delta_t = PARA.INIT_PROP_INTERVAL; % [sec]
            t_start_jd_temp_old = 0;
            flag_t_stop_ok = 0;
            flag_increase_t_start = 1;
            
            while(delta_t > 0.1) % While step-size delta_t is larger than 0.1 sec...
                
                % Get the new scan start and stop time:
                t_start_jd_temp_old = t_start_jd_temp;
                if flag_increase_t_start
                    t_start_jd_temp = t_start_jd_temp + delta_t / 86400; % [jd]
                    t_stop_jd_temp = t_start_jd_temp + scan_duration_sec / 86400; % [jd]
                end
                flag_increase_t_start = 1;
                
                % Check, if the source is observable at (t_start_jd_temp):
                [flag_t_start_ok, error_code, error_msg] = check_quasar_visibility_from_network(stat_data, station_id_list, source_quasar, t_start_jd_temp);
                if error_code > 0
                    error_msg = ['check_quasar_visibility_from_network: ', error_msg];
                    return; 
                end
                
                % Check, if the source is observable at (t_start_jd_temp + scan_duration_sec / 86400 = t_stop_jd_temp):
                if flag_t_start_ok
                    [flag_t_stop_ok, error_code, error_msg] = check_quasar_visibility_from_network(stat_data, station_id_list, source_quasar, t_stop_jd_temp);
                    if error_code > 0
                        error_msg = ['check_quasar_visibility_from_network: ', error_msg];
                        return; 
                    end
                    
                    % Check, if t_stop exceeds session-end time:
                    if flag_t_stop_ok
                        if t_stop_jd_temp > (PARA.endmjd + 2400000.5)
                            flag_t_stop_ok = 0;
                        end
                    end
                    
                end
                
                % If, the current t_start id OK => start iteration to determine t_start more exactly => decrease step-size delta_t by factor 2
                if flag_t_start_ok
                    
                    % initialize subsequent loop iteration:
                    delta_t = delta_t / 2;
                    t_start_jd_temp = t_start_jd_temp_old;
                    flag_increase_t_start = 0;
                    
                    if flag_t_stop_ok
                        
                        % Assign output variables:
                        flag_t_start_found = 1;
                        t_start_jd = t_start_jd_temp;
                    end
                    
                end
                
            end % while
            
        end
        
        % Calc. un_az and el for t_start => Saved in "obs_data":
        [obs_data, error_code, error_msg] = calc_un_az_begin_of_scan(stat_data, station_id_list, obs_data, PARA, source_quasar, satellite_id, t_start_jd);
        if error_code > 0
            flag_t_start_found = 0;
            t_start_jd = 0;
            error_msg = ['calc_un_az_begin_of_scan: ', error_msg];
            return; 
        end
        
        

    end % if ~isempty(satellite_id) && isempty(source_quasar)
    
return;
