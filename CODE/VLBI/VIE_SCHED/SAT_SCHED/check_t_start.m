% #########################################################################
% #     check_t_start
% #########################################################################
%
% DESCRIPTION
%   Calculation the earliest possible start time of a scan for a complete station network, 
%   taking into account the previous antenna postion and the antenna slew time. The 
%   calculation is done iteratively bat the function "calc_start_of_obs".
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
%   2015-03-27     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
% - calc_start_of_scan
% - check_quasar_visibility_from_network
% - calc_un_az_begin_of_scan
%
%
% INPUT
% - stat_data           - station data structure
% - station_id_list     - List of IDs of the observing stations (referring to the required entry of "stat_data.stat" and "obs_data.stat")
% - obs_data            - Observation data structure (preallocated in "find_observation_times.m")
% - PARA                - Global scheduling parameter strucutre
% - source_quasar       - source (quasars) data structure (ONE sub-set/field of the "source" structure, containing the required data)
% - satellite_id        - ID of the satellite to be observed (referring to an entry of "stat_data.stat(i_stat).sat(i_sat)" or "obs_data.sat(i_sat)")
% - t_start_in          - scan start time to be checked for various conditions [JD]
% - flag_statwise_sat_visibility    - Flag. If =1: satellite visibility is only be checked for the stations in the "station_id_list" (optional)
%
%
% OUTPUT
% - error_code          - Error Code (0 = no erros occured)
% - error_msg           - Error Message (empty, if no errors occured)
% - start_jd_min        - Earliest possible start time of the scan for the current station network [JD]
% - flag_t_start_ok     - Flag: 1 => t_start is OK; 0 => t_start is not OK
% - obs_data            - Observation data structure (preallocated in "find_observation_times.m")
%
% CHANGES:
% - 2016-04-21, A. Hellerschmied: Sat. visibility check is based on station ID list (optional).
% - 2016-07-13: A. Hellerschmied: Return the earliest possible scan start time, even if the start time is not valid
% - 2016-11-12: A. Hellerschmied: Bug-fix: Wrong results, if  "flag_statwise_sat_visibility" and a station had more than 1 obs. periods in "obs_data.sat(satellite_id).stat(station_id).obs_times"
%


function [flag_t_start_ok, start_jd_min, obs_data, error_code, error_msg] = check_t_start(stat_data, station_id_list, obs_data, PARA, source_quasar, satellite_id, t_start_in, varargin)
 
    % ##### Init #####
    error_code = 0;
    error_msg = '';
    flag_t_start_ok = 0;
    start_jd_min = 0;
    flag_observable = 0;
    
    % ##### Select, if satellite visibility should only be checked for the stations in the "station_id_list" #####
    % ...Required for sub-netting
    if nargin == 8
        flag_statwise_sat_visibility = varargin{1};
    else
        flag_statwise_sat_visibility = 0;
    end
    
    

    % ##### Calculate the earliest possible scan start time (taking into account slew times, etc...) #####
    
    [start_jd_min_temp, obs_data, error_code, error_msg] = calc_start_of_scan(stat_data, station_id_list, obs_data, PARA, source_quasar, satellite_id); % Calc. un_az
    if error_code > 0
        error_msg = ['calc_start_of_scan: ', error_msg];
        return; 
    end
    
%     start_jd_min_temp = t_start_in;
    
    % ##### Check, if entered t_start is larger than the earliest possible scan start time #####
    if start_jd_min_temp > t_start_in
        start_jd_min = start_jd_min_temp;
        return; 
    end
    
    
    % ##########################
    % ##### Satellite scan #####
    % ##########################
    if ~isempty(satellite_id) && isempty(source_quasar)
        
        % Check, if the satellite is visible at t_start:
        if flag_statwise_sat_visibility
            
            % Init. flag list:
            flag_observable_stat_list = false(length(station_id_list), 1);
        
            % From station network defined in "station_id_list"

            % Loop over all sations in "station_id_list"
            for i_stat = 1 : length(station_id_list)
                
                % Get station ID:
                station_id = station_id_list(i_stat);

                if isempty(obs_data.sat(satellite_id).stat(station_id).obs_times)
                    flag_observable_stat_list(i_stat) = 0;
                end
                
                for i_obs = 1 : size(obs_data.sat(satellite_id).stat(station_id).obs_times, 1)
                    if (t_start_in >= obs_data.sat(satellite_id).stat(station_id).obs_times(i_obs, 1)) && (t_start_in <= obs_data.sat(satellite_id).stat(station_id).obs_times(i_obs, 2))
                        flag_observable_stat_list(i_stat) = true;
                        break;
                    end
                end
            end
            if sum(flag_observable_stat_list) == length(station_id_list)
                flag_observable = 1;
            else
                for i_stat = 1 : length(flag_observable_stat_list)
                    if ~flag_observable_stat_list(i_stat)
                        stat_id = station_id_list(i_stat);
                        fprintf(' => Satellite not observable from station: %s\n', stat_data.stat(stat_id).name);
                    end
                end
            end
        else
            
            % From the complete sat. network:
            for i_obs = 1 : size(obs_data.sat(satellite_id).obs_times, 1)
                if ((t_start_in >= obs_data.sat(satellite_id).obs_times(i_obs, 1)) && (t_start_in <= obs_data.sat(satellite_id).obs_times(i_obs, 2)))
                    flag_observable = 1;
                    break;
                end
            end
            
        end % if flag_statwise_sat_visibility
        
        
        
        if flag_observable
            
            % Calc un_az for t_start:
            [obs_data, error_code, error_msg] = calc_un_az_begin_of_scan(stat_data, station_id_list, obs_data, PARA, source_quasar, satellite_id, t_start_in); %start_jd_min_temp);
            if error_code > 0
                error_msg = ['calc_un_az_begin_of_scan: ', error_msg];
                return; 
            end
            
            % Assign output variables:
            flag_t_start_ok = 1;
            start_jd_min = start_jd_min_temp;
        end
        
        
        % (Check, if t_start_in is within an already scheduled scan)
        % !!! => This is only important for INSERTING a scan, not for APPENDING!
        
        
        
    % ##########################
    % ##### Quasar scan    #####
    % ##########################
    elseif isempty(satellite_id) && ~isempty(source_quasar)

        % Check, if the source is up and observable from the network:
        [flag_t_start_ok, error_code, error_msg] = check_quasar_visibility_from_network(stat_data, station_id_list, source_quasar, t_start_in);
        if error_code > 0
            error_msg = ['check_quasar_visibility_from_network: ', error_msg];
            return; 
        end
        
        if flag_t_start_ok
            
            % Calc un_az for t_start:
            [obs_data, error_code, error_msg] = calc_un_az_begin_of_scan(stat_data, station_id_list, obs_data, PARA, source_quasar, satellite_id, t_start_in); %start_jd_min_temp);
            if error_code > 0
                error_msg = ['calc_un_az_begin_of_scan: ', error_msg];
                return; 
            end
            
            % Assign output variables:
            start_jd_min = start_jd_min_temp;
        end
    end
    
return;
