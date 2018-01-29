% #########################################################################
% #     check_t_end
% #########################################################################
%
% DESCRIPTION
%   This function checks, if the satellite/quasar is observable at the time t_end_in.
%   In the case of quasars, the function also checks the observability in a defined 
%   interval (PARA.PARA.CHECK_SCAN_INT) during the scan (from t_start_in to t_end_in).
%   Furthermore, the antenna axes limits are checks in the defined interval (PARA.CHECK_SCAN_INT)
%   and the unambiguous antenna azimuth is calculated for the scan end time.
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
%   2015-04-01     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
% - check_quasar_visibility_from_network
% - check_axis_limits
%
%
% INPUT
% - stat_data                       - station data structure
% - station_id_list                 - List of IDs of the observing stations (referring to the required entry of "stat_data.stat" and "obs_data.stat")
% - obs_data                        - Observation data structure (preallocated in "find_observation_times.m")
% - PARA                            - Global scheduling parameter strucutre
% - source_quasar                   - source (quasars) data structure (ONE sub-set/field of the "source" structure, containing the required data)
% - satellite_id                    - ID of the satellite to be observed (referring to an entry of "stat_data.stat(i_stat).sat(i_sat)" or "obs_data.sat(i_sat)")
% - t_end_in                        - scan end time to be checked for various conditions [JD]
% - t_start_in                      - scan start time to be checked for various conditions [JD]
% - flag_statwise_sat_visibility    - Flag. If =1: satellite visibility is only be checked for the stations in the "station_id_list" (optional)
%
%
% OUTPUT
% - error_code          - Error Code (0 = no erros occured)
% - error_msg           - Error Message (empty, if no errors occured)
% - flag_t_start_ok     - Flag: 1 => t_start is OK; 0 => t_start is not OK
% - obs_data            - Observation data structure (preallocated in "find_observation_times.m"), "obs_data.stat(station_id).end_of_new_obs_temp.un_az" updated
%
% CHANGES:
% - 2015-07-01, A. Hellerschmied: Use interval defined by PARA.CHECK_SCAN_INT to check various scan conditions
% - 2015-11-16, A. Hellerschmied: Possibility to check satellite visibility only for the stations in "station_id_list".
% - 2016-04-21, A. Hellerschmied: Bug-fix (for last update from 2016-04-16)
% - 2016-05-11, A. Hellerschmied: Bug-fix: Check of scan end time, if a satellite is visible within more than one time windows.
% - 2016-06-29, A. Hellerschmied: - Option for debug output added (PARA.DEBUG_FLAG) 
%                                 - Issue with calculating the correct un_az angle was fixed:
%                                    => Calc un_az step by step in the PARA.CHECK_SCAN_INT interval
%                                    => Take the un_az value of the previous step as start value for the next step, a.s.o..
% - 2016-07-13, A. Hellerschmied: Msg. to CW in case an axis limit is violated.
% - 2016-11-02, A. Hellerschmied: Fixed a problem at checking the axis limits
% - 2016-12-22, A. Hellerschmied: Header updated
%

function [flag_t_end_ok, obs_data, error_code, error_msg] = check_t_end(stat_data, station_id_list, obs_data, PARA, source_quasar, satellite_id, t_end_in, t_start_in, varargin)
 
    % ##### Init #####
    error_code = 0;
    error_msg = '';
    flag_t_end_ok = 0;
    flag_observable = 0;
    
    scan_duration_sec = (t_end_in - t_start_in) * 86400;
    number_of_epochs = floor(scan_duration_sec / PARA.CHECK_SCAN_INT);
    number_of_stat = length(station_id_list);
    
    % (Check, if t_end_in is within an already scheduled scan)
    % !!! => This is only important for INSERTING a scan, not for APPENDING!
    
    
    % ##### Select, if satellite visibility should only be checked for the stations in the "station_id_list" #####
    % ...Required for sub-netting
    if nargin == 9
        flag_statwise_sat_visibility = varargin{1};
    else
        flag_statwise_sat_visibility = 0;
    end
    
    
    
    % ##########################
    % ##### Satellite scan #####
    % ##########################
    if ~isempty(satellite_id) && isempty(source_quasar)
        
        % Check, if the satellite is visible between t_sart_in and t_end_in:
        if ~flag_statwise_sat_visibility
            for i_obs = 1 : size(obs_data.sat(satellite_id).obs_times, 1)
                if ((t_end_in >= obs_data.sat(satellite_id).obs_times(i_obs, 1)) && (t_start_in <= obs_data.sat(satellite_id).obs_times(i_obs, 2)))
                    flag_observable = 1;
                    break;
                end
            end
        else
            flag_observable = 1;
            % ### Loop over all stations in "station_id_list" ###
            for station_id = station_id_list
            
                % Check, if the satellite is visible permanently between t_sart_in and t_end_in:
                flag_observable_from_stat = false;
                if ~isempty(obs_data.sat(satellite_id).stat(station_id).obs_times)
                    for i_obs = 1 : size(obs_data.sat(satellite_id).stat(station_id).obs_times, 1)
                        if ((t_start_in >= obs_data.sat(satellite_id).stat(station_id).obs_times(i_obs, 1)) && (t_start_in <= obs_data.sat(satellite_id).stat(station_id).obs_times(i_obs, 2)) && (t_end_in >= obs_data.sat(satellite_id).stat(station_id).obs_times(i_obs, 1)) && (t_end_in <= obs_data.sat(satellite_id).stat(station_id).obs_times(i_obs, 2)))
                            flag_observable_from_stat = true;
                            break;
                        end
                    end
                end
                if ~flag_observable_from_stat
                    flag_observable = 0;
                end
                
                if ~flag_observable
                    break;
                end
            end
            
        end
        
        % Assign output variables:
        if flag_observable
            flag_t_end_ok = 1;
        else
            return; % Source is not observable...
        end
        
        
        
    % ##########################
    % ##### Quasar scan    #####
    % ##########################
    elseif isempty(satellite_id) && ~isempty(source_quasar)

        % Check, if the source is up and observable from the network during the defined scan duration (from t_start_in to t_end_in; visibility checked in the interval PARA.CHECK_SCAN_INT) and at t_end_in:
        
        t_epoch = t_start_in;
        
        % Check observability during scan: 
        for i_epoch = 1 : number_of_epochs
            [flag_t_end_ok, error_code, error_msg] = check_quasar_visibility_from_network(stat_data, station_id_list, source_quasar, t_epoch);
            if error_code > 0
                error_msg = ['check_quasar_visibility_from_network: ', error_msg];
                return; 
            end
            if ~flag_t_end_ok % Source is not observable...
               return;
            end
            
            t_epoch = t_epoch + PARA.CHECK_SCAN_INT / 86400;
        end
        
        % Check at t_end_in:
        [flag_t_end_ok, error_code, error_msg] = check_quasar_visibility_from_network(stat_data, station_id_list, source_quasar, t_end_in);
        if error_code > 0
            error_msg = ['check_quasar_visibility_from_network: ', error_msg];
            return; 
        end

    end
    
          
            
    % ##### Check antenna axes limts during the scan and calculate un_az #####
    flag_limits_ok_list = zeros(number_of_stat, 1);
    un_az_list = zeros(number_of_stat, 1);
    el_list = zeros(number_of_stat, 1);
    
   
    for i_stat = 1 : number_of_stat
        
        % Loop init.:
        station_id          = station_id_list(i_stat);
        t_epoch             = t_start_in  + PARA.CHECK_SCAN_INT / 86400;
        t_epoch_last_step   = t_start_in;
        un_az_last_step     = obs_data.stat(station_id).begin_of_new_obs_temp.un_az;
        
        if PARA.DEBUG_FLAG
            disp(['-------------- ',num2str(station_id), ' -------------------']);
        end

        for i_epoch = 1 : number_of_epochs
            % Check timeperiods from "t_epoch_last_step" to "t_epoch" step by step (interval = PARA.CHECK_SCAN_INT)
            
            if PARA.DEBUG_FLAG
                disp(['*********** ',num2str(i_epoch), ' ***********']);
            end
            
            [flag_limits_ok_list(i_stat), un_az_last_step, el_list(i_stat), error_code, error_msg] = check_axis_limits(stat_data, station_id, PARA, source_quasar, satellite_id, t_epoch, un_az_last_step, t_epoch_last_step);
            if error_code > 0
                error_msg = ['check_axis_limits: ', error_msg];
                return; 
            end
            if ~flag_limits_ok_list(i_stat) % Axis limit is violated...
               flag_t_end_ok = 0;
               fprintf(' => Axis limit is violated at station: %s; un_az = %f deg; epoch: %s; lim11: %3.1f; lim12: %3.1f\n', stat_data.stat(i_stat).name, un_az_list(i_stat)*180/pi, jd2datestr(t_epoch), stat_data.stat(i_stat).lim11*180/pi, stat_data.stat(i_stat).lim12*180/pi)
               return;
            end
            t_epoch_last_step = t_epoch; 
            t_epoch = t_epoch + PARA.CHECK_SCAN_INT / 86400;
        end
        
        if PARA.DEBUG_FLAG
            disp('*********** t_end_jd ***********');
        end
            
        % Check at t_end_in and calc un_az:
        [flag_limits_ok_list(i_stat), un_az_list(i_stat), el_list(i_stat), error_code, error_msg] = check_axis_limits(stat_data, station_id, PARA, source_quasar, satellite_id, t_end_in, un_az_last_step, t_epoch_last_step);
        if error_code > 0
            error_msg = ['check_axis_limits: ', error_msg];
            return; 
        end
        if ~flag_limits_ok_list(i_stat) % Axis limit is violated...
           flag_t_end_ok = 0;
           fprintf(' => Axis limit is violated at station: %s; un_az = %f deg; epoch: %s; lim11: %3.1f; lim12: %3.1f\n', stat_data.stat(i_stat).name, un_az_list(i_stat)*180/pi, jd2datestr(t_end_in), stat_data.stat(i_stat).lim11*180/pi, stat_data.stat(i_stat).lim12*180/pi)
           return;
        end
        
    end
    
    flag_t_end_ok = (sum(flag_limits_ok_list) == number_of_stat) && flag_t_end_ok;
    
    
    % Save un_az at "t_end_in" for all stations (only, if "t_end_in" is OK!):
    if flag_t_end_ok
        for i_stat = 1 : number_of_stat
            station_id = station_id_list(i_stat);
            obs_data.stat(station_id).end_of_new_obs_temp.un_az = un_az_list(i_stat);  % [rad]
            obs_data.stat(station_id).end_of_new_obs_temp.el = el_list(i_stat);  % [rad]
            obs_data.stat(station_id).end_of_new_obs_temp.jd = t_end_in;  % [JD]
        end
    end

    
return;
