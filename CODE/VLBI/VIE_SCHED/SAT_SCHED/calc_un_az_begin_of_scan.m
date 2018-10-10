% #########################################################################
% #     calc_un_az_begin_of_scan
% #########################################################################
%
% DESCRIPTION
%   This function calculates the unambiguous antenna azimuth ("un_az") for an observation 
%   (satellite or quasar) at the defined time "t_epoch_jd" for a defined statio network.
%   "un_az" is stored in the structure "obs_data". 
%   The elevation angle ans eoch [JD] is also saved in "obs_data.stat.end_of_last_obs_temp...".
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
%   2015-04-08     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
% - tle_2_topo
% - calc_slew_time
% - zazel_s
%
% INPUT
% - stat_data           - station data structure
% - station_id_list     - List of IDs of the observing stations (referring to the required entry of "stat_data.stat" and "obs_data.stat")
% - obs_data            - Observation data structure (preallocated in "find_observation_times.m")
% - PARA                - Global scheduling parameter strucutre
% - source_quasar       - source (quasars) data structure (ONE sub-set/field of the "source" structure, containing the required data)
% - satellite_id        - ID of the satellite to be observed (referring to an entry of "stat_data.stat(i_stat).sat(i_sat)" or "obs_data.sat(i_sat)")
% - t_epoch_jd          - Time epoch for the calculations [JD]
%
%
% OUTPUT
% - error_code          - Error Code (0 = no erros occured)
% - error_msg           - Error Message (empty, if no errors occured)
% - obs_data            - Observation data structure (preallocated in "find_observation_times.m"), "obs_data.stat(station_id).begin_of_new_obs_temp.un_az" updated
%
% CHANGES:
%

function [obs_data, error_code, error_msg] = calc_un_az_begin_of_scan(stat_data, station_id_list, obs_data, PARA, source_quasar, satellite_id, t_epoch_jd)
 
    % ##### Init #####
    error_code = 0;
    error_msg = '';

    number_of_stations = length(station_id_list);

  
    % ##########################
    % ##### Satellite scan #####
    % ##########################
    if ~isempty(satellite_id) && isempty(source_quasar)
        
        
        for i_stat = 1 : number_of_stations
            
            % Get station data:
            station_id = station_id_list(i_stat);
            
            % Calculate topocentric vie directions:
            [az, el, range, az_rate, el_rate, r_rate, tra, dc, tra_rate, tdec_rate, ha, local_hour_angle_rate, error_code, error_msg] = tle_2_topo(stat_data, station_id, PARA, satellite_id, t_epoch_jd);
            if error_code > 0
                error_msg = ['tle_2_topo: ', error_msg];
                return; 
            end
            az = az * pi/180;
            el = el * pi/180;
            ha = ha * pi/180;
            dc = dc * pi/180;
            % Calculate the antenna slew time and un_az:
            [slew_time, un_az] = calc_slew_time(stat_data.stat(station_id), obs_data.stat(station_id).end_of_last_obs, az, el, ha, dc);
            
            % Save un_az:
            obs_data.stat(station_id).begin_of_new_obs_temp.un_az = un_az;
            obs_data.stat(station_id).begin_of_new_obs_temp.el = el;
            obs_data.stat(station_id).begin_of_new_obs_temp.jd = t_epoch_jd;
        
        end
        
        
    % ##########################
    % ##### Quasar scan    #####
    % ##########################
    elseif isempty(satellite_id) && ~isempty(source_quasar)
        
        for i_stat = 1 : number_of_stations
            
            % Get station data:
            station_id = station_id_list(i_stat);
            stat_lon = stat_data.stat(station_id).location.ellipsoid.long;
            stat_lat = stat_data.stat(station_id).location.ellipsoid.lat;
            
            % Calculate topocentric vie directions:
            [az, el, ha, dc] = zazel_s(t_epoch_jd - 2400000.5, stat_lon * pi/180, stat_lat * pi/180, source_quasar.ra, source_quasar.de); 
            
            % Calculate the antenna slew time and un_az:
            [slew_time, un_az] = calc_slew_time(stat_data.stat(station_id), obs_data.stat(station_id).end_of_last_obs, az, el, ha, dc);
            
            % Save un_az:
            obs_data.stat(station_id).begin_of_new_obs_temp.un_az = un_az;
            obs_data.stat(station_id).begin_of_new_obs_temp.el = el;
            obs_data.stat(station_id).begin_of_new_obs_temp.jd = t_epoch_jd;
        end
        
    end
    
return;
