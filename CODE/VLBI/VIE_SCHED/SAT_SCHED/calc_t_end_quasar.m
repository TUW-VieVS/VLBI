% #########################################################################
% #     calc_t_end_quasar
% #########################################################################
%
% DESCRIPTION
%   This function 
%
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
% - stat_data           - station data structure
% - station_id_list     - List of IDs of the observing stations (referring to the required entry of "stat_data.stat" and "obs_data.stat")
% - obs_data            - Observation data structure (preallocated in "find_observation_times.m")
% - PARA                - Global scheduling parameter strucutre
% - source_quasar       - source (quasars) data structure (ONE sub-set/field of the "source" structure, containing the required data)
% - t_end_in            - scan end time to be checked for various conditions [JD]
% - t_start_in          - scan start time to be checked for various conditions [JD]
%
%
% OUTPUT
% - error_code          - Error Code (0 = no erros occured)
% - error_msg           - Error Message (empty, if no errors occured)
% - flag_t_start_ok     - Flag: 1 => t_start is OK; 0 => t_start is not OK
% - obs_data            - Observation data structure (preallocated in "find_observation_times.m"), "obs_data.stat(station_id).end_of_new_obs_temp.un_az" updated
%
% CHANGES:
%

function [t_end_jd, obs_data, error_code, error_msg] = calc_t_end_quasar(stat_data, station_id_list, obs_data, PARA, source_quasar, t_end_in, t_start_in)
 
    % ##### Init #####
    error_code = 0;
    error_msg = '';
    t_end_jd = 0;
    
    
    
    % scantmp.sta(ista).endmjd = scantmp.startmjd + scantmp.sta(ista).duration / 86400.0; % Scan end time = scan start + duration
   
 


    
    
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
