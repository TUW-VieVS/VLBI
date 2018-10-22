% #########################################################################
% #     check_quasar_visibility_from_network
% #########################################################################
%
% DESCRIPTION
% This function checks, if the current quasar is visible simultaneously from the current station network.
% Considered conditions:
%  - Axis limits
%  - Cut-off elevation
%  - Horizontal mask (if available)
%  - Sun distance 
%   
%
% CREATED  
%   2015-03-27     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
% - check_quasar_visibility_from_station
%
%
% INPUT
% - stat_data           - station data structure
% - station_id_list     - List of IDs of the observing stations (referring to the required entry of "stat_data.stat" and "obs_data.stat")
% - source_quasar       - source (quasars) data structure (ONE sub-set/field of the "source" structure, containing the required data)
% - t_obs_jd            - Treated time epoch [JD]
%
%
% OUTPUT
% - error_code          - Error Code (0 = no erros occured)
% - error_msg           - Error Message (empty, if no errors occured)
% - flag_observable     - Flag (1 => observable; 0 => Not observable)
%
% CHANGES:
%

function [flag_observable, error_code, error_msg] = check_quasar_visibility_from_network(stat_data, station_id_list, source_quasar, t_obs_jd)
 
    % ##### Init #####
    error_code = 0;
    error_msg = '';
    flag_observable = 0;
    
    
    % ##### Check observability for each station #####
    num_current_stations = length(station_id_list);
    flag_observable_list = zeros(num_current_stations, 1);
    
    for i_stat = 1 : num_current_stations
        [flag_observable_list(i_stat), error_code, error_msg] = check_quasar_visibility_from_station(stat_data.stat(i_stat), source_quasar, t_obs_jd);
        if error_code > 0
            error_msg = ['check_quasar_visibility_from_station: ', error_msg];
            return; 
        end
    end
    
    
    % ##### Check network #####
    if sum(flag_observable_list) == num_current_stations
        flag_observable = 1;
    end

return;
