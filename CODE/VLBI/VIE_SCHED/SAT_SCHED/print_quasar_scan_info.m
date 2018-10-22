% #########################################################################
% #     print_quasar_scan_info
% #########################################################################
%
% DESCRIPTION
%   This function finds the closest n quasars (in terms of separation angle)
%   to the antenna positions at "t_epoch_jd". 
%   Calculations are done for the station network defined in "station_id_list".
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
%
%
% INPUT
%  
%   - stat_data                     : station data structure
%   - station_id_list               : List of IDs of the observing stations (referring to the required entry of "stat_data.stat" and "obs_data.stat")
%   - source                        : Source structure
%   - obs_data                      : Observation data structure (preallocated in "find_observation_times.m")
%   - flag_observable_quasars_list  : Flag list, to define observable quasars at "t_epoch_jd"
%   - number_quasars_output         : Humber of the closest quasars, which should me returned
%   - t_epoch_jd                    : Calculation epoch
%
%
% OUTPUT
%   - closest_quasars_list          : List of closest quasars (quasar IDs)
%   - error_code                    : Error Code (0 = no erros occured)
%   - error_msg                     : Error Message (empty, if no errors occured)
%
% CHANGES:
%

function [error_code, error_msg] = print_quasar_scan_info(stat_data, station_id_list, source, obs_data, flag_observable_quasars_list, number_quasars_output, t_epoch_jd)

    % ##### Init #####
    error_code = 0;
    error_msg = '';
    
    
    
return;