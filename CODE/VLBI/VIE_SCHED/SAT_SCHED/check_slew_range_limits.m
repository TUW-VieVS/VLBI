% #########################################################################
% #     check_slew_range_limits
% #########################################################################
%
% DESCRIPTION
%   This function checks the antenna axis limits (slew range limits).
%
%   Supported antenna mout types:
%    - Az/el
%    - Ha/Dec
%    - XY/EW
%
%
% CREATED  
%   2015-08-18     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
% - azel2xyew
%
%
% INPUT
% - stat_data           - station data structure
% - station_id          - ID of the observing stations (referring to a entry in "stat_data.stat")
% - az_rad              - Azimuth angle [rad]
% - el_rad              - Elevation angle [rad]
% - ha_rad              - Local hour angle [rad]
% - dec_rad             - Declination [rad]
%
%
% OUTPUT
% - error_code              - Error Code (0 = no erros occured)
% - error_msg               - Error Message (empty, if no errors occured)
% - flag_axis_limits_ok     - Flag: 1 => axis limts are kept; 0 => axis limts are NOT kept
%
% CHANGES:
% - 2015-08-18, Andreas Hellerschmied: Bug-fix: Check slew limits for HADC antennas.
%


function [flag_axis_limits_ok, error_code, error_msg] = check_slew_range_limits(stat_data, station_id, az_rad, el_rad, ha_rad, dec_rad)
 
    % ##### Init #####
    error_code = 0;
    error_msg = '';
    flag_axis_limits_ok = 0;
    
    
    % ##### Check limits: #####
    if (strncmp(stat_data.stat(station_id).axis_type, 'AZEL', 4))
        
        while (az_rad < stat_data.stat(station_id).lim11)
            az_rad = az_rad + 2 * pi;
        end
        while (az_rad > stat_data.stat(station_id).lim12)
            az_rad = az_rad - 2 * pi;
        end
        
        if ((az_rad > stat_data.stat(station_id).lim11) & (az_rad < stat_data.stat(station_id).lim12) & (el_rad > stat_data.stat(station_id).lim21) & (el_rad < stat_data.stat(station_id).lim22))
            flag_axis_limits_ok = true;
        else
            flag_axis_limits_ok = false;
        end
        
    elseif (strncmp(stat_data.stat(station_id).axis_type, 'HADC', 4))
        
        while (ha_rad < stat_data.stat(station_id).lim11)
            ha_rad = ha_rad + 2 * pi;
        end
        while (ha_rad > stat_data.stat(station_id).lim12)
            ha_rad = ha_rad - 2 * pi;
        end
        
        if ((ha_rad > stat_data.stat(station_id).lim11) & (ha_rad < stat_data.stat(station_id).lim12) & (dec_rad > stat_data.stat(station_id).lim21) & (dec_rad < stat_data.stat(station_id).lim22))
            flag_axis_limits_ok = true;
        else 
            flag_axis_limits_ok = false;
% TEST
% fprintf(1, '    - LHA=%5.3f DEC=%5.3f\n', ha_rad*180/pi, dec_rad*180/pi);
        end
        
    elseif (strncmp(stat_data.stat(station_id).axis_type, 'XYEW', 4))
        [x30, y30, error_code, error_msg] = azel2xyew(az_rad, el_rad);
        if error_code ~= 0
            error_msg = ['azel2xyew:', error_msg];
            return;
        end
        if ((x30 > stat_data.stat(station_id).lim11) & (x30 < stat_data.stat(station_id).lim12) & (y30 > stat_data.stat(station_id).lim21) & (y30 < stat_data.stat(station_id).lim22))
            flag_axis_limits_ok = true;
        else 
            flag_axis_limits_ok = false;
        end
    end   

return;
