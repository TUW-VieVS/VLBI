% #########################################################################
% #     check_quasar_visibility_from_station
% #########################################################################
%
% DESCRIPTION
% This function checks, if the current quasar is visible from the given station.
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
% - zazel_s
% - ssunpo
% - azel2xyew
%
%
% INPUT
% - stat                - station data structure (for one station: stat = stat_data.stat(i_stat))
% - source_quasar       - source (quasars) data structure (ONE sub-set/field of the "source" structure, containing the required data)
% - t_obs_jd            - Treated time [JD]
%
%
% OUTPUT
% - error_code          - Error Code (0 = no erros occured)
% - error_msg           - Error Message (empty, if no errors occured)
% - flag_observable     - Flag (1 => observable; 0 => Not observable)
%
% CHANGES:
% - 2015-07-01, A. Hellerschmied: azel2xyew.m used to calculate x and y angles (for XYEW mount)
%

function [flag_observable, error_code, error_msg] = check_quasar_visibility_from_station(stat, source_quasar, t_obs_jd)
 
    % ##### Init #####
    error_code = 0;
    error_msg = '';
    flag_observable = true;
    flag_sun_dist = false;
    

    % ##### Calculate Az/El, Ha/Dec for the current station/source constellation #####
    
    % Station position (conversion deg => rad):
    lon = stat.location.ellipsoid.long * pi / 180;
    lat = stat.location.ellipsoid.lat * pi / 180;
    [az, el, ha, dc] = zazel_s(t_obs_jd - 2400000.5, lon, lat, source_quasar.ra, source_quasar.de);
    
    
    % ##### Check, if the source is observable #####
    % Modified from function "zlup" from Jing Sun 2010 (VIE_SCHED)

    % axis limit
    if (strncmp(stat.axis_type, 'AZEL', 4))
        % Az-check is only userful, if the full azimuth range of the antenna is less than 360 deg!?!?!
        unaz = az;
        while (unaz < stat.lim11)
            unaz = unaz + 2 * pi;
        end
        while (unaz > stat.lim12)
            unaz = unaz - 2 * pi;
        end
        if ((unaz > stat.lim11) & (unaz < stat.lim12))
            lupaz = true;
        else
            lupaz = false;
        end
        if ((el > stat.lim21) & (el < stat.lim22))
            lupel = true;
        else 
            lupel = false;
        end
        lup = lupaz & lupel;
        
    elseif (strncmp(stat.axis_type, 'HADC', 4))
        if ((ha > stat.lim11) & (ha < stat.lim12) & (dc > stat.lim21) & (dc < stat.lim22))
            lup = true;
        else 
            lup = false;
        end
        
    elseif (strncmp(stat.axis_type, 'XYEW', 4))
        [x30, y30, error_code, error_msg] = azel2xyew(az, el, 0);
        if error_code > 0
            error_msg = ['azel2xyew:', error_msg];
            return;
        end
%         cel = cos(el);
%         sel = sin(el);
%         caz = cos(az);
%         saz = sin(az);
%         x30 = -atan2(cel*caz,sel);
%         y30 = asin(cel*saz);
        if ((x30 > stat.lim11) & (x30 < stat.lim12) & (y30 > stat.lim21) & (y30 < stat.lim22))
            lup = true;
        else 
            lup = false;
        end
    end      

    % ##### cut-off elevation #####
    if (el > (stat.min_elevation * pi / 180)) 
        lupc = true;
    else
        lupc = false;
    end
    lup = lup & lupc; 

    % ##### hotizontal mask #####
    hmasknum = stat.horizontal_mask_num;
    hmask    = stat.horizontal_mask;
    if ((lup == true) & (hmasknum > 0))
        if (mod(hmasknum,2) ~= 0)
            hmasktype = 1;   % step functions
        else
            hmasktype = 2;   % line segments
        end
        for i = 1 : floor(hmasknum/2)   
            ib = i * 2;
            if ((az >= hmask(ib-1)) & (az <= hmask(ib+1)))
                if (hmasktype == 1)
                    elmask = hmask(ib);
                elseif (hmasktype == 2)
                    elmask = ((hmask(ib+2) - hmask(ib)) / (hmask(ib+1) - hmask(ib-1))) * (az - hmask(ib-1)) + hmask(ib);
                end
                break;
            end
        end
        lupm = (el > elmask);
        lup = lup & lupm;
    end
    
    
    % ##### Check Sun distance #####
    min_sun_dist = stat.min_sun_dist * pi/180;
    
    % Calculate sun position [rad]:
    [sunra, sunde] = ssunpo(t_obs_jd - 2400000.5);
    
    % Calculate topocentric view angles [rad]:
    [az_sun, el_sun, ha_sun, dc_sun] = zazel_s(t_obs_jd - 2400000.5, lon, lat, sunra, sunde);
    
    % Calculate the separation angle [rad]:
    [sun_dist] = calc_separation_angle(az, el, az_sun, el_sun); 

    % Check the condition:
    if sun_dist > min_sun_dist
        flag_sun_dist = true;
    end

    
    % ##### Assign output variable #####
    flag_observable = lup & flag_sun_dist;


return;
