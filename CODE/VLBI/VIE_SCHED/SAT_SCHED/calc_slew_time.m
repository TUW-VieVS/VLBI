% #########################################################################
% #     calc_slew_time
% #########################################################################
%
% DESCRIPTION
%   Calculation of the antenna slew time between two consecutive
%   observations to quasars/satellites for one particular station. 
%   Adapded from: sslew.m from Jing Sun (2010)
%
% CREATED  
%   2015-03-06     Andreas Hellerschmied
%
% REFERENCES
% - Jing Sun (2013), VLBI Scheduling strategies with respect to VLBI2010,
%   Geowissenschaftliche Mitteilungen, Heft Nr. 92, ISSN 1811-8380.
%
%
% COUPLING
% - sdfaz   - Calculate the azimuth difference between the NOW and the NEW source 
%             positions, taking into account cable wrap.
% - ssva2t  - Calculate the slew time considering acceleration and deceleration [sec]
%
%
% INPUT
% - stat                - station data structure (for one station: stat = stat_data.stat(i_stat))
% - end_of_last_obs     - Structure, cobntaining the antenna position (az, el, ha, dc, un_az) at the end of the last observation (jd) [rad]
% - az                  - Azimuth of the new target [rad]
% - el                  - Elevartion of the new target [rad]
% - ha                  - Hour angle of the new target [rad]
% - dc                  - Declination of the new target [rad]
%
%
% OUTPUT
% - slew_time           - calculated antenna slew time [sec]
% - un_az               - unambiguous NOW source azimuth [rad]
%
% CHANGES:
% - 2015-06-29: A. Hellerschmied: azel2xyew.m used to calc. X/Y EW angles
% - 2016-11-02, A. Hellerschmied: Added possibity to select an alternative cable wrap section ("flag_take_alternative_cable_wrap")

%
%

function [slew_time, un_az] = calc_slew_time(stat, end_of_last_obs, az, el, ha, dc) 

    if (strcmp(stat.axis_type(1:4), 'AZEL'))
        % az slewing time
        if isfield(stat, 'flag_alt_cable_wrap')
            flag_take_alternative_cable_wrap = stat.flag_alt_cable_wrap;
        else
            flag_take_alternative_cable_wrap = false;
        end
        
        [un_az, dfaz] = sdfaz(stat.lim11, stat.lim12, end_of_last_obs.un_az, az, flag_take_alternative_cable_wrap);  % [unaznew, dfaz] = sdfaz(lim11, lim12, unaznow, aznew)
        [staz] = ssva2t(dfaz*180/pi, stat.max_axis1_rate, stat.max_axis1_acc);
        staz = staz + stat.c1;   % [sec]    
        % el slewing time
        dfel = abs(end_of_last_obs.el - el) * 180 / pi;
        [stel] = ssva2t(dfel, stat.max_axis2_rate, stat.max_axis2_acc);
        stel = stel +  stat.c2;   % [sec]
        % slewing time
        slew_time = ceil(max(staz, stel));

    elseif (strcmp(stat.axis_type(1:4), 'HADC'))
        % ha slewing time
        hanow  = end_of_last_obs.ha;
        hanew  = ha;
        dfha = abs(hanew - hanow);
        [stha] = ssva2t(dfha*180/pi, stat.max_axis1_rate, stat.max_axis1_acc);  % Calculate the slew time considering acceleration and deceleration [sec]. 
        stha = stha + stat.c1;   % [sec] 
        % dc slewing time
        dcnow  = end_of_last_obs.dc;
        dcnew  = dc;
        dfdc = abs(dcnew - dcnow);
        [stdc] = ssva2t(dfdc*180/pi, stat.max_axis2_rate, stat.max_axis2_acc);
        stdc = stdc + stat.c2;   % [sec] 
        % slewing time
        slew_time = ceil(max(stha, stdc));
        un_az = az; 

    elseif (strcmp(stat.axis_type(1:4), 'XYEW'))
        
        % xnow, ynow
        [xnow, ynow] = azel2xyew(end_of_last_obs.un_az, end_of_last_obs.el, 0);
%         cel = cos(end_of_last_obs.el);
%         sel = sin(end_of_last_obs.el);
%         caz = cos(end_of_last_obs.un_az);
%         saz = sin(end_of_last_obs.un_az);
%         xnow = -atan2(cel*caz,sel);
%         ynow = asin(cel*saz);
        
        % xnew, ynew
        [xnew, ynew] = azel2xyew(az, el, 0);
%         cel = cos(el);
%         sel = sin(el);
%         caz = cos(az);
%         saz = sin(az);
%         xnew = -atan2(cel*caz,sel);
%         ynew = asin(cel*saz);
        
        % x mount slewing time
        dfx = abs(xnew - xnow);
        [stx] = ssva2t(dfx*180/pi, stat.max_axis1_rate, stat.max_axis1_acc);
        stx = stx + stat.c1;   % [sec] 
        % y mount slewing time
        dfy = abs(ynew - ynow);
        [sty] = ssva2t(dfy*180/pi, stat.max_axis2_rate, stat.max_axis2_acc);
        sty = sty + stat.c2;   % [sec] 
        % slewing time
        slew_time = ceil(max(stx, sty));
        un_az = az; 
    end

return
