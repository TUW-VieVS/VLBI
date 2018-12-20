% -------------------------------------------------------------------------
%
%                              function find_orbit_events.m
%
%   This function calculates the exact point in time of the following orbit
%   events:
%           - Rise : Above cut-off elevation or horizon mask
%           - Set  : Beneath cut-off elevation or horizon mask
%           - Peak : Max. elevation (local maxima)
%           - Exceedance of sun distance limit
%           - Keep sun distance limit.
%           - Exceedance of axes 1 (Az, Ha, XY) rate limit
%           - Keep axes 1 (Az, Ha, XY) rate limit.
%           - Exceedance of axes 2 (El, Dec, EW) rate limit
%           - Keep axes 2 (El, Dec, EW) rate limit.
%           - Exceedance of antenna slew range limits
%           - Keep antenna slew range limits
%
%   The derived data is added to the "stat_data" structure.
%   Calculations are done by an iterative approach.
%
%   Author: 
%       Andreas Hellerschmied (heller182@gmx.at), 27.10.2013
%   
%   changes       :
%   - 2014-01.19: A. Hellerschmied: Support of different antenna axis types (AzEl, HaDec)
%   - 2014-01.19: A. Hellerschmied: Error treatment added.
%   - 2015-02-16: A. Hellerschmied: Remanded (before: "TLE_find_orbit_events.m")
%   - 2015-07-20: A. Hellerschmied: New TLE/SGP4 orbit propagation function used (tle_propagation.m), new input structure: PARA
%   - 2015-07-29: A. Hellerschmied: Horizintal mask support added (for the determination of rise/set).
%   - 2015-07-31: A. Hellerschmied: Support of XYEW and HADC mounted antennas for slew rate checks added.
%   - 2015-08-18: A. Hellerschmied: Check for antenna slew range limits added as an additional observation condition.
%   - 2016-04-11: A. Hellerschmied: Bug-fix: Typo.
%   - 2016-04-13: A. Hellerschmied: Bug-fix: Numerical problems when calculatione elevation with eci2topo.m => El is rounded to the number of digits defined in "el_digits".
%   - 2016-04-25: A. Hellerschmied: Bug-fixes: - Typo.
%                                              - Absolute value of slew rates are now checked.
%           
%
%   inputs        :
%   - stat_data   : station data structure
%   - PARA        : Scheduling parameter structure
%     
%
%   outputs       :
%   - stat_data: station data structure, containing additional
%                    overpass information.
%   - error_code    : Error Code (0 = no erros occured)
%   - error_msg     : Error Message (empty, if no errors occured)
%    
%
%   locals        :
% 
%
%   coupling      :
%       - tle_propagation.m : Orbit propagation
%       - eci2topo.m        : Calc. of topocentric Coordinates and their rates.
%       - invjday.m         : Cal. of Year, Day, Hour, Minutes, Sec. of JD.
%       - calc_sundist.m    : Angular distance between Source and Sun.
%       - azel2xyew         : Conversion of AZEL pointing data to XYEW
%   
%
%   references    :
%
%-------------------------------------------------------------------------


function [stat_data, error_code, error_msg] = find_orbit_events(stat_data, PARA)

    % Init:
    error_code = 0;
    error_msg = '';
    
    
    % Elevation [deg] from eci2topo.m is rounded for the determination of rise and set. 
    % Define the significat digits here:
    el_digits = 1000;


    for i_stat = 1 : stat_data.number_of_stations    % Stations
        
        stat_lon = stat_data.stat(i_stat).location.ellipsoid.long;
        stat_lat = stat_data.stat(i_stat).location.ellipsoid.lat;
        stat_alt = stat_data.stat(i_stat).location.ellipsoid.altitude;
        
        
        for i_sat = 1 : stat_data.number_of_sat   % Satellites
            
            % Loop init.:
            i_overpass = 0;
            stat_data.stat(i_stat).sat(i_sat).number_of_overpasses = 0;
            i_peak = 0;
            flag_above_min_elevation = 0;
            stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_exceed_min_sundist = 0;
            stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_keep_min_sundist = 0;
            
            stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_keep_max_axis1_rate = 0;
            stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_exceed_max_axis1_rate = 0;
            
            stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_keep_max_axis2_rate = 0;
            stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_exceed_max_axis2_rate = 0;
            
            stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_keep_slew_range_limits = 0;
            stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_exceed_slew_range_limits = 0;
 
 
            % First element of satellite ephemeris data is already above
            % min. elevation? 
            if (stat_data.stat(i_stat).sat(i_sat).epoch(1).above_min_elevation == 1)
                flag_above_min_elevation = 1;
                i_overpass = 1;
                i_peak = 0;
                
                stat_data.stat(i_stat).sat(i_sat).number_of_overpasses = i_overpass;
                
                % preallocation
                stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass) = struct('rise',...
                    [], 'set', [], 'peak', [], 'number_of_peaks', []);
                stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise = struct('year', [], 'mon', [], 'day', [],...
                    'h', [], 'min', [], 'sec', [], 'jd', [], 'az', [], 'el', [], 'range', []);
                stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set = struct('year', [], 'mon', [], 'day', [],...
                    'h', [], 'min', [], 'sec', [], 'jd', [], 'az', [], 'el', [], 'range', []);
                stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak = struct('year', [], 'mon', [], 'day', [],...
                    'h', [], 'min', [], 'sec', [], 'jd', [], 'az', [], 'el', [], 'range', []);

                
                stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).number_of_peaks = 0;
                
            end;
            
            for i_epoch = 2 : stat_data.number_of_epochs  % Propagation epochs
                
                    
                if ( (i_epoch >= 3) && (flag_above_min_elevation) )

                    % ----- Overpass -----
                    
                    % Find peak
                    if (  ( (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).el)     < ...
                            (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).el)       ) && ...
                          ( (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch-1).el)     > ...
                            (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 2).el)       )    )
                        
                        i_peak = i_peak + 1;

                        % Start parameter for the iterative rise search:
                        t_start             = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 2).jd;  % !
                        t_stop              = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).jd;      % !
                        delta_t             = stat_data.stat(i_stat).prop_setup.delta_t_min;            % !
                        sat_name_to_prop = stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_header_line;

                        el_0 = 0;
                        el_1 = 0;
                        el_2 = 0;
                        az_0 = 0;
                        az_1 = 0;

                        delta_t = delta_t / 20;

                        % Iterative rise search approach:
                        while (delta_t >= (0.01 / 60))   % /60 = conversion from sec to min!

                            [temp_sat_data, error_code, error_msg] = tle_propagation(t_start, t_stop, delta_t, PARA, sat_name_to_prop);
                            if error_code > 0
                                error_msg = ['tle_propagation:', error_msg];
                                return; 
                            end

                            number_temp_epochs = length(temp_sat_data.sat(1).prop);

                            for i_temp_epoch = 1 : number_temp_epochs

                                el_2 = el_1;
                                el_1 = el_0;
                                az_1 = az_0;

                                [az_0, el_0, range_s_obs, az_rate, el_rate, r_rate] = eci2topo(temp_sat_data.sat(1).prop(i_temp_epoch).jd, temp_sat_data.sat(1).prop(i_temp_epoch).r, temp_sat_data.sat(1).prop(i_temp_epoch).v, stat_lon, stat_lat, stat_alt);
                                if (i_temp_epoch >= 3)
                                    if ( (el_0 < el_1)  && (el_1 > el_2) )
                                        delta_t = delta_t / 20;
                                        t_start = temp_sat_data.sat(1).prop(i_temp_epoch - 2).jd;
                                        t_stop  = temp_sat_data.sat(1).prop(i_temp_epoch).jd;
                                        break;
                                    end;
                                end;

                            end; % i_temp_epoch = 1 : number_temp_epochs

                        end; % while (delta_t >= (0.1 / 60))

                        t_peak = t_start + (t_stop - t_start) / 2;

                        [year, mon, day, hr, min, sec] = invjday(t_peak);

                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak(i_peak).year    = year;                            
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak(i_peak).mon     = mon;
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak(i_peak).day     = day;
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak(i_peak).h       = hr;
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak(i_peak).min     = min;
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak(i_peak).sec     = sec;
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak(i_peak).jd      = t_peak;
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak(i_peak).az      = az_1;
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak(i_peak).el      = el_1;
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak(i_peak).range   = range_s_obs;

                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).number_of_peaks = i_peak;
                        
                    end; % if (...... "find peak")

                end; % if ( (i_epoch >= 3) && (flag_above_min_elevation) )


                
                
                % Find rise
                if ( (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).above_min_elevation == 1)     && ...
                     (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).above_min_elevation == 0)         ) 
                 
                    flag_above_min_elevation = 1;
                    i_overpass = i_overpass + 1;
                    stat_data.stat(i_stat).sat(i_sat).number_of_overpasses = i_overpass;
                    i_peak = 0;
                    
                    
                    % preallocation
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass) = struct('rise',...
                        [], 'set', [], 'peak', [], 'number_of_peaks', []);
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise = struct('year', [], 'mon', [], 'day', [],...
                        'h', [], 'min', [], 'sec', [], 'jd', [], 'az', [], 'el', [], 'range', []);
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set = struct('year', [], 'mon', [], 'day', [],...
                        'h', [], 'min', [], 'sec', [], 'jd', [], 'az', [], 'el', [], 'range', []);
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak = struct('year', [], 'mon', [], 'day', [],...
                        'h', [], 'min', [], 'sec', [], 'jd', [], 'az', [], 'el', [], 'range', []);
                    
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).number_of_peaks = 0;

                    % Start parameter for the iterative rise search:
                    t_start             = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).jd;  % !
                    t_stop              = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).jd;      % !
                    delta_t             = stat_data.stat(i_stat).prop_setup.delta_t_min;            % !
                    sat_name_to_prop = stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_header_line;

                    % Iterative rise search approach:
                    while (delta_t >= (0.01 / 60))   % /60 = conversion from sec to min!

                        delta_t = delta_t / 10;

                        [temp_sat_data, error_code, error_msg] = tle_propagation(t_start, t_stop, delta_t, PARA, sat_name_to_prop);
                        if error_code > 0
                            error_msg = ['tle_propagation:', error_msg];
                            return; 
                        end
                        
                         number_temp_epochs = length(temp_sat_data.sat(1).prop);

                        for i_temp_epoch = 1 : number_temp_epochs

                            [az, el, range_s_obs, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate] = eci2topo(temp_sat_data.sat(1).prop(i_temp_epoch).jd, temp_sat_data.sat(1).prop(i_temp_epoch).r, temp_sat_data.sat(1).prop(i_temp_epoch).v, stat_lon, stat_lat, stat_alt);
                            el = round(el * el_digits) / el_digits;
                            [flag_above_min_el] = check_h_mask(stat_data, i_stat, az*pi/180, el*pi/180);
                            if (el >  stat_data.stat(i_stat).min_elevation) && (flag_above_min_el)
                                t_start = temp_sat_data.sat(1).prop(i_temp_epoch - 1).jd;
                                t_stop  = temp_sat_data.sat(1).prop(i_temp_epoch).jd;
                                break;
                            end;

                        end; % i_temp_epoch = 1 : number_temp_epochs

                    end; % while (delta_t >= (0.1 / 60))

                    t_rise = t_stop;

                    [year, mon, day, hr, min, sec] = invjday(t_rise);

                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise.year    = year;                            
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise.mon     = mon;
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise.day     = day;
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise.h       = hr;
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise.min     = min;
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise.sec     = sec;
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise.jd      = t_rise;
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise.az      = az;
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise.el      = el;
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise.range   = range_s_obs;

                end; % if (...... "find rise"
                    
                
                
                % Find Set
                if ( (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).above_min_elevation == 0)     && ...
                     (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).above_min_elevation == 1)         ) 
                 
                    flag_above_min_elevation = 0;
                    
                    % Start parameter for the iterative rise search:
                    t_start             = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).jd;  % !
                    t_stop              = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).jd;      % !
                    delta_t             = stat_data.stat(i_stat).prop_setup.delta_t_min;            % !
                    sat_name_to_prop = stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_header_line;
                    
                    az_0 = 0;
                    el_0 = 0;
                    az_1 = 0;
                    el_1 = 0;
                    
                    % Iterative set search approach:
                    while (delta_t >= (0.01 / 60))   % /60 = conversion from sec to min!
                        
                        delta_t = delta_t / 10;

                        [temp_sat_data, error_code, error_msg] = tle_propagation(t_start, t_stop, delta_t, PARA, sat_name_to_prop);
                        if error_code > 0
                            error_msg = ['tle_propagation:', error_msg];
                            return; 
                        end
                        
                        number_temp_epochs = length(temp_sat_data.sat(1).prop);
                        
                        for i_temp_epoch = 1 : number_temp_epochs
                            
                            az_1 = az_0;
                            el_1 = el_0;

                            [az_0, el_0, range_s_obs, az_rate, el_rate, r_rate] = eci2topo(temp_sat_data.sat(1).prop(i_temp_epoch).jd, temp_sat_data.sat(1).prop(i_temp_epoch).r, temp_sat_data.sat(1).prop(i_temp_epoch).v, stat_lon, stat_lat, stat_alt);
                            el_0 = round(el_0 * el_digits) / el_digits;
                            [flag_above_min_el] = check_h_mask(stat_data, i_stat, az_0*pi/180, el_0*pi/180);
                            
                            if (el_0 <  stat_data.stat(i_stat).min_elevation) || (flag_above_min_el == 0)
                                t_start = temp_sat_data.sat(1).prop(i_temp_epoch - 1).jd;
                                t_stop  = temp_sat_data.sat(1).prop(i_temp_epoch).jd;
                                break;
                            end;

                        end; % i_temp_epoch = 1 : number_temp_epochs
                        
                    end; % while (delta_t >= (0.1 / 60))
                    
                    t_set = t_start;
                    
                    [year, mon, day, hr, min, sec] = invjday(t_set);

                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set.year    = year;                            
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set.mon     = mon;
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set.day     = day;
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set.h       = hr;
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set.min     = min;
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set.sec     = sec;
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set.jd      = t_set;
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set.az      = az_1;
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set.el      = el_1;
                    stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set.range   = range_s_obs;
               
                end; % if (...... "find set"
                
                
                
                %% ----- Sun distance -----
                
                % Find "keep min. sun distance":
                if ( (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_min_sun_dist == 0)     && ...
                     (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).exceed_min_sun_dist == 1)         ) 
                 
                    
                    % Start parameter for the iterative rise search:
                    sun_dist_residual_limit = 0.001; %[deg]
                    sun_dist = 180;
                    t_start             = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).jd;  % !
                    t_stop              = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).jd;      % !
                    delta_t             = stat_data.stat(i_stat).prop_setup.delta_t_min;            % !
                    sat_name_to_prop = stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_header_line;
                    
                    
                    % Iterative rise search approach:
                    while ( (delta_t >= (0.0001 / 60)) && (abs(sun_dist - stat_data.stat(i_stat).min_sun_dist) > sun_dist_residual_limit) )   % /60 = conversion from sec to min!
                        
                        delta_t = delta_t / 10;

                        [temp_sat_data, error_code, error_msg] = tle_propagation(t_start, t_stop, delta_t, PARA, sat_name_to_prop);
                        if error_code > 0
                            error_msg = ['tle_propagation:', error_msg];
                            return; 
                        end
                        
                        number_temp_epochs = length(temp_sat_data.sat(1).prop);
                        
                        for i_temp_epoch = 1 : number_temp_epochs

                            [az, el, range_s_obs, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate] = eci2topo(temp_sat_data.sat(1).prop(i_temp_epoch).jd, temp_sat_data.sat(1).prop(i_temp_epoch).r, temp_sat_data.sat(1).prop(i_temp_epoch).v, stat_lon, stat_lat, stat_alt);
                            [sun_dist] = calc_sundist(temp_sat_data.sat(1).prop(i_temp_epoch).jd, stat_lon, stat_lat, az, el);
                            
                            % Break Condition
                            if (sun_dist > stat_data.stat(i_stat).min_sun_dist) % el_0 <  stat_data.stat(i_stat).min_elevation
                                if (i_temp_epoch == 1)
                                    t_start = temp_sat_data.sat(1).prop(i_temp_epoch).jd - delta_t / (24*60);
                                    t_stop  = temp_sat_data.sat(1).prop(i_temp_epoch).jd;
                                else
                                    t_start = temp_sat_data.sat(1).prop(i_temp_epoch - 1).jd;
                                    t_stop  = temp_sat_data.sat(1).prop(i_temp_epoch).jd;
                                end
                                break;
                            end;

                        end; % i_temp_epoch = 1 : number_temp_epochs
                        
                    end; % while (delta_t >= (0.1 / 60))
                    
                    t_keep_min_sun_dist = t_stop;
                    [temp_sat_data, error_code, error_msg] = tle_propagation(t_keep_min_sun_dist, t_keep_min_sun_dist, delta_t, PARA, sat_name_to_prop);
                    if error_code > 0
                        error_msg = ['tle_propagation:', error_msg];
                        return; 
                    end
                    
                    [az, el, range_s_obs, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate] = eci2topo(temp_sat_data.sat(1).prop(1).jd, temp_sat_data.sat(1).prop(1).r, temp_sat_data.sat(1).prop(1).v, stat_lon, stat_lat, stat_alt);
                    [sun_dist] = calc_sundist(temp_sat_data.sat(1).prop(1).jd, stat_lon, stat_lat, az, el);
                            
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_keep_min_sundist = stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_keep_min_sundist + 1;
                    
                    [year, mon, day, hr, min, sec] = invjday(t_keep_min_sun_dist);

                    stat_data.stat(i_stat).sat(i_sat).sun_dist.keep_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_keep_min_sundist).year     = year;                            
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.keep_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_keep_min_sundist).mon      = mon;
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.keep_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_keep_min_sundist).day      = day;
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.keep_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_keep_min_sundist).h        = hr;
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.keep_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_keep_min_sundist).min      = min;
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.keep_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_keep_min_sundist).sec      = sec;
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.keep_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_keep_min_sundist).jd       = t_keep_min_sun_dist;
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.keep_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_keep_min_sundist).az       = az;
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.keep_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_keep_min_sundist).el       = el;
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.keep_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_keep_min_sundist).range    = range_s_obs;
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.keep_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_keep_min_sundist).sun_dist = sun_dist;

               
                end; % if (...... Find "keep min. sun distance"
                
                

                % Find "exceed min. sun distance":
                if ( (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_min_sun_dist == 1)     && ...
                     (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).exceed_min_sun_dist == 0)         ) 
                 
                    
                    % Start parameter for the iterative rise search:
                    sun_dist_residual_limit = 0.001; %[deg]
                    sun_dist = 180;
                    t_start             = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).jd;  % !
                    t_stop              = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).jd;      % !
                    delta_t             = stat_data.stat(i_stat).prop_setup.delta_t_min;            % !
                    sat_name_to_prop = stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_header_line;
                    
                    
                    % Iterative rise search approach:
                    while ( (delta_t >= (0.0001 / 60)) && (abs(sun_dist - stat_data.stat(i_stat).min_sun_dist) > sun_dist_residual_limit) )   % /60 = conversion from sec to min!
                        
                        delta_t = delta_t / 10;
                        [temp_sat_data, error_code, error_msg] = tle_propagation(t_start, t_stop, delta_t, PARA, sat_name_to_prop);
                        if error_code > 0
                            error_msg = ['tle_propagation:', error_msg];
                            return; 
                        end
                        
                        number_temp_epochs = length(temp_sat_data.sat(1).prop);
                        
                        for i_temp_epoch = 1 : number_temp_epochs

                            [az, el, range_s_obs, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate] = eci2topo(temp_sat_data.sat(1).prop(i_temp_epoch).jd, temp_sat_data.sat(1).prop(i_temp_epoch).r, temp_sat_data.sat(1).prop(i_temp_epoch).v, stat_lon, stat_lat, stat_alt);

                            [sun_dist] = calc_sundist(temp_sat_data.sat(1).prop(i_temp_epoch).jd, stat_lon, stat_lat, az, el);
                            
                            % Break Condition
                            if (sun_dist < stat_data.stat(i_stat).min_sun_dist) % el_0 <  stat_data.stat(i_stat).min_elevation
                                if (i_temp_epoch == 1)
                                    t_start = temp_sat_data.sat(1).prop(i_temp_epoch).jd - delta_t / (24*60);
                                    t_stop  = temp_sat_data.sat(1).prop(i_temp_epoch).jd;
                                else
                                    t_start = temp_sat_data.sat(1).prop(i_temp_epoch - 1).jd;
                                    t_stop  = temp_sat_data.sat(1).prop(i_temp_epoch).jd;
                                end
                                break;
                            end;

                        end; % i_temp_epoch = 1 : number_temp_epochs
                        
                    end; % while (delta_t >= (0.1 / 60))
                    
                    t_exceed_min_sun_dist = t_start;
                    [temp_sat_data, error_code, error_msg] = tle_propagation(t_exceed_min_sun_dist, t_exceed_min_sun_dist, delta_t, PARA, sat_name_to_prop);
                    if error_code > 0
                        error_msg = ['tle_propagation:', error_msg];
                        return; 
                    end
                    
                    [az, el, range_s_obs, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate] = eci2topo(temp_sat_data.sat(1).prop(1).jd, temp_sat_data.sat(1).prop(1).r, temp_sat_data.sat(1).prop(1).v, stat_lon, stat_lat, stat_alt);
                    [sun_dist] = calc_sundist(temp_sat_data.sat(1).prop(1).jd, stat_lon, stat_lat, az, el);
                            
                    
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_exceed_min_sundist = stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_exceed_min_sundist + 1;
                    
                    [year, mon, day, hr, min, sec] = invjday(t_exceed_min_sun_dist);

                    stat_data.stat(i_stat).sat(i_sat).sun_dist.exceed_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_exceed_min_sundist).year     = year;                            
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.exceed_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_exceed_min_sundist).mon      = mon;
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.exceed_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_exceed_min_sundist).day      = day;
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.exceed_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_exceed_min_sundist).h        = hr;
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.exceed_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_exceed_min_sundist).min      = min;
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.exceed_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_exceed_min_sundist).sec      = sec;
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.exceed_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_exceed_min_sundist).jd       = t_exceed_min_sun_dist;
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.exceed_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_exceed_min_sundist).az       = az;
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.exceed_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_exceed_min_sundist).el       = el;
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.exceed_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_exceed_min_sundist).range    = range_s_obs;
                    stat_data.stat(i_stat).sat(i_sat).sun_dist.exceed_min_sundist(stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_exceed_min_sundist).sun_dist = sun_dist;
               
                end; % if (...... Find "exceed min. sun distance"
                

                %% ----- Antenna Axis 1 rate (AZ, HA) -----
                
                % Find "keep_max_axis1_rate":
                if ( (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_max_axis1_rate == 0)     && ...
                     (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).exceed_max_axis1_rate == 1)         ) 
                 
                    % Start parameter for the iterative rise search:
                    az_rate_residual_limit = 0.01; %[deg/s]
                    az_rate = 0.0;
                    t_start             = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).jd;  % !
                    t_stop              = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).jd;      % !
                    delta_t             = stat_data.stat(i_stat).prop_setup.delta_t_min;            % !
                    sat_name_to_prop = stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_header_line;
                    
                    % Iterative search approach:
                    while ( (delta_t >= (0.0001 / 60)) && (abs(az_rate - stat_data.stat(i_stat).max_axis1_rate) > az_rate_residual_limit) )   % /60 = conversion from sec to min!
                        
                        delta_t = delta_t / 10;

                        [temp_sat_data, error_code, error_msg] = tle_propagation(t_start, t_stop, delta_t, PARA, sat_name_to_prop);
                        if error_code > 0
                            error_msg = ['tle_propagation:', error_msg];
                            return; 
                        end
                        
                        number_temp_epochs = length(temp_sat_data.sat(1).prop);
                        
                        for i_temp_epoch = 1 : number_temp_epochs
                            
                            
                            % ### Break Condition for the iteration ###
                            
                            % Calc. slew rates for the specific antenna mount type:
                            [az, el, range_s_obs, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate] = eci2topo(temp_sat_data.sat(1).prop(i_temp_epoch).jd, temp_sat_data.sat(1).prop(i_temp_epoch).r, temp_sat_data.sat(1).prop(i_temp_epoch).v, stat_lon, stat_lat, stat_alt);

                            switch(stat_data.stat(i_stat).axis_type)
                                case 'AZEL'
                                    slew_rate_deg = az_rate;
                                case 'HADC'
                                    slew_rate_deg = local_hour_angle_rate;
                                case 'XYEW'
                                    % Conversion AzEl => XYew (in [deg] and [deg/sec]):
                                    [x, y, error_code, error_msg, x_rate, y_rate] = azel2xyew(az, el, 1, az_rate, el_rate);
                                    if error_code ~= 0
                                        error_msg = ['azel2xyew:', error_msg];
                                        return;
                                    end
                                    slew_rate_deg = x_rate;
                                otherwise
                                    error_code = 2;
                                    error_msg = 'Unknown antenna mount type.';
                                    return;
                            end % switch(stat_data.stat(i_stat).axis_type)

                            if (abs(slew_rate_deg) < stat_data.stat(i_stat).max_axis1_rate)
                                
                                if (i_temp_epoch == 1)
                                    t_start = temp_sat_data.sat(1).prop(i_temp_epoch).jd - delta_t / (24*60);
                                    t_stop  = temp_sat_data.sat(1).prop(i_temp_epoch).jd;
                                else
                                    t_start = temp_sat_data.sat(1).prop(i_temp_epoch - 1).jd;
                                    t_stop  = temp_sat_data.sat(1).prop(i_temp_epoch).jd;
                                end
                                break;
                            end; % if (slew_rate_deg < stat_data.stat(i_stat).max_axis1_rate)
                            
                        end; % i_temp_epoch = 1 : number_temp_epochs
                        
                    end; % while (delta_t >= (0.1 / 60))
                    
                    t_keep_max_axis1_rate = t_stop;
                    [temp_sat_data, error_code, error_msg] = tle_propagation(t_keep_max_axis1_rate, t_keep_max_axis1_rate, delta_t, PARA, sat_name_to_prop);
                    if error_code > 0
                        error_msg = ['tle_propagation:', error_msg];
                        return; 
                    end   
                        
                    [az, el, range_s_obs, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate] = eci2topo(temp_sat_data.sat(1).prop(1).jd, temp_sat_data.sat(1).prop(1).r, temp_sat_data.sat(1).prop(1).v, stat_lon, stat_lat, stat_alt);   
                    
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_keep_max_axis1_rate = stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_keep_max_axis1_rate + 1;
                    
                    [year, mon, day, hr, min, sec] = invjday(t_keep_max_axis1_rate);

                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.keep_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_keep_max_axis1_rate).year     = year;                            
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.keep_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_keep_max_axis1_rate).mon      = mon;
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.keep_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_keep_max_axis1_rate).day      = day;
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.keep_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_keep_max_axis1_rate).h        = hr;
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.keep_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_keep_max_axis1_rate).min      = min;
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.keep_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_keep_max_axis1_rate).sec      = sec;
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.keep_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_keep_max_axis1_rate).jd       = t_keep_max_axis1_rate;
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.keep_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_keep_max_axis1_rate).az       = az;
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.keep_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_keep_max_axis1_rate).el       = el;
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.keep_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_keep_max_axis1_rate).range    = range_s_obs;
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.keep_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_keep_max_axis1_rate).az_rate   = az_rate;

                end; % if (...... Find "keep_max_axis1_rate"
                
                
                % Find "exceed_max_axis1_rate":
                if ( (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_max_axis1_rate == 1)     && ...
                     (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).exceed_max_axis1_rate == 0)         ) 
                 
                    % Start parameter for the iterative rise search:
                    az_rate_residual_limit = 0.01; %[deg/s]
                    az_rate = 0.0;
                    t_start             = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).jd;  % !
                    t_stop              = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).jd;      % !
                    delta_t             = stat_data.stat(i_stat).prop_setup.delta_t_min;            % !
                    sat_name_to_prop = stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_header_line;
                    
                    % Iterative search approach:
                    while ( (delta_t >= (0.0001 / 60)) && (abs(az_rate - stat_data.stat(i_stat).max_axis1_rate) > az_rate_residual_limit) )   % /60 = conversion from sec to min!
                        
                        delta_t = delta_t / 10;

                        [temp_sat_data, error_code, error_msg] = tle_propagation(t_start, t_stop, delta_t, PARA, sat_name_to_prop);
                        if error_code > 0
                            error_msg = ['tle_propagation:', error_msg];
                            return; 
                        end
                        
                        number_temp_epochs = length(temp_sat_data.sat(1).prop);
                        
                        for i_temp_epoch = 1 : number_temp_epochs

                            % ### Break Condition for the iteration ###
                            
                            % Calc. slew rates for the specific antenna mount type:
                            [az, el, range_s_obs, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate] = eci2topo(temp_sat_data.sat(1).prop(i_temp_epoch).jd, temp_sat_data.sat(1).prop(i_temp_epoch).r, temp_sat_data.sat(1).prop(i_temp_epoch).v, stat_lon, stat_lat, stat_alt);
                            
                            switch(stat_data.stat(i_stat).axis_type)
                                case 'AZEL'
                                    slew_rate_deg = az_rate;
                                case 'HADC'
                                    slew_rate_deg = local_hour_angle_rate;
                                case 'XYEW'
                                    % Conversion AzEl => XYew (in [deg] and [deg/sec]):
                                    [x, y, error_code, error_msg, x_rate, y_rate] = azel2xyew(az, el, 1, az_rate, el_rate);
                                    if error_code ~= 0
                                        error_msg = ['azel2xyew:', error_msg];
                                        return;
                                    end
                                    slew_rate_deg = x_rate;
                                otherwise
                                    error_code = 2;
                                    error_msg = 'Unknown antenna mount type.';
                                    return;
                            end % switch(stat_data.stat(i_stat).axis_type)
                            
                            if (abs(slew_rate_deg) > stat_data.stat(i_stat).max_axis1_rate)
                                if (i_temp_epoch == 1)
                                    t_start = temp_sat_data.sat(1).prop(i_temp_epoch).jd - delta_t / (24*60);
                                    t_stop  = temp_sat_data.sat(1).prop(i_temp_epoch).jd;                                  
                                else
                                    t_start = temp_sat_data.sat(1).prop(i_temp_epoch - 1).jd;
                                    t_stop  = temp_sat_data.sat(1).prop(i_temp_epoch).jd;
                                end
                                break;
                            end;

                        end; % i_temp_epoch = 1 : number_temp_epochs
                        
                    end; % while (delta_t >= (0.1 / 60))
                    
                    t_exceed_max_axis1_rate = t_start;
                    [temp_sat_data, error_code, error_msg] = tle_propagation(t_exceed_max_axis1_rate, t_exceed_max_axis1_rate, delta_t, PARA, sat_name_to_prop);
                        if error_code > 0
                            error_msg = ['tle_propagation:', error_msg];
                            return; 
                        end
                        
                    [az, el, range_s_obs, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate] = eci2topo(temp_sat_data.sat(1).prop(1).jd, temp_sat_data.sat(1).prop(1).r, temp_sat_data.sat(1).prop(1).v, stat_lon, stat_lat, stat_alt);     
                    
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_exceed_max_axis1_rate = stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_exceed_max_axis1_rate + 1;
                    
                    [year, mon, day, hr, min, sec] = invjday(t_exceed_max_axis1_rate);

                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.exceed_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_exceed_max_axis1_rate).year     = year;                            
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.exceed_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_exceed_max_axis1_rate).mon      = mon;
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.exceed_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_exceed_max_axis1_rate).day      = day;
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.exceed_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_exceed_max_axis1_rate).h        = hr;
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.exceed_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_exceed_max_axis1_rate).min      = min;
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.exceed_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_exceed_max_axis1_rate).sec      = sec;
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.exceed_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_exceed_max_axis1_rate).jd       = t_exceed_max_axis1_rate;
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.exceed_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_exceed_max_axis1_rate).az       = az;
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.exceed_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_exceed_max_axis1_rate).el       = el;
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.exceed_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_exceed_max_axis1_rate).range    = range_s_obs;
                    stat_data.stat(i_stat).sat(i_sat).axis1_rate.exceed_max_axis1_rate(stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_exceed_max_axis1_rate).axis1_rate  = az_rate;
                    
                end; % if (...... Find "exceed_max_axis1_rate"
                
                
                
                
                
                %% ----- Antenna Axis 2 rate (EL, DC) -----    
                
                % Find "keep_max_axis2_rate":
                if ( (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_max_axis2_rate == 0)     && ...
                     (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).exceed_max_axis2_rate == 1)         ) 
                 
                    % Start parameter for the iterative rise search:
                    el_rate_residual_limit = 0.01; %[deg/s]
                    el_rate = 0.0;
                    t_start             = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).jd;  % !
                    t_stop              = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).jd;      % !
                    delta_t             = stat_data.stat(i_stat).prop_setup.delta_t_min;            % !
                    sat_name_to_prop = stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_header_line;
                    
                    % Iterative rise search approach:
                    while ( (delta_t >= (0.0001 / 60)) && (abs(el_rate - stat_data.stat(i_stat).max_axis2_rate) > el_rate_residual_limit) )   % /60 = conversion from sec to min!
                        
                        delta_t = delta_t / 10;

                        [temp_sat_data, error_code, error_msg] = tle_propagation(t_start, t_stop, delta_t, PARA, sat_name_to_prop);
                        if error_code > 0
                            error_msg = ['tle_propagation:', error_msg];
                            return; 
                        end
                        
                        number_temp_epochs = length(temp_sat_data.sat(1).prop);
                        
                        for i_temp_epoch = 1 : number_temp_epochs

                            % ### Break Condition for the iteration ###
                            
                            % Calc. slew rates for the specific antenna mount type:
                            [az, el, range_s_obs, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate] = eci2topo(temp_sat_data.sat(1).prop(i_temp_epoch).jd, temp_sat_data.sat(1).prop(i_temp_epoch).r, temp_sat_data.sat(1).prop(i_temp_epoch).v, stat_lon, stat_lat, stat_alt);
                            
                            switch(stat_data.stat(i_stat).axis_type)
                                case 'AZEL'
                                    slew_rate_deg = el_rate;
                                case 'HADC'
                                    slew_rate_deg = tdec_rate;
                                case 'XYEW'
                                    % Conversion AzEl => XYew (in [deg] and [deg/sec]):
                                    [x, y, error_code, error_msg, x_rate, y_rate] = azel2xyew(az, el, 1, az_rate, el_rate);
                                    if error_code ~= 0
                                        error_msg = ['azel2xyew:', error_msg];
                                        return;
                                    end
                                    slew_rate_deg = y_rate;
                                otherwise
                                    error_code = 2;
                                    error_msg = 'Unknown antenna mount type.';
                                    return;
                            end % switch(stat_data.stat(i_stat).axis_type)

                            if (abs(slew_rate_deg) < stat_data.stat(i_stat).max_axis2_rate)
                                if (i_temp_epoch == 1)
                                    t_start = temp_sat_data.sat(1).prop(i_temp_epoch).jd - delta_t / (24*60);
                                    t_stop  = temp_sat_data.sat(1).prop(i_temp_epoch).jd;
                                else
                                    t_start = temp_sat_data.sat(1).prop(i_temp_epoch - 1).jd;
                                    t_stop  = temp_sat_data.sat(1).prop(i_temp_epoch).jd;
                                end
                                break;
                            end;

                        end; % i_temp_epoch = 1 : number_temp_epochs
                        
                    end; % while (delta_t >= (0.1 / 60))
                    
                    t_keep_max_axis2_rate = t_stop;
                    [temp_sat_data, error_code, error_msg] = tle_propagation(t_keep_max_axis2_rate, t_keep_max_axis2_rate, delta_t, PARA, sat_name_to_prop);
                    if error_code > 0
                        error_msg = ['tle_propagation:', error_msg];
                        return; 
                    end

                    [az, el, range_s_obs, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate] = eci2topo(temp_sat_data.sat(1).prop(1).jd, temp_sat_data.sat(1).prop(1).r, temp_sat_data.sat(1).prop(1).v, stat_lon, stat_lat, stat_alt);   
                    
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_keep_max_axis2_rate = stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_keep_max_axis2_rate + 1;
                    
                    [year, mon, day, hr, min, sec] = invjday(t_keep_max_axis2_rate);

                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.keep_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_keep_max_axis2_rate).year     = year;                            
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.keep_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_keep_max_axis2_rate).mon      = mon;
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.keep_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_keep_max_axis2_rate).day      = day;
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.keep_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_keep_max_axis2_rate).h        = hr;
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.keep_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_keep_max_axis2_rate).min      = min;
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.keep_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_keep_max_axis2_rate).sec      = sec;
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.keep_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_keep_max_axis2_rate).jd       = t_keep_max_axis2_rate;
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.keep_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_keep_max_axis2_rate).az       = az;
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.keep_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_keep_max_axis2_rate).el       = el;
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.keep_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_keep_max_axis2_rate).range    = range_s_obs;
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.keep_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_keep_max_axis2_rate).el_rate   = el_rate;

                end; % if (...... Find "keep_max_axis2_rate"
                
                
                % Find "exceed_max_axis2_rate":
                if ( (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_max_axis2_rate == 1)     && ...
                     (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).exceed_max_axis2_rate == 0)         ) 
                 
                    % Start parameter for the iterative rise search:
                    el_rate_residual_limit = 0.01; %[deg/s]
                    el_rate = 0.0;
                    t_start             = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).jd;  % !
                    t_stop              = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).jd;      % !
                    delta_t             = stat_data.stat(i_stat).prop_setup.delta_t_min;            % !
                    sat_name_to_prop = stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_header_line;
                    
                    % Iterative rise search approach:
                    while ( (delta_t >= (0.0001 / 60)) && (abs(el_rate - stat_data.stat(i_stat).max_axis2_rate) > el_rate_residual_limit) )   % /60 = conversion from sec to min!
                        
                        delta_t = delta_t / 10;

                        [temp_sat_data, error_code, error_msg] = tle_propagation(t_start, t_stop, delta_t, PARA, sat_name_to_prop);
                        if error_code > 0
                            error_msg = ['tle_propagation:', error_msg];
                            return; 
                        end
                        
                        number_temp_epochs = length(temp_sat_data.sat(1).prop);
                        
                        for i_temp_epoch = 1 : number_temp_epochs

                            
                            % ### Break Condition for the iteration ###
                            
                            % Calc. slew rates for the specific antenna mount type:
                            [az, el, range_s_obs, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate] = eci2topo(temp_sat_data.sat(1).prop(i_temp_epoch).jd, temp_sat_data.sat(1).prop(i_temp_epoch).r, temp_sat_data.sat(1).prop(i_temp_epoch).v, stat_lon, stat_lat, stat_alt);
                            
                            switch(stat_data.stat(i_stat).axis_type)
                                case 'AZEL'
                                    slew_rate_deg = el_rate;
                                case 'HADC'
                                    slew_rate_deg = tdec_rate;
                                case 'XYEW'
                                    % Conversion AzEl => XYew (in [deg] and [deg/sec]):
                                    [x, y, error_code, error_msg, x_rate, y_rate] = azel2xyew(az, el, 1, az_rate, el_rate);
                                    if error_code ~= 0
                                        error_msg = ['azel2xyew:', error_msg];
                                        return;
                                    end
                                    slew_rate_deg = y_rate;
                                otherwise
                                    error_code = 2;
                                    error_msg = 'Unknown antenna mount type.';
                                    return;
                            end % switch(stat_data.stat(i_stat).axis_type)

                            if (abs(slew_rate_deg) > stat_data.stat(i_stat).max_axis2_rate)
                                if (i_temp_epoch == 1)
                                    t_start = temp_sat_data.sat(1).prop(i_temp_epoch).jd - delta_t / (24*60);
                                    t_stop  = temp_sat_data.sat(1).prop(i_temp_epoch).jd;
                                else
                                    t_start = temp_sat_data.sat(1).prop(i_temp_epoch - 1).jd;
                                    t_stop  = temp_sat_data.sat(1).prop(i_temp_epoch).jd;
                                end
                                break;
                            end;

                        end; % i_temp_epoch = 1 : number_temp_epochs
                        
                    end; % while (delta_t >= (0.1 / 60))
                    
                    t_exceed_max_axis2_rate = t_start;
                    [temp_sat_data, error_code, error_msg] = tle_propagation(t_exceed_max_axis2_rate, t_exceed_max_axis2_rate, delta_t, PARA, sat_name_to_prop);
                    if error_code > 0
                        error_msg = ['tle_propagation:', error_msg];
                        return; 
                    end

                    [az, el, range_s_obs, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate] = eci2topo(temp_sat_data.sat(1).prop(1).jd, temp_sat_data.sat(1).prop(1).r, temp_sat_data.sat(1).prop(1).v, stat_lon, stat_lat, stat_alt);     
                    
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_exceed_max_axis2_rate = stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_exceed_max_axis2_rate + 1;
                    
                    [year, mon, day, hr, min, sec] = invjday(t_exceed_max_axis2_rate);

                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.exceed_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_exceed_max_axis2_rate).year     = year;                            
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.exceed_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_exceed_max_axis2_rate).mon      = mon;
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.exceed_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_exceed_max_axis2_rate).day      = day;
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.exceed_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_exceed_max_axis2_rate).h        = hr;
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.exceed_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_exceed_max_axis2_rate).min      = min;
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.exceed_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_exceed_max_axis2_rate).sec      = sec;
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.exceed_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_exceed_max_axis2_rate).jd       = t_exceed_max_axis2_rate;
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.exceed_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_exceed_max_axis2_rate).az       = az;
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.exceed_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_exceed_max_axis2_rate).el       = el;
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.exceed_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_exceed_max_axis2_rate).range    = range_s_obs;
                    stat_data.stat(i_stat).sat(i_sat).axis2_rate.exceed_max_axis2_rate(stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_exceed_max_axis2_rate).axis2_rate  = el_rate;
                    
                end; % if (...... Find "exceed_max_axis2_rate"
                
                
                
                %% ----- Slew range limits -----
                
                % Find "keep_slew_range_limits":
                if ( (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_axis_limits == 0)     && ...
                     (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).exceed_axis_limits == 1)         ) 
                 
                    % Start parameter for the iterative rise search:
                    t_start             = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).jd;  % !
                    t_stop              = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).jd;      % !
                    delta_t             = stat_data.stat(i_stat).prop_setup.delta_t_min;            % delta_t in [min] !
                    sat_name_to_prop = stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_header_line;
                    
                    % Iterative search approach:
                    while (delta_t >= (0.01 / 60))  % /60 = conversion from sec to min!
                        
                        delta_t = delta_t / 10;

                        [temp_sat_data, error_code, error_msg] = tle_propagation(t_start, t_stop, delta_t, PARA, sat_name_to_prop);
                        if error_code > 0
                            error_msg = ['tle_propagation:', error_msg];
                            return; 
                        end
                        
                        number_temp_epochs = length(temp_sat_data.sat(1).prop);
                        
                        for i_temp_epoch = 1 : number_temp_epochs

                            [az, el, range_s_obs, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate] = eci2topo(temp_sat_data.sat(1).prop(i_temp_epoch).jd, temp_sat_data.sat(1).prop(i_temp_epoch).r, temp_sat_data.sat(1).prop(i_temp_epoch).v, stat_lon, stat_lat, stat_alt);
                            
                            [flag_axis_limits_ok, error_code, error_msg] = check_slew_range_limits(stat_data, i_stat, az*pi/180, el*pi/180, local_hour_angle*pi/180, tdec*pi/180);
                            if error_code ~= 0
                                error_msg = ['check_slew_range_limits:', error_msg];
                                return;
                            end
                            
                            % Break Condition
                            if (flag_axis_limits_ok)
                                if (i_temp_epoch == 1)
                                    t_start = temp_sat_data.sat(1).prop(i_temp_epoch).jd - delta_t / (24*60);
                                    t_stop  = temp_sat_data.sat(1).prop(i_temp_epoch).jd;
                                else
                                    t_start = temp_sat_data.sat(1).prop(i_temp_epoch - 1).jd;
                                    t_stop  = temp_sat_data.sat(1).prop(i_temp_epoch).jd;
                                end
                                break;
                            end
                                                        
                        end % i_temp_epoch = 1 : number_temp_epochs
                        
                    end %while (delta_t >= (0.01 / 60))

                    t_keep_slew_range_limits = t_stop;
                    
                    [temp_sat_data, error_code, error_msg] = tle_propagation(t_keep_slew_range_limits, t_keep_slew_range_limits, delta_t, PARA, sat_name_to_prop);
                    if error_code > 0
                        error_msg = ['tle_propagation:', error_msg];
                        return; 
                    end
                    
                    [az, el, range_s_obs, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate] = eci2topo(temp_sat_data.sat(1).prop(1).jd, temp_sat_data.sat(1).prop(1).r, temp_sat_data.sat(1).prop(1).v, stat_lon, stat_lat, stat_alt);
                            
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_keep_slew_range_limits = stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_keep_slew_range_limits + 1;
                    
                    [year, mon, day, hr, min, sec] = invjday(t_keep_slew_range_limits);

                    % Save data:
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.keep_slew_range_limits(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_keep_slew_range_limits).year     = year;                            
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.keep_slew_range_limits(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_keep_slew_range_limits).mon      = mon;
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.keep_slew_range_limits(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_keep_slew_range_limits).day      = day;
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.keep_slew_range_limits(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_keep_slew_range_limits).h        = hr;
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.keep_slew_range_limits(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_keep_slew_range_limits).min      = min;
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.keep_slew_range_limits(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_keep_slew_range_limits).sec      = sec;
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.keep_slew_range_limits(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_keep_slew_range_limits).jd       = t_keep_slew_range_limits;
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.keep_slew_range_limits(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_keep_slew_range_limits).az       = az;
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.keep_slew_range_limits(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_keep_slew_range_limits).el       = el;
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.keep_slew_range_limits(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_keep_slew_range_limits).range    = range_s_obs;

               
                end % if (...... % Find "keep_slew_range_limits"
                

                % Find "exceed_slew_range_limits":
                if ( (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_axis_limits == 1)     && ...
                     (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).exceed_axis_limits == 0)         ) 
                    
                    % Start parameter for the iterative rise search:
                    t_start             = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch - 1).jd;  % !
                    t_stop              = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).jd;      % !
                    delta_t             = stat_data.stat(i_stat).prop_setup.delta_t_min;            % delta_t in [min] !
                    sat_name_to_prop = stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_header_line;
                    
                    % Iterative search approach:
                    while (delta_t >= (0.01 / 60))  % /60 = conversion from sec to min!
                        
                        delta_t = delta_t / 10;
                        [temp_sat_data, error_code, error_msg] = tle_propagation(t_start, t_stop, delta_t, PARA, sat_name_to_prop);
                        if error_code > 0
                            error_msg = ['tle_propagation:', error_msg];
                            return; 
                        end
                        
                        number_temp_epochs = length(temp_sat_data.sat(1).prop);
                        
                        for i_temp_epoch = 1 : number_temp_epochs
                            
                            [az, el, range_s_obs, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate] = eci2topo(temp_sat_data.sat(1).prop(i_temp_epoch).jd, temp_sat_data.sat(1).prop(i_temp_epoch).r, temp_sat_data.sat(1).prop(i_temp_epoch).v, stat_lon, stat_lat, stat_alt);
                            
                            [flag_axis_limits_ok, error_code, error_msg] = check_slew_range_limits(stat_data, i_stat, az*pi/180, el*pi/180, local_hour_angle*pi/180, tdec*pi/180);
                            if error_code ~= 0
                                error_msg = ['check_slew_range_limits:', error_msg];
                                return;
                            end
                            
                            % Break Condition
                            if (flag_axis_limits_ok == 0)
                                if (i_temp_epoch == 1)
                                    t_start = temp_sat_data.sat(1).prop(i_temp_epoch).jd - delta_t / (24*60);
                                    t_stop  = temp_sat_data.sat(1).prop(i_temp_epoch).jd;
                                else
                                    t_start = temp_sat_data.sat(1).prop(i_temp_epoch - 1).jd;
                                    t_stop  = temp_sat_data.sat(1).prop(i_temp_epoch).jd;
                                end
                                break;
                            end

                        end % i_temp_epoch = 1 : number_temp_epochs
                        
                    end; % while (delta_t >= (0.1 / 60))
                    
                    
                    t_exceed_slew_range_limits = t_start;
                    
                    [temp_sat_data, error_code, error_msg] = tle_propagation(t_exceed_slew_range_limits, t_exceed_slew_range_limits, delta_t, PARA, sat_name_to_prop);
                    if error_code > 0
                        error_msg = ['tle_propagation:', error_msg];
                        return; 
                    end
                    
                    [az, el, range_s_obs, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate] = eci2topo(temp_sat_data.sat(1).prop(1).jd, temp_sat_data.sat(1).prop(1).r, temp_sat_data.sat(1).prop(1).v, stat_lon, stat_lat, stat_alt);                            
                    
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_exceed_slew_range_limits = stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_exceed_slew_range_limits + 1;
                    
                    [year, mon, day, hr, min, sec] = invjday(t_exceed_slew_range_limits);

                    % Save data:
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.exceed_slew_range_limits(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_exceed_slew_range_limits).year     = year;                            
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.exceed_slew_range_limits(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_exceed_slew_range_limits).mon      = mon;
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.exceed_slew_range_limits(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_exceed_slew_range_limits).day      = day;
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.exceed_slew_range_limits(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_exceed_slew_range_limits).h        = hr;
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.exceed_slew_range_limits(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_exceed_slew_range_limits).min      = min;
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.exceed_slew_range_limits(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_exceed_slew_range_limits).sec      = sec;
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.exceed_slew_range_limits(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_exceed_slew_range_limits).jd       = t_exceed_slew_range_limits;
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.exceed_slew_range_limits(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_exceed_slew_range_limits).az       = az;
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.exceed_slew_range_limits(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_exceed_slew_range_limits).el       = el;
                    stat_data.stat(i_stat).sat(i_sat).slew_range_limits.exceed_slew_range_limits(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_exceed_slew_range_limits).range    = range_s_obs;
               
                end % if (...... Find "exceed min. sun distance"

                
            end    % i_epoch = 1 : number_of_epochs
            
        end    % for i_sat = 1 : number_of_sat

    end    % for i_stat = 1 : stat_data.number_of_stations

return;



