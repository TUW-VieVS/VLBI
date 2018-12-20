% -------------------------------------------------------------------------
%
%                              calc_pointing_data.m
%
%   This function calculates Antenna pointing data in terms of azimut and
%   elevation angls, topocentric right ascension and declination, local 
%   hour angle and declination, their
%   according rates and the sun distance, for an arbtrary number of 
%   staionss and satellites.
%
%       Initial checks for:
%           - visibility (above min elevation, horizontal mask)
%           - slew rates (axes 1 {az, ha} and axes 2 {el, dec})
%           - sun dist
%           - 
%
%       Note:   Support of 'XYSN' mounted Antennas is not implemented
%               yet!
%
%
%   Author: 
%       Andreas Hellerschmied, 22.9.2013
%   
%   changes       :
%   - 27.10.2013 : A. Hellerschmied, Calculation of Az/El - rates and 
%                      sun dist. implemented. 
%   - 2014-01.19: A. Hellerschmied: Support of different antenna axis types
%   - 2014-01.19: A. Hellerschmied: Error treatment added.
%            Remark: Output of coupled procedures is currently not checked!
%   - 2015-02-16: A. Hellerschmied: Remanded (before:
%            "TLE_calc_pointing_data.m")
%   - 2015-07-29: A. Hellerschmied: Horizintal mask support added.
%   - 2015-07-30: A. Hellerschmied: Bug-fix, error in units.
%   - 2015-07-30: A. Hellerschmied: Check for antenna sxis limits added (function "check_slew_range_limits")
%   - 2016-04-25: A. Hellerschmied: Bug-fix: Check of slew-rate limits => Now the absolute value is checked! 
%           
%
%   inputs        :
%       - stat_data: station data structure
%       - sat_data: satellite data structure
%     
%
%   outputs       :
%   - stat_data: station data structure, containing additional
%                    pointing data.
%   - error_code    : Error Code (0 = no erros occured)
%   - error_msg     : Error Message (empty, if no errors occured)
%    
%
%   locals        :
% 
%
%   coupling      :
%       - eci2topo.m
%       - calc_sundist.m
%       - azel2xyew.m
%       - check_h_mask.m
%       - check_slew_range_limits
%   
%
%   references    :
%
%-------------------------------------------------------------------------


function [stat_data, error_code, error_msg] = calc_pointing_data(sat_data, stat_data)

    % Init
    error_code = 0;
    error_msg = '';

    stat_data.number_of_sat = length(sat_data.sat);
    stat_data.number_of_epochs = length(sat_data.sat(1).prop);
    
%     % ##### Open Waitbar #####
%     total_number_of_epochs = stat_data.number_of_sat * stat_data.number_of_sat * stat_data.number_of_epochs;
%     epoch_count = 0;
%     h_wb = waitbar(0, '1' ,'Name', 'Antenna pointing data'); 
    tic
    
    for i_stat = 1 : stat_data.number_of_stations    % Stations
        
        stat_lon = stat_data.stat(i_stat).location.ellipsoid.long;
        stat_lat = stat_data.stat(i_stat).location.ellipsoid.lat;
        stat_alt = stat_data.stat(i_stat).location.ellipsoid.altitude;
        
        stat_data.stat(i_stat).prop_setup.delta_t_min = sat_data.prop_setup.delta_t_min;
        
        
        for i_sat = 1 : stat_data.number_of_sat   % Satellites
            
            stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_header_line = sat_data.sat(i_sat).TLE_header_line;
            stat_data.stat(i_stat).sat(i_sat).TLE_data.sat_number      = sat_data.sat(i_sat).sat_number;
            stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_filepath    = sat_data.prop_setup.tle_filepath;
            stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_filename    = sat_data.prop_setup.tle_filename;
            stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_epoch.year  = sat_data.sat(i_sat).epoch.year; 
            stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_epoch.mon   = sat_data.sat(i_sat).epoch.mon;
            stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_epoch.day   = sat_data.sat(i_sat).epoch.day;
            stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_epoch.h     = sat_data.sat(i_sat).epoch.h;
            stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_epoch.min   = sat_data.sat(i_sat).epoch.min;
            stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_epoch.sec   = sat_data.sat(i_sat).epoch.sec;
            stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_epoch.jd    = sat_data.sat(i_sat).epoch.jd;

            
            for i_epoch = 1 : stat_data.number_of_epochs  % Propagation epochs
                
%                 % Update waitbar
%                 epoch_count = epoch_count + 1;
%                 waitbar(epoch_count/total_number_of_epochs, h_wb, sprintf('Epoch %1.0f of %1.0f processed.',epoch_count, total_number_of_epochs))
                
                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).year   = sat_data.sat(i_sat).prop(i_epoch).year;
                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).mon    = sat_data.sat(i_sat).prop(i_epoch).mon;
                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).day    = sat_data.sat(i_sat).prop(i_epoch).day;
                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).h      = sat_data.sat(i_sat).prop(i_epoch).h;
                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).min    = sat_data.sat(i_sat).prop(i_epoch).min;
                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).sec    = sat_data.sat(i_sat).prop(i_epoch).sec;
                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).jd     = sat_data.sat(i_sat).prop(i_epoch).jd;
        
                % Calculation of Az, El, topo. Ra/Dec, LHA, range and the according rates for
                % the given epoch and station - satellite - constellation:
                % results of the function "eci2topo" are in [deg] or [deg/sec]!!!
                [Az, El, range_s_obs, az_rate, el_rate, r_vel, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate] = eci2topo(sat_data.sat(i_sat).prop(i_epoch).jd, sat_data.sat(i_sat).prop(i_epoch).r, sat_data.sat(i_sat).prop(i_epoch).v, stat_lon, stat_lat, stat_alt);
                
                
                % Convert RA from [-180, +180 deg] to [0, 360 deg]:
                if (tra < 0)
                    tra = tra + 360;
                end
                
                % Convert LHA from [-180, +180 deg] to [0, +360 deg]:
                while(local_hour_angle < 0)
                    local_hour_angle = local_hour_angle + 360;
                end
                while(local_hour_angle > 360)
                    local_hour_angle = local_hour_angle - 360;
                end
                

                % #### Above min. Elevation and horizontal mask? ####
                [flag_above_min_el] = check_h_mask(stat_data, i_stat, Az*pi/180, El*pi/180);
                if (El > stat_data.stat(i_stat).min_elevation) && (flag_above_min_el)
                    stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).above_min_elevation = 1;
                else
                    stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).above_min_elevation = 0; 
                end;
                
                
                % #### Exceed max Axis 1 (AZ, HA) and Axis 2 (EL, DEC) slew rate of antenna? ####
                if (strncmp(stat_data.stat(i_stat).axis_type, 'AZEL', 4))
                    
                    % AZ-Axis
                    if (abs(az_rate) > stat_data.stat(i_stat).max_axis1_rate)
                        stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_max_axis1_rate = 1;
                    else
                        stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_max_axis1_rate = 0; 
                    end;
                    
                    % EL-Axis
                    if (abs(el_rate) > stat_data.stat(i_stat).max_axis2_rate)
                        stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_max_axis2_rate = 1;
                    else
                        stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_max_axis2_rate = 0; 
                    end;
                    
                elseif (strncmp(stat_data.stat(i_stat).axis_type, 'HADC', 4))
                    
                    % HA-Axis
                    if (abs(local_hour_angle_rate) > stat_data.stat(i_stat).max_axis1_rate)
                        stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_max_axis1_rate = 1;
                    else
                        stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_max_axis1_rate = 0; 
                    end;
                    
                    % DC-Axis
                    if (abs(tdec_rate) > stat_data.stat(i_stat).max_axis2_rate)
                        stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_max_axis2_rate = 1;
                    else
                        stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_max_axis2_rate = 0; 
                    end;
                      
                elseif (strncmp(stat_data.stat(i_stat).axis_type, 'XYNS', 4))
                    % Support of 'XYSN' mounted Antennas is not implemented
                    % yet! Still to do.
                    error_code = 1;
                    error_msg = ['"XYSN" mounted Antenna Axis are not supported yet.', '. Station: ' stat_data.stat(i_stat).name, '.'];
                    return;
                    
                elseif (strncmp(stat_data.stat(i_stat).axis_type, 'XYEW', 4)) 
                    
                    % Conversion AzEl => XYew (in [deg] and [deg/sec]):
                    [x, y, error_code, error_msg, x_rate, y_rate] = azel2xyew(Az, El, 1, az_rate, el_rate);
                    if error_code ~= 0
                        error_msg = ['azel2xyew:', error_msg];
                        return;
                    end
                    
                    % X-Axis
                    if (abs(x_rate) > stat_data.stat(i_stat).max_axis1_rate)
                        stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_max_axis1_rate = 1;
                    else
                        stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_max_axis1_rate = 0; 
                    end;
                    
                    % Y-Axis
                    if (abs(y_rate) > stat_data.stat(i_stat).max_axis2_rate)
                        stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_max_axis2_rate = 1;
                    else
                        stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_max_axis2_rate = 0; 
                    end;
                    
                else
                    error_code = 1;
                    error_msg = ['Unknown Antenna Axis Type: ', stat_data.stat(i_stat).axis_type, '. At station: ' stat_data.stat(i_stat).name, '.'];
                    return;
                end
                
                
                % ##### Axis limits #####
% Test:
% fprintf(1, '%s: LHA=%5.3f DEC=%5.3f\n', jd2datestr(sat_data.sat(i_sat).prop(i_epoch).jd), local_hour_angle, tdec);

                [flag_axis_limits_ok, error_code, error_msg] = check_slew_range_limits(stat_data, i_stat, Az*pi/180, El*pi/180, local_hour_angle*pi/180, tdec*pi/180);
                if error_code ~= 0
                    error_msg = ['check_slew_range_limits:', error_msg];
                    return;
                end
                if flag_axis_limits_ok
                    stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_axis_limits = 0; 
                else
                    stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_axis_limits = 1; 
                end

                
                % ##### Sun distance #####
                % Calculation of the sun distance [deg] from the satellite,
                % observed at the given station
                [sun_dist] = calc_sundist(sat_data.sat(i_sat).prop(i_epoch).jd, stat_lon, stat_lat, Az, El);
                
                % Keep min. Sun distance? 
                if (sun_dist > stat_data.stat(i_stat).min_sun_dist)
                    stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_min_sun_dist = 0;
                else
                    stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_min_sun_dist = 1; 
                end; 
                
                % Conversions:
                local_hour_angle = local_hour_angle / 15;           % deg => hour [0, 24 h]
                local_hour_angle_rate = local_hour_angle_rate / 15; % deg / sec => hour / sec
                tra = tra / 15;             % deg => hour [0, 24 h]
                tra_rate = tra_rate / 15;   % deg / sec => hour / sec
                
                % Write data to "stat_data" struct:
                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).az = Az;
                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).el = El;
                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).range  = range_s_obs;
                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).sun_dist = sun_dist;
                
                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).az_rate = az_rate;
                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).el_rate = el_rate;
                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).r_vel = r_vel;
                
                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).tra = tra;
                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).tdec = tdec;
                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).tra_rate = tra_rate;
                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).tdec_rate = tdec_rate;
                
                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).local_hour_angle = local_hour_angle;
                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).local_hour_angle_rate = local_hour_angle_rate;
                
            end;    % i_epoch = 1 : number_of_epochs
            
            % Test:
            % disp(sum_rate_lha);
            % disp(sum_tra_rate);
            
        end;    % for i_sat = 1 : number_of_sat

    end;    % for i_stat = 1 : stat_data.number_of_stations
    
%     % ##### Close waitbar (if it is still there...) #####
%     try
%         close(h_wb)
%     catch
%     end
    
    disp(toc);

return;



