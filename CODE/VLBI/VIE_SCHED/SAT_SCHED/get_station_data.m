% -------------------------------------------------------------------------
%
%                              get_station_data.m
%
%   Get station data and preallocate station data
%   structure (stat_data). Merge station data for satellite and/or quasar
%   observations into one struct, if available.
%
%   Author: 
%       Andreas Hellerschmied, 18.10.2013
%   
%   changes       : 
%   - 2014-01.19: A. Hellerschmied: Support of different antenna axis types
%   - 2014-01.19: A. Hellerschmied: Error treatment added.
%   - 2015-02-16: A. Hellerschmied: Remanded from "TLE_get_station_data.m"
%         to "get_station_data.m"; TRR_label referrs to resource of TRF
%         data.
%   - 2015-06-29: A. Hellerschmied: Fields for cable wrap section limits added
%   - 2015-07-30: A. Hellerschmied: Field "exceed_axis_limits" added
%   - 2016-11-10: A. Hellerschmied: Acceleration values are taken from station structures
%           
%
%   inputs        :
%       - PARA      : VieVS parameter structure
%       - station   : VieVS station structure
%       - INFILE    : File path structure
%
%   outputs       :
%   - stat_data     : station data structure
%   - error_code    : Error Code (0 = no erros occured)
%   - error_msg     : Error Message (empty, if no errors occured)
%    
%
%   locals        :
% 
%
%   coupling      :
%   - xyz2ell.m : Calculation of elipsoid coordinates (geodetic lat.)
%   
%
%   references    :
%
%-------------------------------------------------------------------------



function [stat_data, error_code, error_msg] = get_station_data(station_sat, station_quasar, PARA, INFILE) 

    % Init
    rad2deg = 180/pi;
    error_code = 0;                 % 0 = No error
    error_msg = '';                 % '' = No error
    
    

    % preallocating station data structure:
    stat_data = struct('stat', [], 'number_of_stations', [], 'number_of_stations_quasar', [], 'error', [],... 
                        'error_message', [], 'number_of_sat', [], 'number_of_epochs', []);
    
    stat_data.stat = struct('label', [], 'name', [], 'location', [], 'min_elevation', [],...
                        'max_axis2_rate', [], 'max_axis1_rate', [],  'max_axis1_acc', [], 'max_axis2_acc', [], 'lim11', [], 'lim12', [], 'lim21', [], 'lim22', [], 'c1', [], 'c2', [], 'horizontal_mask', [], 'horizontal_mask_num', [],  ...
                        'sefd_para', [], 'min_snr', [], 'axis_offset', [], 'az_n1', [], 'az_n2', [], 'az_cw1', [], 'az_cw2', [], 'az_ccw1', [], 'az_ccw2', [],...
                        'error', [], 'error_message', [], 'prop_setup', [], 'axis_type', [], 'min_sun_dist', [], 'obs_type', [], 'twin_partner', []); % obs_type = 'sat_only' or 'quasar_only' or 'sat_and_quasar' or 'twin_sat' or 'twin_quasar'
                    
    stat_data.stat.location = struct('ellipsoid', [], 'TRF', []);
    
    stat_data.stat.location.ellipsoid = struct('long', [], 'lat', [],'altitude', [], 'ref_ellipsoid_label', [] );
    
    stat_data.stat.location.TRF = struct('x', [], 'y', [], 'z', [], 'vx', [], 'vy', [], 'vz', [], 'epoch', [], 'start', [], 'end',  [], 'TRF_label', []);
    
    stat_data.stat.sat = struct('TLE_data', [], 'epoch', [], 'overpass', [], 'number_of_overpasses', [], 'sun_dist', [], 'az_rate', [], 'el_rate', [], ... 
                        'number_of_sun_dist_events', [], 'number_of_az_rate_events', [], 'number_of_el_rate_events', [], 'slew_range_limits', []);
    
    stat_data.stat.prop_setup =struct('delta_t_min', []);
    
    stat_data.stat.sat.TLE_data = struct('TLE_header_line', [], 'sat_number', [], 'TLE_filepath', [], 'TLE_filename', [], 'TLE_epoch', []);
    
    stat_data.stat.sat.TLE_data.TLE_epoch = struct('year', [], 'mon', [], 'day', [], 'h', [], 'min', [], 'sec', [], 'jd', []);
    
    
    stat_data.stat.sat.epoch = struct('year', [], 'mon', [], 'day', [], 'h', [], 'min', [], 'sec', [], 'jd', [],...
                        'az', [], 'el', [], 'range', [], 'above_min_elevation', [], 'az_rate', [], 'el_rate', [], 'range_rate', [],...
                        'exceed_max_axis1_rate', [], 'exceed_max_axis2_rate', [], 'az_accel', [], 'el_accel', [], ....
                        'exceed_max_axis1_accel', [], 'exceed_max_axis2_accel', [], 'sun_dist', [], 'exceed_min_sun_dist', [], 'tra', [], 'tdec', [],...
                        'tra_rate', [], 'tdec_rate', [], 'local_hour_angle', [], 'local_hour_angle_rate', [], 'exceed_axis_limits', []);
                    
    stat_data.stat.sat.overpass = struct('rise', [], 'set', [], 'peak', [], 'number_of_peaks', []);
    
    stat_data.stat.sat.overpass.rise = struct('year', [], 'mon', [], 'day', [], 'h', [], 'min', [], 'sec', [], 'jd', [],...
                        'az', [], 'el', [], 'range', []);
                    
    stat_data.stat.sat.overpass.set = struct('year', [], 'mon', [], 'day', [], 'h', [], 'min', [], 'sec', [], 'jd', [],...
                        'az', [], 'el', [], 'range', []);
                    
    stat_data.stat.sat.overpass.peak = struct('year', [], 'mon', [], 'day', [], 'h', [], 'min', [], 'sec', [], 'jd', [],...
                        'az', [], 'el', [], 'range', []);
                    
                    
    stat_data.stat.sat.sun_dist = struct('keep_min_sundist', [], 'exceed_min_sundist', [], 'number_of_keep_min_sundist', [], 'number_of_exceed_min_sundist', []);
                              
    stat_data.stat.sat.sun_dist.keep_min_sundist = struct('year', [], 'mon', [], 'day', [], 'h', [], 'min', [], 'sec', [], 'jd', [],...
                        'az', [], 'el', [], 'range', [], 'sun_dist', []);
                    
    stat_data.stat.sat.sun_dist.exceed_min_sundist = struct('year', [], 'mon', [], 'day', [], 'h', [], 'min', [], 'sec', [], 'jd', [],...
                        'az', [], 'el', [], 'range', [], 'sun_dist', []);
                    
                    
                    
                    
    stat_data.stat.sat.axis1_rate = struct('keep_max_axis1_rate', [], 'exceed_max_axis1_rate', [], 'number_of_keep_max_axis1_rate', [], 'number_of_exceed_max_axis1_rate', []);
                              
    stat_data.stat.sat.axis1_rate.keep_max_axis1_rate = struct('year', [], 'mon', [], 'day', [], 'h', [], 'min', [], 'sec', [], 'jd', [],...
                        'az', [], 'el', [], 'range', [], 'az_rate', []);
                    
    stat_data.stat.sat.axis1_rate.exceed_max_axis1_rate = struct('year', [], 'mon', [], 'day', [], 'h', [], 'min', [], 'sec', [], 'jd', [],...
                        'az', [], 'el', [], 'range', [], 'axis1_rate', []);
                    
                    
                    
    stat_data.stat.sat.axis2_rate = struct('keep_max_axis2_rate', [], 'exceed_max_axis2_rate', [], 'number_of_keep_max_axis2_rate', [], 'number_of_exceed_max_axis2_rate', []);
                              
    stat_data.stat.sat.axis2_rate.keep_max_axis2_rate = struct('year', [], 'mon', [], 'day', [], 'h', [], 'min', [], 'sec', [], 'jd', [],...
                        'az', [], 'el', [], 'range', [], 'el_rate', []);
                    
    stat_data.stat.sat.axis2_rate.exceed_max_axis2_rate = struct('year', [], 'mon', [], 'day', [], 'h', [], 'min', [], 'sec', [], 'jd', [],...
                        'az', [], 'el', [], 'range', [], 'axis2_rate', []);
                    
                    
                    
    stat_data.stat.sat.slew_range_limits = struct('keep_slew_range_limits', [], 'exceed_slew_range_limits', [], 'number_of_keep_slew_range_limits', [], 'number_of_exceed_slew_range_limits', []);
                              
    stat_data.stat.sat.slew_range_limits.keep_slew_range_limits = struct('year', [], 'mon', [], 'day', [], 'h', [], 'min', [], 'sec', [], 'jd', [],...
                        'az', [], 'el', [], 'range', [], 'el_rate', []);
                    
    stat_data.stat.sat.slew_range_limits.exceed_slew_range_limits = struct('year', [], 'mon', [], 'day', [], 'h', [], 'min', [], 'sec', [], 'jd', [],...
                        'az', [], 'el', [], 'range', [], 'axis2_rate', []);
                    
                    
                     
                    
    % Init:
    stat_data.error = 0;
    stat_data.error_message = '';
    
    % Only if separate staton networks for satellites/quasars are defined
%     if PARA.USE_SEPARATE_STAT_NWW_FOR_SATS

        % ##### Find stations observing satellites AND quasars #####

        % in station_quasar
        both_q = zeros(1,length(station_quasar));
        if ~isempty(station_quasar)
            for i_stat = 1 : length(station_sat)
                both_q = both_q | strcmp({station_quasar(:).name}, station_sat(i_stat).name);
            end
        end

        % in station_sat
        temp_station_name = {station_quasar(both_q).name};
        both_s = zeros(1,length(station_sat));
        for i_stat = 1 : length(station_quasar(both_q))
            both_s = both_s | strcmp(temp_station_name{i_stat}, {station_sat(:).name});
        end

        % ##### Delete all stations which observe satellites AND quasars from the "quasar" station network #####
        station_quasar = station_quasar(~both_q);

%     end
    
    
    % ##### Save data in "stat_data" structure #####
    
    % ## Satellit network ##
    number_of_stations_sat = length(station_sat);
    for i_stat = 1 : number_of_stations_sat
         
        try
            [lat,lon,h] = xyz2ell((station_sat(i_stat).trf_xyz)');
        catch error
            error_code = 1;
            error_msg = 'Calculation of ellipsoid station coordinates failed (xyz2ell.m)';  
            return;
         end

         try
             stat_data.stat(i_stat).label = station_sat(i_stat).po;
             stat_data.stat(i_stat).name = station_sat(i_stat).name;
             stat_data.stat(i_stat).location.ellipsoid.long = lon *rad2deg;
             stat_data.stat(i_stat).location.ellipsoid.lat = lat *rad2deg;
             stat_data.stat(i_stat).location.ellipsoid.altitude = h;
             stat_data.stat(i_stat).location.ellipsoid.ref_ellipsoid_label = '';                % Informaion currently not available
             stat_data.stat(i_stat).location.TRF.x = station_sat(i_stat).trf_xyz(1);
             stat_data.stat(i_stat).location.TRF.y = station_sat(i_stat).trf_xyz(2);
             stat_data.stat(i_stat).location.TRF.z = station_sat(i_stat).trf_xyz(3);
             stat_data.stat(i_stat).location.TRF.vx = station_sat(i_stat).trf_vel_xyz(1);
             stat_data.stat(i_stat).location.TRF.vy = station_sat(i_stat).trf_vel_xyz(2);
             stat_data.stat(i_stat).location.TRF.vz = station_sat(i_stat).trf_vel_xyz(3);
             stat_data.stat(i_stat).location.TRF.epoch = station_sat(i_stat).trf_epoch;
             stat_data.stat(i_stat).location.TRF.start = station_sat(i_stat).trf_start;
             stat_data.stat(i_stat).location.TRF.end = station_sat(i_stat).trf_end;
             stat_data.stat(i_stat).location.TRF.TRF_label = station_sat(i_stat).trf_source;
             stat_data.stat(i_stat).min_elevation = PARA.MIN_CUTEL * rad2deg;
             stat_data.stat(i_stat).max_axis2_rate = station_sat(i_stat).rate2;                 % axis 2 = Elevation (EL), DECLINATION (DC), (NS) [deg/sec]
             stat_data.stat(i_stat).max_axis1_rate = station_sat(i_stat).rate1;                 % axis 1 = Azimuth (AZ), Hour Angle (HA), (XY) [deg/sec]
             stat_data.stat(i_stat).max_axis2_acc = station_sat(i_stat).acc1;                   % If it is not defined in acceleration.cat the globally in param.txt is taken [deg/sec^2]
             stat_data.stat(i_stat).max_axis1_acc = station_sat(i_stat).acc2;                   % If it is not defined in acceleration.cat the globally in param.txt is taken [deg/sec^2]
             stat_data.stat(i_stat).lim11 = station_sat(i_stat).lim11;                          % Axis 1, limit 1 (Limit from antenna.cat incl. marge from param.txt)[rad]
             stat_data.stat(i_stat).lim12 = station_sat(i_stat).lim12;                          % Axis 1, limit 2 (Limit from antenna.cat incl. marge from param.txt)[rad]
             stat_data.stat(i_stat).lim21 = station_sat(i_stat).lim21;                          % Axis 2, limit 1 (Limit from antenna.cat incl. marge from param.txt)[rad]
             stat_data.stat(i_stat).lim22 = station_sat(i_stat).lim22;                          % Axis 2, limit 2 (Limit from antenna.cat incl. marge from param.txt)[rad]
             stat_data.stat(i_stat).c1 = station_sat(i_stat).c1;                                % Slew time constant [sec]
             stat_data.stat(i_stat).c2 = station_sat(i_stat).c2;                                % Slew time constant [sec]
             stat_data.stat(i_stat).horizontal_mask = station_sat(i_stat).hmask;                % Horizontal mask
             stat_data.stat(i_stat).horizontal_mask_num = station_sat(i_stat).hmasknum;         % Number of elements the horiziontal mask consists of
             stat_data.stat(i_stat).min_sun_dist =  PARA.MIN_SUNDIST*rad2deg;                   % Minimal separation angle between source and sun, set in the GUI
             stat_data.stat(i_stat).axis_type = station_sat(i_stat).axis;                       % Antenna axis type
             stat_data.stat(i_stat).obs_type = 'sat_only';
             stat_data.stat(i_stat).twin_partner = '';                                          % Is empty, if the antenna does not operate in twin/sibling mode
             stat_data.stat(i_stat).sefd_para = station_sat(i_stat).sefdpara;                   % SEFD parameters from "equip.cat"
             stat_data.stat(i_stat).min_snr = station_sat(i_stat).minsnr;                       % Minimum SNR for X and S band, set in the GUI
             stat_data.stat(i_stat).axis_offset = station_sat(i_stat).offset;                   % Axis offset [m]
             stat_data.stat(i_stat).az_n1 = station_sat(i_stat).azn1;                           % Cable wrap neutral pointing sector limit 1 [rad]
             stat_data.stat(i_stat).az_n2 = station_sat(i_stat).azn2;                           % Cable wrap neutral pointing sector limit 2 [rad]
             stat_data.stat(i_stat).az_cw1 = station_sat(i_stat).azc1;                          % Cable wrap clockwise pointing sector limit 1 [rad]
             stat_data.stat(i_stat).az_cw2 = station_sat(i_stat).azc2;                          % Cable wrap clockwise pointing sector limit 2 [rad]
             stat_data.stat(i_stat).az_ccw1 = station_sat(i_stat).azw1;                          % Cable wrap counter-clockwise pointing sector limit 1 [rad]
             stat_data.stat(i_stat).az_ccw2 = station_sat(i_stat).azw2;                          % Cable wrap counter-clockwise pointing sector limit 2 [rad]
             
             
             
             if PARA.USE_SEPARATE_STAT_NWW_FOR_SATS
                 if both_s(i_stat) == 1
                     stat_data.stat(i_stat).obs_type = 'sat_and_quasar';
                 else
                    stat_data.stat(i_stat).obs_type = 'sat_only';                                      
                 end
             else
                stat_data.stat(i_stat).obs_type = 'sat_and_quasar';
             end
             stat_data.stat(i_stat).twin_partner = '';                                          % Is empty, if the antenna does not operate in twin/sibling mode
         catch error
            error_code = 1;
            error_msg = 'Station Parameter not available';  
            return;
         end
    end
    
    stat_data.number_of_stations_quasar = length(station_quasar);
    stat_data.number_of_stations = stat_data.number_of_stations_quasar + number_of_stations_sat;
    
    % Only if separate staton networks for satellites/quasars are defined
    if PARA.USE_SEPARATE_STAT_NWW_FOR_SATS
        
        % ## Quasar network ##
        

        for i_stat = 1 : stat_data.number_of_stations_quasar

            i_stat_2 = i_stat + number_of_stations_sat;

            try
                [lat,lon,h] = xyz2ell((station_quasar(i_stat).trf_xyz)');
            catch error
                error_code = 1;
                error_msg = 'Calculation of ellipsoid station coordinates failed (xyz2ell.m)';  
                return;
             end

             try
                 stat_data.stat(i_stat_2).label = station_quasar(i_stat).po;
                 stat_data.stat(i_stat_2).name = station_quasar(i_stat).name;
                 stat_data.stat(i_stat_2).location.ellipsoid.long = lon *rad2deg;
                 stat_data.stat(i_stat_2).location.ellipsoid.lat = lat *rad2deg;
                 stat_data.stat(i_stat_2).location.ellipsoid.altitude = h;
                 stat_data.stat(i_stat_2).location.ellipsoid.ref_ellipsoid_label = '';                     % Informaion currently not available
                 stat_data.stat(i_stat_2).location.TRF.x = station_quasar(i_stat).trf_xyz(1);
                 stat_data.stat(i_stat_2).location.TRF.y = station_quasar(i_stat).trf_xyz(2);
                 stat_data.stat(i_stat_2).location.TRF.z = station_quasar(i_stat).trf_xyz(3);
                 stat_data.stat(i_stat_2).location.TRF.vx = station_quasar(i_stat).trf_vel_xyz(1);
                 stat_data.stat(i_stat_2).location.TRF.vy = station_quasar(i_stat).trf_vel_xyz(2);
                 stat_data.stat(i_stat_2).location.TRF.vz = station_quasar(i_stat).trf_vel_xyz(3);
                 stat_data.stat(i_stat_2).location.TRF.epoch = station_quasar(i_stat).trf_epoch;
                 stat_data.stat(i_stat_2).location.TRF.start = station_quasar(i_stat).trf_start;
                 stat_data.stat(i_stat_2).location.TRF.end = station_quasar(i_stat).trf_end;
                 stat_data.stat(i_stat_2).location.TRF.TRF_label = station_quasar(i_stat).trf_source;
                 stat_data.stat(i_stat_2).min_elevation = PARA.MIN_CUTEL * rad2deg;
                 stat_data.stat(i_stat_2).max_axis2_rate = station_quasar(i_stat).rate2;                    % axis 2 = Elevation (EL), DECLINATION (DC), (NS) [deg/sec]
                 stat_data.stat(i_stat_2).max_axis1_rate = station_quasar(i_stat).rate1;                    % axis 1 = Azimuth (AZ), Hour Angle (HA), (XY) [deg/sec]
                 stat_data.stat(i_stat_2).max_axis2_acc = station_quasar(i_stat).acc1;                      % If it is not defined in acceleration.cat the globally in param.txt is taken [deg/sec^2]
                 stat_data.stat(i_stat_2).max_axis1_acc = station_quasar(i_stat).acc2;                      % If it is not defined in acceleration.cat the globally in param.txt is taken [deg/sec^2]
                 stat_data.stat(i_stat_2).lim11 = station_quasar(i_stat).lim11;                             % Axis 1, limit 1 (Limit from antenna.cat incl. marge from param.txt)[rad]
                 stat_data.stat(i_stat_2).lim12 = station_quasar(i_stat).lim12;                             % Axis 1, limit 2 (Limit from antenna.cat incl. marge from param.txt)[rad]
                 stat_data.stat(i_stat_2).lim21 = station_quasar(i_stat).lim21;                             % Axis 2, limit 1 (Limit from antenna.cat incl. marge from param.txt)[rad]
                 stat_data.stat(i_stat_2).lim22 = station_quasar(i_stat).lim22;                             % Axis 2, limit 2 (Limit from antenna.cat incl. marge from param.txt)[rad]
                 stat_data.stat(i_stat_2).c1 = station_quasar(i_stat).c1;                                   % Slew time constant [sec]
                 stat_data.stat(i_stat_2).c2 = station_quasar(i_stat).c2;                                   % Slew time constant [sec]
                 stat_data.stat(i_stat_2).horizontal_mask = station_quasar(i_stat).hmask;                   % Horizontal mask
                 stat_data.stat(i_stat_2).horizontal_mask_num = station_quasar(i_stat).hmasknum;            % Number of elements the horiziontal mask consists of
                 stat_data.stat(i_stat_2).min_sun_dist =  PARA.MIN_SUNDIST*rad2deg;                         % Minimal separation angle between source and sun, set in the GUI
                 stat_data.stat(i_stat_2).axis_type =  station_quasar(i_stat).axis;                         % Antenna axis type
                 stat_data.stat(i_stat_2).obs_type = 'quasar_only';
                 stat_data.stat(i_stat_2).twin_partner = '';                                                % Is empty, if the antenna does not operate in twin/sibling mode
                 stat_data.stat(i_stat_2).sefd_para =  station_quasar(i_stat).sefdpara;                     % SEFD parameters from "equip.cat"
                 stat_data.stat(i_stat_2).min_snr =  station_quasar(i_stat).minsnr;                         % Minimum SNR for X and S band, set in the GUI
                 stat_data.stat(i_stat_2).axis_offset = station_quasar(i_stat).offset;                         % Axis offset [m]
                 stat_data.stat(i_stat_2).az_n1 = station_quasar(i_stat).azn1;                                 % Cable wrap neutral pointing sector limit 1 [rad]
                 stat_data.stat(i_stat_2).az_n2 = station_quasar(i_stat).azn2;                                 % Cable wrap neutral pointing sector limit 2 [rad]
                 stat_data.stat(i_stat_2).az_cw1 = station_quasar(i_stat).azc1;                                % Cable wrap clockwise pointing sector limit 1 [rad]
                 stat_data.stat(i_stat_2).az_cw2 = station_quasar(i_stat).azc2;                                % Cable wrap clockwise pointing sector limit 2 [rad]
                 stat_data.stat(i_stat_2).az_ccw1 = station_quasar(i_stat).azw1;                               % Cable wrap counter-clockwise pointing sector limit 1 [rad]
                 stat_data.stat(i_stat_2).az_ccw2 = station_quasar(i_stat).azw2;                               % Cable wrap counter-clockwise pointing sector limit 2 [rad]


             catch error
                error_code = 1;
                error_msg = 'Station Parameter not available';  
                return;
             end
        end
    

        % ##### Find twin/sibling telescopes #####
   
    
        % Get info from file (twin.txt):
        fid = fopen(INFILE.twin, 'r');
        if (fid == -1)
             error_code = 1; 
             error_msg = ['Can not open the file: ',INFILE.twin];
            return;
        end
        twin_data = textscan(fid, '%s %s', 'CommentStyle', '*');
        fclose(fid);

        number_of_twins = length(twin_data{1,1});
        twin_q = zeros(1,length(station_quasar));
        twin_s = zeros(1,length(station_sat));

        for i_twin = 1 : number_of_twins % loop over all twin entries

            for i_col = 1 : 2 % Search in first and second column

                twin_q = strcmp(twin_data{1,i_col}{i_twin}, {station_quasar(:).name});

                if sum(twin_q == 1) % Station was found in twin data (first column)

                    % Compare second column of twin_data with satellite station network: 
                    twin_s  = strcmp(twin_data{1, (3-i_col)}{i_twin}, {station_sat(:).name});

                    if sum(twin_s == 1)
                        % Exclude stations, which are included in both, the quasar and the satellite network: 
                        temp_s = twin_s | both_s;
                        twin_s = twin_s & temp_s;

                        if sum(twin_q == 1)
                            % ########## Found twins #####################

                            % Set observation type and "twin pratner" in output structure:
                            stat_data.stat(logical([twin_s, zeros(1, stat_data.number_of_stations_quasar)])).obs_type         = 'twin_sat';
                            stat_data.stat(logical([twin_s, zeros(1, stat_data.number_of_stations_quasar)])).twin_partner     = stat_data.stat(logical([zeros(1, number_of_stations_sat), twin_q])).name;

                            stat_data.stat(logical([zeros(1, number_of_stations_sat), twin_q])).obs_type      = 'twin_quasar';
                            stat_data.stat(logical([zeros(1, number_of_stations_sat), twin_q])).twin_partner  = stat_data.stat(logical([twin_s, zeros(1, stat_data.number_of_stations_quasar)])).name;

                            % ###############################
                        end

                    elseif sum(twin_s > 1) % Error, is station is included more often than once in twin.txt:
                        error_code = 3; 
                        error_msg = ['A Station is included in more than one twon/sibling telescopes in the file: ',INFILE.twin];
                        return;
                    end

                elseif sum(twin_q > 1) % Error, is station is included more often than once in twin.txt:
                    error_code = 2; 
                    error_msg = ['A Station is included in more than one twon/sibling telescopes in the file: ',INFILE.twin];
                    return;
                end  
            end 
        end
        
    end
        
return
         
    

