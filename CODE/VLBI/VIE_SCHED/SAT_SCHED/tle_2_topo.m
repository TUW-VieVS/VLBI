% #########################################################################
% #     tle_2_topo
% #########################################################################
%
% DESCRIPTION
%   This function uses TLE datasets and the SGP4 propagator model to calculate 
%   topocentric satellite coordiante and the according change rates.
%
%
% CREATED  
%   2015-04-07     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
% - tle_propagation
% - eci2topo
%
%
% INPUT
% - stat_data           - station data structure
% - station_id          - ID of the observing station (referring to the required entry of "stat_data.stat" and "obs_data.stat")
% - PARA                - Global scheduling parameter strucutre
% - satellite_id        - ID of the satellite to be observed (referring to an entry of "stat_data.stat.sat" )
% - t_epoch_jd          - Time epoch for the calculations
%
%
% OUTPUT
%   - error_code                - Error Code (0 = no erros occured)
%   - error_msg                 - Error Message (empty, if no errors occured)
%   - az                        - Azimuth (Station -> Source) [deg]
%   - el                        - Elevation (Station -> Source) [deg]
%   - range                     - Range (Station -> Source) [m]
%   - az_rate                   - Azimuth rate (Station -> Source) [deg/sec]
%   - el_rate                   - Elevation rate (Station -> Source) [deg/sec]
%   - r_rate                    - Range rate (Station -> Source) [m/sec]
%   - tra                       - Topocentric right ascension [deg]
%   - tdec                      - Topocentric declination [deg]
%   - tra_rate                  - Topocentric right ascension rate [deg/sec]
%   - tdec_rate                 - Topocentric declination rate [deg/sec]
%   - local_hour_angle          - Local Hour Angle [deg]; [180 deg, +180 deg]
%   - local_hour_angle_rate     - Local Hour Angle change rate [deg/sec]
%
% CHANGES:
%   - 2015-07-20: A. Hellerschmied: New TLE/SGP4 orbit propagation function used (tle_propagation.m)
%   - 2016-12-22, A. Hellerschmied: - PARA.INIT_PROP_INTERVAL in [sec] instead of [min]
%

function [az, el, range, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate, error_code, error_msg] = tle_2_topo(stat_data, station_id, PARA, satellite_id, t_epoch_jd)
    
    % ##### Init #####
    error_code = 0;
    error_msg = '';
    
    % General initialization for SGP4 orbit propagation:    
%     grav_const          = 72;
%     write_file          = 0;
%     path_out            = 'dummy';
%     filename_out        = 'dummy';
%     verification_mode   = 0;
%     path_tle            = PARA.TLE_FILEPATH;
%     filename_tle        = PARA.TLE_FILENAME;
    delta_t             = PARA.INIT_PROP_INTERVAL/60;
    
    sat_name_to_prop    = stat_data.stat(station_id).sat(satellite_id).TLE_data.TLE_header_line;

    
    % Get station data
    stat_lon = stat_data.stat(station_id).location.ellipsoid.long;
    stat_lat = stat_data.stat(station_id).location.ellipsoid.lat;
    stat_alt = stat_data.stat(station_id).location.ellipsoid.altitude;

    % SGP4 orbit propagation:
%     [temp_sat_data] = TLE_propagation(path_tle, filename_tle, t_epoch_jd, t_epoch_jd, delta_t, grav_const, write_file, path_out, filename_out, verification_mode, sat_name_to_prop);
    [temp_sat_data, error_code, error_msg] = tle_propagation(t_epoch_jd, t_epoch_jd, delta_t, PARA, sat_name_to_prop);
    if error_code > 0
        error_msg = ['tle_propagation:', error_msg];
        return; 
    end

    % Calculation of topocentric coordinates and change rates:
    [az, el, range, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate] = eci2topo(t_epoch_jd, temp_sat_data.sat(1).prop(1).r, temp_sat_data.sat(1).prop(1).v, stat_lon, stat_lat, stat_alt);

    
return;
