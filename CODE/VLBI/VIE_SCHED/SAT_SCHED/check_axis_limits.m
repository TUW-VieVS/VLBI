% #########################################################################
% #     check_axis_limits
% #########################################################################
%
% DESCRIPTION
%   This function checks if the antenna axes limits are kept at the epoch defined
%   at t_epoch_jd at one station, i.e. if the antenna "slews through" the axis 
%   limit during a scan.
% 
%   Works if the angular distance in azimuth of source positions between start (t_start_jd) and end time (t_epoch_jd) is not too large!  
%       => Set "treshold_angle_rad" reasonably! ...this parameter is required to recognise the case, when the antenna moves across the 0/360 deg mark in azimuth!
%                                               ...this is critical for calculating un_az!
%              => This issue is addressed in check_t_end.m, where this function is called from: 
%                    - Calc un_az step by step in the PARA.CHECK_SCAN_INT interval
%                    - Take the un_az value of the previous step as start value for the next step, a.s.o..
%
%
%   The ambiguous antenna azimuth (un_az) is calculated and tracked rigorously!
%
%   To specify an satellite observation:
%       - satellite_id...is not empty
%       - source_quasar...is empty
%
%   To specify an quasar observation:
%       - satellite_id...is empty
%       - source_quasar...is not empty
%
%
% CREATED  
%   2015-04-07     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
% - tle_2_topo
% - zazel_s
% - azel2xyew
%
%
% INPUT
% - un_az_t_start_jd    - Unambiguous antenna-azimuthat the scan start time (t_start_jd) [rad]
% - t_epoch_jd          - Time epoch for the calculations [JD]
% - stat_data           - station data structure
% - station_id          - ID of the observing stations (referring to the required entry of "stat_data.stat" and "obs_data.stat")
% - obs_data            - Observation data structure (preallocated in "find_observation_times.m")
% - PARA                - Global scheduling parameter strucutre
% - source_quasar       - source (quasars) data structure (ONE sub-set/field of the "source" structure, containing the required data)
% - satellite_id        - ID of the satellite to be observed (referring to an entry of "stat_data.stat(i_stat).sat(i_sat)" or "obs_data.sat(i_sat)")
% - t_start_jd          - Scan start time [JD]
%
%
% OUTPUT
% - error_code          - Error Code (0 = no erros occured)
% - error_msg           - Error Message (empty, if no errors occured)
% - flag_limits_ok      - Flag: 1 => axis 1 limts are kept; 0 => axis 1 limts are NOT kept
% - un_az               - Unambiguous azimuth at t_epoch_jd [rad]
% - el                  - Elevation angle at t_epoch_jd [rad]
%
% CHANGES:
% - 2015-06-29: A. Hellerschmied: azel2xyew.m used to calc. X/Y EW angles
% - 2016-06-29: A. Hellerschmied: - Calculation of un_az corrected. un_az was wrong, if the source crossed the 0/360 deg azimuth during an observation! 
%                                   => Now it should work at least, if the satellite does not cross the whole sky during an observation => See comment in the function description above!
%                                 - Option for debug output added (PARA.DEBUG_FLAG) 
% - 2016-09-30: A. Hellerschmied: Bug fixed: un_az was not assigned for HADC and XYEW antennas.
%                                              
%

function [flag_limits_ok, un_az, el, error_code, error_msg] = check_axis_limits(stat_data, station_id, PARA, source_quasar, satellite_id, t_epoch_jd, un_az_t_start_jd, t_start_jd)
 
    % ##### Init #####
    error_code = 0;
    error_msg = '';
    flag_limits_ok = 0;
    
    % ##### Options #####
    treshold_angle_rad = pi;

    % ##########################
    % ##### Satellite scan #####
    % ##########################
    if ~isempty(satellite_id) && isempty(source_quasar)
        
        if (strncmp(stat_data.stat(station_id).axis_type, 'AZEL', 4))
            
            % ##### To calculate the differenze between azimuth at t_epoch_jd and the azimuth at the scan start time (t_start_jd) [rad] #####
            [az_t_start, el_t_start, range, az_rate, el_rate, r_rate, tra, dc, tra_rate, tdec_rate, ha, local_hour_angle_rate, error_code, error_msg] = tle_2_topo(stat_data, station_id, PARA, satellite_id, t_start_jd);
            if error_code ~= 0
                error_msg = ['tle_2_topo:', error_msg];
                return;
            end
            az_t_start = az_t_start * pi/180;
%             delta_az = un_az_t_start_jd - az;
        end
        
        
        
        % ##### Calc. axis 1 angle at t_epoch_jd #####
        [az, el, range, az_rate, el_rate, r_rate, tra, dc, tra_rate, tdec_rate, ha, local_hour_angle_rate, error_code, error_msg] = tle_2_topo(stat_data, station_id, PARA, satellite_id, t_epoch_jd);
        if error_code ~= 0
            error_msg = ['tle_2_topo:', error_msg];
            return;
        end
        az = az * pi/180;
        el = el * pi/180;
        ha = ha * pi/180;
        dc = dc * pi/180;
        

    % ##########################
    % ##### Quasar scan    #####
    % ##########################
    elseif isempty(satellite_id) && ~isempty(source_quasar)

        lon = stat_data.stat(station_id).location.ellipsoid.long * pi / 180;
        lat = stat_data.stat(station_id).location.ellipsoid.lat * pi / 180;
        
        if (strncmp(stat_data.stat(station_id).axis_type, 'AZEL', 4))
            % ##### To calculate the differenze between azimuth at t_epoch_jd and the azimuth at the scan start time (t_start_jd) [rad] #####
            [az_t_start, el, ha, dc] = zazel_s(t_start_jd - 2400000.5, lon, lat, source_quasar.ra, source_quasar.de);
%             delta_az = un_az_t_start_jd - az;
        end
        
        
        % ##### Calc. axis 1 angle at t_epoch_jd #####
        [az, el, ha, dc] = zazel_s(t_epoch_jd - 2400000.5, lon, lat, source_quasar.ra, source_quasar.de);

    end
    
    
    % ##### Check axis limits ######
    
    
    % axis limit
    if (strncmp(stat_data.stat(station_id).axis_type, 'AZEL', 4))
        
        flag_test = 0;
        % Calc unambiguous azimuth at t_epoch_jd:
        d_az = az - az_t_start;                 % Difference in az between teh epochs t_epoch_jd and t_start_jd
        if (d_az > treshold_angle_rad)          % Antenna azimuth moved across 0/360 deg, ccw direction
            d_az = d_az - 2*pi;
            flag_test = 1;
        elseif (d_az < (treshold_angle_rad)*-1) % Antenna azimuth moved across 0/360 deg, ccw direction
            d_az = d_az + 2*pi;
            flag_test = 2;
        end

        un_az = un_az_t_start_jd + d_az;
        
        if PARA.DEBUG_FLAG
            disp(['############### ',num2str(station_id), ' ##########################']);
            disp(['   => un_az in:  ', num2str(un_az_t_start_jd*180/pi)]);
            disp(['   => un_az new: ', num2str(un_az*180/pi)]);
            disp(['   => az in:     ', num2str(az_t_start*180/pi)]);
            disp(['   => az new:    ', num2str(az*180/pi)]);
            disp(['   => d_az:      ', num2str(d_az*180/pi)]);
			disp(['   => dt [sec]:  ', num2str((t_epoch_jd - t_start_jd)*3600*24)]);
            if flag_test == 1
                disp('   => CCW over 0/360 deg !!!!!!');
            elseif flag_test == 2
                disp('   => CW over 0/360 deg !!!!!!');
            end
            disp('############################################');
        end
    
        if ((un_az > stat_data.stat(station_id).lim11) & (un_az < stat_data.stat(station_id).lim12))
            lupaz = true;
        else
            lupaz = false;
        end
        if ((el > stat_data.stat(station_id).lim21) & (el < stat_data.stat(station_id).lim22))
            lupel = true;
        else 
            lupel = false;
        end
        lup = lupaz & lupel;
        
    elseif (strncmp(stat_data.stat(station_id).axis_type, 'HADC', 4))
        if ((ha > stat_data.stat(station_id).lim11) & (ha < stat_data.stat(station_id).lim12) & (dc > stat_data.stat(station_id).lim21) & (dc < stat_data.stat(station_id).lim22))
            lup = true;
        else 
            lup = false;
        end
        un_az = az;
        
    elseif (strncmp(stat_data.stat(station_id).axis_type, 'XYEW', 4))
        [x30, y30, error_code, error_msg] = azel2xyew(az, el);
        if error_code ~= 0
            error_msg = ['azel2xyew:', error_msg];
            return;
        end
        
        if ((x30 > stat_data.stat(station_id).lim11) & (x30 < stat_data.stat(station_id).lim12) & (y30 > stat_data.stat(station_id).lim21) & (y30 < stat_data.stat(station_id).lim22))
            lup = true;
        else 
            lup = false;
        end
        un_az = az;
    end   

    flag_limits_ok = lup;
    
return;
