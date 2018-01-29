% #########################################################################
% #     calc_start_of_obs
% #########################################################################
%
% DESCRIPTION
%   Calculation the earliest possible start time of an observation for one single station, 
%   taking into account the previous antenna postion and the antenna slew time. The 
%   calculation is done iteratively (the start time is determined iteratively => The final antenna poiting angles change, etc...).
%   Observations to quasars and to satellites are considered.
%   Satellite positions are calculated with the SGP4 model + TLE data.
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
%   2015-03-17     Andreas Hellerschmied
%
% REFERENCES
% - Jing Sun (2013), VLBI Scheduling strategies with respect to VLBI2010,
%   Geowissenschaftliche Mitteilungen, Heft Nr. 92, ISSN 1811-8380.
%
%
% COUPLING
% - calc_slew_time
% - tle_propagation
% - eci2topo
% - zazel_s
%
%
% INPUT
% - stat_data           - station data structure
% - station_id          - ID of the observing station (referring to the required entry of "stat_data.stat" and "obs_data.stat")
% - obs_data            - Observation data structure (preallocated in "find_observation_times.m")
% - PARA                - Global scheduling parameter strucutre
% - source_quasar       - source (quasars) data structure (ONE sub-set/field of the "source" structure, containing the required data)
% - satellite_id        - ID of the satellite to be observed (referring to an entry of "stat_data.stat.sat" )
% - t_min_jd            - Earliest possible observation start time. Used, if there hasn't been any observations before (optional)
%
%
% OUTPUT
% - error_code          - Error Code (0 = no erros occured)
% - error_msg           - Error Message (empty, if no errors occured)
% - start_jd            - Earliest possible start time of the observation from the current station [JD]
% - slew_time           - Antenna slew time of the current station [sec]
% - obs_data            - Observation data structure (preallocated in "find_observation_times.m") - unambiguous NOW source azimuth [rad] updated
%
% CHANGES:
% - 2015-07-08, A. Hellerschmied: Possibility to add preob time constant (PARA.PREOB_NEW_SOURCE) after slewing/before next scan added.
% - 2015-07-20: A. Hellerschmied: New TLE/SGP4 orbit propagation function used (tle_propagation.m)
% - 2015-07-20: A. Hellerschmied: Bug-fix
% - 2015-11-18: A. Hellerschmied: Bug-fix: "delta_start_time_treshold" increased fom "1" to "1.01". With "1.0" the iteration loop sometimes was a dead end..
% - 2015-11-19: A. Hellerschmied: Check for the max. number of iterations for the calcualtion of the obs start time added. In case, "error_code" is set to "5"
% - 2016-05-11: A. Hellerschmied: Input an earliest possible observation start time (optional)
% - 2016-08-02: A. Hellerschmied: Max. number of iteration ("max_num_of_iterations") increased from 5 to 10.
% - 2016-11-25: A. Hellerschmied: Option to set PARA.SOURCE, PARA.TAPETM, PARA.IDLE, PARA.CALIBRATION to zero (PARA.DISABLE_SLEW_CONST)

function [start_jd, slew_time_sec, obs_data, error_code, error_msg] = calc_start_of_obs(stat_data, station_id, obs_data, PARA, source_quasar, satellite_id, varargin)

    % ##### Init #####
    error_code = 0;
    error_msg = '';
    preob_time_const_sec = 0;
    
    PARA_tmp = PARA;
    
    if isfield(PARA_tmp, 'DISABLE_SLEW_CONST')
        if PARA_tmp.DISABLE_SLEW_CONST
            PARA_tmp.SOURCE         = 0;
            PARA_tmp.TAPETM         = 0;
            PARA_tmp.IDLE           = 0;
            PARA_tmp.CALIBRATION    = 0;
        end
    end

    
    % ##### Handle input #####
    switch(nargin)
        case 6 % Without optional input
            t_min_jd = [];
        case 7 % t_min_jd
            t_min_jd = varargin{1};
        otherwise
            error_code = 1;
            error_msg = 'Invalid number of input arguments!';
            return;
    end
    
    % Iteration setup:
    delta_start_time_treshold = 1.01; % [sec] % Do not choose this value smaller than 1 sec, because the slew time is calculated with a resolution of 1 sec (calc_slew_time)!
    delta_start_time = delta_start_time_treshold + 1; % Init. for first iteration (to enter the iteration loop)
    iteration_counter = 0;
    max_num_of_iterations = 10; % Maximim number of iterations => After "max_num_of_iterations" iterations "flag_iteration_error" is set!
    
    % Observation start epoch/position initialisation:
    if obs_data.number_of_scans == 0    % first scan
        % For the first scan, the slew time is not critical, because the
        % session can be started a sufficient amount of time before the
        % first observation is scheduled at the stations.
        % ==> The observation start time (start_jd) is set to the session
        % start time in the GUI!
        if isempty(t_min_jd)
            start_jd = PARA_tmp.startmjd + 2400000.5;
        else
            start_jd = t_min_jd;
        end
        start_jd = ceil(start_jd * 86400) / 86400; % round up to integer seconds...
        slew_time_sec = 0;
        return;
    else
        start_jd        = obs_data.end_of_last_scan_jd;  % Init. the calculation epoch for iteration
        start_jd_new    = 0;                             % Required to calculate the difference in start time between 2 iteration steps 
    end
    
    
    % ##### Get source data #####
    if isempty(satellite_id) && ~isempty(source_quasar) % Quasar observation
        % ##########################
        % ##### Quasar scan    #####
        % ##########################
        obs_type = 1;
        % Get Ra/Dec for observed source
        ra = source_quasar.ra; 
        de = source_quasar.de;
        
        % Get preob time constant, which is added after slewing:
        %  => Check, if the station participated in the last scan (add preob time only in this case!)
        if obs_data.stat(station_id).end_of_last_obs.jd == obs_data.end_of_last_scan_jd
            % => Check, if the station oberved in this session before:
            if isfield(obs_data.stat(station_id).end_of_last_obs, 'sat_id') && isfield(obs_data.stat(station_id).end_of_last_obs, 'quasar_id')
                %  => Add, if a satellite was observed by this station in the last scan:
                if ~isempty(obs_data.stat(station_id).end_of_last_obs.sat_id) && isempty(obs_data.stat(station_id).end_of_last_obs.quasar_id)
                    preob_time_const_sec = PARA_tmp.PREOB_NEW_SOURCE;
                end
            end
        end
        
    elseif ~isempty(satellite_id) && isempty(source_quasar) % Satellite observation
        % ##########################
        % ##### Satellite scan #####
        % ##########################
        obs_type = 2;
        
        % Get preob time constant, which is added after slewing:
        %  => Check, if the station participated in the last scan (add preob time only in this case!)
        if obs_data.stat(station_id).end_of_last_obs.jd == obs_data.end_of_last_scan_jd
            % => Check, if the station oberved in this session before:
            if isfield(obs_data.stat(station_id).end_of_last_obs, 'sat_id') && isfield(obs_data.stat(station_id).end_of_last_obs, 'quasar_id')
                %  => Add, if a different satellite or any quasar was observed by this station in the last scan:
                if isempty(obs_data.stat(station_id).end_of_last_obs.sat_id) && ~isempty(obs_data.stat(station_id).end_of_last_obs.quasar_id) % last obs.: quasar
                    preob_time_const_sec = PARA_tmp.PREOB_NEW_SOURCE;
                elseif ~isempty(obs_data.stat(station_id).end_of_last_obs.sat_id) && isempty(obs_data.stat(station_id).end_of_last_obs.quasar_id) % last obs.: satellite
                    % check, if a different satellite was observed before:
                    if satellite_id ~= obs_data.stat(station_id).end_of_last_obs.sat_id
                        preob_time_const_sec = PARA_tmp.PREOB_NEW_SOURCE;
                    end
                end
            end
        end
        

        sat_name_to_prop    = stat_data.stat(station_id).sat(satellite_id).TLE_data.TLE_header_line;
        delta_t             = stat_data.stat(station_id).prop_setup.delta_t_min; 

    else % Error case
        error_code = 1;
        error_msg = 'Source data error.';
        return;
    end
    
    
    % ##### Get station data #####
    stat_lon = stat_data.stat(station_id).location.ellipsoid.long;
    stat_lat = stat_data.stat(station_id).location.ellipsoid.lat;
    stat_alt = stat_data.stat(station_id).location.ellipsoid.altitude;
    
    
    % ##### Iteration loop #####
    while(delta_start_time > delta_start_time_treshold) 

        iteration_counter = iteration_counter + 1;
        %         fprintf('---------- calc_start_of_obs: iteration_counter: %d -------------- \n', iteration_counter);  
        
        % #### Check the max number of allowed iterations ####
        if iteration_counter > max_num_of_iterations
            % Error init.:
            start_jd = 0;
            slew_time_sec = 0;
            error_code = 5;
            error_msg = ['Max. number of iterations (',num2str(max_num_of_iterations),') reached. Station: ', stat_data.stat(station_id).name];
            return; % while(delta_start_time > delta_start_time_treshold)
        end
        
        % ##### Calculation of the source position in a topocentric system (az/el, ha/dc) #####
        
        % ##########################
        % ##### Quasar scan    #####
        % ##########################
        if obs_type == 1
            % Calc. Az/El and Ha/Dec of current source at given jd
            [az, el, ha, dc] = zazel_s(start_jd - 2400000.5, stat_lon * pi/180, stat_lat * pi/180, ra, de); 
            
        % ##########################
        % ##### Satellite scan #####
        % ##########################
        elseif obs_type == 2
            
            % ### calculate antenna pointing data (az, el) for the given epoch (jd_current_epoch_marker) [JD] ###
            
            % SGP4 orbit propagation:
%             [temp_sat_data] = TLE_propagation(path_tle, filename_tle, start_jd, start_jd, delta_t, grav_const, write_file, path_out, filename_out, verification_mode, sat_name_to_prop);
            [temp_sat_data, error_code, error_msg] = tle_propagation(start_jd, start_jd, delta_t, PARA_tmp, sat_name_to_prop);
            if error_code > 0
                error_msg = ['tle_propagation:', error_msg];
                return; 
            end

            % Calculation of topocentric view angles:
            [az, el, range_s_obs, az_rate, el_rate, r_rate, tra, dc, tra_rate, tdec_rate, ha, local_hour_angle_rate] = eci2topo(start_jd, temp_sat_data.sat(1).prop(1).r, temp_sat_data.sat(1).prop(1).v, stat_lon, stat_lat, stat_alt);
            az = az * pi/180;
            el = el * pi/180;
            ha = ha * pi/180;
            dc = dc * pi/180;
        end
        

        % #### Calculation of the antenna slew time ####
        [slew_time_sec, un_az] = calc_slew_time(stat_data.stat(station_id), obs_data.stat(station_id).end_of_last_obs, az, el, ha, dc);
%         fprintf('st: %5.1f; az: %7.2f, un_az: %7.2f\n', slew_time_sec, az*180/pi, un_az*180/pi);
        

        % #### Calculation observation start time ####
        start_jd_new = obs_data.stat(station_id).end_of_last_obs.jd + slew_time_sec/86400.0 + PARA_tmp.SOURCE/86400 + PARA_tmp.TAPETM/86400 + PARA_tmp.IDLE/86400 + PARA_tmp.CALIBRATION/86400 + preob_time_const_sec/86400; % Calc observation start time [JD]
        delta_start_time = abs(start_jd_new - start_jd) * 86400;
%         fprintf(' delta_start_time = %f, it = %d\n', delta_start_time, iteration_counter);
        start_jd = ceil(start_jd_new * 86400) / 86400; % round up to integer seconds...
%         fprintf('  => start_jd: %s; dt: %5.1f sec\n', jd2datestr(start_jd), delta_start_time);
        
    end
    
return;
