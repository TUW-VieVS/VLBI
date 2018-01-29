% #########################################################################
% #     save_scan
% #########################################################################
%
% DESCRIPTION
%   This function saves the new scan to the required structures:
%    - sched_data
%    - obs_data
%
%   To specify an satellite observation:
%       - satellite_id...is not empty
%       - quasar_id...is empty
%
%   To specify an quasar observation:
%       - satellite_id...is empty
%       - quasar_id...is not empty
%
% CREATED  
%   2015-03-31     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
% - zazel_s
% - tle_2_topo
% - check_scan_duration_and_calc_repos_epochs
%
%
% INPUT
% - stat_data           - station data structure
% - station_id_list     - List of IDs of the observing stations (referring to the required entry of "stat_data.stat" and "obs_data.stat")
% - obs_data            - Observation data structure (preallocated in "find_observation_times.m")
% - sched_data          - Scheduling data structure
% - source              - source (quasars) data structure
% - quasar_id           - ID of the observed quasar (ref. to entry of "source" structure)
% - satellite_id        - ID of the satellite to be observed (referring to an entry of "stat_data.stat(i_stat).sat(i_sat)" or "obs_data.sat(i_sat)")
% - PARA                - Global scheduling parameter structure
%
%
% OUTPUT
% - error_code          - Error Code (0 = no erros occured)
% - error_msg           - Error Message (empty, if no errors occured)
% - obs_data            - Observation data structure (preallocated in "find_observation_times.m")
% - sched_data          - Scheduling data structure
%
% CHANGES:
% - 2015-11-06, A. Hellerschmied: Minor bug-fix.
% - 2016-04-28, A. Hellerschmied: Save updated "scan_data" and "obs_data" structure to the LEVEL5 sub-dir.
% - 2016-05-11, A. Hellerschmied: Sky-coverage is now saved in obs_data struct, bug-fixes
% - 2016-10-28, A. Hellerschmied: First make all checks, than assign the output structures!

function [sched_data, obs_data, error_code, error_msg] = save_scan(stat_data, station_id_list, sched_data, obs_data, t_start_jd, t_end_jd, source, quasar_id, satellite_id, PARA)
 
    % ##### Init #####
    error_code = 0;
    error_msg = '';    
    
    obs_data_tmp    = obs_data;
    sched_data_tmp  = sched_data;
    
    % ##### Save data to "obs_data" structure #####
    
    obs_data_tmp.number_of_scans = obs_data_tmp.number_of_scans + 1;
    obs_data_tmp.end_of_last_scan_jd = t_end_jd;
    
    

    % ##### Save data to "sched_data" structure #####
    
    sched_data_tmp.number_of_scans = sched_data_tmp.number_of_scans + 1;
    i_scan = sched_data_tmp.number_of_scans;
    
    % Nominal session starte/end time:
    if isempty(sched_data_tmp.t_nominal_start_jd)
        sched_data_tmp.t_nominal_start_jd = t_start_jd;
    elseif sched_data_tmp.t_nominal_start_jd > t_start_jd
        sched_data_tmp.t_nominal_start_jd = t_start_jd;
    end
    if isempty(sched_data_tmp.t_nominal_end_jd)
        sched_data_tmp.t_nominal_end_jd = t_end_jd;
    elseif sched_data_tmp.t_nominal_end_jd < t_end_jd
        sched_data_tmp.t_nominal_end_jd = t_end_jd;
    end

    sched_data_tmp.scan(i_scan).number_of_stat = length(station_id_list);
    sched_data_tmp.scan(i_scan).repos_int_sec = PARA.REPOS_INT;
    sched_data_tmp.scan(i_scan).flag_record_scan = 1;                       % !!!! Hard coded here! If it is required, change that later on.... 
    sched_data_tmp.scan(i_scan).t_start_jd = t_start_jd;
    sched_data_tmp.scan(i_scan).t_end_jd = t_end_jd;
    
    
    

    % ##### Observation type dependend parameters #####
    
    % ##########################
    % ##### Satellite scan #####
    % ##########################
    if ~isempty(satellite_id) && isempty(quasar_id)

        
        % #### Save data to "sched_data" structure ####
        
        % Preallocate "epoch" data:
        [flag_scan_duration_ok, alt_t_end_jd, repos_epochs_jd, error_code, error_msg] = check_scan_duration_and_calc_repos_epochs(PARA, t_start_jd, t_end_jd);
        if error_code > 0
            error_msg = ['check_scan_duration_and_calc_repos_epochs: ', error_msg];
            return;
        else
            if ~flag_scan_duration_ok % the scan duration is not a multiple of the antenna repos. interval...
                % ERROR
                error_code = 1;
                error_msg = 'The scan duration is not a multiple of the antenna repos. interval (taking into account the threshold).';
                return;
            end
        end
        sched_data_tmp.scan(i_scan).number_of_epochs = length(repos_epochs_jd);

        sched_data_tmp.scan(i_scan).sat_id          = satellite_id;
        sched_data_tmp.scan(i_scan).sat_name        = stat_data.stat(1).sat(satellite_id).TLE_data.TLE_header_line;
        sched_data_tmp.scan(i_scan).sat_number      = stat_data.stat(1).sat(satellite_id).TLE_data.sat_number;
        sched_data_tmp.scan(i_scan).quasar_id       = [];
        sched_data_tmp.scan(i_scan).quasar_name     = [];
        sched_data_tmp.scan(i_scan).obs_type        = 'sat';
        
        obs_data_tmp.sat(satellite_id).number_of_scans  = obs_data_tmp.sat(satellite_id).number_of_scans + 1;
        if t_end_jd > obs_data_tmp.sat(satellite_id).last_obs_jd
            obs_data_tmp.sat(satellite_id).last_obs_jd  = t_end_jd;
        end
        

        % Topocentric pointing angles at scan start and end:
        for i_stat = 1 : length(station_id_list)
            
            station_id = station_id_list(i_stat);
            
            % % Preallocate "epoch" data for stepwise tracking:
            for i_epoch = 1 : sched_data_tmp.scan(i_scan).number_of_epochs
                sched_data_tmp.scan(i_scan).stat(i_stat).epoch(i_epoch).jd = repos_epochs_jd(i_epoch);
            end
            sched_data_tmp.scan(i_scan).stat(i_stat).stat_id = station_id;

            
            % ### calculate antenna pointing data (az, el, ha, dc) ###
            
            % t_start_jd:
            % Check t_start_jd (Just to be sure....):
            if (t_start_jd ~= obs_data_tmp.stat(station_id).begin_of_new_obs_temp.jd)
                t_diff_sec = abs(t_start_jd - obs_data_tmp.stat(station_id).begin_of_new_obs_temp.jd) * 86400;
                fprintf(1, 'WARNING: Scan start time in "obs_data" and "t_start_jd is not equal! Difference = %1.1f sec."', t_diff_sec);
                input_str = input(' Continue? (y = yes, other characters = no => ERROR mesg.) ', 's');
                switch(input_str)
                	case 'y'
                         % continue...
                    otherwise
                        error_code = 1;
                        error_msg = '';
                        return;
                end
            end
            % SGP4 orbit propagation:
            [az, el, range, az_rate, el_rate, r_rate, tra, dc, tra_rate, tdec_rate, ha, local_hour_angle_rate, error_code, error_msg] = tle_2_topo(stat_data, station_id, PARA, satellite_id, t_start_jd);
            sched_data_tmp.scan(i_scan).stat(i_stat).start.az = az * pi/180;
            sched_data_tmp.scan(i_scan).stat(i_stat).start.el = el * pi/180;
            sched_data_tmp.scan(i_scan).stat(i_stat).start.ha = ha * pi/180;
            sched_data_tmp.scan(i_scan).stat(i_stat).start.dc = dc * pi/180;
            sched_data_tmp.scan(i_scan).stat(i_stat).start.un_az = obs_data_tmp.stat(station_id).begin_of_new_obs_temp.un_az;
            sched_data_tmp.scan(i_scan).stat(i_stat).start.jd = t_start_jd;
            
            % Calc. sky coverage number at scan start for this station and source:
            [sky_cov_num] = zcoverage(az * pi/180, el * pi/180);
            
            sched_data.scan(i_scan).stat(i_stat).sky_cov_num  = sky_cov_num;
            
            % t_end_jd:
            % Check t_end_jd (Just to be sure....):
            if (t_end_jd ~= obs_data_tmp.stat(station_id).end_of_new_obs_temp.jd)
                t_diff_sec = abs(t_end_jd - obs_data_tmp.stat(station_id).end_of_new_obs_temp.jd) * 86400;
                fprintf(1, 'WARNING: Scan end time in "obs_data" and "t_end_jd is not equal! Difference = %1.1f sec."', t_diff_sec);
                input_str = input(' Continue? (y = yes, other characters = no => ERROR mesg.) ', 's');
                switch(input_str)
                	case 'y'
                         % continue...
                    otherwise
                        error_code = 1;
                        error_msg = '';
                        return;
                end
            end
            % SGP4 orbit propagation:
            [az, el, range, az_rate, el_rate, r_rate, tra, dc, tra_rate, tdec_rate, ha, local_hour_angle_rate, error_code, error_msg] = tle_2_topo(stat_data, station_id, PARA, satellite_id, t_end_jd);
            sched_data_tmp.scan(i_scan).stat(i_stat).end.az = az * pi/180;
            sched_data_tmp.scan(i_scan).stat(i_stat).end.el = el * pi/180;
            sched_data_tmp.scan(i_scan).stat(i_stat).end.ha = ha * pi/180;
            sched_data_tmp.scan(i_scan).stat(i_stat).end.dc = dc * pi/180;
            sched_data_tmp.scan(i_scan).stat(i_stat).end.un_az = obs_data_tmp.stat(station_id).end_of_new_obs_temp.un_az;
            sched_data_tmp.scan(i_scan).stat(i_stat).end.jd = t_end_jd;
            
            % #### Save data to "obs_data" structure ####
            obs_data_tmp.stat(station_id).end_of_last_obs.jd = t_end_jd;             % [JD]
            obs_data_tmp.stat(station_id).end_of_last_obs.az = az * pi/180;          % [rad]
            obs_data_tmp.stat(station_id).end_of_last_obs.el = el * pi/180;          % [rad]
            obs_data_tmp.stat(station_id).end_of_last_obs.ha = ha * pi/180;          % [rad]
            obs_data_tmp.stat(station_id).end_of_last_obs.dc = dc * pi/180;          % [rad]
            obs_data_tmp.stat(station_id).end_of_last_obs.un_az = obs_data_tmp.stat(station_id).end_of_new_obs_temp.un_az; % [rad]
            obs_data_tmp.stat(station_id).end_of_last_obs.sat_id = satellite_id;
            obs_data_tmp.stat(station_id).end_of_last_obs.quasar_id = [];
            
            % Update sky-coverage 
            if PARA.CALC_SKYCOV_SAT
                obs_data_tmp.stat(station_id).sky_cov_sat(sky_cov_num)  = t_start_jd;
            else
                obs_data_tmp.stat(station_id).sky_cov(sky_cov_num)      = t_start_jd;
            end
            
            
            % ##### Reset temp. "obs_data" #####
            obs_data_tmp.stat(station_id).begin_of_new_obs_temp.un_az = [];
            obs_data_tmp.stat(station_id).begin_of_new_obs_temp.el = [];
            obs_data_tmp.stat(station_id).begin_of_new_obs_temp.jd = [];
            obs_data_tmp.stat(station_id).end_of_new_obs_temp.un_az = [];
            obs_data_tmp.stat(station_id).end_of_new_obs_temp.el = [];
            obs_data_tmp.stat(station_id).end_of_new_obs_temp.jd = [];
        end

        
    % ##########################
    % ##### Quasar scan    #####
    % ##########################
    elseif isempty(satellite_id) && ~isempty(quasar_id)
        
   
        % #### Save data to "sched_data" structure ####
        
        % Preallocate "epoch" data:
        sched_data_tmp.scan(i_scan).number_of_epochs = 1;
        
        sched_data_tmp.scan(i_scan).sat_id                  = [];
        sched_data_tmp.scan(i_scan).sat_name                = [];
        sched_data_tmp.scan(i_scan).sat_number              = [];
        sched_data_tmp.scan(i_scan).quasar_id               = quasar_id;
        sched_data_tmp.scan(i_scan).quasar_name             = source(quasar_id).name;
        sched_data_tmp.scan(i_scan).obs_type                = 'quasar';
        
        obs_data_tmp.quasars(quasar_id).number_of_scans     = obs_data_tmp.quasars(quasar_id).number_of_scans + 1;
        if t_end_jd > obs_data_tmp.quasars(quasar_id).last_obs_jd
            obs_data_tmp.quasars(quasar_id).last_obs_jd         = t_end_jd;
        end
        
        
        % Topocentric pointing angles at scan start and end:
        for i_stat = 1 : length(station_id_list)
            
            station_id = station_id_list(i_stat);
            
            % Preallocate "epoch" data for stepwise tracking
            sched_data_tmp.scan(i_scan).stat(i_stat).epoch(1).jd = t_start_jd;
            sched_data_tmp.scan(i_scan).stat(i_stat).stat_id = station_id;

            % t_start_jd:
            % Check t_start_jd (Just to be sure....):
            if (t_start_jd ~= obs_data_tmp.stat(station_id).begin_of_new_obs_temp.jd)
                t_diff_sec = abs(t_start_jd - obs_data_tmp.stat(station_id).begin_of_new_obs_temp.jd) * 86400;
                fprintf(1, 'WARNING: Scan start time in "obs_data" and "t_start_jd is not equal! Difference = %1.1f sec."', t_diff_sec);
                input_str = input(' Continue? (y = yes, other characters = no => ERROR mesg.) ', 's');
                switch(input_str)
                	case 'y'
                         % continue...
                    otherwise
                        error_code = 1;
                        error_msg = '';
                        return;
                end
            end
            [az, el, ha, dc] = zazel_s(t_start_jd - 2400000.5, stat_data.stat(station_id).location.ellipsoid.long * pi/180, stat_data.stat(station_id).location.ellipsoid.lat * pi/180, source(quasar_id).ra, source(quasar_id).de);
            sched_data_tmp.scan(i_scan).stat(i_stat).start.az = az;
            sched_data_tmp.scan(i_scan).stat(i_stat).start.el = el;
            sched_data_tmp.scan(i_scan).stat(i_stat).start.ha = ha;
            sched_data_tmp.scan(i_scan).stat(i_stat).start.dc = dc;
            sched_data_tmp.scan(i_scan).stat(i_stat).start.un_az = obs_data_tmp.stat(station_id).begin_of_new_obs_temp.un_az;
            sched_data_tmp.scan(i_scan).stat(i_stat).start.jd = t_start_jd;
            
            % Calc. sky coverage number at scan start for this station and source:
            [sky_cov_num] = zcoverage(az, el);
            
            sched_data.scan(i_scan).stat(i_stat).sky_cov_num  = sky_cov_num;
                        
            % t_end_jd:
            % Check t_end_jd (Just to be sure....):
            if (t_end_jd ~= obs_data_tmp.stat(station_id).end_of_new_obs_temp.jd)
                t_diff_sec = abs(t_end_jd - obs_data_tmp.stat(station_id).end_of_new_obs_temp.jd) * 86400;
                fprintf(1, 'WARNING: Scan end time in "obs_data" and "t_end_jd is not equal! Difference = %1.1f sec."', t_diff_sec);
                input_str = input(' Continue? (y = yes, other characters = no => ERROR mesg.) ', 's');
                switch(input_str)
                	case 'y'
                         % continue...
                    otherwise
                        error_code = 1;
                        error_msg = '';
                        return;
                end
            end
            [az, el, ha, dc] = zazel_s(t_end_jd - 2400000.5, stat_data.stat(station_id).location.ellipsoid.long * pi/180, stat_data.stat(station_id).location.ellipsoid.lat * pi/180, source(quasar_id).ra, source(quasar_id).de);
            sched_data_tmp.scan(i_scan).stat(i_stat).end.az = az;
            sched_data_tmp.scan(i_scan).stat(i_stat).end.el = el;
            sched_data_tmp.scan(i_scan).stat(i_stat).end.ha = ha;
            sched_data_tmp.scan(i_scan).stat(i_stat).end.dc = dc;
            sched_data_tmp.scan(i_scan).stat(i_stat).end.un_az = obs_data_tmp.stat(station_id).end_of_new_obs_temp.un_az;
            sched_data_tmp.scan(i_scan).stat(i_stat).end.jd = t_end_jd;
            
            % #### Save data to "obs_data" structure ####
            obs_data_tmp.stat(station_id).end_of_last_obs.jd        = t_end_jd;    % [JD]
            obs_data_tmp.stat(station_id).end_of_last_obs.az        = az;          % [rad]
            obs_data_tmp.stat(station_id).end_of_last_obs.el        = el;          % [rad]
            obs_data_tmp.stat(station_id).end_of_last_obs.ha        = ha;          % [rad]
            obs_data_tmp.stat(station_id).end_of_last_obs.dc        = dc;          % [rad]
            obs_data_tmp.stat(station_id).end_of_last_obs.un_az     = obs_data_tmp.stat(station_id).end_of_new_obs_temp.un_az; % [rad]
            obs_data_tmp.stat(station_id).end_of_last_obs.sat_id    = [];
            obs_data_tmp.stat(station_id).end_of_last_obs.quasar_id = quasar_id;
            obs_data_tmp.stat(station_id).sky_cov(sky_cov_num)      = t_start_jd; % Update sky-coverage
                        
            % ##### Reset temp. "obs_data" #####
            obs_data_tmp.stat(station_id).begin_of_new_obs_temp.un_az = [];
            obs_data_tmp.stat(station_id).begin_of_new_obs_temp.el = [];
            obs_data_tmp.stat(station_id).begin_of_new_obs_temp.jd = [];
            obs_data_tmp.stat(station_id).end_of_new_obs_temp.un_az = [];
            obs_data_tmp.stat(station_id).end_of_new_obs_temp.el = [];
            obs_data_tmp.stat(station_id).end_of_new_obs_temp.jd = [];
            
        end % for i_stat = 1 : length(station_id_list)
        
    else % ERROR
        error_code = 1;
        error_msg = 'Observation type undefined.';
        return;
    end
    
    % #### Assign output data here #####
    % => If an error appears earlier in this script, the original input structures will be returned!
    obs_data    = obs_data_tmp;
    sched_data  = sched_data_tmp;
    
    % ##### Save updated mat-file #####
    save([PARA.pfolder 'sched_data.mat'], 'sched_data');
    save([PARA.pfolder 'obs_data.mat'], 'obs_data');
        

   
return;
