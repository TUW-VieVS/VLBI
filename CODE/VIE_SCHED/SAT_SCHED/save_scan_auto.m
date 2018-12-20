% #########################################################################
% #     save_scan_auto
% #########################################################################
%
% DESCRIPTION
%   This function saves the new scan to the required structures (during the automatic scheduling!):
%    - sched_data
%    - obs_data
%
%
% CREATED  
%   2015-11-17     Andreas Hellerschmied
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
% - sched_data          - Scheduling data structure
% - obs_data            - Observation data structure (preallocated in "find_observation_times.m")
% - source              - source (quasars) data structure
% - PARA                - Global scheduling parameter structure
% - scan_cons           - Scan configurations structure
% - scon_id             - Scan configurations ID
%
%
% OUTPUT
% - error_code          - Error Code (0 = no erros occured)
% - error_msg           - Error Message (empty, if no errors occured)
% - obs_data            - Observation data structure (preallocated in "find_observation_times.m")
% - sched_data          - Scheduling data structure
%
% CHANGES:
% - 2016-04-28, A. Hellerschmied: Save updated "scan_data" and "obs_data" structure to the LEVEL5 sub-dir.
% - 2016-05-11, A. Hellerschmied: Handling of sky coverage for satellite scans individually added.
% - 2016-06-21, A. Hellerschmied: Indexing error fixed. Caused a problem, if the station network for satellite observations was changed.
% - 2016-10-28, A. Hellerschmied: Minor revision

function [sched_data, obs_data, error_code, error_msg] = save_scan_auto(stat_data, sched_data, obs_data, source, PARA, scan_cons, scon_id)
 
    % ##### Init #####
    error_code = 0;
    error_msg = ''; 
    
    if isempty(obs_data.end_of_last_scan_jd)
        obs_data.end_of_last_scan_jd = 0;
    end
        
    
    % ##### Loop over all sub-scans #####
    for i_scan = 1 : scan_cons(scon_id).number_of_scans
        
        % ##### Save data to "obs_data" structure #####
        obs_data.number_of_scans = obs_data.number_of_scans + 1;
        
        % End time of last scan (latest scan end of all sub-scans!)
        if (scan_cons(scon_id).scan(i_scan).scan_end_jd > obs_data.end_of_last_scan_jd)
            obs_data.end_of_last_scan_jd = scan_cons(scon_id).scan(i_scan).scan_end_jd;
        end

        
        % ##### Save data to "sched_data" structure #####
        sched_data.number_of_scans = sched_data.number_of_scans + 1;
        
        % Nominal session starte/end time:
        if isempty(sched_data.t_nominal_start_jd)
            sched_data.t_nominal_start_jd = scan_cons(scon_id).scan(i_scan).scan_start_jd;
        elseif sched_data.t_nominal_start_jd > scan_cons(scon_id).scan(i_scan).scan_start_jd
            sched_data.t_nominal_start_jd = scan_cons(scon_id).scan(i_scan).scan_start_jd;
        end
        if isempty(sched_data.t_nominal_end_jd)
            sched_data.t_nominal_end_jd = scan_cons(scon_id).scan(i_scan).scan_end_jd;
        elseif sched_data.t_nominal_end_jd < scan_cons(scon_id).scan(i_scan).scan_end_jd
            sched_data.t_nominal_end_jd = scan_cons(scon_id).scan(i_scan).scan_end_jd;
        end

        sched_data.scan(sched_data.number_of_scans).number_of_stat      = length(scan_cons(scon_id).scan(i_scan).list_stat_id_obs);
        sched_data.scan(sched_data.number_of_scans).repos_int_sec       = PARA.REPOS_INT;
        sched_data.scan(sched_data.number_of_scans).flag_record_scan    = 1;                       % !!!! Hard coded here! If it is required, change that later on.... 
        sched_data.scan(sched_data.number_of_scans).t_start_jd          = scan_cons(scon_id).scan(i_scan).scan_start_jd;
        sched_data.scan(sched_data.number_of_scans).t_end_jd            = scan_cons(scon_id).scan(i_scan).scan_end_jd;




        % ##### Observation type dependend parameters #####

        % ##########################
        % ##### Satellite scan #####
        % ##########################
        switch(scan_cons(scon_id).scan(i_scan).scan_type)
            
            case 's'

                % #### Save data to "sched_data" structure ####

                % Preallocate "epoch" data:
                [flag_scan_duration_ok, alt_t_end_jd, repos_epochs_jd, error_code, error_msg] = check_scan_duration_and_calc_repos_epochs(PARA, scan_cons(scon_id).scan(i_scan).scan_start_jd, scan_cons(scon_id).scan(i_scan).scan_end_jd);
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
                sched_data.scan(sched_data.number_of_scans).number_of_epochs = length(repos_epochs_jd);

                sched_data.scan(sched_data.number_of_scans).sat_id      = scan_cons(scon_id).scan(i_scan).sat_id;
                sched_data.scan(sched_data.number_of_scans).sat_name    = stat_data.stat(1).sat(scan_cons(scon_id).scan(i_scan).sat_id).TLE_data.TLE_header_line;
                sched_data.scan(sched_data.number_of_scans).sat_number  = stat_data.stat(1).sat(scan_cons(scon_id).scan(i_scan).sat_id).TLE_data.sat_number;
                sched_data.scan(sched_data.number_of_scans).quasar_id   = [];
                sched_data.scan(sched_data.number_of_scans).quasar_name = [];
                sched_data.scan(sched_data.number_of_scans).obs_type    = 'sat';
                
                obs_data.sat(scan_cons(scon_id).scan(i_scan).sat_id).number_of_scans      = obs_data.sat(scan_cons(scon_id).scan(i_scan).sat_id).number_of_scans + 1;


                % Topocentric pointing angles at scan start and end:
                for i_stat = 1 : length(scan_cons(scon_id).scan(i_scan).list_stat_id_obs)

                    station_id = scan_cons(scon_id).scan(i_scan).list_stat_id_obs(i_stat);

                    % % Preallocate "epoch" data for stepwise tracking:
                    for i_epoch = 1 : sched_data.scan(sched_data.number_of_scans).number_of_epochs
                        sched_data.scan(sched_data.number_of_scans).stat(i_stat).epoch(i_epoch).jd = repos_epochs_jd(i_epoch);
                    end
                    
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).stat_id        = station_id;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).duration_sec   = scan_cons(scon_id).scan(i_scan).stat(i_stat).duration_sec;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).slew_time_sec  = scan_cons(scon_id).scan(i_scan).stat(i_stat).slew_time_sec;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).sky_cov_num    = scan_cons(scon_id).scan(i_scan).stat(i_stat).sky_cov_num;


                    % ### calculate antenna pointing data (az, el, ha, dc) ###
                    
                    % SGP4 orbit propagation:
                    [az, el, range, az_rate, el_rate, r_rate, tra, dc, tra_rate, tdec_rate, ha, local_hour_angle_rate, error_code, error_msg] = tle_2_topo(stat_data, station_id, PARA, scan_cons(scon_id).scan(i_scan).sat_id, scan_cons(scon_id).scan(i_scan).scan_start_jd);
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).start.az       = az * pi/180;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).start.el       = el * pi/180;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).start.ha       = ha * pi/180;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).start.dc       = dc * pi/180;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).start.un_az    = scan_cons(scon_id).scan(i_scan).stat(i_stat).begin_of_new_obs_temp.un_az;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).start.jd       = scan_cons(scon_id).scan(i_scan).scan_start_jd;

                    % SGP4 orbit propagation (obs. end time) - individual end time for each station!
                    [az, el, range, az_rate, el_rate, r_rate, tra, dc, tra_rate, tdec_rate, ha, local_hour_angle_rate, error_code, error_msg] = tle_2_topo(stat_data, station_id, PARA, scan_cons(scon_id).scan(i_scan).sat_id, scan_cons(scon_id).scan(i_scan).stat(i_stat).end_of_new_obs_temp.jd);
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).end.az     = az * pi/180;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).end.el     = el * pi/180;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).end.ha     = ha * pi/180;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).end.dc     = dc * pi/180;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).end.un_az  = scan_cons(scon_id).scan(i_scan).stat(i_stat).end_of_new_obs_temp.un_az;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).end.jd     = scan_cons(scon_id).scan(i_scan).stat(i_stat).end_of_new_obs_temp.jd;

                    % #### Save data to "obs_data" structure ####
                    obs_data.stat(station_id).end_of_last_obs.jd        = scan_cons(scon_id).scan(i_scan).stat(i_stat).end_of_new_obs_temp.jd;             % [JD]
                    obs_data.stat(station_id).end_of_last_obs.az        = az * pi/180;          % [rad]
                    obs_data.stat(station_id).end_of_last_obs.el        = el * pi/180;          % [rad]
                    obs_data.stat(station_id).end_of_last_obs.ha        = ha * pi/180;          % [rad]
                    obs_data.stat(station_id).end_of_last_obs.dc        = dc * pi/180;          % [rad]
                    obs_data.stat(station_id).end_of_last_obs.un_az     = scan_cons(scon_id).scan(i_scan).stat(i_stat).end_of_new_obs_temp.un_az; % [rad]
                    obs_data.stat(station_id).end_of_last_obs.sat_id    = scan_cons(scon_id).scan(i_scan).sat_id;
                    obs_data.stat(station_id).end_of_last_obs.quasar_id = [];
                    
                    % Update sky-coverage 
                    if PARA.CALC_SKYCOV_SAT 
                        obs_data.stat(station_id).sky_cov_sat(scan_cons(scon_id).scan(i_scan).stat(i_stat).sky_cov_num) = scan_cons(scon_id).scan(i_scan).scan_start_jd; 
                    else
                        obs_data.stat(station_id).sky_cov(scan_cons(scon_id).scan(i_scan).stat(i_stat).sky_cov_num)     = scan_cons(scon_id).scan(i_scan).scan_start_jd; 
                    end
                    
                    if scan_cons(scon_id).scan(i_scan).stat(i_stat).end_of_new_obs_temp.jd > obs_data.sat(scan_cons(scon_id).scan(i_scan).sat_id).last_obs_jd
                        obs_data.sat(scan_cons(scon_id).scan(i_scan).sat_id).last_obs_jd = scan_cons(scon_id).scan(i_scan).stat(i_stat).end_of_new_obs_temp.jd;
                    end

                    % ##### Reset temp. "obs_data" #####
                    obs_data.stat(station_id).begin_of_new_obs_temp.un_az   = [];
                    obs_data.stat(station_id).begin_of_new_obs_temp.el      = [];
                    obs_data.stat(station_id).begin_of_new_obs_temp.jd      = [];
                    obs_data.stat(station_id).end_of_new_obs_temp.un_az     = [];
                    obs_data.stat(station_id).end_of_new_obs_temp.el        = [];
                    obs_data.stat(station_id).end_of_new_obs_temp.jd        = [];

                end


        % ##########################
        % ##### Quasar scan    #####
        % ##########################
            case 'q'


                % #### Save data to "sched_data" structure ####

                % Preallocate "epoch" data:
                sched_data.scan(sched_data.number_of_scans).number_of_epochs = 1;


                sched_data.scan(sched_data.number_of_scans).sat_id      = [];
                sched_data.scan(sched_data.number_of_scans).sat_name    = [];
                sched_data.scan(sched_data.number_of_scans).sat_number  = [];
                sched_data.scan(sched_data.number_of_scans).quasar_id   = scan_cons(scon_id).scan(i_scan).quasar_id;
                sched_data.scan(sched_data.number_of_scans).quasar_name = source(scan_cons(scon_id).scan(i_scan).quasar_id).name;
                sched_data.scan(sched_data.number_of_scans).obs_type    = 'quasar';
                
                obs_data.quasars(scan_cons(scon_id).scan(i_scan).quasar_id).number_of_scans      = obs_data.quasars(scan_cons(scon_id).scan(i_scan).quasar_id).number_of_scans + 1;


                % Topocentric pointing angles at scan start and end:
                for i_stat = 1 : length(scan_cons(scon_id).scan(i_scan).list_stat_id_obs)

                    station_id = scan_cons(scon_id).scan(i_scan).list_stat_id_obs(i_stat);

                    % Preallocate "epoch" data for stepwise tracking
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).epoch(1).jd = scan_cons(scon_id).scan(i_scan).scan_start_jd;
                    
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).stat_id        = station_id;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).duration_sec   = scan_cons(scon_id).scan(i_scan).stat(i_stat).duration_sec;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).slew_time_sec  = scan_cons(scon_id).scan(i_scan).stat(i_stat).slew_time_sec;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).sky_cov_num    = scan_cons(scon_id).scan(i_scan).stat(i_stat).sky_cov_num;

                    [az, el, ha, dc] = zazel_s(scan_cons(scon_id).scan(i_scan).stat(i_stat).begin_of_new_obs_temp.jd - 2400000.5, stat_data.stat(station_id).location.ellipsoid.long * pi/180, stat_data.stat(station_id).location.ellipsoid.lat * pi/180, source(scan_cons(scon_id).scan(i_scan).quasar_id).ra, source(scan_cons(scon_id).scan(i_scan).quasar_id).de);
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).start.az       = az;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).start.el       = el;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).start.ha       = ha;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).start.dc       = dc;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).start.un_az    = scan_cons(scon_id).scan(i_scan).stat(i_stat).begin_of_new_obs_temp.un_az;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).start.jd       = scan_cons(scon_id).scan(i_scan).stat(i_stat).begin_of_new_obs_temp.jd;
                    
                    [az, el, ha, dc] = zazel_s(scan_cons(scon_id).scan(i_scan).stat(i_stat).end_of_new_obs_temp.jd - 2400000.5, stat_data.stat(station_id).location.ellipsoid.long * pi/180, stat_data.stat(station_id).location.ellipsoid.lat * pi/180, source(scan_cons(scon_id).scan(i_scan).quasar_id).ra, source(scan_cons(scon_id).scan(i_scan).quasar_id).de);
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).end.az     = az;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).end.el     = el;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).end.ha     = ha;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).end.dc     = dc;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).end.un_az  = scan_cons(scon_id).scan(i_scan).stat(i_stat).end_of_new_obs_temp.un_az;
                    sched_data.scan(sched_data.number_of_scans).stat(i_stat).end.jd     = scan_cons(scon_id).scan(i_scan).stat(i_stat).end_of_new_obs_temp.jd;

                    % #### Save data to "obs_data" structure ####
                    obs_data.stat(station_id).end_of_last_obs.jd        = scan_cons(scon_id).scan(i_scan).stat(i_stat).end_of_new_obs_temp.jd; % [JD]
                    obs_data.stat(station_id).end_of_last_obs.az        = az;          % [rad]
                    obs_data.stat(station_id).end_of_last_obs.el        = el;          % [rad]
                    obs_data.stat(station_id).end_of_last_obs.ha        = ha;          % [rad]
                    obs_data.stat(station_id).end_of_last_obs.dc        = dc;          % [rad]
                    obs_data.stat(station_id).end_of_last_obs.un_az     = scan_cons(scon_id).scan(i_scan).stat(i_stat).end_of_new_obs_temp.un_az; % [rad]
                    obs_data.stat(station_id).end_of_last_obs.sat_id    = [];
                    obs_data.stat(station_id).end_of_last_obs.quasar_id = scan_cons(scon_id).scan(i_scan).quasar_id;
                    obs_data.stat(station_id).sky_cov(scan_cons(scon_id).scan(i_scan).stat(i_stat).sky_cov_num) = scan_cons(scon_id).scan(i_scan).scan_start_jd; % Update sky-coverage 
                    
                    if scan_cons(scon_id).scan(i_scan).stat(i_stat).end_of_new_obs_temp.jd > obs_data.quasars(scan_cons(scon_id).scan(i_scan).quasar_id).last_obs_jd
                        obs_data.quasars(scan_cons(scon_id).scan(i_scan).quasar_id).last_obs_jd = scan_cons(scon_id).scan(i_scan).stat(i_stat).end_of_new_obs_temp.jd;
                    end
                    

                    % ##### Reset temp. "obs_data" #####
                    obs_data.stat(station_id).begin_of_new_obs_temp.un_az   = [];
                    obs_data.stat(station_id).begin_of_new_obs_temp.el      = [];
                    obs_data.stat(station_id).begin_of_new_obs_temp.jd      = [];
                    obs_data.stat(station_id).end_of_new_obs_temp.un_az     = [];
                    obs_data.stat(station_id).end_of_new_obs_temp.el        = [];
                    obs_data.stat(station_id).end_of_new_obs_temp.jd        = [];

                end % for i_stat = 1 : length(station_id_list)

            otherwise % ERROR
                error_code = 1;
                error_msg = 'Observation type undefined.';
                return;
        end % switch(scan_cons(scon_id).scan(i_scan).scan_type)

    end % for i_scan = 1 : scan_cons(scon_id).number_of_scans
    
    % ##### Save updated mat-file #####
    save([PARA.pfolder 'sched_data.mat'], 'sched_data');
    save([PARA.pfolder 'obs_data.mat'], 'obs_data');
   
return;
