% #########################################################################
% #     auto_sched
% #########################################################################
%
% DESCRIPTION
%   Automatic station based scheduling for quasars and satellites
%   
%
% CREATED  
%   2015-11-??     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
%   - calc_scan_configurations
%   - sort_scan_cons
%   - calc_start_of_obs
%   - calc_scan_duration_baseline
%   - calc_un_az_begin_of_scan
%   - save_scan_auto
%
%
% INPUT
%   - 
%
%   - stat_data                     : station data structure
%   - source                        : Source structure
%   - PARA                          : Global scheduling parameter strucutre
%   - sched_data                    : Scheduling data structure (for VieVS satellite scheduling)
%   - obs_data                      : Observation data structure (preallocated in "find_observation_times.m")
%   - t_start_jd                    : Start time of the session
%   - t_end_jd                      : End time of the session
%   - quasar_network_id_list        : Station IDs of quasar network
%   - sat_network_id_list           : Station IDs of satellite network
%
%
% OUTPUT
%   - sched_data                    : Scheduling data structure (for VieVS satellite scheduling) - updated
%   - obs_data                      : Observation data structure (preallocated in "find_observation_times.m") - updated
%   - error_code                    : Error Code (0 = no erros occured)
%   - error_msg                     : Error Message (empty, if no errors occured)
%
% CHANGES:
% - 2016-05-04: A. Hellerschmied: Fixed a bug at sorting the start_jd list to determine the scan start times.
% - 2016-05-11: A. Hellerschmied: Bug fixes
% - 2016-10-27: A. Hellerschmied: Only initialize "obs_data.sat" and "obs_data.quasars", if they haven't been initialized before.
%

function [sched_data, obs_data, error_code, error_msg] = auto_sched(stat_data, source, PARA, sched_data, obs_data, t_start_jd, t_end_jd, quasar_network_id_list, sat_network_id_list)

    tic;

    % ##### Init #####
    error_code = 0;
    error_msg = '';
    flag_exit_auto_sched = false;
    number_of_scheduled_scans = 0;
    number_of_quasars  = length(source);
    number_of_sat      = length(obs_data.sat);
    
    % ##### Preallocation #####
    % obs_data...if it hasn't been done before
    if (length(obs_data.quasars) == 1) && isempty(obs_data.quasars(1).number_of_scans)
        for i = 1 : number_of_quasars
           obs_data.quasars(i).last_obs_jd      = 0;
           obs_data.quasars(i).number_of_scans  = 0;
        end
    end
    if (length(obs_data.sat) == 1) && isempty(obs_data.sat(1).number_of_scans)
        for i = 1 : number_of_sat
           obs_data.sat(i).last_obs_jd          = 0;
           obs_data.sat(i).number_of_scans      = 0;
        end
    end
    
    
    % ##### Options #####
    % Use the earliest possible scan start time (one baseline = 2 stations), or the earliest possible observation time of an individual station as reference time for the PARA.MAX_WAIT condition:
    % 1 => scan
    % 2 => observation
    switch_scan_obs_as_ref          = 2;
    block_duration_treshold         = 30;   % [sec]
    
    % ### Scan type determination ###
    init_scan_type                  = PARA.INIT_SCAN_TYPE;
    if ~strcmp(init_scan_type, 'q') && ~strcmp(init_scan_type, 's')
        error_code = 2;
        error_msg = 'PARA.INIT_SCAN_TYPE not valid!';
        return; 
    end
    sat_block_duration              = (PARA.SAT_BLOCK_DUR_MIN * 60);    % [sec]
    quasar_block_duration           = (PARA.Q_BLOCK_DUR_MIN * 60);      % [sec]
    sat_scan_duration_sec           = PARA.SAT_SCAN_DUR_SEC;            % [sec]
    
    
    % ##### Loop init #####
    switch(init_scan_type)
        case 'q'
            scan_type = 'q';
            stat_id_list = quasar_network_id_list;
        case 's'
            scan_type = 's';
            stat_id_list = sat_network_id_list;
    end
    block_start_time_jd = t_start_jd;
    
    % ##### Main Loop #####
    while ~flag_exit_auto_sched
        
        
        % ###################################    
        % ##### 1.) Determine scan type #####
        % ################################### 
        
        % ### Switch scan type depending on scan block duration: ###
        
        % Calc. duration of current scan block (sum) [sec]
        if number_of_scheduled_scans == 0
            block_duration = 0;   
        else
            block_duration = (obs_data.end_of_last_scan_jd - block_start_time_jd) * 86400;
        end

        switch(scan_type)
            
            case 'q'
                if block_duration >= (quasar_block_duration - block_duration_treshold);
                    scan_type           = 's';
                    stat_id_list        = sat_network_id_list;
                    block_start_time_jd = obs_data.end_of_last_scan_jd;
                end
                
            case 's'
                if block_duration >= (sat_block_duration - block_duration_treshold);
                    scan_type           = 'q';
                    stat_id_list        = quasar_network_id_list;
                    block_start_time_jd = obs_data.end_of_last_scan_jd;
                end
        end % switch(scan_type)
        
        
        
        % ####################################################
        % ##### 2.) Get all possible scan configurations #####
        % ####################################################
        
        % Determine all possible scans (configurations) for the current station network (sat./quasar):
        
        % Define the time epoch for the inital calcualtion of possible scan configurations:
        if number_of_scheduled_scans == 0 % First scan is goning to be scheduled...
            t_epoch_jd = t_start_jd;
        else
            t_epoch_jd = obs_data.end_of_last_scan_jd + PARA.SOURCE/86400 + PARA.TAPETM/86400 + PARA.IDLE/86400 + PARA.CALIBRATION/86400; % + estimated slew time ???
        end
        
        [scan_cons, error_code, error_msg] = calc_scan_configurations(stat_data, source, PARA, t_epoch_jd, stat_id_list, scan_type, obs_data);
        if error_code > 0
            error_msg = ['calc_scan_configurations: ', error_msg];
            return; 
        end
        
        
        % ##### 3.) Initially sort all scan configurations #####
        [sort_index, scan_cons, obs_data, error_code, error_msg] = sort_scan_cons(scan_cons, stat_data, source, PARA, t_epoch_jd, scan_type, obs_data, stat_id_list);
        if error_code > 0
            error_msg = ['sort_scan_cons: ', error_msg];
            return; 
        end
        
        % ##### Loop over the scan configurations #####
        % Take into account the best "PARA.SORTNUM" scan cons. 
        % Begin with the best-rated scan con.
        %
        % Within the following loop, again various schduling conditions are checked:
        %  - PARA.MAX_WAIT
        %  - PARA.MIN_CUTEL (already done in the function "calc_scan_configurations")
        %  - PARA.MAXSLEWTIME
        %  - PARA.MAX_SCAN
        %  - PARA.MIN_SCAN
        %  - 
        % => These conditions influence the data in "scan_cons", e.g. stations may drop out,...
        % => Finally, if the requirements are not met any more, a whole scan con may be thrown out. 
        
        scon_count = 0;
        
        for i_sort_index = length(sort_index) : -1 : 1
            
            % ##### Check, if max. number of scan cons, which should be considered, is already reached: #####
            if scon_count == PARA.SORTNUM
                break;
            end
            
            % Loop init.:
            i_scon = sort_index(i_sort_index); % index of currently treated scan con.
            flag_next_scan_con = false;


            % #### Loop over all scans within this scan con. ####
            for i_scan = 1 : scan_cons(i_scon).number_of_scans
                
                % Loop init.:
                number_of_stat_in_scan = length(scan_cons(i_scon).scan(i_scan).list_stat_id_obs);
                del_index = true(1, number_of_stat_in_scan);
                start_jd_list = zeros(1, number_of_stat_in_scan);
                
                switch(scan_type)
                    case 'q' % Quasar scan
                        source_quasar = source(scan_cons(i_scon).scan(i_scan).quasar_id);
                        satellite_id = [];
                    case 's' % Satellite scan
                        source_quasar = [];
                        satellite_id = scan_cons(i_scon).scan(i_scan).sat_id;
                end % switch(scan_type) 
                
                
                % #####################################
                % ##### 4.) Cal. scan start times #####
                % #####################################
                
                % ##### Loop over all stations within this scan #####
                for i_stat = 1 : number_of_stat_in_scan

                    stat_id = scan_cons(i_scon).scan(i_scan).list_stat_id_obs(i_stat);
% TEST TEST TEST  
% fprintf('====>> i_scan: %d, i_stat: %d, i_scon: %d\n', i_scan, i_stat, i_scon);
% TEST TEST TEST  
                    [start_jd, slew_time_sec, obs_data, error_code, error_msg] = calc_start_of_obs(stat_data, stat_id, obs_data, PARA, source_quasar, satellite_id, t_start_jd);
                    if error_code  == 5 % Max. number of iterations for the calculation of the obs start time exceeded!
                        fprintf(' => ERROR: calc_start_of_obs: "%s"; i_scan=%d; i_scon=%d\n', error_msg, i_scan, i_scon);
                        % Remove station from station ID list
                        del_index(i_stat) = false;
                        continue;
                    elseif error_code > 0
                        error_msg = ['calc_start_of_obs: ', error_msg];
                        return; 
                    end
% TEST TEST TEST        
% fprintf(' => slew time (%s): %3d sec (start time: %s)\n',stat_data.stat(stat_id).name ,slew_time_sec, jd2datestr(obs_data.end_of_last_scan_jd));
% TEST TEST TEST  
                    
                    % #### Check max. slew time (PARA.MAXSLEWTIME) ####
                    if slew_time_sec > PARA.MAXSLEWTIME
                        % Remove station from station ID list
                        del_index(i_stat) = false;
                        continue;
                    end
                    
                    % #### Save data ####
                    scan_cons(i_scon).scan(i_scan).stat(i_stat).begin_of_new_obs_temp.jd = start_jd;
                    scan_cons(i_scon).scan(i_scan).stat(i_stat).slew_time_sec = slew_time_sec;
                    start_jd_list(i_stat) = start_jd;

                end % for i_stat = 1 : number_of_stat_in_scan
                
                % ##### Sort scan start times #####
                start_jd_list = sort(start_jd_list);
                
                start_jd_list_temp = start_jd_list(del_index); % Only consider stations which are not already thrown out due to the slew-time condition.
                
                
                % ### Check, if there are enough stations left (1) ###
                if length(start_jd_list_temp) >= max(2, PARA.MIN_STANUM)
                    
%                     % ### Get earliest start time [jd] within this scan ###
%                     start_jd_list_temp = sort(start_jd_list_temp);
                    
                    % Option: Set earliest possible obsertvation-, or scan-start time as reference.
                    if switch_scan_obs_as_ref == 1
                        t_ref_jd = start_jd_list_temp(2); % Earliest SCAN START TIME (one scan requires two stations = one baseline)
                    elseif switch_scan_obs_as_ref == 2
                        t_ref_jd = start_jd_list_temp(1); % Earliest OBSERVATION START TIME
                    end
                
                    % #### Check if stations arrive late (PARA.MAX_WAIT) within this scan => set flags in "del_index" ####
                    for i_stat = 1 : number_of_stat_in_scan
                        if del_index(i_stat) == true % Only consider stations which are not already thrown out!
                            if (((scan_cons(i_scon).scan(i_scan).stat(i_stat).begin_of_new_obs_temp.jd - t_ref_jd) * 86400) > PARA.MAX_WAIT)
                                del_index(i_stat) = false; % Update "del_index"
                            end
                        end
                    end % for i_stat = 1 : length(scan_cons(i_scon).scan(i_scan).list_stat_id_obs)
                    
                else % Not enough stations left!
                    flag_next_scan_con = true;
                    break;
                end % if length(start_jd_list_temp) >= max(2, PARA.MIN_STANUM)
                
                start_jd_list_temp = start_jd_list(del_index); % Only consider stations which are not already thrown out due to the slew-time and max. wait conditions.
                
                % ### Check, if there are enough stations left (2) ###
                if length(start_jd_list_temp) >= max(2, PARA.MIN_STANUM)
                
                     % #### Calc. start time of this scan ####
                     scan_cons(i_scon).scan(i_scan).scan_start_jd = start_jd_list_temp(end); % Latest start time of the individual stations within this scan
                
                else % Not enough stations left!
                    flag_next_scan_con = true;
                    break;
                end % if length(start_jd_list_temp) >= max(2, PARA.MIN_STANUM)
                
                % #### Apply the "del_index" to remove the flagged stations (arriving too late or slew-time too large) ####
                % Delete bad stations and their IDs from the station ID-list
                scan_cons(i_scon).scan(i_scan).list_stat_id_obs = scan_cons(i_scon).scan(i_scan).list_stat_id_obs(del_index);
                scan_cons(i_scon).scan(i_scan).stat(~del_index) = [];
                
                start_jd_list_2 = start_jd_list(del_index);
                
                
                % #####################################
                % ##### 5.) Calc scan durations   #####
                % #####################################
                
                % Init.:
                flag_bad_stat_in_scan = true;

                % Init a new list for flagging bad stations ("del_index")
                % The indizes of "del_index" and "initial_stat_id_list" have to match!!
                initial_stat_id_list = scan_cons(i_scon).scan(i_scan).list_stat_id_obs;
                initial_number_of_stat_in_scan = length(initial_stat_id_list);
                del_index = true(1, initial_number_of_stat_in_scan);

                % #### Iteration loop ####
                % To throw out bas stations successively...
                % Whenever a station is thrown out of the scan, the iteration loop starts again from beginning.
                % Done within this loop:
                % - Select scan start time from "start_jd_list_2"
                % - Calc./check the observation duration
                % - Calc./check the obs. end time
                % - Check various conditions for a valid scan with the "check_t_end" function
                % - Save data for this scan to "obs_data.stat.begin/end_of_new_obs_temp.un_az/el/jd"

                while flag_bad_stat_in_scan

                    % #### Check the min. number of stations in this scan ####
                    if sum(del_index) < max(2, PARA.MIN_STANUM)
                        flag_next_scan_con = true;
                        break; % while flag_bad_stat_in_scan
                    end


                    % ### Init.: ###
                    i_baseline = 0;
                    flag_bad_stat_in_scan = false;

                    temp_start_jd_list_2 = start_jd_list_2(del_index);

                    % #### Remove bad stations from initial ID list with each iteration step ####
                    temp_stat_id_list = initial_stat_id_list(del_index);
                    number_of_stat_in_scan = sum(del_index);
                    stat_duration_sum = zeros(1, number_of_stat_in_scan);
                    number_of_baselines = (number_of_stat_in_scan * (number_of_stat_in_scan - 1)/2);
                    scan_duration_data = zeros((number_of_baselines), 4); % One line per basline; columns: | stat ID 1 | stat ID 2 | duration [sec] | flag: duration OK? |

                    % #### Calc. start time of this scan & initialize "obs_data.stat.begin_of_new_obs_temp" structure for this start time ####
                    % Calc for scan start time:
                    % - un_az
                    % - el
                    %  => Save these values to: "obs_data.stat.begin_of_new_obs_temp"
                    t_scan_start = temp_start_jd_list_2(end);
                    scan_cons(i_scon).scan(i_scan).scan_start_jd = t_scan_start;

                    % Update "obs_data.stat.begin_of_new_obs_temp" for scan start time:
                    [obs_data, error_code, error_msg] = calc_un_az_begin_of_scan(stat_data, temp_stat_id_list, obs_data, PARA, source_quasar, satellite_id, t_scan_start);
                    if error_code > 0
                        error_msg = ['calc_un_az_begin_of_scan: ', error_msg];
                        return; 
                    end
                    
                    switch(scan_type)
                        % ##########################
                        % ##### Quasar scan    #####
                        % ##########################
                        
                        % - 1.) Calc. scan duration for each scan
                        % - 2.) Check if the treshold for the max. scan duration is exceeded
                        % - 3.) Calc. the obs. end time for each station
                        % - 4.) Calc the scan end time (= max. obs. time of participating stations)
                        
                        case 'q' % Quasar scan

                            % #### Loop over all possible (remaining) baselines ####
                            for i_stat_1 = 1 : (number_of_stat_in_scan - 1)
                                for i_stat_2 = (i_stat_1 + 1) : number_of_stat_in_scan

                                    % Loop init.:
                                    i_baseline = i_baseline + 1;
                                    flag_duration_OK = false;

                                    % #### Calc. observation duration for the baseline for the current scan-start epoch ####

                                    % Get station IDs:
                                    stat_id_1 = temp_stat_id_list(i_stat_1);
                                    stat_id_2 = temp_stat_id_list(i_stat_2);

                                    % Calculate observation angle (elevation):
                                    [az, el_1, ha, dc] = zazel_s(t_scan_start - 2400000.5, stat_data.stat(stat_id_1).location.ellipsoid.long * pi/180, stat_data.stat(stat_id_1).location.ellipsoid.lat * pi/180, source_quasar.ra, source_quasar.de);
                                    [az, el_2, ha, dc] = zazel_s(t_scan_start - 2400000.5, stat_data.stat(stat_id_2).location.ellipsoid.long * pi/180, stat_data.stat(stat_id_2).location.ellipsoid.lat * pi/180, source_quasar.ra, source_quasar.de);

                                    [scan_duration_sec, error_code, error_msg] = calc_scan_duration_baseline(stat_data, stat_id_1, stat_id_2, el_1, el_2, PARA, source_quasar, t_scan_start);

                                    % ### Check the scan duration ###
                                    if scan_duration_sec > PARA.MAX_SCAN
                                        stat_duration_sum(i_stat_1) = stat_duration_sum(i_stat_1) + scan_duration_sec;
                                        stat_duration_sum(i_stat_2) = stat_duration_sum(i_stat_2) + scan_duration_sec;
                                    else
                                        flag_duration_OK = true;
                                    end

                                    scan_duration_data(i_baseline, 1) = stat_id_1;
                                    scan_duration_data(i_baseline, 2) = stat_id_2;
                                    scan_duration_data(i_baseline, 3) = scan_duration_sec;
                                    scan_duration_data(i_baseline, 4) = flag_duration_OK;

                                end % for i_stat_2 = (i_stat_1 + 1) : number_of_stat_in_scan
                            end % for i_stat_1 = 1 : (number_of_stat_in_scan - 1)

                            % #### Check, if there are bad stations ####
                            if sum(scan_duration_data(:, 4)) < number_of_baselines % There is at least one bad station...
                                % Remove the station with the largest obs-duration sum (add it to the flag-list "del_index"):
                                [max_duration_sum, i_bad_stat] = max(stat_duration_sum);
                                bad_stat_id = temp_stat_id_list(i_bad_stat);
                                del_index(initial_stat_id_list == bad_stat_id) = false;
                                flag_bad_stat_in_scan = true; % => Another iteration!
                                continue; % => Another iteration! (while flag_bad_stat_in_scan)
                            end

                            % #### Calc. the observation duration and the observation end time for each station ####
                            %  - Using "scan_duration_data"
                            %  - Save the values to "scan_cons" structure
                            %  - Also calc the lastes obs. end time of all stations => "scan end jd"
                            scan_end_jd = 0;
                            for i_bl = 1 : number_of_baselines
                                stat_id_1 = scan_duration_data(i_bl, 1);
                                stat_id_2 = scan_duration_data(i_bl, 2);
                                duration  = scan_duration_data(i_bl, 3);
                                i_stat_1 = find(temp_stat_id_list == stat_id_1);
                                i_stat_2 = find(temp_stat_id_list == stat_id_2);
                                % Station one of BL:
                                if scan_cons(i_scon).scan(i_scan).stat(i_stat_1).duration_sec < duration
                                    scan_cons(i_scon).scan(i_scan).stat(i_stat_1).duration_sec = duration;
                                    scan_cons(i_scon).scan(i_scan).stat(i_stat_1).end_of_new_obs_temp.jd = scan_cons(i_scon).scan(i_scan).scan_start_jd + duration / 86400;
                                    if scan_cons(i_scon).scan(i_scan).stat(i_stat_1).end_of_new_obs_temp.jd > scan_end_jd
                                        scan_end_jd = scan_cons(i_scon).scan(i_scan).stat(i_stat_1).end_of_new_obs_temp.jd;
                                    end
                                end
                                % Station two of BL:
                                if scan_cons(i_scon).scan(i_scan).stat(i_stat_2).duration_sec < duration
                                    scan_cons(i_scon).scan(i_scan).stat(i_stat_2).duration_sec = duration;
                                    scan_cons(i_scon).scan(i_scan).stat(i_stat_2).end_of_new_obs_temp.jd = scan_cons(i_scon).scan(i_scan).scan_start_jd + duration / 86400;
                                    if scan_cons(i_scon).scan(i_scan).stat(i_stat_2).end_of_new_obs_temp.jd > scan_end_jd
                                        scan_end_jd = scan_cons(i_scon).scan(i_scan).stat(i_stat_2).end_of_new_obs_temp.jd;
                                    end
                                end
                            end % for i_bl = 1 : number_of_baselines


                        % ##########################
                        % ##### Satellite scan #####
                        % ##########################
                        
                        % - 1.) Take the obs. duration defined for satellites ("sat_scan_duration_sec")
                        % - 3.) Calc. the obs. end time for each station (= scan start time + duration)
                        % - 4.) Calc the scan end time (= max. obs. time of participating stations)
                        
                        case 's' % Satellite scan
                            
                            % - 2.)
                            scan_end_jd = scan_cons(i_scon).scan(i_scan).scan_start_jd + sat_scan_duration_sec / 86400;
                            
                            % #### Loop over all stations ####
                            for i_stat = 1 : number_of_stat_in_scan
                                
                                % - 1.)
                                scan_cons(i_scon).scan(i_scan).stat(i_stat).duration_sec = sat_scan_duration_sec;

                                % - 3.)
                                scan_cons(i_scon).scan(i_scan).stat(i_stat).end_of_new_obs_temp.jd = scan_end_jd;
                                
                            end % for i_stat = 1 : number_of_stat_in_scan
                            
                    end % switch(scan_type)

                    scan_cons(i_scon).scan(i_scan).scan_end_jd = scan_end_jd;


                    % ##### Check end time of the scan (of individual stations) #####
                    % Furthermore, initialization of:
                    %     - obs_data.stat.end_of_new_obs_temp.un_az
                    %     - obs_data.stat.end_of_new_obs_temp.el
                    %     - obs_data.stat.end_of_new_obs_temp.jd

                    % ### Loop over all stations in this scan ###
                    flag_t_end_ok = 1;
                    for i_stat = 1 : number_of_stat_in_scan
                        t_end_in        = scan_cons(i_scon).scan(i_scan).stat(i_stat).end_of_new_obs_temp.jd;
                        t_start_in      = scan_cons(i_scon).scan(i_scan).stat(i_stat).begin_of_new_obs_temp.jd;
                        stat_id         = temp_stat_id_list(i_stat);          % Call "check_t_end" for all station individually to be able to consider differing observation end times!
                        % Check_t_end checks observation conditons for the whole scan duration!
                        [flag_t_end_ok, obs_data, error_code, error_msg] = check_t_end(stat_data, stat_id, obs_data, PARA, source_quasar, satellite_id, t_end_in, t_start_in, 1);
                        if error_code > 0
                            error_msg = ['check_t_end: ', error_msg];
                            return; 
                        end
                        if ~flag_t_end_ok
                            % Remove current station => add it to "del_index" and start the next iteration round! 
                            del_index(initial_stat_id_list == stat_id) = false;
                            flag_bad_stat_in_scan = true; % => Another iteration!
                            break; % for i_stat = 1 : number_of_stat_in_scan
                        end
                         
                    end % for i_stat = 1 : number_of_stat_in_scan
                    
                    if ~flag_t_end_ok
                        continue; % => Another iteration! (while flag_bad_stat_in_scan)
                    end

                end % while flag_bad_stat_in_scan
                
                if flag_next_scan_con
                    break; % for i_scan = 1 : scan_cons(i_scon).number_of_scans
                end

                % ##### Update scan_cons #####

                % ### Apply del_index ###
                % Delete bad stations and their IDs from the station ID-list
                scan_cons(i_scon).scan(i_scan).list_stat_id_obs = scan_cons(i_scon).scan(i_scan).list_stat_id_obs(del_index);
                scan_cons(i_scon).scan(i_scan).stat(~del_index) = [];
                
                % ### Update begin/end of scan fields ###
                % Loop over all stations:
                for i_stat = 1 : length(scan_cons(i_scon).scan(i_scan).list_stat_id_obs)
                    stat_id = scan_cons(i_scon).scan(i_scan).list_stat_id_obs(i_stat);
                    scan_cons(i_scon).scan(i_scan).stat(i_stat).begin_of_new_obs_temp.un_az     = obs_data.stat(stat_id).begin_of_new_obs_temp.un_az;
                    scan_cons(i_scon).scan(i_scan).stat(i_stat).begin_of_new_obs_temp.el        = obs_data.stat(stat_id).begin_of_new_obs_temp.el;
                    % scan_cons(i_scon).scan(i_scan).stat(i_stat).begin_of_new_obs_temp.jd % ==> Already assigned!!!
                    scan_cons(i_scon).scan(i_scan).stat(i_stat).end_of_new_obs_temp.un_az       = obs_data.stat(stat_id).end_of_new_obs_temp.un_az;
                    scan_cons(i_scon).scan(i_scan).stat(i_stat).end_of_new_obs_temp.el          = obs_data.stat(stat_id).end_of_new_obs_temp.el;
                    % scan_cons(i_scon).scan(i_scan).stat(i_stat).end_of_new_obs_temp.jd % ==> Already assigned!!!
                end
                

            end % for i_scan = 1 : scan_cons(i_scon).number_of_scans
            
            if flag_next_scan_con
                continue;
            end
            
            
            % ##### Save scan con, if it meets all requirements ##### 
            scon_count = scon_count + 1;
            scan_cons_temp(scon_count) = scan_cons(i_scon);
            
% TEST TEST TEST     
% fprintf(1, 'sort_id = %4d, numb. of scons = %4d\n', i_sort_index, scon_count );
% TEST TEST TEST     
            
        end % for i_sort_index = length(sort_index) : -1 : 1 => Loop over all scan cons.
        
        
% TEST TEST TEST     
% fprintf(1, '==> Excluded scons = %4d\n', (length(sort_index) - i_sort_index) - scon_count - 1);
% TEST TEST TEST
        
        % ##### Check if at least one Scan configuration is included in "scan_cons_temp" #####
        if scon_count == 0
            fprintf('  => No valid scan configurations available!\n');
            break; % while ~flag_exit_auto_sched
        end


        % ##### Save updated scan cons. structure #####
        scan_cons = scan_cons_temp;
        clear scan_cons_temp;
        
        
        % ##### 6.) Final sorting #####
        t_epoch_jd = []; % Take the start of the individual scans as reference epoch.
        [sort_index, scan_cons, obs_data, error_code, error_msg] = sort_scan_cons(scan_cons, stat_data, source, PARA, t_epoch_jd, scan_type, obs_data, stat_id_list);
        if error_code > 0
            error_msg = ['sort_scan_cons: ', error_msg];
            return; 
        end
        
        % ##### 7.) Select best scan configuration #####
        best_scon_id = sort_index(end);
        
        
        % ##### Exit condition for the automatic scheduling #####
        % Get latest scan end time within this constellation:
        latest_scan_end_time_jd = 0;
        for i_scan = 1 : scan_cons(best_scon_id).number_of_scans
            if scan_cons(best_scon_id).scan(i_scan).scan_end_jd > latest_scan_end_time_jd
               latest_scan_end_time_jd = scan_cons(best_scon_id).scan(i_scan).scan_end_jd;
            end
        end
        if latest_scan_end_time_jd >= t_end_jd
            % flag_exit_auto_sched = true;
            break; % while ~flag_exit_auto_sched
        end
        
  
        % ##### Save scan: #####
        number_of_scheduled_scans = number_of_scheduled_scans + 1;
       
% TEST TEST TEST
% result(number_of_scheduled_scans) = scan_cons(best_scon_id);
% TEST TEST TEST
        
        % Update:
        % - sched_data
        % - obs_data
        % Preallocation:
        % - "epoch" data for calc. of antenna pointing angles  
        
        % #### Prepare data from "scan_cons" to be used by "save_scan" ####
        [sched_data, obs_data, error_code, error_msg] = save_scan_auto(stat_data, sched_data, obs_data, source, PARA, scan_cons, best_scon_id);
        if error_code > 0
            error_msg = ['save_scan_auto: ', error_msg];
            return; 
        end
        
        
        % ##### Print status msg. to CW #####
        fprintf('---------- Scan-Configuration %3d; w = %7.4f--------------------------------------------------\n', number_of_scheduled_scans, scan_cons(best_scon_id).weight);
        for i_scan = 1 : scan_cons(best_scon_id).number_of_scans
            switch(scan_type)
                case 'q' % Quasar scan
                    source_name_str = source(scan_cons(best_scon_id).scan(i_scan).quasar_id).name;
                    scan_type_str = 'quasar';
                case 's' % Satellite scan
                    source_name_str = obs_data.sat(scan_cons(best_scon_id).scan(i_scan).sat_id).name(1:24);
                    scan_type_str = 'satellite';
            end % switch(scan_type)
        fprintf(' - Sub-Scan %2d, start: %s (%s: %s)\n', i_scan, jd2datestr(scan_cons(best_scon_id).scan(i_scan).scan_start_jd),  scan_type_str, source_name_str);
        for i_stat = 1 : length(scan_cons(best_scon_id).scan(i_scan).list_stat_id_obs)
        fprintf('    - %8s: %19s - %19s; slew: %3.0fs; dur: %3.0fs; sky_cov_num: %2d \n', stat_data.stat(scan_cons(best_scon_id).scan(i_scan).list_stat_id_obs(i_stat)).name, jd2datestr(scan_cons(best_scon_id).scan(i_scan).stat(i_stat).begin_of_new_obs_temp.jd), jd2datestr(scan_cons(best_scon_id).scan(i_scan).stat(i_stat).end_of_new_obs_temp.jd), scan_cons(best_scon_id).scan(i_scan).stat(i_stat).slew_time_sec, scan_cons(best_scon_id).scan(i_scan).stat(i_stat).duration_sec, scan_cons(best_scon_id).scan(i_scan).stat(i_stat).sky_cov_num);
        end
        end
        fprintf('--------------------------------------------------------------------------------------------- ----\n');

        
% TEST TEST TEST   
% fprintf('\n  ############# Number of scheduled scans = %d ############\n\n\n', number_of_scheduled_scans);
% if number_of_scheduled_scans >= 100
%    disp('stop'); 
% end
% TEST TEST TEST   


    end % while ~flag_exit_auto_sched
    
    elapsed_time_sec = toc;
    
    fprintf('-------------------------------------------------------------------------------------------------\n');
    fprintf('  => Automatic scheduling finished (elapsed time: %7.1f sec)\n', elapsed_time_sec); 
    fprintf('  => The last scan ends at: %s \n', jd2datestr(obs_data.end_of_last_scan_jd));
    fprintf('  => Defined sessions end:  %s \n', jd2datestr(t_end_jd));
    fprintf('-------------------------------------------------------------------------------------------------\n');

end % function


% ##### To Do ##### 
% - Integrate satellite support! (Which things have to be calculated/checked sepatetly for sats/quasars?) - done!!!!
% - Maybe it is better/easier to initialize "del_index"(...etc...) only once at the very beginning of checking a scan within scan_cons. - yes, but.... too much work...
%    - At the moment, it is initialized again, after calculating obs. start times...
% - Save scans => done!!!
% - Update Sky coverage!!!! Where exactly? - done!
% - slew time checken!  - done!
% - Gewichte checken! - done!





