% #########################################################################
% #     calc_scan_configurations
% #########################################################################
%
% DESCRIPTION
%   This function calculates all possible scan configurations for the defined station network, sources (satellites/quasars), time and observation conditions.
%
% CREATED  
%   2015-10-27     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING (external function calls)
%   - check_quasar_visibility_from_station
%   - spmcb2
%   - tle_propagation
%   - rv2radec
%
%
% SUB-ROUTINES defined within this file
%   - check_flux_at_baselines
%   - extend_scan_cons_struct
%
%
% INPUT
%   - stat_data                     : station data structure
%   - source                        : Source structure
%   - PARA                          : Global scheduling parameter strucutre
%   - t_epoch_jd                    : Scan start time
%   - stat_id_list                  : Station IDs of the observing station network
%   - scan_type                     : Scan type (sat. or quasar)
%   - obs_data                      : Observation data structure (preallocated in "find_observation_times.m")

%
%
% OUTPUT
%   - scan_cons                     : Structure containing all possible scan-configurations for the defined time and observation constellation
%   - error_code                    : Error Code (0 = no erros occured)
%   - error_msg                     : Error Message (empty, if no errors occured)
%
% CHANGES:
% - 2016-05-09: A. Hellerschmied: Separate source repetition interval treatment for satellites added (PARA.MIN_SRCRP_SAT).
% - 2016-12-22, A. Hellerschmied: - PARA.INIT_PROP_INTERVAL in [sec] instead of [min]
%

function [scan_cons, error_code, error_msg] = calc_scan_configurations(stat_data, source, PARA, t_epoch_jd, stat_id_list, scan_type, obs_data)


    % ##### Init #####
    error_code = 0;
    error_msg = '';
    
    number_of_stations = length(stat_id_list);
    number_of_quasars  = length(source);
    number_of_sat      = length(obs_data.sat);
    
    i_scan_cons = 0;
    
    % max_number_of_scan_cons = number_of_quasars; % + ....
    max_number_of_scans     = 2;
    
    % Initialize obs_data here, if no scan has been schedules so far?
    if obs_data.number_of_scans == 0
        obs_data.end_of_last_scan_jd = 99999;
    end
    
    
%     % ##### Preallocation #####
%     % obs_data:
%     
%     obs_data.quasars(number_of_quasars).last_obs_jd         = 0;
%     obs_data.quasars(number_of_quasars).number_of_scans     = 0;
%     for i = 1 : number_of_quasars
%        obs_data.quasars(i).last_obs_jd      = 0;
%        obs_data.quasars(i).number_of_scans  = 0;
%     end
%     
%     obs_data.sat(number_of_sat).last_obs_jd                 = 0;
%     obs_data.sat(number_of_sat).number_of_scans             = 0;
%     for i = 1 : number_of_sat
%        obs_data.sat(i).last_obs_jd          = 0;
%        obs_data.sat(i).number_of_scans      = 0;
%     end
    
    
    % scan_cons ("scan configurations")
    % scan_cons = struct('number_of_scans', 0, 'scan', repmat(struct('scan_type', '','quasar_id', [],'sat_id', [], 'scan_start_jd', [], 'scan_end_jd', [], 'flag_list_stat_obs', false(1, number_of_stations), 'stat', repmat(struct('begin_of_new_obs_temp', struct('un_az', [], 'el', [], 'jd', []) ,'end_of_new_obs_temp', struct('un_az', [], 'el', [], 'jd', [])), 1, number_of_stations)), 1, max_number_of_scans));
    % scan_cons = struct('number_of_scans', 0, 'scan', repmat(struct('scan_type', '','quasar_id', [],'sat_id', [], 'scan_start_jd', [], 'scan_end_jd', [], 'list_stat_id_obs', [], 'stat', repmat(struct('begin_of_new_obs_temp', struct('un_az', [], 'el', [], 'jd', 0) ,'end_of_new_obs_temp', struct('un_az', [], 'el', [], 'jd', 0), 'duration_sec', 0, 'slew_time_sec', []), 1, number_of_stations)), 1, max_number_of_scans));
    scan_cons = [];
    % ##### Options #####
    
    
    
    switch(scan_type)
        
        % #############################
        % ####     Quasar scan     ####
        % #############################
        case 'q'
            
            % ##### Prelesection of sources: Check the visibility of all sources at all stations #####
            
            % ### Preallocations ###
            quasars = repmat(struct('flag_list_stat_obs', false(1, number_of_stations)),1 , number_of_quasars);
            
            % ### Loop over all sources ###
            for i_q = 1 : number_of_quasars
                                
                % ### Loop over all stations ###
                for i_stat = 1 : number_of_stations
                    % Get the station ID:
                    stat_id = stat_id_list(i_stat);
                    % Check the source-visibility at scan start:
                    [quasars(i_q).flag_list_stat_obs(i_stat), error_code, error_msg] = check_quasar_visibility_from_station(stat_data.stat(stat_id), source(i_q), t_epoch_jd);
                    if error_code > 0
                        error_msg = ['check_quasar_visibility_from_station: ', error_msg];
                        return; 
                    end
                    
                end % for i_stat = 1 : number_of_stations
                
            end % for i_q = 1 : length(source))
        
            
            
            % ##### Get the possible scan-configuration with one source #####
            
            % #### Loop over all sources ####
            for i_q = 1 : number_of_quasars
                
                obs_stat_id_list = stat_id_list(quasars(i_q).flag_list_stat_obs);
                nsta = size(obs_stat_id_list, 2);
                
                % ### Check if a min. number of stations is able to observe the source ###
                % Min. 2 stations = one baseline
                if nsta < max(2, PARA.MIN_STANUM)
                    continue;
                end

                % ### Check if the source has been observed repeatedly within a defined interval ###
                if ((obs_data.end_of_last_scan_jd - obs_data.quasars(i_q).last_obs_jd) < PARA.MIN_SRCRP/1440)
                   continue; 
                end

                % ### Check if the projected flux is larger than the defined threshold at all baselines ###
                [flag_flux_ok] = check_flux_at_baselines(source, stat_data, obs_stat_id_list, i_q, PARA, t_epoch_jd);
                if ~flag_flux_ok
                    continue; % Check next source
                end
                
                % ### If all conditions are met, initialize the configuration in the scan_con structure ###
                i_scan_cons = i_scan_cons + 1;
                
                
                % Extend the scan_cons structre by another ...
                [scan_cons] = extend_scan_cons_struct(scan_cons, max_number_of_scans, number_of_stations);
                
                scan_cons(i_scan_cons).scan(1).scan_type = scan_type;
                scan_cons(i_scan_cons).scan(1).quasar_id = i_q;
                
                scan_cons(i_scan_cons).number_of_scans = 1;
                
                % Stations which are able to observe this source:
                scan_cons(i_scan_cons).scan(1).list_stat_id_obs = stat_id_list(quasars(i_q).flag_list_stat_obs);
                
% ADD FURTHER PARAMETERS HERE, IF REQUIRED!!!
% Don't forget to preallocate!!!
                
            end % for i_q = 1 : number_of_quasars
            
            
            % ##### Get the possible scan-configuration with two source/subnets #####
            % Subnetting is only possible, if there are enough stations in the defined network:
            %  => ("min. number of stations per subnet (max(2, PARA.MIN_STANUM))" * "number of subnets (2)" >= "number of stations")
            if (max(2, PARA.MIN_STANUM) * 2) <= number_of_stations

                % ##### Nested Loops over all sources #####
                for i_q_1 = 1 : (number_of_quasars - 1)
                    for i_q_2 = (i_q_1 + 1) : number_of_quasars

                        obs_stat_id_list_1 = stat_id_list(quasars(i_q_1).flag_list_stat_obs);
                        nsta_1 = size(obs_stat_id_list_1, 2);
                        obs_stat_id_list_2 = stat_id_list(quasars(i_q_2).flag_list_stat_obs);
                        nsta_2 = size(obs_stat_id_list_2, 2);


                        % ### Check if a min. number of stations is able to observe each source ###
                        % (a) Min. 2 stations per source, or at least the number of stations defined in PARA.MIN_STANUM
                        % At this stage it is not considered, if 2 simultanous scans to two sources with a min. number of stations each (2 or PARA.MIN_STANUM), are actually possible (one station can only observe one source at a time)!!! 
                        %   => This is done later on, after calculating all possible sub-configurations!!! 
                        %   => But nevertheless it is a prerequisit for 2 simultaneous scans, that the min. sub-net size is met!
                        if (nsta_1 < max(2, PARA.MIN_STANUM)) || (nsta_2 < max(2, PARA.MIN_STANUM))
                            continue;
                        end


                        % ### Check if one of the two sources has been observed repeatedly within a defined interval ###
                        if ((obs_data.end_of_last_scan_jd - obs_data.quasars(i_q_1).last_obs_jd) < PARA.MIN_SRCRP/1440) || ((obs_data.end_of_last_scan_jd - obs_data.quasars(i_q_2).last_obs_jd) < PARA.MIN_SRCRP/1440)
                           continue; 
                        end


                        % ### Check the angular distance between the two sources ###
                        ra1 = source(i_q_1).ra;
                        de1 = source(i_q_1).de;
                        ra2 = source(i_q_2).ra;
                        de2 = source(i_q_2).de;
                        arg = cos(de1) * cos(de2) * cos(ra1 - ra2) + sin(de1) * sin(de2);
                        srcd = acos(arg);
                        if (srcd < PARA.MIN_SRC2ANG)
                            continue;
                        end 


                        % ### Calculate all possible sub-configurations ###
                        [number_of_sub_cons, sub_cons] = spmcb2(nsta_1, obs_stat_id_list_1, nsta_2, obs_stat_id_list_2);

                        % ### Sort out ###
                        % Throw out sub-sonfigurations...
                        % Scans (two scans)
                        % - with less than 2 (or PARA.MIN_STANUM) stations => min sub-net size is not met!
                        % - 

                        % sub-cons:
                        % - where only one (of two) is possible => configurations with only one scan to a single source are already convered!
                        % - where no scan is possible

                        % ### Loop over all sub-configurations ###
                        for i_sub = 1 : number_of_sub_cons

                            % Loop init.:
                            sub_stat_id_list_1 = sub_cons(i_sub).downstasn1;
                            sub_nsta_1 = sub_cons(i_sub).downstanum1;
                            sub_stat_id_list_2 = sub_cons(i_sub).downstasn2;
                            sub_nsta_2 =sub_cons(i_sub).downstanum2;


                            % ### Min. sub-net site (source/subnet 1)? ###
                            if sub_nsta_1 < max(2, PARA.MIN_STANUM)
                                continue;
                            end

                            % ### Min. sub-net size (source/subnet 2)? ###
                            if sub_nsta_2 < max(2, PARA.MIN_STANUM)
                                continue;
                            end

                            % ### Check if the projected flux is larger than the defined threshold at all baselines (source/subnet 1) ###
                            [flag_flux_ok] = check_flux_at_baselines(source, stat_data, sub_stat_id_list_1, i_q_1, PARA, t_epoch_jd);
                            if ~flag_flux_ok
                                continue; % Check next source
                            end

                            % ### Check if the projected flux is larger than the defined threshold at all baselines (source/subnet 2) ###
                            [flag_flux_ok] = check_flux_at_baselines(source, stat_data, sub_stat_id_list_2, i_q_2, PARA, t_epoch_jd);
                            if ~flag_flux_ok
                                continue; % Check next source
                            end


                            % ##### If all requirements are met => Save this sub-configuration to "scan_cons" #####

                            % Extend the scan_cons structre by another ...
                            [scan_cons] = extend_scan_cons_struct(scan_cons, max_number_of_scans, number_of_stations);

                            scan_cons(i_scan_cons).number_of_scans = 2;

                            % ### Scan/source 1 ###
                            scan_cons(i_scan_cons).scan(1).scan_type = scan_type; 
                            scan_cons(i_scan_cons).scan(1).quasar_id = i_q_1;
                            % Stations which are able to observe this source:
%                             scan_cons(i_scan_cons).scan(1).flag_list_stat_obs(sub_stat_id_list_1) = true; 
                            scan_cons(i_scan_cons).scan(1).list_stat_id_obs = sub_stat_id_list_1;

                            % ### Scan/source 2 ###
                            scan_cons(i_scan_cons).scan(2).scan_type = scan_type;
                            scan_cons(i_scan_cons).scan(2).quasar_id = i_q_2;
                            % Stations which are able to observe this source:
%                             scan_cons(i_scan_cons).scan(2).flag_list_stat_obs(sub_stat_id_list_2) = true;
                            scan_cons(i_scan_cons).scan(2).list_stat_id_obs = sub_stat_id_list_2;

                        end % for i_sub = 1 : number_of_sub_cons

                    end % for i_q_2 = (i_q_1 + 1) : number_of_quasars

                end % for i_q_1 = 1 : (number_of_quasars - 1)
                
            end % if (max(2, PARA.MIN_STANUM) * 2) >= number_of_stations
    
            
            
    
        % ################################
        % ####     Satellite scan     ####
        % ################################
        case 's'
            
            % ##### Check which satellites are actually observable from which stations at the defined epoch #####
            
            % ### Preallocations ###
            sats = repmat(struct('flag_list_stat_obs', false(1, number_of_stations)),1 , obs_data.number_of_sat);
            
            % #### Loop over all satellites ####
            for i_sat = 1 : obs_data.number_of_sat
                
                % ### Loop over all stations ###
                for i_stat = 1 : number_of_stations
                     % Get the station ID:
                    stat_id = stat_id_list(i_stat);

                    % #### Loop over all possible observation periods ####
                    for i_obs_period = 1 : size(obs_data.sat(i_sat).stat(stat_id).obs_times, 1)
                        % Check visibility:
                        if (t_epoch_jd >= obs_data.sat(i_sat).stat(stat_id).obs_times(i_obs_period, 1)) && (t_epoch_jd <= obs_data.sat(i_sat).stat(stat_id).obs_times(i_obs_period, 2))
                            sats(i_sat).flag_list_stat_obs(i_stat) = true;
                        end

                    end % for i_obs_period = 1 : size(obs_data.sat(i_sat).obs_times, 1)
                end % for i_stat = 1 : number_of_stations
            end % for i_sat = 1 : number_of_satellites
            
            % ##### Get the possible scan-configuration with one satellite #####
            % #### Loop over all satellites ####
            for i_sat = 1 : obs_data.number_of_sat
                
                % Get stations which are able to observe the current satellite (i_sat):
                obs_stat_id_list = stat_id_list(sats(i_sat).flag_list_stat_obs);
                nsta = size(obs_stat_id_list, 2);
                
                % ### Check if a min. number of stations is able to observe the source ###
                % Min. 2 stations = one baseline
                if nsta < max(2, PARA.MIN_STANUM)
                    continue;
                end

                % ### Check if the source has been observed repeatedly within a defined interval ###
                if ((obs_data.end_of_last_scan_jd - obs_data.sat(i_sat).last_obs_jd) < PARA.MIN_SRCRP_SAT/1440)
                   continue;    
                end

                % ### Check if the projected flux is larger than the defined threshold at all baselines ###
                % Cannot be considered for satellites so far!
                
                % ### If all conditions are met, initialize the configuration in the scan_con structure ###
                i_scan_cons = i_scan_cons + 1;
                
                % Extend the scan_cons structre by another ...
                [scan_cons] = extend_scan_cons_struct(scan_cons, max_number_of_scans, number_of_stations);
                
                scan_cons(i_scan_cons).scan(1).scan_type = scan_type;
                scan_cons(i_scan_cons).scan(1).sat_id = i_sat;
                
                scan_cons(i_scan_cons).number_of_scans = 1;
                
                % Stations which are able to observe this source:
                scan_cons(i_scan_cons).scan(1).list_stat_id_obs = stat_id_list(sats(i_sat).flag_list_stat_obs);
                
% ADD FURTHER PARAMETERS HERE, IF REQUIRED!!!
% Don't forget to preallocate!!!
                
            end % for i_sat = 1 : obs_data.number_of_sat
            
                        
            % ##### Get the possible scan-configuration with two source/subnets #####
            % Subnetting is only possible, if there are enough stations in the defined network:
            %  => ("min. number of stations per subnet (max(2, PARA.MIN_STANUM))" * "number of subnets (2)" >= "number of stations")
            if (max(2, PARA.MIN_STANUM) * 2) <= number_of_stations

                % ##### Nested Loops over all sources #####
                for i_sat_1 = 1 : (obs_data.number_of_sat - 1)
                    for i_sat_2 = (i_sat_1 + 1) : obs_data.number_of_sat

                        obs_stat_id_list_1 = stat_id_list(sats(i_sat_1).flag_list_stat_obs);
                        nsta_1 = size(obs_stat_id_list_1, 2);
                        obs_stat_id_list_2 = stat_id_list(sats(i_sat_2).flag_list_stat_obs);
                        nsta_2 = size(obs_stat_id_list_2, 2);


                        % ### Check if a min. number of stations is able to observe each source ###
                        % (a) Min. 2 stations per source, or at least the number of stations defined in PARA.MIN_STANUM
                        % At this stage it is not considered, if 2 simultanous scans to two sources with a min. number of stations each (2 or PARA.MIN_STANUM), are actually possible (one station can only observe one source at a time)!!! 
                        %   => This is done later on, after calculating all possible sub-configurations!!! 
                        %   => But nevertheless it is a prerequisit for 2 simultaneous scans, that the min. sub-net size is met!
                        if (nsta_1 < max(2, PARA.MIN_STANUM)) || (nsta_2 < max(2, PARA.MIN_STANUM))
                            continue;
                        end


                        % ### Check if one of the two sources has been observed repeatedly within a defined interval ###
                        if ((obs_data.end_of_last_scan_jd - obs_data.sat(i_sat_1).last_obs_jd) < PARA.MIN_SRCRP_SAT/1440) || ((obs_data.end_of_last_scan_jd - obs_data.sat(i_sat_2).last_obs_jd) < PARA.MIN_SRCRP_SAT/1440)
                           continue; 
                        end
                       

                        % ### Check the angular distance between the two sources ###

                        % Satellite 1:
                        [temp_sat_data, error_code, error_msg] = tle_propagation(t_epoch_jd, t_epoch_jd, PARA.INIT_PROP_INTERVAL/60, PARA, obs_data.sat(i_sat_1).name);
                        if error_code > 0
                            error_msg = ['tle_propagation: ', error_msg];
                            return; 
                        end
                        [radius_1_km, ra_1_rad, de_1_rad] = rv2radec(temp_sat_data.sat(1).prop(1).r, temp_sat_data.sat(1).prop(1).v);
                        
                        % Satellite 2:
                        [temp_sat_data, error_code, error_msg] = tle_propagation(t_epoch_jd, t_epoch_jd, PARA.INIT_PROP_INTERVAL/60, PARA, obs_data.sat(i_sat_2).name);
                        if error_code > 0
                            error_msg = ['tle_propagation: ', error_msg];
                            return; 
                        end
                        [radius_2_km, ra_2_rad, de_2_rad] = rv2radec(temp_sat_data.sat(1).prop(1).r, temp_sat_data.sat(1).prop(1).v);
                        
                        % Calc. separation angle [rad]:
                        arg = cos(de_1_rad) * cos(de_2_rad) * cos(ra_1_rad - ra_2_rad) + sin(de_1_rad) * sin(de_2_rad);
                        sep_angle_rad = acos(arg);
                        
                        if (sep_angle_rad < PARA.MIN_SRC2ANG)
                            continue;
                        end 


                        % ### Calculate all possible sub-configurations ###
                        [number_of_sub_cons, sub_cons] = spmcb2(nsta_1, obs_stat_id_list_1, nsta_2, obs_stat_id_list_2);

                        % ### Sort out ###
                        % Throw out sub-sonfigurations...
                        % Scans (two scans)
                        % - with less than 2 (or PARA.MIN_STANUM) stations => min sub-net size is not met!
                        % - 

                        % sub-cons:
                        % - where only one (of two) is possible => configurations with only one scan to a single source are already convered!
                        % - where no scan is possible

                        % ### Loop over all sub-configurations ###
                        for i_sub = 1 : number_of_sub_cons

                            % Loop init.:
                            sub_stat_id_list_1 = sub_cons(i_sub).downstasn1;
                            sub_nsta_1 = sub_cons(i_sub).downstanum1;
                            sub_stat_id_list_2 = sub_cons(i_sub).downstasn2;
                            sub_nsta_2 =sub_cons(i_sub).downstanum2;


                            % ### Min. sub-net site (source/subnet 1)? ###
                            if sub_nsta_1 < max(2, PARA.MIN_STANUM)
                                continue;
                            end

                            % ### Min. sub-net size (source/subnet 2)? ###
                            if sub_nsta_2 < max(2, PARA.MIN_STANUM)
                                continue;
                            end

                            % ### Check if the projected flux is larger than the defined threshold at all baselines (source/subnet 1) ###
                            % This is not possible so far for satellites


                            % ##### If all requirements are met => Save this sub-configuration to "scan_cons" #####

                            % Extend the scan_cons structre by another ...
                            [scan_cons] = extend_scan_cons_struct(scan_cons, max_number_of_scans, number_of_stations);

                            scan_cons(i_scan_cons).number_of_scans = 2;

                            % ### Scan/source 1 ###
                            scan_cons(i_scan_cons).scan(1).scan_type = scan_type; 
                            scan_cons(i_scan_cons).scan(1).sat_id = i_sat_1;
                            % Stations which are able to observe this source:
                            scan_cons(i_scan_cons).scan(1).list_stat_id_obs = sub_stat_id_list_1;

                            % ### Scan/source 2 ###
                            scan_cons(i_scan_cons).scan(2).scan_type = scan_type;
                            scan_cons(i_scan_cons).scan(2).sat_id = i_sat_2;
                            % Stations which are able to observe this source:
                            scan_cons(i_scan_cons).scan(2).list_stat_id_obs = sub_stat_id_list_2;

                        end % for i_sub = 1 : number_of_sub_cons

                    end % i_sat_2 = (i_sat_1 + 1) : obs_data.number_of_sat

                end % for i_sat_1 = 1 : (obs_data.number_of_sat - 1)
                
            end % if (max(2, PARA.MIN_STANUM) * 2) >= number_of_stations
            
            
            
    end % switch(scan_type)
    
end % function


%% #################### Sub routines #########################

function [flag_flux_ok] = check_flux_at_baselines(source, stat_data, obs_stat_id_list, quasar_id, PARA, t_epoch_jd)

    % Init.:
    flag_flux_ok = true;
    nsta = size(obs_stat_id_list, 2);
    flag_flux_break = false;
    fluxpara(1 : PARA.MAX_BANDNUM, 1 : PARA.MAX_FLUXPARA) = source(quasar_id).fluxpara(1 : PARA.MAX_BANDNUM, 1 : PARA.MAX_FLUXPARA);

    % Loop over all baselines:
    for i_bl_1 = 1 : (nsta - 1)
        for i_bl_2 = (i_bl_1 + 1) : nsta   

            % Get station IDs of the baseline:
            stat_id_1 = obs_stat_id_list(i_bl_1);
            stat_id_2 = obs_stat_id_list(i_bl_2);

            % Calc. coordinate differences:
            blx = stat_data.stat(stat_id_1).location.TRF.x - stat_data.stat(stat_id_2).location.TRF.x;
            bly = stat_data.stat(stat_id_1).location.TRF.y - stat_data.stat(stat_id_2).location.TRF.y;
            blz = stat_data.stat(stat_id_1).location.TRF.z - stat_data.stat(stat_id_2).location.TRF.z;

            % Calc. projected flux:
            [obsflux] = sobsflux(t_epoch_jd, blx, bly, blz, source(quasar_id).ra, source(quasar_id).de, fluxpara, PARA);

            % Observed flux large enough?
            if (min(obsflux(1:PARA.MAX_BANDNUM)) < PARA.MIN_FLUX)  
                flag_flux_break = true;
                flag_flux_ok = false;
                break;
            end
        end
        if flag_flux_break
            break;
        end
    end
end

function [scan_cons] = extend_scan_cons_struct(scan_cons, max_number_of_scans, number_of_stations)
    % Update the structure
    if ~isempty(scan_cons)
        scan_cons_index = length(scan_cons) + 1;
    else
        scan_cons_index = 1;
        clear scan_cons;
    end
    scan_cons(scan_cons_index) = struct('weight', [], 'number_of_scans', 0, 'scan', repmat(struct('scan_type', '','quasar_id', [],'sat_id', [], 'scan_start_jd', [], 'scan_end_jd', [], 'list_stat_id_obs', [], 'stat', repmat(struct('begin_of_new_obs_temp', struct('un_az', [], 'el', [], 'jd', 0) ,'end_of_new_obs_temp', struct('un_az', [], 'el', [], 'jd', 0), 'duration_sec', 0, 'slew_time_sec', [], 'sky_cov_num', []), 1, number_of_stations)), 1, max_number_of_scans));
end


