% #########################################################################
% #     get_station_statistics
% #########################################################################
%
% DESCRIPTION
%   This function calculates some statistical values for each station. 
%   E.g. it finds all observed sources of each station in the "sched_data" structure.
%
%   It adds the following fields to:
%    - sched_data.stat(i_stat).statistics.sat_name_list                : names of observed satellites (string, 24 char long) (unique!)
%                                         .sat_norad_id_vector          : vector (col) with NORAD ID of observed satellites (unique!)
%                                         .quasar_name_list             : names of all observed quasars (unique!)
%
% CREATED  
%   2015-06-19     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
% - 
%
%
% INPUT
% - sched_data          - Scheduling data structure
% - t_min_jd            - Lower limit of time window (optional)
% - t_max_jd            - Upper limit of time window (optional)
%
%
% OUTPUT
% - error_code          - Error Code (0 = no erros occured)
% - error_msg           - Error Message (empty, if no errors occured)
% - sched_data          - Scheduling data structure
%
% CHANGES:
% - 2015-11-26, A. Hellerschmied: Additional statistics added (num_of_sat_scans, num_of_quasar_scans, num_of_all_scans, sat_id_vector, quasar_id_vector)
% - 2016-04-28, A. Hellerschmied: Added possibility to define time limits for calculating the statistics (optional input argument)
% - 2016-11-15, A. Hellerschmied: Additionally get statistics for whole session


function [sched_data, error_code, error_msg] = get_station_statistics(sched_data, varargin)

    % Init.:
    error_code = 0;
    error_msg = '';  
    
    sat_norad_id_vector_session = [];
    sat_name_list_session       = {};
    quasar_id_vector_session    = [];
    quasar_name_list_session    = {};
    sat_id_vector_session       = [];
    
    % ##### Treat optinal input arguments #####
    switch(nargin)
        
        case 1
            flag_check_time_window = false;
            
        case 3
            flag_check_time_window = true;
            t_min_jd   = varargin{1};
            t_max_jd     = varargin{2};
            
        otherwise
            error_code = 1;
            error_msg = 'Invalid number of input arguments';
            return;
    end

    % ##### Loop over all stations #####
    for i_stat = 1 : length(sched_data.stat)
        
        % loop init.:
        sat_name_list = {};                 % List of the names of all satellites which are observed by this station 
        quasar_name_list = {};              % List of the names of all satellites which are observed by this station 
        sat_norad_id_vector = [];           % Vector of the NORAD IDs of all satellites which are observed by this station 
        sat_id_vector = [];                 % Vector of IDs of all satellites which are observed by this station
        quasar_id_vector = [];               % Vector of IDs of all quasars which are observed by this station
        sat_count = 0;
        quasar_count = 0;
        num_of_sat_scans = 0;
        num_of_quasar_scans = 0;


        % ##### Finde all different sources which are observed in this session by the specific antenna #####

        % ##### loop over all scans #####
        for i_scan = 1 : sched_data.number_of_scans
            if flag_check_time_window
                if ((sched_data.scan(i_scan).t_start_jd < t_min_jd) || (sched_data.scan(i_scan).t_end_jd > t_max_jd))
                    continue; 
                end
            end

            % Loop init.:
            flag_stat_joins_scan = 0;

            % #### Select scans, where the station participates ####
            % loop over all stations of the scan:
            for i_stat_2 = 1 : sched_data.scan(i_scan).number_of_stat
                if sched_data.scan(i_scan).stat(i_stat_2).stat_id == i_stat
                    % Set flag, if the station is included in the scan
                    flag_stat_joins_scan = 1;
                    break;
                end
            end % for i_stat_2 = 1 : sched_data.scan(i_scan).number_of_stat


            % ##### If the station joins this scan... #####
            if flag_stat_joins_scan


                % #### Finde all observed source ####

                % ### Distinguish between observation-type ###
                switch(sched_data.scan(i_scan).obs_type)

                    % satellite scan:
                    case 'sat'                    
                        
                        sat_count = sat_count + 1;
                        
                        % save satellite name to list:
                        sat_name_list{sat_count} = sched_data.scan(i_scan).sat_name;

                        % save NORAD ID to vector:
                        sat_norad_id_vector = [sat_norad_id_vector; sched_data.scan(i_scan).sat_number];
                        
                        % save satellite ID
                        sat_id_vector = [sat_id_vector; sched_data.scan(i_scan).sat_id];

                    % quasar scan:
                    case 'quasar'
                        
                        quasar_count = quasar_count + 1;
                        
                        % save quasar name to list:
                        quasar_name_list{quasar_count} = sched_data.scan(i_scan).quasar_name;
                        
                        % save quasar ID
                        quasar_id_vector = [quasar_id_vector; sched_data.scan(i_scan).quasar_id];
                        
                end % switch(sched_data.scan(i_scan).obs_type)

            end % if flag_stat_joins_scan

        end % for i_scan = 1 : sched_data.number_of_scans


        % ##### Find unique names of observed satellites #####
        if ~isempty(sat_name_list)
            % Set flag:
%             flag_sat_scan_included = 1;

            % Total number of satellite scans:
            num_of_sat_scans = length(sat_id_vector);
            % Throw out multiple entries in the lists:
            [sat_name_list, m, n] = unique(sat_name_list);
            sat_norad_id_vector = sat_norad_id_vector(m,:);
            sat_id_vector = sat_id_vector(m,:);
        end
        
        % ##### Find unique names of observed quasars #####
        if ~isempty(quasar_name_list)
            % Set flag:
%             flag_quasar_scan_inlcuded = 1;

            % Total number of satellite scans:
            num_of_quasar_scans = length(quasar_id_vector);
            
            % Throw out multiple entries in the lists:
            [quasar_name_list, m, n] = unique(quasar_name_list);
            quasar_id_vector = quasar_id_vector(m,:);
            
        end
        
        
        % #### Assign output variables ####
        sched_data.stat(i_stat).statistics.sat_name_list               = sat_name_list;
        sched_data.stat(i_stat).statistics.sat_norad_id_vector         = sat_norad_id_vector;
        sched_data.stat(i_stat).statistics.quasar_name_list            = quasar_name_list;
        sched_data.stat(i_stat).statistics.quasar_id_vector            = quasar_id_vector;
        sched_data.stat(i_stat).statistics.sat_id_vector               = sat_id_vector;
        sched_data.stat(i_stat).statistics.num_of_sat_scans            = num_of_sat_scans;
        sched_data.stat(i_stat).statistics.num_of_quasar_scans         = num_of_quasar_scans;
        sched_data.stat(i_stat).statistics.num_of_all_scans            = num_of_quasar_scans + num_of_sat_scans;
        
        % #### Collect infos for complete session ####
        sat_norad_id_vector_session     = [sat_norad_id_vector_session; sched_data.stat(i_stat).statistics.sat_norad_id_vector];
        sat_name_list_session           = [sat_name_list_session, sched_data.stat(i_stat).statistics.sat_name_list];
        sat_id_vector_session           = [sat_id_vector_session; sched_data.stat(i_stat).statistics.sat_id_vector];
        quasar_id_vector_session        = [quasar_id_vector_session; sched_data.stat(i_stat).statistics.quasar_id_vector];
        quasar_name_list_session        = [quasar_name_list_session, sched_data.stat(i_stat).statistics.quasar_name_list];

    end % for i_stat = 1 : length(sched_data.stat)
    
    % ##### Get session statistics #####
    [sat_name_list_session, m, n]       = unique(sat_name_list_session);
    sat_norad_id_vector_session         = sat_norad_id_vector_session(m,:);
    sat_id_vector_session               = sat_id_vector_session(m);
    [quasar_id_vector_session, m, n]    = unique(quasar_id_vector_session);
    quasar_name_list_session            = quasar_name_list_session(m);
    number_of_quasars_session           = length(quasar_id_vector_session);
    number_of_sat_session               = length(sat_id_vector_session);
    
    % Save data:
    sched_data.exper_statistics.sat_name_list           = sat_name_list_session;
    sched_data.exper_statistics.sat_norad_id_vector     = sat_norad_id_vector_session;
    sched_data.exper_statistics.sat_id_vector           = sat_id_vector_session;
    sched_data.exper_statistics.quasar_name_list        = quasar_name_list_session;
    sched_data.exper_statistics.quasar_id_vector        = quasar_id_vector_session;
    sched_data.exper_statistics.number_of_quasars       = number_of_quasars_session;
    sched_data.exper_statistics.number_of_sat           = number_of_sat_session;

return;