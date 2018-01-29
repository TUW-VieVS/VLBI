% -------------------------------------------------------------------------
%
%                   function find_observation_times.m
%
%   This function finds the possible observation periods for a given number
%   of satellites and observing stations. Sun distance, simultanious 
%   visibility from stations (cut-off elevation & hizozon mask), antenna slew 
%   rate limits and slew range limits are checked.
%
%   Author: 
%       Andreas Hellerschmied, 2013.10.27
%   
%   changes       :
%   - 2014-01.19: A. Hellerschmied: Support of different antenna axis types
%   - 2014-01.19: A. Hellerschmied: Error treatment added, but no error
%               checks done yet!
%   - 2015-02-16: A. Hellerschmied: Remanded (before: "TLE_find_observation_times.m")
%   - 2015-07-08: A. Hellerschmied: obs_data preallocation: fields added
%   - 2015-08-18: A. Hellerschmied: check for slew range limits added
%   - 2015-10-28: A. Hellerschmied: Preallocation of "obs_data" updated. 
%   - 2015-10-30: A. Hellerschmied: Calc. possible observation time windows for individual stations. This infos is stored in "obs_data.sat.stat.obs_times". 
%   - 2015-10-30: A. Hellerschmied: Field for sky-coverage number added to obs_data strucuture
%   - 2016-05-09: A. Hellerschmied: Field for sky-coverage number added to obs_data strucuture ("sky_cov_sat")
%   - 2016-08-02: A. Hellerschmied: Round start/stop epochs of available observation windows to integer seconds (with round_jd2integer_sec.m)
%   - 2016-11-08: A. Hellerschmied: Bug-fix: axis 2 slew range limit was not considered correctly in the event matrix 
%           
%
%   inputs        :
%   - stat_data         : station data structure
%   - PARA              : VieVS parameter vector 
%     
%
%   outputs         :
%   - obs_data      : observation data structure
%   - error_code
%   - error_msg
%    
%
%   locals        :
% 
%
%   coupling      :
%   - round_jd2integer_sec
%   
%
%   references    :
%
%-------------------------------------------------------------------------

function [obs_data, error_code, error_msg] = find_observation_times(stat_data, PARA)

    % Init:
    error_code = 0;
    error_msg = '';

    % preallocation:
    obs_data = struct('sat', [], 'stat', [], 'number_of_sat', [], 'number_of_stat', [], 'end_of_last_scan_jd', [], 'number_of_scans', [], 'quasars', []);
    obs_data.sat = struct('name', [], 'state_vector', [], 'event_matrix', [], 'obs_times', [], 'last_obs_jd', [], 'number_of_scans', [], 'stat', []);
    obs_data.sat.stat = struct('obs_times', []);
    obs_data.stat = struct('name', [], 'end_of_last_obs', [], 'begin_of_new_obs_temp', [], 'end_of_new_obs_temp', [], 'sky_cov', [], 'sky_cov_sat', []);
    obs_data.stat.end_of_last_obs = struct('az', [], 'el', [], 'ha', [], 'dc', [], 'jd', [], 'un_az', [], 'quasar_id', [], 'sat_id', []);
    obs_data.stat.begin_of_new_obs_temp = struct('un_az', [], 'el', [], 'jd', []);
    obs_data.stat.end_of_new_obs_temp = struct('un_az', [], 'el', [], 'jd', []);
    obs_data.end_of_last_scan_jd = [];
    obs_data.quasars = struct('last_obs_jd', [], 'number_of_scans', []);
    obs_data.number_of_scans = 0;

    % Init.:
    for i_stat = 1 : stat_data.number_of_stations
        
       obs_data.stat(i_stat).name = stat_data.stat(i_stat).name;
       
       obs_data.stat(i_stat).end_of_last_obs.az = [];
       obs_data.stat(i_stat).end_of_last_obs.el = [];
       obs_data.stat(i_stat).end_of_last_obs.ha = [];
       obs_data.stat(i_stat).end_of_last_obs.dc = [];
       obs_data.stat(i_stat).end_of_last_obs.un_az = [];
       obs_data.stat(i_stat).end_of_last_obs.jd = PARA.startmjd + 2400000.5;
       
       % temp:
       obs_data.stat(i_stat).begin_of_new_obs_temp.un_az = [];
       obs_data.stat(i_stat).end_of_new_obs_temp.un_az = [];
       obs_data.stat(i_stat).begin_of_new_obs_temp.el = [];
       obs_data.stat(i_stat).begin_of_new_obs_temp.jd = [];
       obs_data.stat(i_stat).end_of_new_obs_temp.el = [];
       obs_data.stat(i_stat).end_of_new_obs_temp.jd = [];
       
       % Sky coverage numbers
       obs_data.stat(i_stat).sky_cov     = zeros(13, 1);
       obs_data.stat(i_stat).sky_cov_sat = zeros(13, 1);
       
       switch(stat_data.stat(i_stat).axis_type)
           
           case 'AZEL'
               obs_data.stat(i_stat).end_of_last_obs.az = stat_data.stat(i_stat).lim11 + (stat_data.stat(i_stat).lim12 - stat_data.stat(i_stat).lim11) / 2;
               obs_data.stat(i_stat).end_of_last_obs.el = stat_data.stat(i_stat).lim21 + (stat_data.stat(i_stat).lim22 - stat_data.stat(i_stat).lim21) / 2;
               obs_data.stat(i_stat).end_of_last_obs.un_az = stat_data.stat(i_stat).lim11 + (stat_data.stat(i_stat).lim12 - stat_data.stat(i_stat).lim11) / 2;
               
%                if strcmp(stat_data.stat(i_stat).name, 'WETTDBBC') || strcmp(stat_data.stat(i_stat).name, 'WETTZELL')
%                    obs_data.stat(i_stat).end_of_last_obs.un_az = 630/180*pi ;
%                end
%                if strcmp(stat_data.stat(i_stat).name, 'ONSALA60')
%                    obs_data.stat(i_stat).end_of_last_obs.un_az = 630/180*pi ;
%                end
               
           case 'HADC'
               obs_data.stat(i_stat).end_of_last_obs.ha = stat_data.stat(i_stat).lim11 + (stat_data.stat(i_stat).lim12 - stat_data.stat(i_stat).lim11) / 2;
               obs_data.stat(i_stat).end_of_last_obs.dc = stat_data.stat(i_stat).lim21 + (stat_data.stat(i_stat).lim22 - stat_data.stat(i_stat).lim21) / 2;
               obs_data.stat(i_stat).end_of_last_obs.el = 45 * pi/180; % Dummy!!!!...
               obs_data.stat(i_stat).end_of_last_obs.un_az = 0; % Dummy!!!!...
                              
           case 'XYEW'
               obs_data.stat(i_stat).end_of_last_obs.az = stat_data.stat(i_stat).lim11 + (stat_data.stat(i_stat).lim12 - stat_data.stat(i_stat).lim11) / 2;
               obs_data.stat(i_stat).end_of_last_obs.el = stat_data.stat(i_stat).lim21 + (stat_data.stat(i_stat).lim22 - stat_data.stat(i_stat).lim21) / 2;
               obs_data.stat(i_stat).end_of_last_obs.un_az = stat_data.stat(i_stat).lim11 + (stat_data.stat(i_stat).lim12 - stat_data.stat(i_stat).lim11) / 2; % Not important for this antenna type...
               
                      
       end
       
    end
    
    obs_data.number_of_sat = stat_data.number_of_sat;
    
    
    
    
    % ##### Find observation times #####
    
    % number_of_columns = 1 + ((stat_data.number_of_stations - stat_data.number_of_stations_quasar) * 8);
    number_of_columns = 1 + ((stat_data.number_of_stations - stat_data.number_of_stations_quasar) * 10);
    time_tag_count = 0;

     for i_sat = 1 : stat_data.number_of_sat
        
        % loop init.:
        i_row_obs_times = 0;
        flag_obs_possible = 0;
        i_row_obs_times_stat = zeros(1, stat_data.number_of_stations - stat_data.number_of_stations_quasar);
        flag_obs_possible_stat = false(1, stat_data.number_of_stations - stat_data.number_of_stations_quasar);
        obs_data.sat(i_sat).name = stat_data.stat(1).sat(i_sat).TLE_data.TLE_header_line;
        obs_data.sat(i_sat).stat(stat_data.number_of_stations - stat_data.number_of_stations_quasar).obs_times = [];

        
        
        % #### Get Event Matrix for each Satellite ####
        
        for i_stat = 1 : (stat_data.number_of_stations - stat_data.number_of_stations_quasar)
            
            % ---- Check first considered epoch of the dataset for any events ----
            if (    (stat_data.stat(i_stat).sat(i_sat).epoch(1).above_min_elevation == 1)   || ...
                    (stat_data.stat(i_stat).sat(i_sat).epoch(1).exceed_min_sun_dist == 1)   || ...
                    (stat_data.stat(i_stat).sat(i_sat).epoch(1).exceed_max_axis1_rate == 1)     || ...
                    (stat_data.stat(i_stat).sat(i_sat).epoch(1).exceed_max_axis2_rate == 1)     || ...
                    (stat_data.stat(i_stat).sat(i_sat).epoch(1).exceed_axis_limits == 1)             )
                
                temp_line = zeros(1,number_of_columns);
                temp_line(1,1) = stat_data.stat(i_stat).sat(i_sat).epoch(1).jd;
                    
                % Is satellite up at the first considered epoch? Yes => rise!
                if (stat_data.stat(i_stat).sat(i_sat).epoch(1).above_min_elevation == 1)
                    temp_line(1,(1 + (i_stat - 1) * 10) + 1) = 1;
                end

                % Is satellite exceeding min. sun dist. limit at the first considered epoch? Yes => exceed_min_sun_dist!
                if (stat_data.stat(i_stat).sat(i_sat).epoch(1).exceed_min_sun_dist == 1)
                    temp_line(1,(1 + (i_stat - 1) * 10) + 4) = 1;
                end
                
                % Is satellite exceeding max. axis 1 rate limit at the first considered epoch? Yes => exceed_max_axis1_rate!
                if (stat_data.stat(i_stat).sat(i_sat).epoch(1).exceed_max_axis1_rate == 1)
                    temp_line(1,(1 + (i_stat - 1) * 10) + 6) = 1;
                end

                % Is satellite exceeding max. axis 2 rate limit at the first considered epoch? Yes => exceed_max_axis2_rate!
                if (stat_data.stat(i_stat).sat(i_sat).epoch(1).exceed_max_axis2_rate == 1)
                    temp_line(1,(1 + (i_stat - 1) * 10) + 8) = 1;
                end

                % Is satellite exceeding any slew range limit at the first considered epoch? Yes => exceed_slew_range_limits!
                if (stat_data.stat(i_stat).sat(i_sat).epoch(1).exceed_axis_limits == 1)
                    temp_line(1,(1 + (i_stat - 1) * 10) + 10) = 1;
                end
                
                obs_data.sat(i_sat).event_matrix = [obs_data.sat(i_sat).event_matrix; temp_line];
                
            end
            
            
            

            % ---- Check Overpass data ----

            % Is satellite up at the last considered epoch? Yes => set!
            if (stat_data.stat(i_stat).sat(i_sat).epoch(stat_data.number_of_epochs).above_min_elevation == 1)
                temp_line = zeros(1,number_of_columns);
                temp_line(1,1) = stat_data.stat(i_stat).sat(i_sat).epoch(stat_data.number_of_epochs).jd;
                temp_line(1,(1 + (i_stat - 1) * 10) + 2) = 1;
                obs_data.sat(i_sat).event_matrix = [obs_data.sat(i_sat).event_matrix; temp_line];
            end
            
            % Go through overpass data:
            for i_over = 1 : stat_data.stat(i_stat).sat(i_sat).number_of_overpasses
                % rise:
                if ~isempty(stat_data.stat(i_stat).sat(i_sat).overpass(i_over).rise.jd)
                    temp_line = zeros(1,number_of_columns);
                    temp_line(1,1) = stat_data.stat(i_stat).sat(i_sat).overpass(i_over).rise.jd;
                    temp_line(1,(1 + (i_stat - 1) * 10) + 1) = 1;
                    obs_data.sat(i_sat).event_matrix = [obs_data.sat(i_sat).event_matrix; temp_line];
                end
                % set:
                if ~isempty(stat_data.stat(i_stat).sat(i_sat).overpass(i_over).set.jd)
                    temp_line = zeros(1,number_of_columns);
                    temp_line(1,1) = stat_data.stat(i_stat).sat(i_sat).overpass(i_over).set.jd;
                    temp_line(1,(1 + (i_stat - 1) * 10) + 2) = 1;
                    obs_data.sat(i_sat).event_matrix = [obs_data.sat(i_sat).event_matrix; temp_line];
                end
                
            end
            
            
            % ---- Check sun dist. data ---- 

            for i_sun = 1 : stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_keep_min_sundist
                if ~isempty(stat_data.stat(i_stat).sat(i_sat).sun_dist.keep_min_sundist(i_sun).jd)
                    temp_line = zeros(1,number_of_columns);
                    temp_line(1,1) = stat_data.stat(i_stat).sat(i_sat).sun_dist.keep_min_sundist(i_sun).jd;
                    temp_line(1,(1 + (i_stat - 1) * 10) + 3) = 1;
                    obs_data.sat(i_sat).event_matrix = [obs_data.sat(i_sat).event_matrix; temp_line];
                end
            end
            
            for i_sun = 1 : stat_data.stat(i_stat).sat(i_sat).sun_dist.number_of_exceed_min_sundist
                if ~isempty(stat_data.stat(i_stat).sat(i_sat).sun_dist.exceed_min_sundist(i_sun).jd)
                    temp_line = zeros(1,number_of_columns);
                    temp_line(1,1) = stat_data.stat(i_stat).sat(i_sat).sun_dist.exceed_min_sundist(i_sun).jd;
                    temp_line(1,(1 + (i_stat - 1) * 10) + 4) = 1;
                    obs_data.sat(i_sat).event_matrix = [obs_data.sat(i_sat).event_matrix; temp_line];
                end
            end
            
            
            % ---- Check Axis 1 rate data ----
            
            for i_axis1 = 1 : stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_keep_max_axis1_rate
                if ~isempty(stat_data.stat(i_stat).sat(i_sat).axis1_rate.keep_max_axis1_rate(i_axis1).jd)
                    temp_line = zeros(1,number_of_columns);
                    temp_line(1,1) = stat_data.stat(i_stat).sat(i_sat).axis1_rate.keep_max_axis1_rate(i_axis1).jd;
                    temp_line(1,(1 + (i_stat - 1) * 10) + 5) = 1;
                    obs_data.sat(i_sat).event_matrix = [obs_data.sat(i_sat).event_matrix; temp_line];
                end
            end
            
            for i_axis1 = 1 : stat_data.stat(i_stat).sat(i_sat).axis1_rate.number_of_exceed_max_axis1_rate
                if ~isempty(stat_data.stat(i_stat).sat(i_sat).axis1_rate.exceed_max_axis1_rate(i_axis1).jd)
                    temp_line = zeros(1,number_of_columns);
                    temp_line(1,1) = stat_data.stat(i_stat).sat(i_sat).axis1_rate.exceed_max_axis1_rate(i_axis1).jd;
                    temp_line(1,(1 + (i_stat - 1) * 10) + 6) = 1;
                    obs_data.sat(i_sat).event_matrix = [obs_data.sat(i_sat).event_matrix; temp_line];
                end
            end
            
            
            % ---- Check Axis 2 rate data ---- 
            
            for i_axis2 = 1 : stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_keep_max_axis2_rate
                if ~isempty(stat_data.stat(i_stat).sat(i_sat).axis2_rate.keep_max_axis2_rate(i_axis2).jd)
                    temp_line = zeros(1,number_of_columns);
                    temp_line(1,1) = stat_data.stat(i_stat).sat(i_sat).axis2_rate.keep_max_axis2_rate(i_axis2).jd;
                    temp_line(1,(1 + (i_stat - 1) * 10) + 7) = 1;
                    obs_data.sat(i_sat).event_matrix = [obs_data.sat(i_sat).event_matrix; temp_line];
                end
            end
            
            for i_axis2 = 1 : stat_data.stat(i_stat).sat(i_sat).axis2_rate.number_of_exceed_max_axis2_rate
                if ~isempty(stat_data.stat(i_stat).sat(i_sat).axis2_rate.exceed_max_axis2_rate(i_axis2).jd)
                    temp_line = zeros(1,number_of_columns);
                    temp_line(1,1) = stat_data.stat(i_stat).sat(i_sat).axis2_rate.exceed_max_axis2_rate(i_axis2).jd;
                    temp_line(1,(1 + (i_stat - 1) * 10) + 8) = 1;
                    obs_data.sat(i_sat).event_matrix = [obs_data.sat(i_sat).event_matrix; temp_line];
                end
            end
            

            % ---- Check axis limits data ---- 

            for i_axis_lim = 1 : stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_keep_slew_range_limits
                if ~isempty(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.keep_slew_range_limits(i_axis_lim).jd)
                    temp_line = zeros(1,number_of_columns);
                    temp_line(1,1) = stat_data.stat(i_stat).sat(i_sat).slew_range_limits.keep_slew_range_limits(i_axis_lim).jd;
                    temp_line(1,(1 + (i_stat - 1) * 10) + 9) = 1;
                    obs_data.sat(i_sat).event_matrix = [obs_data.sat(i_sat).event_matrix; temp_line];
                end
            end
            
            for i_axis_lim = 1 : stat_data.stat(i_stat).sat(i_sat).slew_range_limits.number_of_exceed_slew_range_limits
                if ~isempty(stat_data.stat(i_stat).sat(i_sat).slew_range_limits.exceed_slew_range_limits(i_axis_lim).jd)
                    temp_line = zeros(1,number_of_columns);
                    temp_line(1,1) = stat_data.stat(i_stat).sat(i_sat).slew_range_limits.exceed_slew_range_limits(i_axis_lim).jd;
                    temp_line(1,(1 + (i_stat - 1) * 10) + 10) = 1;
                    obs_data.sat(i_sat).event_matrix = [obs_data.sat(i_sat).event_matrix; temp_line];
                end
            end
            

        end % for i_stat = 1 : (stat_data.number_of_stations - stat_data.number_of_stations_quasar)
        
        
        % sort event Matrix based on the entries of the first row (JD).
        if ~isempty(obs_data.sat(i_sat).event_matrix)
            obs_data.sat(i_sat).event_matrix = sortrows(obs_data.sat(i_sat).event_matrix, 1);
        end
        
        
        
        
        % #### Get Observation Times for each Satellite ####
        
        rows = size(obs_data.sat(i_sat).event_matrix, 1);
        
        % Set up initial state vector:
        obs_data.sat(i_sat).state_vector = ones(1, (5 * (stat_data.number_of_stations - stat_data.number_of_stations_quasar)));
        
        for i_stat_2 = 1 : (stat_data.number_of_stations - stat_data.number_of_stations_quasar)
            % Satellite below min. elevation:
            obs_data.sat(i_sat).state_vector(1, ((i_stat_2 - 1) * 5) + 1) = 0; 
        end
        
        % Go through each row of the Event Matrix:
        for i_row = 1 : rows
            
            % Go through each element of the state vector row:
            for i_stat_2 = 1 : (stat_data.number_of_stations - stat_data.number_of_stations_quasar)
                
                % Set state vector elements according to the entries in the
                % event matrix:
                
                % -- Overpass: --
                % rise:
                if (obs_data.sat(i_sat).event_matrix(i_row, 1 + ((i_stat_2 - 1) * 10) + 1) == 1)
                    obs_data.sat(i_sat).state_vector(1, ((i_stat_2 - 1) * 5) + 1) = 1;
                end
                % set:
                if (obs_data.sat(i_sat).event_matrix(i_row, 1 + ((i_stat_2 - 1) * 10) + 2) == 1)
                    obs_data.sat(i_sat).state_vector(1, ((i_stat_2 - 1) * 5) + 1) = 0;
                end
                
                % Sun dist.
                % keep min. sun dist.:
                if (obs_data.sat(i_sat).event_matrix(i_row, 1 + ((i_stat_2 - 1) * 10) + 3) == 1)
                    obs_data.sat(i_sat).state_vector(1, ((i_stat_2 - 1) * 5) + 2) = 1;
                end
                % exceed min. sun dist.:
                if (obs_data.sat(i_sat).event_matrix(i_row, 1 + ((i_stat_2 - 1) * 10) + 4) == 1)
                    obs_data.sat(i_sat).state_vector(1, ((i_stat_2 - 1) * 5) + 2) = 0;
                end
                
                % Axis 1 rate:
                % keep max. axis 1 rate limit:
                if (obs_data.sat(i_sat).event_matrix(i_row, 1 + ((i_stat_2 - 1) * 10) + 5) == 1)
                    obs_data.sat(i_sat).state_vector(1, ((i_stat_2 - 1) * 5) + 3) = 1;
                end
                % exceed max. axis 1 rate limit:
                if (obs_data.sat(i_sat).event_matrix(i_row, 1 + ((i_stat_2 - 1) * 10) + 6) == 1)
                    obs_data.sat(i_sat).state_vector(1, ((i_stat_2 - 1) * 5) + 3) = 0;
                end
                
                % Axis 2 rate:
                % keep max. axis 2 rate limit:
                if (obs_data.sat(i_sat).event_matrix(i_row, 1 + ((i_stat_2 - 1) * 10) + 7) == 1)
                    obs_data.sat(i_sat).state_vector(1, ((i_stat_2 - 1) * 5) + 4) = 1;
                end
                % exceed max.axis 2 rate limit:
                if (obs_data.sat(i_sat).event_matrix(i_row, 1 + ((i_stat_2 - 1) * 10) + 8) == 1)
                    obs_data.sat(i_sat).state_vector(1, ((i_stat_2 - 1) * 5) + 4) = 0;
                end
                
                % Antenna slew range limits:
                % keep slew range limits:
                if (obs_data.sat(i_sat).event_matrix(i_row, 1 + ((i_stat_2 - 1) * 10) + 9) == 1)
                    obs_data.sat(i_sat).state_vector(1, ((i_stat_2 - 1) * 5) + 5) = 1;
                end
                % exceed slew range limits:
                if (obs_data.sat(i_sat).event_matrix(i_row, 1 + ((i_stat_2 - 1) * 10) + 10) == 1)
                    obs_data.sat(i_sat).state_vector(1, ((i_stat_2 - 1) * 5) + 5) = 0;
                end
               
                % Check state vector, whether observation is possible at the considered epoch for the individual station:
                if (sum(obs_data.sat(i_sat).state_vector( ((i_stat_2-1)*5 + 1) : (((i_stat_2-1)*5) + 5) )) == 5)
                    
                    if ~flag_obs_possible_stat(i_stat_2)
                        % Observation is possilble from now on => Set obs_times table!
                        i_row_obs_times_stat(i_stat_2) = i_row_obs_times_stat(i_stat_2) + 1;
                        temp_line = zeros(1, 2); % start_jd, end_jd
                        % temp_line(1,1) = obs_data.sat(i_sat).event_matrix(i_row,1); % JD - start

                        temp_line(1,1) = round_jd2integer_sec(obs_data.sat(i_sat).event_matrix(i_row,1), 'up');

                        obs_data.sat(i_sat).stat(i_stat_2).obs_times = [obs_data.sat(i_sat).stat(i_stat_2).obs_times; temp_line];
                        flag_obs_possible_stat(i_stat_2) = 1;
                    end % if ~flag_obs_possible_stat(i_stat_2)

                elseif (flag_obs_possible_stat(i_stat_2))
                    % Observation is not possilble any more => Set obs_times table!
                    % obs_data.sat(i_sat).stat(i_stat_2).obs_times(i_row_obs_times_stat(i_stat_2), 2) = obs_data.sat(i_sat).event_matrix(i_row,1); % JD - end
                    obs_data.sat(i_sat).stat(i_stat_2).obs_times(i_row_obs_times_stat(i_stat_2), 2) = round_jd2integer_sec(obs_data.sat(i_sat).event_matrix(i_row,1), 'down'); % JD - end
                    flag_obs_possible_stat(i_stat_2) = 0;
                end
                 
            end % for i_stat_2 = 1 : (stat_data.number_of_stations - stat_data.number_of_stations_quasar)
            
            
            % Check state vector, whether observation is possible at the considered epoch for the whole (satellite) station network:
            if (sum(obs_data.sat(i_sat).state_vector) == (5 * (stat_data.number_of_stations - stat_data.number_of_stations_quasar)))
                % Observation is possilble from now on => Set obs_times table!
                i_row_obs_times = i_row_obs_times + 1;
                time_tag_count = time_tag_count + 1;
                temp_line = zeros(1,4);
                temp_line(1,1) = round_jd2integer_sec(obs_data.sat(i_sat).event_matrix(i_row,1), 'up'); % JD - start
                temp_line(1,3) = time_tag_count;                            % Time Tag Label
                obs_data.sat(i_sat).obs_times = [obs_data.sat(i_sat).obs_times; temp_line];
                flag_obs_possible = 1;

            elseif (flag_obs_possible)
                % Observation is not possilble any more => Set obs_times table!
                time_tag_count = time_tag_count + 1;
                temp_line = zeros(1,4);
                obs_data.sat(i_sat).obs_times(i_row_obs_times,2) = round_jd2integer_sec(obs_data.sat(i_sat).event_matrix(i_row,1), 'down'); % JD - end
                obs_data.sat(i_sat).obs_times(i_row_obs_times,4) = time_tag_count;                            % Time Tag Label
                flag_obs_possible = 0;
            end
            
        end % for i_row = 1 : rows
        
    end % for i_sat = 1 : stat_data.number_of_sat

return;

