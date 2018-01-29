% #########################################################################
% #     calc_start_of_scan
% #########################################################################
%
% DESCRIPTION
%   Calculation the earliest possible start time of a scan for a complete station network, 
%   taking into account the previous antenna postion and the antenna slew time. The 
%   calculation is done iteratively bat the function "calc_start_of_obs".
%   Observations to quasars and to satellites are considered.
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
%   2015-03-26     Andreas Hellerschmied
%
% REFERENCES
% - Jing Sun (2013), VLBI Scheduling strategies with respect to VLBI2010,
%   Geowissenschaftliche Mitteilungen, Heft Nr. 92, ISSN 1811-8380.
%
%
% COUPLING
% - calc_start_of_obs
%
%
% INPUT
% - stat_data           - station data structure
% - station_id_list     - List of IDs of the observing stations (referring to the required entry of "stat_data.stat" and "obs_data.stat")
% - obs_data            - Observation data structure (preallocated in "find_observation_times.m")
% - PARA                - Global scheduling parameter strucutre
% - source_quasar       - source (quasars) data structure (ONE sub-set/field of the "source" structure, containing the required data)
% - satellite_id        - ID of the satellite to be observed (referring to an entry of "stat_data.stat.sat" )
%
%
% OUTPUT
% - error_code          - Error Code (0 = no erros occured)
% - obs_data            - Observation data structure, "obs_data.stat(station_id).begin_of_new_obs_temp.un_az" updated. (preallocated in "find_observation_times.m")
% - error_msg           - Error Message (empty, if no errors occured)
% - start_jd            - Earliest possible start time of the scan for the current station network [JD]
%
% CHANGES:
%   - 2016-08-02: A. Hellerschmied: start_jd = 0 initialised.
%

function [start_jd, obs_data, error_code, error_msg] = calc_start_of_scan(stat_data, station_id_list, obs_data, PARA, source_quasar, satellite_id)
 
    % ##### Init #####
    error_code = 0;
    error_msg = '';
    start_jd = 0;
    
    % ##### Check, if it is the first scan in this session #####
    if obs_data.number_of_scans == 0    % first scan
        % For the first scan, the slew time is not critical, because the
        % session can be started a sufficient amount of time before the
        % first observation is scheduled at the stations.
        % ==> The observation start time (start_jd) is set to the session
        % start time in the GUI!
        start_jd = PARA.startmjd + 2400000.5;
        return;
    end
    

    % ##### Calculate the earliest possible scan start for the current station network #####
    
    % Calculate start times for each station:
    num_current_stations = length(station_id_list);
    start_jd_list = zeros(num_current_stations, 1);
    
    for i_stat = 1 : num_current_stations
        station_id = station_id_list(i_stat);
        [start_jd_list(i_stat), slew_time, obs_data, error_code, error_msg] = calc_start_of_obs(stat_data, station_id, obs_data, PARA, source_quasar, satellite_id);
        if error_code > 0
            error_msg = ['calc_start_of_obs: ', error_msg];
            return; 
        end
    end

    % Find earliest possible start time by sorting the data:
    start_jd_list = sort(start_jd_list);
    start_jd = start_jd_list(end);
    
    
    % ##### Check, if the calculated scan start of the current (sub-)network is later than the end of the last scan (obs_data.end_of_last_scan_jd) #####
    if (start_jd < obs_data.end_of_last_scan_jd)
        start_jd = obs_data.end_of_last_scan_jd; % Set the start time of the next scan equal to the end time of the last scan
    end
    
    
return;
