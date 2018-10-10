% #########################################################################
% #     calc_scan_duration_network
% #########################################################################
%
% DESCRIPTION
%   This function calculates the duration of a quasar scan for a VLBI station 
%   network. 
%   The elevation angles, required for the calculation of the adjusted SEFD values,
%   are calculated for the epoch defined by "t_epoch_jd".
%
%
% CREATED  
%   2015-03-31     Andreas Hellerschmied
%
% REFERENCES
% - Jing Sun (2013), VLBI Scheduling strategies with respect to VLBI2010,
%   Geowissenschaftliche Mitteilungen, Heft Nr. 92, ISSN 1811-8380.
%
% COUPLING
% - calc_scan_duration_baseline
% - zazel_s
%
%
% INPUT
% - stat_data           - station data structure
% - station_id_list     - List of IDs of the observing stations (referring to the required entry of "stat_data.stat" and "obs_data.stat")
% - PARA                - Global scheduling parameter strucutre
% - source_quasar       - source (quasars) data structure (ONE sub-set/field of the "source" structure, containing the required data)
% - t_epoch_jd          - Time epoch required to calculate the adjusted SEFD values [JD]
%
% OUTPUT
% - error_code          - Error Code (0 = no erros occured)
% - error_msg           - Error Message (empty, if no errors occured)
% - scan_duration_sec   - Calculated scan duration [sec]
%
% CHANGES:
% - 2015-07-01, A. Hellerschmied: Check network size (min, 2 stations are required to calc. scan durations per baseline)
%

function [scan_duration_sec, error_code, error_msg] = calc_scan_duration_network(stat_data, station_id_list, PARA, source_quasar, t_epoch_jd)
 
    % ##### Init #####
    error_code = 0;
    error_msg = '';
    
    % Loop init.:
    i_baseline = 0;
    number_of_stations = length(station_id_list);
    
    % Check, if min. 2 stations:
    if number_of_stations < 2
       error_code = 1;
       error_msg = ['A network of min. 2 stations is required to be able to calculte the scan duration for baselines. Number of stations = ', num2str(number_of_stations)];
       return;
    end
    
    number_of_baselines = nchoosek(number_of_stations,2); % Binomial coeff.
    baseline_data = zeros(number_of_baselines, 3); % | scan duration [sec] | station 1 id | station 2 id |
    
    
    % Loop over all possible baselines:
    for i_stat_1 = 1 : (number_of_stations - 1)
       
        for i_stat_2 = (i_stat_1 + 1) : number_of_stations
            
            i_baseline = i_baseline + 1;
            
            % Get station IDs:
            stat_id_1 = station_id_list(i_stat_1);
            stat_id_2 = station_id_list(i_stat_2);
            
            % Calculate observation angle (elevation):
            [az, el_1, ha, dc] = zazel_s(t_epoch_jd - 2400000.5, stat_data.stat(stat_id_1).location.ellipsoid.long * pi/180, stat_data.stat(stat_id_1).location.ellipsoid.lat * pi/180, source_quasar.ra, source_quasar.de);
            [az, el_2, ha, dc] = zazel_s(t_epoch_jd - 2400000.5, stat_data.stat(stat_id_1).location.ellipsoid.long * pi/180, stat_data.stat(stat_id_1).location.ellipsoid.lat * pi/180, source_quasar.ra, source_quasar.de);
            
            % ##### Calc. the scan duration for the current baseline #####
            [scan_duration_sec, error_code, error_msg] = calc_scan_duration_baseline(stat_data, stat_id_1, stat_id_2, el_1, el_2, PARA, source_quasar, t_epoch_jd);
            
            baseline_data(i_baseline, 1) = scan_duration_sec;
            baseline_data(i_baseline, 2) = i_stat_1;
            baseline_data(i_baseline, 3) = i_stat_2;
        
            % Test:
            % fprintf(1, '%1.0f - %1.0f \n', i_stat_1, i_stat_2);
        end
    end
    
    
    % ##### Find the maximum scan duration of all possible baselines #####
    
    baseline_data = sortrows(baseline_data);
    scan_duration_sec = baseline_data(end, 1);
    

return;
