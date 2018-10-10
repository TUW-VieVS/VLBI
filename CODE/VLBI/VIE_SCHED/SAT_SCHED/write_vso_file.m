% #########################################################################
% #     write_schedule_file
% #########################################################################
%
% DESCRIPTION
%   This function writes the schedule defined in sched_data to a VSO file.
% 
%   One satellite scan schedules in VIE_SCHED can be fragmented into individual observations by defining the 
%   desried observation interval via the input argument "sat_int_sec" [sec].
%   E.g.    Whole satellite track lasts 300 sec => Divide whole track into 30 sub-observations with 10 sec duration (sat_int_sec = 10).
%
%   sat_int_sec = 0 => No sub-observations in VSO file
%   
%   Output format:
% | Delay reference epoch            | Stat. 1 name  | Stat. 2 name  | source name | obs. type |
% | YYYY MM DD hh mm ss.ssssssssssss | <max 8 char.> | <max 8 char.> | <string>    |   s/q     |
%
% CREATED  
%   2017-01-23     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
%
% INPUT
%   - sched_data                : schedule data structure
%   - sat_int_sec               : Interval for satellite sub-observations [sec]
%   - filepath_out              : Filename of output file
%   - filename_out              : Filepath of output file
%   - flag_geocentric_delay     : Filepath of output file (optional)
%   - pre_sec                   : Each scan starts "pre_sec" seconds earlier (optional)
%   - post_sec                  : Add "post_sec" seconds to scan end time (optional)
%
%
% OUTPUT
%   - error_code          : Error Code (0 = no erros occured)
%   - error_msg           : Error Message (empty, if no errors occured)
%
% CHANGES:
%   - YYYY-MM-DD, author name: description
%   - 2017-02-06, Hellerschmied A.: Option added to write VSO files with geocentric baselines (used for generating the IM files for correlation in DiFX)
%
function [error_code, error_msg] = write_vso_file(sched_data, sat_int_sec, filepath_out, filename_out, varargin)

% ##### Init #####
error_code = 0;
error_msg = '';


% ##### Options #####
flag_write_header = true;

% ##### Optional input arguments #####
switch(nargin)
    case 4 % Standard
        flag_geocentric_delay   = false;
        pre_sec                 = 0;
        post_sec                = 0;
    case 7 % Geocentric delays
        flag_geocentric_delay   = varargin{1};
        pre_sec                 = varargin{2};
        post_sec                = varargin{3};
    otherwise
        error_code = 1;
        error_msg = 'Invald input arguments!';
end


% ##### open file #####
fid = fopen([filepath_out, filename_out], 'w');
if (fid == -1)
    error_msg = ['Can not create/open VSO file at the defined filepath: ', filepath_out, filename_out];
    error_code = 1;
    return;
end

% ##### Write header (optional) #####
if flag_write_header
    station_list = {sched_data.stat.name};
    station_list_str = '';
    for i_stat = 1 : length(sched_data.stat)
        station_list_str = [station_list_str, deblank(station_list{i_stat}), ' '];
    end
    
    fprintf(fid, '################################################################################################\n');
    fprintf(fid, '# VSO file for experiment:  %s\n', sched_data.exper_name);
    fprintf(fid, '# Created in VieVS:         %s\n', datestr(clock));
    fprintf(fid, '# Nominal exper. start:     %s\n', jd2datestr(sched_data.t_nominal_start_jd));
    fprintf(fid, '# Nominal exper. end:       %s\n', jd2datestr(sched_data.t_nominal_end_jd));
    fprintf(fid, '# Stations:                 %s\n', station_list_str);
    fprintf(fid, '################################################################################################\n');
    if flag_geocentric_delay
        fprintf(fid, '# geocentricAprioriDelay\n');
        fprintf(fid, '################################################################################################\n');
        fprintf(fid, '# Reference epoch                | Baseline        | Source name               | Obs. type |     | Geocentric Delay \n');
    else
        fprintf(fid, '# noDelay\n');
        fprintf(fid, '################################################################################################\n');
        fprintf(fid, '# Reference epoch                | Baseline        | Source name               | Obs. type |\n');
    end
    
end

% ##### loop over all scans #####
for i_scan = 1 : sched_data.number_of_scans
    
    
    % ### Determine the baselines of the current scan ###
    number_of_stations  = length(sched_data.scan(i_scan).stat);
    if ~flag_geocentric_delay
        number_of_baselines = nchoosek(number_of_stations,2); % Binomial coeff.
        baseline_names = cell(number_of_baselines, 2); % | station 1 name | station 2 name |
        i_baseline          = 0;

        % ### Loop over all possible baselines ###
        for i_stat_1 = 1 : (number_of_stations - 1)
            for i_stat_2 = (i_stat_1 + 1) : number_of_stations
                i_baseline = i_baseline + 1;
                % Write station names:
                baseline_names{i_baseline, 1} = sched_data.stat(sched_data.scan(i_scan).stat(i_stat_1).stat_id).name; 
                baseline_names{i_baseline, 2} = sched_data.stat(sched_data.scan(i_scan).stat(i_stat_2).stat_id).name; 
            end
        end
    else
        number_of_baselines = number_of_stations; % One baseline per station (GEOCENTER => STATION)
        baseline_names = cell(number_of_baselines, 2); % | 'GEOCENTER' | station 2 name |
        
        % ### Loop over all baselines ###
        for i_stat = 1 : number_of_stations
            % Write station names:
            baseline_names{i_stat, 1} = 'GEOCENTR';
            baseline_names{i_stat, 2} = sched_data.stat(sched_data.scan(i_scan).stat(i_stat).stat_id).name; 
        end
    end

    % ### Get observation type, source name and (manipulated) scan start/end times: ###
    switch(sched_data.scan(i_scan).obs_type)
        case 'sat'
            src_name_str = deblank(sched_data.scan(i_scan).sat_name);
            src_type_str = 's';
            % Get start end end epoch of the current scan:
            t_start_jd  = sched_data.scan(i_scan).t_start_jd    - pre_sec / (24*60*60);
            t_end_jd    = sched_data.scan(i_scan).t_end_jd      + post_sec / (24*60*60);

        case 'quasar'
            src_name_str = sched_data.scan(i_scan).quasar_name;
            src_type_str = 'q';
            % Get start end end epoch of the current scan:
            t_start_jd  = sched_data.scan(i_scan).t_start_jd;
            t_end_jd    = sched_data.scan(i_scan).t_end_jd;
    end    
    
    % Scan epoch = scan start:
    [year, mon, day, hr, min, sec] = invjday(t_start_jd); % Conversion

    % Round to integer seconds:
    sec = round(sec);
    
    % pre/obs/post tags:
    if flag_geocentric_delay
        if  t_start_jd < sched_data.scan(i_scan).t_start_jd % pre
           tmp_str = '          pre';
        else % obs
           tmp_str = '          obs';
        end
    else
       tmp_str = '';
    end


    % ### loop over all possible baselines ###
    for i_bl = 1 : number_of_baselines
        % Format:
        % yyyy mm dd hh mm ss.sss <station_1, 8 char> <station_2, 8 char> <source name> <src type: sc/q>

        % station names:
        station_name_1_str = baseline_names{i_bl, 1}; % stat_data.stat(baseline_data(i_bl, 1)).name;
        station_name_2_str = baseline_names{i_bl, 2}; % stat_data.stat(baseline_data(i_bl, 2)).name;
           
        % Write line to file:
        % | YYYY MM DD hh mm ss.ssssssssssss | <max 8 char.> | <max 8 char.> | <string>    |   s/q     |
        fprintf(fid, '%4d %02d %02d %02d %02d %015.12f  %8s %8s  %25s  %s%s\n', year, mon, day, hr, min, sec, station_name_1_str, station_name_2_str, src_name_str, src_type_str, tmp_str);
    end %  for i_bl = 1 : number_of_baselines


    % ##### Write sub-observations for satellite tracks (optional) #####)
    if strcmp(src_type_str, 's') % satellite scan
        if (sat_int_sec ~= 0)

            % Scan duration [integer sec]:
            scan_duration_sec = round((t_end_jd - t_start_jd) * (24*60*60));

            % number of sub-scans (in the antenna reposition interval):
            number_of_sub_scans = floor(scan_duration_sec / sat_int_sec);

            % #### Loop over all sub-scans ####
            for i_sub_scan = 1 : number_of_sub_scans

                t_epoch_jd = t_start_jd + i_sub_scan * sat_int_sec / (24*60*60);
                
                [year, mon, day, hr, min, sec] = invjday(t_epoch_jd); % Conversion

                % Round to integer seconds:
                sec = round(sec);
                
                % pre/obs/post tags:
                if flag_geocentric_delay
                    if (t_epoch_jd >= sched_data.scan(i_scan).t_start_jd) && (t_epoch_jd <= sched_data.scan(i_scan).t_end_jd) % obs
                        tmp_str = '          obs';
                    elseif  t_epoch_jd < sched_data.scan(i_scan).t_start_jd % pre
                        tmp_str = '          pre';
                    else
                        tmp_str = '          post';
                    end
                else
                   tmp_str = '';
                end
                
                % ### loop over all possible baselines ###
                for i_bl = 1 : number_of_baselines
                   
                    % station names:
                    station_name_1_str = baseline_names{i_bl, 1};
                    station_name_2_str = baseline_names{i_bl, 2};

                    % Write line to file:
                    fprintf(fid, '%4d %02d %02d %02d %02d %015.12f  %8s %8s  %25s  %s%s\n', year, mon, day, hr, min, sec, station_name_1_str, station_name_2_str, src_name_str, src_type_str, tmp_str);

                end % i_bl = 1 : number_of_baselines

            end % i_sub_scan = 1 : number_of_sub_scans

        end % if (sat_int_sec ~= 0)

    end % if strcmp(sched_data.scan(i_scan).obs_type, 's')
    
end % for i_scan = 1 : sched_data.number_of_scans



% ##### close  file #####
if ( (exist('fid', 'var') ) && (fid ~= -1) )
    fclose(fid);
end

return;