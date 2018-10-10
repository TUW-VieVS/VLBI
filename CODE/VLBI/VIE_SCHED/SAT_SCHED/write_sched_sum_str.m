% #########################################################################
% #     write_sched_sum_str
% #########################################################################
%
% DESCRIPTION
%   This function writes a schedule summary string containing information of all 
%   scans in the "sched_data" structure.
%
%
% CREATED  
%   2015-08-28     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
% - jd2datestr
%
%
% INPUT
% - sched_data              - schedule data structure
%
%
% OUTPUT
% - error_code              - Error Code (0 = no erros occured)
% - error_msg               - Error Message (empty, if no errors occured)
% - sched_sum_str           - String containing the schedule summary (all scans included in "sched_data" structure)
%
% CHANGES:
%


function [final_sched_sum_str, error_code, error_msg] = write_sched_sum_str(sched_data)
 
    % ##### Init #####
    error_code = 0;
    error_msg = '';
    final_sched_sum_str = ''; % Init. empty output string
    scan_num_str = '';
    sched_sum_str = ''; 
    temp_line_str = '';
    type_str = '';
    source_name_str = '';
    t_scan_start_str = '';
    t_scan_end_str = '';
    duration_str ='';
    
    
    % ##### Inital checks #####
    % Does the schedule contain at least on valid scan?
    if sched_data.number_of_scans < 1
        error_code = 1;
        error_msg = 'The schedule has to cotain at least one scan!';
        return;
    end
    
    t_start_str = jd2datestr(sched_data.t_nominal_start_jd);
    t_end_str = jd2datestr(sched_data.t_nominal_end_jd);
    
    % ##### Creat output string #####
    sched_sum_str = [sched_sum_str, sprintf('###################################################################################\n')];
    sched_sum_str = [sched_sum_str, sprintf('# Schedule summary for session: %s                                              #\n', sched_data.exper_name)];
    sched_sum_str = [sched_sum_str, sprintf('###################################################################################\n')];
    sched_sum_str = [sched_sum_str, sprintf('Nominal start:     %s\n', t_start_str)];
    sched_sum_str = [sched_sum_str, sprintf('Nominal end:       %s\n', t_end_str)];
    % sched_sum_str = [sched_sum_str, sprintf('Used TLE file:     %s\n', sched_data.TLE_filename)];
    sched_sum_str = [sched_sum_str, sprintf('-----------------------------------------------------------------------------------\n')];
    sched_sum_str = [sched_sum_str, sprintf('Scan#  Source name              Type  Start [UT]  End [UT]  Dur [s]  Sations\n')];
    
    % Loop over all scheduled scans:
    for i_scan = 1 : sched_data.number_of_scans
        
        % #### prepare data ####
        
        % Scan number:
        scan_num_str = sprintf('%1.0f%', i_scan);
        
        % observation type and source name:
        switch(sched_data.scan(i_scan).obs_type)
            case 'sat'
                type_str = 's';
                source_name_str = deblank(sched_data.scan(i_scan).sat_name);
            case 'quasar'
                type_str = 'q';
                source_name_str = deblank(sched_data.scan(i_scan).quasar_name);
            otherwise
                error_code = 2;
                error_msg = 'Unknown scan type!';
                return;
        end % switch(sched_data.scan(i_scan).obs_type)
        
        % Check if the soure name exists:
        if isempty(source_name_str)
            error_code = 3;
            error_msg = 'Source name is not available!';
            return;
        end
        
        % Scan start and stop time:
        t_scan_start_str = jd2datestr(sched_data.scan(i_scan).t_start_jd, 'HH:MM:SS');
        t_scan_end_str = jd2datestr(sched_data.scan(i_scan).t_end_jd, 'HH:MM:SS');
        
        % duration [sec]
        duration_str = sprintf('%1.0f', (sched_data.scan(i_scan).t_end_jd - sched_data.scan(i_scan).t_start_jd) * 24*60*60);
        
        % Stations:
        stat_label_str = '';
        for i_stat = 1 : length(sched_data.scan(i_scan).stat)
            stat_id = sched_data.scan(i_scan).stat(i_stat).stat_id;
            stat_label_str = [stat_label_str, ' ', sched_data.stat(stat_id).label];
        end
        
        % Prepare line:
        temp_line_str = sprintf('%-6s %-24s %-5s %-11s %-9s %-7s %s\n', scan_num_str, source_name_str, type_str, t_scan_start_str, t_scan_end_str, duration_str, stat_label_str);
        
        % Write line to output string:
        sched_sum_str = [sched_sum_str, temp_line_str];
    end
    
    sched_sum_str = [sched_sum_str, sprintf('-----------------------------------------------------------------------------------\n')];
    
    
    % ##### Assign output string if no error occured #####
    final_sched_sum_str = sched_sum_str;

return;
