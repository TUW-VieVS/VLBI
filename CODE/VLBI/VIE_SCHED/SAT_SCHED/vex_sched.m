% -------------------------------------------------------------------------
%
%                              vex_sched
%
%   Writes the $SCHED Block to a VEX File.
%
%   Note on sub-netting:    The creation of "combined" VEX files, used as correlator input,
%                           is curently not supported if sub-netting is
%                           applied! The ref. station has to join all
%                           scans in case of sub-netting!
%
%   Author: 
%       2013-11-20 : Andreas Hellerschmied
%   
%   changes       :
%   - 2015-06-19, A. Hellerschmied: Updated for VieVS 2.3
%   - 2015-06-29, A. Hellerschmied: Cable wrap sector can now be written to the VEX file
%   - 2015-07-07, A. Hellerschmied: Bug-fix: station index error for sched_data.scan.stat
%   - 2015-10-21, A. Hellerschmied: Possibility added to write combined VEX files.
%   - 2016-05-12, A. Hellerschmied: Added option to use one static mode definition for all satellites ("flag_use_static_freq_setup")
%   - 2016-10-28, A. Hellerschmied: Get start time of scans and scan durations individually for each station
%   - 2016-11-02, A. Hellerschmied: Option to deselect stepwise satellite tracking added 
%   - 2016-11-15, A. Hellerschmied: Changes to write proper combined VEX file with scans of ALL stations, not only the ref. station!
%   - 2016-11-24, A. Hellerschmied: Bug. fix: Number of repos. int. per station ins now checked correctly!
%           
%
%   inputs        :
%   - fid_vex           : File ID of VEX File
%   - vex_para          : Data from VEX paramer file, structure
%   - sched_data        : scheduling data structure
%   - stat_id_list      : List of IDs of stations which should be considered
%     
%
%   outputs       :
%   - error_code        : Error Code (0 = no erros occured)
%   - error_msg         : Error Message (empty, if no errors occured)
%    
%
%   locals        :
%  
%
%   coupling      :
%   - round_sec_of_date_to_n_dec_figures
%   
%
%   references    :
%   - VEX File Definition/Example, Rev 1.5b1, 30. Jan. 2002, available
%     online at: http://www.vlbi.org/vex/
%   - VEX Parameter Tables, Rev 1.5b1, 29. Jan. 2002, available
%     online at: http://www.vlbi.org/vex/
%
%-------------------------------------------------------------------------

function [error_code, error_msg] = vex_sched(fid_vex, vex_para, sched_data, stat_id_list) 

    % ##### Init #####
    error_code = 0;
    error_msg = '';
    
    
%     flag_use_static_freq_setup = 0;
    
    if isfield(vex_para, 'flag_use_stepwise_sat_tracking')
        flag_stepwise_tracking = logical(str2double(vex_para.flag_use_stepwise_sat_tracking));
    else
        flag_stepwise_tracking = true;
    end
    
    
    % ##### Wite block #####
    if isempty(stat_id_list) % Error
        error_code = 2;
        error_msg = 'Station ID list is empty';  
        return;
    else
        
        % Init.:
        number_of_stat = length(stat_id_list);
        total_amount_of_data_gbyte = zeros(length(sched_data.stat), 1); % Amount of recorded data [GByte], per station!
        flag_stat_dependend_vex_file = true;
        

        % ##### Distinguish between station-dependent and combined VEX files #####
        number_of_stat = length(stat_id_list);
        if number_of_stat == 1 % Station-dependent
            station_id = stat_id_list(1);
            % ##### Get data (for ref. station in case of comb. vex files): #####
            station_name_ref         = sched_data.stat(station_id).name;
            number_of_sat           = length(sched_data.stat(station_id).statistics.sat_name_list);
            sat_num_vector          = sched_data.stat(station_id).statistics.sat_norad_id_vector;
            sat_name_list           = sched_data.stat(station_id).statistics.sat_name_list;
            number_of_quasars       = length(sched_data.stat(station_id).statistics.quasar_name_list);
        else % combined
            % ### Prepare reference station used for the combined VEX file ###
            flag_stat_dependend_vex_file = false;
            if isfield(vex_para, 'comb_vex_file_ref_station')
                ref_station_name = sprintf('%-*s' , 8 , vex_para.comb_vex_file_ref_station);
                % Ref. station ID:
                if strcmp(ref_station_name, 'none    ')
                   station_id = stat_id_list(1);
                else
                    station_id = find(strncmp({sched_data.stat.name}, ref_station_name, 8));
                end
            else
                fprintf( ' No reference station defined for combined VEX files.\n'); 
                station_id = stat_id_list(1);
            end
            fprintf( ' => Obs. mode names in $MODE are taken from station: %s\n', sched_data.stat(station_id).name);
            
            % ##### Get data (label for ref. station in case of comb. vex files): #####
            station_name_ref    = sched_data.stat(station_id).name;
            sat_name_list       = sched_data.exper_statistics.sat_name_list;
            sat_num_vector      = sched_data.exper_statistics.sat_norad_id_vector;
            number_of_quasars   = sched_data.exper_statistics.number_of_quasars;
            number_of_sat       = sched_data.exper_statistics.number_of_sat;
        end
        
        % #### Get obs. mode labels and check, if a static frequ. setup/mode should be used  ####
        if number_of_sat > 0 % If satellites are observed...
            i_mode = 1; % Take first Mode element defined in VEX parameter file for the ref. station
            try
                mode_label_sat = eval(['vex_para.', station_name_ref , '.mode(', num2str(i_mode), ').label;']);
            catch error
                error_code = 1;
                error_msg = 'Mode Label for satellites not available';  
                return;
            end
            
            % ##### Check, if a static frequ. setup/mode should be used #####
            try 
                i_frequ = 1; % satellite obs. setup
                temp_str_1 = eval(['vex_para.', station_name_ref , '.freq(', num2str(i_frequ), ').flag_use_static_freq_setup.str;']);
                flag_use_static_freq_setup = logical(str2double(temp_str_1));
            catch
                flag_use_static_freq_setup = 0;
            end
        end
        if number_of_quasars > 0 % If quasars are observed...
            i_mode = 2; % Take second Mode element defined in VEX parameter file  for the ref. station
            try
                mode_label_quasar = eval(['vex_para.', station_name_ref , '.mode(', num2str(i_mode), ').label;']);
            catch error
                error_code = 1;
                error_msg = 'Mode Label for quasars not available';  
                return;
            end
        end


        % ##### Prepare data #####

        % #### Get information to calculate the amount of data recorded: ####
        % sampling_rate
        try
            temp_str_1 = vex_para.sampling_rate;
            sampling_rate = str2num(temp_str_1);
        catch error
            error_code = 1;
            error_msg = '"sampling_rate" parameter not available.';  
            return;
        end

        % data_resolution
        try
            temp_str_1 = vex_para.data_resolution;
            data_resolution = str2num(temp_str_1);
        catch error
            error_code = 1;
            error_msg = '"data_resolution" parameter not available.';  
            return;
        end

        % number_of_channels
        try
            temp_str_1 = vex_para.number_of_channels;
            number_of_channels = str2num(temp_str_1);
        catch error
            error_code = 1;
            error_msg = '"number_of_channels" parameter not available.';  
            return;
        end

        recording_rate = (sampling_rate * data_resolution * number_of_channels); % [bit/sec] 

        % #### Define number of decimal places of seconds for the source_label: ####
        switch (vex_para.source_notation)
            case '1' % <satellite_label>_<hhmmss>
                dec_places = 0;
            case '2' % <satellite_label>_<hhmmss.ss>
                dec_places = 2;
            case '3' % <satellite_label>_<doyhhmmss>
                dec_places = 0;
            case '4' % <hhmmss>
                dec_places = 0;
            case '5' % <doyhhmmss>
                dec_places = 0;
            otherwise
                error_code = 1;
                error_msg = 'Source notation info not available.';  
                return;
        end
        if (dec_places == 0)
            digits = 2;
        else
            digits = dec_places + 3;
        end

        % ##### Write SCHED block #####
        fprintf(fid_vex, '$SCHED;\n');
        fprintf(fid_vex, '*\n');

        
        % ##### loop over all scans #####
        for i_scan = 1 : sched_data.number_of_scans

            % Loop init.:
            flag_any_stat_in_list_joins_scan = 0;

            % #### Check, if any station in "stat_id_list" joins this scan ####
            % loop over all stations of the scan:
            for i_stat_2 = 1 : sched_data.scan(i_scan).number_of_stat
                if sum(sched_data.scan(i_scan).stat(i_stat_2).stat_id == stat_id_list)
                    % Set flag, if the station is included in the scan
                    flag_any_stat_in_list_joins_scan = 1;
                    % Get (ref.) station ID for struct "sched_data.scan.stat"
                    break;
                end
            end % for i_stat_2 = 1 : sched_data.scan(i_scan).number_of_stat


            % ##### If the (ref.) station joins this scan... #####
            if flag_any_stat_in_list_joins_scan

                % ### Distinguish between observation-type ###
                switch(sched_data.scan(i_scan).obs_type)

                    % satellite scan:
                    case 'sat'                    
                        flag_sat_scan       = 1;
                        flag_quasar_scan    = 0;

                    % quasar scan:
                    case 'quasar'
                        flag_sat_scan       = 0;
                        flag_quasar_scan    = 1;

                end % switch(sched_data.scan(i_scan).obs_type)

            else % No station in "stat_id_list" joins this scan
                continue;
            end % if flag_stat_joins_scan


            % ##### Write scan definition #####

            % ######################################################################
            % # Satellite scan
            % ######################################################################
            if flag_sat_scan

                % Satellite Label:
                sat_label = sched_data.scan(i_scan).sat_name;
                sat_label = sat_label(~isspace(sat_label)); % Remove white space

                % Get satelite NORAD ID (as char string):
                sat_number_str = num2str(sched_data.scan(i_scan).sat_number);

                % Get scan duration:
                scan_duration_sec = (sched_data.scan(i_scan).t_end_jd - sched_data.scan(i_scan).t_start_jd) * 24*60*60;
%                 scan_duration_sec = (sched_data.scan(i_scan).stat(stat_id_scan).end.jd - sched_data.scan(i_scan).stat(stat_id_scan).start.jd) * 24*60*60; % Individually for each station

                % Write scan header for this satellite scan:
                fprintf(fid_vex, '* ---- Scan %d of %d: %s (NORAD ID: %s) ----\n', i_scan, sched_data.number_of_scans, sat_label, sat_number_str);
                fprintf(fid_vex, '*      Scan duration: %1.2f sec\n', scan_duration_sec);
                fprintf(fid_vex, '*\n');

                % Create Mode Label for satellite observation:
                if ~flag_use_static_freq_setup
                    mode_label = [sat_number_str, '.', mode_label_sat]; % <NORAD_ID>.<MODE_LABEL>
                else
                    mode_label = ['sat.', mode_label_sat]; % <NORAD_ID>.<MODE_LABEL>
                end

                % ##### loop over all antenna reposition epochs within this scan (stepwise tracking) #####
                if flag_stepwise_tracking
                    num_of_repos_int = length(sched_data.scan(i_scan).stat(1).epoch) - 1; % The number of repos. int. has to be the same for all stations in a scan => define a new scan, if a station drops out!
                else
                    num_of_repos_int = 1;
                end
                for i_epoch = 1 : num_of_repos_int

                    % Loop init.:
                    temp_str = '';
                    temp_str_1 = '';

                    % Get scan start epoch:
                    [year, mon, day, hr, min, sec] = invjday(sched_data.scan(i_scan).stat(1).epoch(i_epoch).jd); % Start epochs have to be the same for all stations in a scan 
%                     [year, mon, day, hr, min, sec] = invjday(sched_data.scan(i_scan).t_start_jd); % All staions start at the same time
                    [days_of_year] = tdays(year, mon, day);
                    % Get source label (ref. to $SOURCE block):
                    switch (vex_para.source_notation)
                        case '1' % <satellite_label>_<hhmmss>
                            [year, days_of_year, hr, min, sec] = round_sec_of_date_to_n_dec_figures(year, days_of_year, hr, min, sec, dec_places);
                            source_label = [sat_label, '.', sprintf('%02.0f',hr), sprintf('%02.0f',min), sprintf(['%0',num2str(digits) ,'.',num2str(dec_places) ,'f'],sec)];
                        case '2' % <satellite_label>_<hhmmss.ss>
                            [year, days_of_year, hr, min, sec] = round_sec_of_date_to_n_dec_figures(year, days_of_year, hr, min, sec, dec_places);
                            source_label = [sat_label, '.', sprintf('%02.0f',hr), sprintf('%02.0f',min), sprintf(['%0',num2str(digits) ,'.',num2str(dec_places) ,'f'],sec)];
                        case '3' % <satellite_label>_<doyhhmmss>
                            [year, days_of_year, hr, min, sec] = round_sec_of_date_to_n_dec_figures(year, days_of_year, hr, min, sec, dec_places);
                            source_label = [sat_label, '.', sprintf('%s',days_of_year), sprintf('%02.0f',hr), sprintf('%02.0f',min), sprintf(['%0',num2str(digits) ,'.',num2str(dec_places) ,'f'],sec)];
                        case '4' % <hhmmss>
                            [year, days_of_year, hr, min, sec] = round_sec_of_date_to_n_dec_figures(year, days_of_year, hr, min, sec, dec_places);
                            source_label = [sprintf('%02.0f',hr), sprintf('%02.0f',min), sprintf(['%0',num2str(digits) ,'.',num2str(dec_places) ,'f'],sec)];
                        case '5' % <doyhhmmss>
                            [year, days_of_year, hr, min, sec] = round_sec_of_date_to_n_dec_figures(year, days_of_year, hr, min, sec, dec_places);
                            source_label = [sprintf('%s',days_of_year), sprintf('%02.0f',hr), sprintf('%02.0f',min), sprintf(['%0',num2str(digits) ,'.',num2str(dec_places) ,'f'],sec)];  
                    end

                    % --- scan definition: ---
                    temp_str = sprintf('%04d', i_scan);
                    scan_label = [temp_str, '.', source_label];
                    fprintf(fid_vex, 'scan %s;\n', scan_label);

                    % --- start: ---
                    [year, days_of_year, hr, min, sec] = round_sec_of_date_to_n_dec_figures(year, days_of_year, hr, min, sec, dec_places);
                    temp_str = [sprintf('%04.0fy',year), sprintf('%sd',days_of_year), sprintf('%02.0fh',hr), sprintf('%02.0fm',min), sprintf(['%0',num2str(digits) ,'.',num2str(dec_places) ,'fs'],sec)];
                    fprintf(fid_vex, '    start = %s;\n', temp_str);

                    % --- mode: ---
                    fprintf(fid_vex, '    mode = %s;\n', mode_label);

                    % --- source: ---
                    fprintf(fid_vex, '    source = %s.%s;\n', sat_number_str, source_label);


                    % --- station: ---

                    % #### Loop over all stations defined in this scan ####
                    for i_stat = 1 : length(sched_data.scan(i_scan).stat)
                        
                        % Check, if the number of repos. int is the same for all stations in scan:
                        if flag_stepwise_tracking
                            if length(sched_data.scan(i_scan).stat(i_stat).epoch) <= num_of_repos_int
                                % ERROR 
                                % The number of repos. int. has to be the same for all stations in a scan 
                                % => define a new scan, if a station drops out!
                                error_code = 1;
                                error_msg = sprintf('Not enough repo. interval (%d) defined for station %s in scan %d!', length(sched_data.scan(i_scan).stat(i_stat).epoch), sched_data.stat(station_id).name, i_scan);
                                return;
                            end
                        end
                        
                        % Check, if the current station ("i_stat" in "sched_data.scan(i_scan).stat") joins this scan:
                        if sum(sched_data.scan(i_scan).stat(i_stat).stat_id == stat_id_list) == 1
                            
                            station_id = sched_data.scan(i_scan).stat(i_stat).stat_id; % ref. to "sched_data.stat(station_id)"
                            stat_label = sched_data.stat(station_id).label;
                            
                            % Get observation duration [sec]:
                            if flag_stepwise_tracking
                                obs_duration_sec = (sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch + 1).jd - sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).jd) * 24*60*60;
                            else
                                obs_duration_sec = (sched_data.scan(i_scan).stat(i_stat).end.jd - sched_data.scan(i_scan).stat(i_stat).start.jd) * 24*60*60;
                            end

                            % Amount of recorded data [GByte]:
                            recorded_data_gbyte                 = (recording_rate * obs_duration_sec) /(8 * 1024); % [GByte / this obs.]
                            total_amount_of_data_gbyte(station_id)  = total_amount_of_data_gbyte(station_id) + recorded_data_gbyte;
                            
                            % Init:
                            pointsectr_id_str = ''; % empty string

                            % Determine Cable wrap pointing sector (only for AZEL antenna mount type and "write_cable_wrap_sector" flag in VEX parameter file = 1)
                            if (strncmpi(sched_data.stat(station_id).axis_type_1, 'az', 2) && strncmpi(sched_data.stat(station_id).axis_type_2, 'el', 2)) && strncmpi(vex_para.write_cable_wrap_sector, '1', 1)
                                % Get pointsectr_id_str
                                try
                                un_az_rad = sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).un_az;
                                catch
                                   error_code = 5;
                                   error_msg = sprintf('Tracking data is still missing => Calculate it first! (i_stat: %d; i_scan: %d; i_epoch: %d)', i_stat, i_scan, i_epoch) ;
                                   return;
                                end
                                if (un_az_rad >= sched_data.stat(station_id).az_n1) && (un_az_rad <= sched_data.stat(station_id).az_n2) % neutral
                                    pointsectr_id_str = '&n';
                                elseif (un_az_rad >= sched_data.stat(station_id).az_cw1) && (un_az_rad <= sched_data.stat(station_id).az_cw2) % clockwise
                                    pointsectr_id_str = '&cw';
                                elseif (un_az_rad >= sched_data.stat(station_id).az_ccw1) && (un_az_rad <= sched_data.stat(station_id).az_ccw2) % counter-clockwise
                                    pointsectr_id_str = '&ccw';       
                                end
                            end

                            % Station id 
                            temp_str = stat_label;

                            % data_start
                            try
                                temp_str_1 = vex_para.data_start;
                                temp_str = [temp_str, ' : ', temp_str_1];
                            catch error
                                error_code = 1;
                                error_msg = '"data_start" parameter not available.';  
                                return;
                            end

                            % data_stop
                            try
                                temp_str_1 = sprintf('%4.1f sec', obs_duration_sec);
                                temp_str = [temp_str, ' : ', temp_str_1];
                            catch error
                                error_code = 1;
                                error_msg = '"data_stop" parameter not available.';  
                                return;
                            end

                            % start_pos
                            % recording_rate [bit/sec] = (sampling rate [Hz] * resolution [bit] *
                            % number of channels)

                            temp_str_1 = sprintf('%7.3f GB', total_amount_of_data_gbyte(station_id));
                            temp_str = [temp_str, ' : ', temp_str_1];

                            % pass
                            try
                                 temp_str_1 = vex_para.pass;
                            catch error
                                temp_str_1 = ' ';
                            end
                            temp_str = [temp_str, ' : ', temp_str_1];

                            % pointsectr
                            if ~isempty(pointsectr_id_str) % If AzEl antenna AND poinsector should be written to VEX file
                                temp_str_1 = sprintf('%4s', pointsectr_id_str);
                            else
                                try
                                     temp_str_1 = vex_para.pointsectr;
                                catch error
                                    temp_str_1 = '    ';
                                end
                            end
                            temp_str = [temp_str, ' : ', temp_str_1];

                            % drive
                            try
                                if (sched_data.scan(i_scan).flag_record_scan == 0)
                                    temp_str_1 = '0';
                                else
                                    temp_str_1 = '1';
                                end
                                % temp_str_1 = vex_para.drive;
                            catch error
                                error_code = 1;
                                error_msg = '"drive" parameter not available.';  
                                return;
                            end
                            temp_str = [temp_str, ' : ', temp_str_1];

                            fprintf(fid_vex, '    station = %s;\n', temp_str);
                            
                        end % if sum(sched_data.scan(i_scan).stat(i_stat).stat_id == stat_id_list) == 1
                        
                    end % for i_stat = 1 : length(sched_data.scan(i_scan).stat)

                    
                    fprintf(fid_vex, 'endscan;\n');
                    fprintf(fid_vex, '*\n');

                end % for i_epoch = 1 : num_of_repos_int

            end % if flag_sat_scan



            % ######################################################################
            % # Quasars scan
            % ######################################################################
            if flag_quasar_scan

                temp_str = '';
                temp_str_1 = '';

                % Get scan duration:
%                 scan_duration_sec = sched_data.scan(i_scan).stat(i_stat).duration_sec; % Individually for each station
                scan_duration_sec = (sched_data.scan(i_scan).t_end_jd - sched_data.scan(i_scan).t_start_jd) * 24*60*60;

                % Get scan start time:
                [year, mon, day, hr, min, sec] = invjday (sched_data.scan(i_scan).t_start_jd); % Same for each station
                [days_of_year] = tdays(year, mon, day);

                
                % Get source label (ref. to $SOURCE block):
                source_label = sched_data.scan(i_scan).quasar_name;

                % Write scan header for this quasar scan:
                fprintf(fid_vex, '* ---- Scan %d of %d: quasar %s ----\n', i_scan, sched_data.number_of_scans, source_label);
                fprintf(fid_vex, '*      Scan duration: %1.2f sec\n', scan_duration_sec);
                fprintf(fid_vex, '*\n');


                % --- scan definition: ---
                temp_str = sprintf('%04d', i_scan);
                scan_label = [temp_str, '.', source_label];
                fprintf(fid_vex, 'scan %s;\n', scan_label);

                % --- start: ---
                [year, days_of_year, hr, min, sec] = round_sec_of_date_to_n_dec_figures(year, days_of_year, hr, min, sec, dec_places);
                temp_str = [sprintf('%04.0fy',year), sprintf('%sd',days_of_year), sprintf('%02.0fh',hr), sprintf('%02.0fm',min), sprintf(['%0',num2str(digits) ,'.',num2str(dec_places) ,'fs'],sec)];
                fprintf(fid_vex, '    start = %s;\n', temp_str);

                % --- mode: ---
                fprintf(fid_vex, '    mode = %s;\n', mode_label_quasar);

                % --- source: ---
                fprintf(fid_vex, '    source = %s;\n', source_label);

                
                % #### Loop over all stations in this scan ####
                for i_stat = 1 : length(sched_data.scan(i_scan).stat)
                    
                    
                    % Check, if the current station ("i_stat" in "sched_data.scan(i_scan).stat") joins this scan:
                    if sum(sched_data.scan(i_scan).stat(i_stat).stat_id == stat_id_list) == 1

                        station_id = sched_data.scan(i_scan).stat(i_stat).stat_id; % ref. to "sched_data.stat(station_id)"
                        stat_label = sched_data.stat(station_id).label;
                        
%                         scan_duration_stat_sec = sched_data.scan(i_scan).stat(i_stat).duration_sec;
                        scan_duration_stat_sec = (sched_data.scan(i_scan).stat(i_stat).end.jd - sched_data.scan(i_scan).stat(i_stat).start.jd) * 24*60*60;
                        
                        % Amount of recorded data [GByte]:
                        recorded_data_gbyte = (recording_rate * scan_duration_stat_sec) /(8 * 1024); % [GByte / this obs.]
                        total_amount_of_data_gbyte(station_id) = total_amount_of_data_gbyte(station_id) + recorded_data_gbyte;

                        % Init:
                        pointsectr_id_str = ''; % empty string

                        % --- station: ---

                        % Determine Cable wrap pointing sector (only for AZEL antenna mount type and "write_cable_wrap_sector" flag in VEX parameter file = 1)
                        if (strncmpi(sched_data.stat(station_id).axis_type_1, 'az', 2) && strncmpi(sched_data.stat(station_id).axis_type_2, 'el', 2)) && strncmpi(vex_para.write_cable_wrap_sector, '1', 1)
                            % Get pointsectr_id_str
                            un_az_rad = sched_data.scan(i_scan).stat(i_stat).epoch(1).un_az;
                            try
                                if (un_az_rad >= sched_data.stat(station_id).az_n1) && (un_az_rad <= sched_data.stat(station_id).az_n2) % neutral
                                    pointsectr_id_str = '&n';
                                elseif (un_az_rad >= sched_data.stat(station_id).az_cw1) && (un_az_rad <= sched_data.stat(station_id).az_cw2) % clockwise
                                    pointsectr_id_str = '&cw';
                                elseif (un_az_rad >= sched_data.stat(station_id).az_ccw1) && (un_az_rad <= sched_data.stat(station_id).az_ccw2) % counter-clockwise
                                    pointsectr_id_str = '&ccw';       
                                end
                            catch
                                disp('tst');
                                keyboard
                            end
                        end

                        % Station id
                        temp_str = stat_label;

                        % data_start
                        try
                            temp_str_1 = vex_para.data_start;
                            temp_str = [temp_str, ' : ', temp_str_1];
                        catch error
                            error_code = 1;
                            error_msg = '"data_start" parameter not available.';  
                            return;
                        end

                        % data_stop
                        try
                            temp_str_1 = sprintf('%4.1f sec', scan_duration_stat_sec);
                            temp_str = [temp_str, ' : ', temp_str_1];
                        catch error
                            error_code = 1;
                            error_msg = '"data_stop" parameter not available.';  
                            return;
                        end

                        % start_pos
                        % recording_rate [bit/sec] = (sampling rate [Hz] * resolution [bit] *
                        % number of channels)

                        temp_str_1 = sprintf('%7.3f GB', total_amount_of_data_gbyte(station_id));
                        temp_str = [temp_str, ' : ', temp_str_1];

                        % pass
                        try
                             temp_str_1 = vex_para.pass;
                        catch error
                            temp_str_1 = ' ';
                        end
                        temp_str = [temp_str, ' : ', temp_str_1];

                        % pointsectr
                        if ~isempty(pointsectr_id_str) % If AzEl antenna AND poinsector should be written to VEX file
                            temp_str_1 = sprintf('%4s', pointsectr_id_str);
                        else
                            try
                                 temp_str_1 = vex_para.pointsectr;
                            catch error
                                temp_str_1 = '    ';
                            end
                        end
                        temp_str = [temp_str, ' : ', temp_str_1];


                        % drive
                        try
                            if (sched_data.scan(i_scan).flag_record_scan == 0)
                                temp_str_1 = '0';
                            else
                                temp_str_1 = '1';
                            end
                            % temp_str_1 = vex_para.drive;
                        catch error
                            error_code = 1;
                            error_msg = '"drive" parameter not available.';  
                            return;
                        end
                        temp_str = [temp_str, ' : ', temp_str_1];

                        fprintf(fid_vex, '    station = %s;\n', temp_str);
                        
                    end % if sum(sched_data.scan(i_scan).stat(i_stat).stat_id == stat_id_list) == 1
                     
                end % for i_stat = 1 : length(sched_data.scan(i_scan).stat)

                fprintf(fid_vex, 'endscan;\n');
                fprintf(fid_vex, '*\n');

            end % if flag_quasar_scan

        end % for i_scan = 1 : sched_data.number_of_scans



        % ##### Additional notes, if available: ##### 
        try
            fprintf(fid_vex, '*   %s\n', vex_para.sched_note);
        catch error

        end

        fprintf(fid_vex, '* ---------------------------------------------------\n');
        fprintf(fid_vex, '*\n');
        
    end % if isempty(stat_id_list)  
    
return;

