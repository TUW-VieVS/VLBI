% -------------------------------------------------------------------------
%
%                              vex_source
%
%   Writes the $SOURCE Block to a VEX File.
%       - 
%
%   Author: 
%       2013-11-14 : Andreas Hellerschmied (heller182@gmx.at)
%   
%   changes       :
%   - 2014-01-14 : Andreas Hellerschmied: 
%       -   Topo. Radec Output is now always
%           positiv and within the interval [0 h, 24 h].
%       -   Resolved problem with output of negative declination angle 
%           values. 
%   - 2015-06-22, A. Hellerschmied: Updated for VieVS 2.3
%   - 2015-07-07, A. Hellerschmied: Bug-fix: station index error for sched_data.scan.stat
%   - 2016-11-02, A. Hellerschmied: Option deselect stepwise satellite tracking added 
%   - 2016-11-15, A. Hellerschmied: Changes to write proper combined VEX file with scans of ALL stations, not only the ref. station!
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
%   
%
%   references    :
%   - VEX File Definition/Example, Rev 1.5b1, 30. Jan. 2002, available
%     online at: http://www.vlbi.org/vex/
%   - VEX Parameter Tables, Rev 1.5b1, 29. Jan. 2002, available
%     online at: http://www.vlbi.org/vex/
%
%-------------------------------------------------------------------------

% function [error_code, error_msg] = vex_source(fid_vex, vex_para, sched_data, i_stat)
function [error_code, error_msg] = vex_source(fid_vex, vex_para, sched_data, stat_id_list)

    % Init
    error_code = 0;
    error_msg = '';
    
    all_src_labels_str_array = {};
    
    if isfield(vex_para, 'flag_use_stepwise_sat_tracking')
        flag_stepwise_tracking = logical(str2double(vex_para.flag_use_stepwise_sat_tracking));
    else
        flag_stepwise_tracking = true;
    end
    
%     % Total number of observation targets for this station:
%     number_of_sat           = length(sched_data.stat(i_stat).statistics.sat_name_list);
%     number_of_quasars       = length(sched_data.stat(i_stat).statistics.quasar_name_list);
%     quasar_name_list        = sched_data.stat(i_stat).statistics.quasar_name_list;
    
    % Define number of decimal places of seconds for the source_label:
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

    fprintf(fid_vex, '$SOURCE;\n');
    fprintf(fid_vex, '*\n');
    
    % ##### loop over all scans #####
    for i_scan = 1 : sched_data.number_of_scans
        
        % #### Loop over all stations defined in this scan ####
        for i_stat = 1 : length(sched_data.scan(i_scan).stat)

            % #### Check, if the station participates the scan ####
            if sum(sched_data.scan(i_scan).stat(i_stat).stat_id == stat_id_list) == 1
                % ### Distinguish between observation-type ###
                switch(sched_data.scan(i_scan).obs_type)
                    case 'sat'                    
                        flag_sat_scan       = 1;
                        flag_quasar_scan    = 0;
                    case 'quasar'
                        flag_sat_scan       = 0;
                        flag_quasar_scan    = 1;
                end % switch(sched_data.scan(i_scan).obs_type)
            else % Station does not join this scan
                continue;
            end 


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

                % Write scan header for this satellite scan:
                fprintf(fid_vex, '* ---- Scan %d of %d: %s (NORAD ID: %s) ----\n', i_scan, sched_data.number_of_scans, sat_label, sat_number_str);
                fprintf(fid_vex, '*      Scan duration: %1.2f sec\n', scan_duration_sec);
                fprintf(fid_vex, '*\n');


                % ##### loop over all antenna reposition epochs within this scan (stepwise tracking) #####
                if flag_stepwise_tracking
                    num_of_repos_int = length(sched_data.scan(i_scan).stat(i_stat).epoch) - 1;
                else
                    num_of_repos_int = 1;
                end
                
                for i_epoch = 1 : num_of_repos_int
                    
                    % Loop init.:
                    temp_str = '';
                    temp_str_1 = '';

                    % Get obs. start epoch:
                    [year, mon, day, hr, min, sec] = invjday(sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).jd); % Conversion
%                     [year, mon, day, hr, min, sec] = invjday(sched_data.scan(i_scan).stat(1).epoch(i_epoch).jd); % Start epochs have to be the same for all stations in a scan
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
                    try
                        fprintf(fid_vex, 'def %s.%s;\n', sat_number_str, source_label);
                    catch error
                        error_code = 1;
                        error_msg = 'Source label not available.';  
                        return;
                    end

                    % source_name:
                    try
                        fprintf(fid_vex, '    source_name = %s.%s;\n', sat_number_str, source_label);
                    catch error
                        error_code = 1;
                        error_msg = 'Source label not available.';  
                        return;
                    end

                    % ra:
                    try
                        ra_deg = sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).ra;
                    catch error
                        error_code = 1;
                        error_msg = 'Ra. not available.';  
                        return; 
                    end

                    while(ra_deg < 0)
                        ra_deg = ra_deg + 360.0;
                    end
                    while(ra_deg > 360.0)
                        ra_deg = ra_deg - 360.0;
                    end

                    hms = degrees2dms(ra_deg / 15); % Conversion to hour, min, sec
                    fprintf(fid_vex, '    ra = %02.0fh%02.0fm%09.6fs;', hms(1), hms(2), hms(3));


                    % dec:
                    try
                        % Get dec and conversion to hour, min, sec:
                        hms = degrees2dms(sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).dec);
                    catch error
                        error_code = 1;
                        error_msg = 'Declination not available.';  
                        return;
                    end

                    if (sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).dec < 0) % negative declination value
                        fprintf(fid_vex, ' dec = -%02.0fd%02.0f''%09.6f";', abs(hms(1)), abs(hms(2)), abs(hms(3)));
                    else % positiv declination value
                        fprintf(fid_vex, ' dec = %02.0fd%02.0f''%09.6f";', hms(1), hms(2), hms(3));
                    end

                    % ref_coord_frame:
                    try
                        fprintf(fid_vex, ' ref_coord_frame = %s;\n', vex_para.ref_coord_frame);
                    catch error
                        error_code = 1;
                        error_msg = 'Ref. Coordinate Frame info. not available.';  
                        return;
                    end

                    fprintf(fid_vex, 'enddef;\n');
                    fprintf(fid_vex, '*\n');

                end % for i_epoch = 1 : (length(sched_data.scan(i_scan).stat(i_stat).epoch) - 1)

            end % if flag_sat_scan



            % ######################################################################
            % # Quasars scans
            % ######################################################################
            if flag_quasar_scan

                % init.:
                temp_str = '';
                temp_str_1 = '';

                % Get source label (ref. to $SCHED block):
                source_label = sched_data.scan(i_scan).quasar_name;
                
                % Do not write the source def. to the VEX file, if the label was already written to the VEX file!
                if ~logical(sum(ismember(all_src_labels_str_array, source_label)))

                    % Maintain list of already defined sources:
                    all_src_labels_str_array = union(all_src_labels_str_array, source_label);
                    
                    % Get coordinates:
                    try
                        ra_deg      = sched_data.scan(i_scan).stat(i_stat).epoch.ra;
                        dec_deg     = sched_data.scan(i_scan).stat(i_stat).epoch.dec;
                    catch error
                        error_code = 1;
                        error_msg = ['Source coordinates are not available not available: ', sched_data.scan(i_scan).quasar_name];  
                        return; 
                    end

                    % Write scan header for this quasar scan:
                    fprintf(fid_vex, '* ---- Quasar %s ----\n', source_label);

                    % def. source:
                    try
                        fprintf(fid_vex, 'def %s;\n', source_label);
                    catch error
                        error_code = 1;
                        error_msg = 'Source label not available.';  
                        return;
                    end

                    % source_name:
                    try
                        fprintf(fid_vex, '    source_name = %s;\n', source_label);
                    catch error
                        error_code = 1;
                        error_msg = 'Source label not available.';  
                        return;
                    end

                    % ra:
                    while(ra_deg < 0)
                        ra_deg = ra_deg + 360.0;
                    end
                    while(ra_deg > 360.0)
                        ra_deg = ra_deg - 360.0;
                    end
                    hms = degrees2dms(ra_deg / 15); % Conversion to hour, min, sec
                    fprintf(fid_vex, '    ra = %02.0fh%02.0fm%09.6fs;', hms(1), hms(2), hms(3));

                    % dec:
                    % Convert dec to hour, min, sec:
                    hms = degrees2dms(dec_deg);
                    if (dec_deg < 0) % negative declination value
                        fprintf(fid_vex, ' dec = -%02.0fd%02.0f''%09.6f";', abs(hms(1)), abs(hms(2)), abs(hms(3)));
                    else % positiv declination value
                        fprintf(fid_vex, ' dec = %02.0fd%02.0f''%09.6f";', hms(1), hms(2), hms(3));
                    end

                    % ref_coord_frame:
                    try
                        fprintf(fid_vex, ' ref_coord_frame = %s;\n', vex_para.ref_coord_frame);
                    catch error
                        error_code = 1;
                        error_msg = 'Ref. Coordinate Frame info. not available.';  
                        return;
                    end

                    fprintf(fid_vex, 'enddef;\n');
                    fprintf(fid_vex, '*\n');
                    
                end % if ~logical(sum(ismember(all_src_labels_str_array, source_label)))

            end % if flag_quasar_scan
            
            
            % ##### Write source definitions only for the first station in a scan!!!! #####
            break; % ...the for look over stations in scan! 
            
        end % for i_stat = 1 : length(sched_data.scan(i_scan).stat)
        
    end % for i_scan = 1 : sched_data.number_of_scans

    
    % ##### Additional notes, if available #####
    try
        fprintf(fid_vex, '*   %s\n', vex_para.source_note);
    catch error

    end
    
    fprintf(fid_vex, '* ---------------------------------------------------\n');
    fprintf(fid_vex, '*\n');
    
return;
