% -------------------------------------------------------------------------
%
%                              vex_mode
%
%   Writes the $MODE Block to a VEX File.
% 
%   Regarding combined VEX files, which can be used as correlator input:
%   - Modes (and mode labels!) from a defined reference station are taken initially and the
%     references valid for other stations are simply added!
%       
%
%   Author: 
%       2013-11-20 : Andreas Hellerschmied (andreas.hellerschmied@geo.tuwien.ac.at)
%   
%   changes       :
%   - 2015-06-19, A. Hellerschmied: Updated for VieVS 2.3
%   - 2015-10-21, A. Hellerschmied: Possibility added to write combined VEX files.
%   - 2015-11-16, A. Hellerschmied: Option added to use a static freq. setup for satellite scans
%   - 2016-05-12, A. Hellerschmied: Added option to use one static mode definition for all satellites ("flag_use_static_freq_setup")
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

function [error_code, error_msg] = vex_mode(fid_vex, vex_para, sched_data, stat_id_list)

    % ##### Init #####
    error_code = 0;
    error_msg = '';
    temp_str_array = {};
    flag_use_static_freq_setup = 0;
    all_labels_str_array = {};
    
    
    % ##### Wite block #####
    if isempty(stat_id_list) % Error
        error_code = 2;
        error_msg = 'Station ID list is empty';  
        return;
    else
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

        fprintf(fid_vex, '$MODE;\n');
        fprintf(fid_vex, '*\n');


        % ##### Check if all required modes are defined in the VEX parameter file: #####
        number_of_mode = eval(['length(vex_para.', station_name_ref , '.mode);']);
        % Mode 1: <station name>.mode(1).<...>  => for satellites
        % Mode 2: <station name>.mode(2).<...>  => for quasars
        
        % If satellites and quasars are observed in this session by this station => two obs. modes have to be defined!
        if (number_of_sat > 1) && (number_of_quasars > 1)
            if number_of_mode ~= 2
                error_code = 1;
                error_msg = 'Two obs. modes have to be defined in the VEX parameter file (for satellite (i_mode=1) and quasar (i_mode=2) obs.). ';
                return;
            end
        end


        % ######################################################################
        % # Satellites
        % ######################################################################
        % Write one mode definition for each observed satellite (to take into account varying frequency configirations):

        for i_sat = 1 : number_of_sat
            
            if flag_use_static_freq_setup && (i_sat == 2)
                break;
            end

            i_mode = 1; % => for satellites

            % Create satellite label:
            sat_name = sat_name_list{i_sat};
            sat_label = sat_name;
            sat_label = sat_label(~isspace(sat_label)); % Remove white space
            
            % Check, if a static frequ. setup/mode should be used:
            try 
                i_frequ = 1; % satellite obs. setup
                temp_str_1 = eval(['vex_para.', station_name_ref , '.freq(', num2str(i_frequ), ').flag_use_static_freq_setup.str;']);
                flag_use_static_freq_setup = logical(str2double(temp_str_1));
            catch
                flag_use_static_freq_setup = 0;
            end
            
            % Get satelite NORAD ID (as char string):
            if ~flag_use_static_freq_setup
                sat_number_str = num2str(sat_num_vector(i_sat));
            end

            % init:
            temp_str = '';
            temp_str_1 = '';

            % --- Definition: ---
            try
                temp_str_1 = eval(['vex_para.', station_name_ref , '.mode(', num2str(i_mode), ').label;']);
                
                if ~flag_use_static_freq_setup
                    temp_str = [sat_number_str, '.', temp_str_1];
                else
                    temp_str = ['sat.', temp_str_1];
                end
                
% #####################################
% CHECK LABEL HERE FOR MULTIPLE ENTRIES                
% #####################################

                % Do not write the current def. to the VEX file, if the label was already written to the VEX file!
                if logical(sum(ismember(all_labels_str_array, temp_str)))
                    continue;
                end
                all_labels_str_array = union(all_labels_str_array, temp_str);                
                
                fprintf(fid_vex, 'def %s;\n', temp_str);
            catch error
                error_code = 1;
                error_msg = 'Mode Label not available';  
                return;
            end

            % --- Satellite Description: ---
            if ~flag_use_static_freq_setup
                try
                    fprintf(fid_vex, '* Satellite name: %s, NORAD ID: %s\n', sat_label, sat_number_str);
                catch error

                end 
            end


            % --- mode_note: ---
            try
                number_of_notes = eval(['length(vex_para.', station_name_ref , '.mode(', num2str(i_mode), ').mode_note)']);

                for i_note = 1 : number_of_notes
                    temp_str = eval(['vex_para.', station_name_ref , '.mode(', num2str(i_mode), ').mode_note(', num2str(i_note), ').str']);
                    fprintf(fid_vex, '* %s\n', temp_str);
                end

            catch error

            end


            % --- ref_procedures: ---
            temp_str_array = {};
            i_temp = 0;
            % ### Loop over all stations ###
            for i_stat = 1 : number_of_stat
                stat_id      = stat_id_list(i_stat);
                station_name = sched_data.stat(stat_id).name;
                stat_label   = sched_data.stat(stat_id).label;
                try
                    number_of_ref_procedures = eval(['length(vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_procedures);']);
                catch error
                    number_of_ref_procedures = 0;
                end
                for i_index = 1 : number_of_ref_procedures
                    temp_str = eval(['vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_procedures(', num2str(i_index), ').str;']);
                    %fprintf(fid_vex, '    ref $PROCEDURES = %s;\n', temp_str);
                    i_temp = i_temp + 1;
                    temp_str_array{i_temp,1} = temp_str;
                    temp_str_array{i_temp,2} = stat_label;
                end % for i_index = 1 : number_of_ref_procedures
            end % for i_stat = 1 : number_of_stat
            print_references(fid_vex, temp_str_array, 'PROCEDURES');
            


            % --- ref_freq: ---
            temp_str_array = {};
            i_temp = 0;
            % ### Loop over all stations ###
            for i_stat = 1 : number_of_stat
                stat_id      = stat_id_list(i_stat);
                station_name = sched_data.stat(stat_id).name;
                stat_label   = sched_data.stat(stat_id).label;
                try
                    number_of_ref_freq = eval(['length(vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_freq);']);
                catch error
                    error_code = 1;
                    error_msg = '"ref_freq" parameter not available';  
                    return;
                end
                for i_index = 1 : number_of_ref_freq
                    temp_str_1 = eval(['vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_freq(', num2str(i_index), ').str;']);
                    if ~flag_use_static_freq_setup
                        temp_str = [sat_number_str, '.', temp_str_1];
                    else
                        temp_str = ['sat.', temp_str_1];
                    end
                    % fprintf(fid_vex, '    ref $FREQ = %s : %s;\n', temp_str, stat_label);
                    i_temp = i_temp + 1;
                    temp_str_array{i_temp,1} = temp_str;
                    temp_str_array{i_temp,2} = stat_label;
                end % for i_index = 1 : number_of_ref_frequ
            end % for i_stat = 1 : number_of_stat
            print_references(fid_vex, temp_str_array, 'FREQ');

            

            % --- ref_if: ---
            temp_str_array = {};
            i_temp = 0;
            % ### Loop over all stations ###
            for i_stat = 1 : number_of_stat
                stat_id      = stat_id_list(i_stat);
                station_name = sched_data.stat(stat_id).name;
                stat_label   = sched_data.stat(stat_id).label;
                try
                    number_of_ref_frequ = eval(['length(vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_if);']);
                catch error
                    error_code = 1;
                    error_msg = '"ref_if" parameter not available';  
                    return;
                end
                for i_index = 1 : number_of_ref_frequ
                    temp_str = eval(['vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_if(', num2str(i_index), ').str;']);
                    % fprintf(fid_vex, '    ref $IF = %s : %s;\n', temp_str, stat_label);
                    i_temp = i_temp + 1;
                    temp_str_array{i_temp,1} = temp_str;
                    temp_str_array{i_temp,2} = stat_label;
                end % for i_index = 1 : number_of_ref_if
            end % for i_stat = 1 : number_of_stat
            print_references(fid_vex, temp_str_array, 'IF');


            % --- ref_bbc: ---
            temp_str_array = {};
            i_temp = 0;
            % ### Loop over all stations ###
            for i_stat = 1 : number_of_stat
                stat_id      = stat_id_list(i_stat);
                station_name = sched_data.stat(stat_id).name;
                stat_label   = sched_data.stat(stat_id).label;
                try
                    number_of_ref_bbc = eval(['length(vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_bbc);']);
                catch error
                    error_code = 1;
                    error_msg = '"ref_bbc" parameter not available';  
                    return;
                end
                for i_index = 1 : number_of_ref_bbc

                    temp_str = eval(['vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_bbc(', num2str(i_index), ').str;']);
                    % fprintf(fid_vex, '    ref $BBC = %s : %s;\n', temp_str, stat_label);
                    i_temp = i_temp + 1;
                    temp_str_array{i_temp,1} = temp_str;
                    temp_str_array{i_temp,2} = stat_label;
                end % for i_index = 1 : number_of_ref_bbc
            end % for i_stat = 1 : number_of_stat
            print_references(fid_vex, temp_str_array, 'BBC');
            
            

            % --- ref_tracks: ---
            temp_str_array = {};
            i_temp = 0;
            % ### Loop over all stations ###
            for i_stat = 1 : number_of_stat
                stat_id      = stat_id_list(i_stat);
                station_name = sched_data.stat(stat_id).name;
                stat_label   = sched_data.stat(stat_id).label;
                try
                    number_of_ref_tracks = eval(['length(vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_tracks);']);
                catch error
                    error_code = 1;
                    error_msg = '"ref_tracks" parameter not available';  
                    return;
                end
                for i_index = 1 : number_of_ref_tracks
                    temp_str = eval(['vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_tracks(', num2str(i_index), ').str;']);
                    % fprintf(fid_vex, '    ref $TRACKS = %s : %s;\n', temp_str, stat_label);
                    i_temp = i_temp + 1;
                    temp_str_array{i_temp,1} = temp_str;
                    temp_str_array{i_temp,2} = stat_label;
                end % for i_index = 1 : number_of_ref_tracks
            end % for i_stat = 1 : number_of_stat
            print_references(fid_vex, temp_str_array, 'TRACKS');

            
            % --- ref_roll: ---
            temp_str_array = {};
            i_temp = 0;
            % ### Loop over all stations ###
            for i_stat = 1 : number_of_stat
                stat_id      = stat_id_list(i_stat);
                station_name = sched_data.stat(stat_id).name;
                stat_label   = sched_data.stat(stat_id).label;
                try
                    number_of_ref_roll = eval(['length(vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_roll);']);
                catch error
                    error_code = 1;
                    error_msg = '"ref_roll" parameter not available';  
                    return;
                end
                for i_index = 1 : number_of_ref_roll
                    temp_str = eval(['vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_roll(', num2str(i_index), ').str;']);
                    % fprintf(fid_vex, '    ref $ROLL = %s : %s;\n', temp_str, stat_label);
                    i_temp = i_temp + 1;
                    temp_str_array{i_temp,1} = temp_str;
                    temp_str_array{i_temp,2} = stat_label;
                end % for i_index = 1 : number_of_ref_roll
            end % for i_stat = 1 : number_of_stat
            print_references(fid_vex, temp_str_array, 'ROLL');
            

            % --- ref_phase_cal_detect: ---
            temp_str_array = {};
            i_temp = 0;
            % ### Loop over all stations ###
            for i_stat = 1 : number_of_stat
                stat_id      = stat_id_list(i_stat);
                station_name = sched_data.stat(stat_id).name;
                stat_label   = sched_data.stat(stat_id).label;
                try
                    number_of_ref_phase_cal_detect = eval(['length(vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_phase_cal_detect);']);
                catch error
                    error_code = 1;
                    error_msg = '"ref_phase_cal_detect" parameter not available';  
                    return;
                end
                for i_index = 1 : number_of_ref_phase_cal_detect
                    temp_str = eval(['vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_phase_cal_detect(', num2str(i_index), ').str;']);
                    % fprintf(fid_vex, '    ref $PHASE_CAL_DETECT = %s : %s;\n', temp_str, stat_label);
                    i_temp = i_temp + 1;
                    temp_str_array{i_temp,1} = temp_str;
                    temp_str_array{i_temp,2} = stat_label;
                end % for i_index = 1 : number_of_ref_phase_cal_detect
            end % for i_stat = 1 : number_of_stat
            print_references(fid_vex, temp_str_array, 'PHASE_CAL_DETECT');
            

            fprintf(fid_vex, 'enddef;\n');
            fprintf(fid_vex, '*\n');



        end % for i_sat = 1 : number_of_sat



        % ######################################################################
        % # Quasars
        % ######################################################################
        % Write one mode definition for all quasar observations:

        if number_of_quasars > 0 % If quasars are observed...

            i_mode = 2; % => for quasars

            % init:
            temp_str = '';
            temp_str_1 = '';
            flag_write_def = true;

            % --- Definition: ---
            try
                temp_str = eval(['vex_para.', station_name_ref , '.mode(', num2str(i_mode), ').label;']);

% #####################################
% CHECK LABEL HERE FOR MULTIPLE ENTRIES                
% #####################################

                % Do not write the current def. to the VEX file, if the label was already written to the VEX file!
                if logical(sum(ismember(all_labels_str_array, temp_str)))
                    flag_write_def = false;
                else
                    all_labels_str_array = union(all_labels_str_array, temp_str);
                    fprintf(fid_vex, 'def %s;\n', temp_str);
                end
            catch error
                error_code = 1;
                error_msg = 'Mode Label not available';  
                return;
            end

            if flag_write_def
            
                % --- Comment: ---
                fprintf(fid_vex, '* Mode for quasar observations.\n');


                % --- mode_note: ---
                try
                    number_of_notes = eval(['length(vex_para.', station_name_ref , '.mode(', num2str(i_mode), ').mode_note)']);

                    for i_note = 1 : number_of_notes
                        temp_str = eval(['vex_para.', station_name_ref , '.mode(', num2str(i_mode), ').mode_note(', num2str(i_note), ').str']);
                        fprintf(fid_vex, '* %s\n', temp_str);
                    end

                catch error

                end


                % --- ref_procedures: ---
                temp_str_array = {};
                i_temp = 0;
                % ### Loop over all stations ###
                for i_stat = 1 : number_of_stat
                    stat_id      = stat_id_list(i_stat);
                    station_name = sched_data.stat(stat_id).name;
                    stat_label   = sched_data.stat(stat_id).label;
                    try
                        number_of_ref_procedures = eval(['length(vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_procedures);']);
                    catch error
                        number_of_ref_procedures = 0;
                    end
                    for i_index = 1 : number_of_ref_procedures
                        temp_str = eval(['vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_procedures(', num2str(i_index), ').str;']);
                        % fprintf(fid_vex, '    ref $PROCEDURES = %s;\n', temp_str);
                        i_temp = i_temp + 1;
                        temp_str_array{i_temp,1} = temp_str;
                        temp_str_array{i_temp,2} = stat_label;
                    end % for i_index = 1 : number_of_ref_procedures
                end % for i_stat = 1 : number_of_stat
                print_references(fid_vex, temp_str_array, 'PROCEDURES');


                % --- ref_freq: ---
                temp_str_array = {};
                i_temp = 0;
                % ### Loop over all stations ###
                for i_stat = 1 : number_of_stat
                    stat_id      = stat_id_list(i_stat);
                    station_name = sched_data.stat(stat_id).name;
                    stat_label   = sched_data.stat(stat_id).label;
                    try
                        number_of_ref_freq = eval(['length(vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_freq);']);
                    catch error
                        error_code = 1;
                        error_msg = '"ref_freq" parameter not available';  
                        return;
                    end
                    for i_index = 1 : number_of_ref_freq
                        temp_str = eval(['vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_freq(', num2str(i_index), ').str;']);
                        % fprintf(fid_vex, '    ref $FREQ = %s : %s;\n', temp_str, stat_label);
                        i_temp = i_temp + 1;
                        temp_str_array{i_temp,1} = temp_str;
                        temp_str_array{i_temp,2} = stat_label;
                    end % for i_index = 1 : number_of_ref_frequ
                end % for i_stat = 1 : number_of_stat
                print_references(fid_vex, temp_str_array, 'FREQ');


                % --- ref_if: ---
                temp_str_array = {};
                i_temp = 0;
                % ### Loop over all stations ###
                for i_stat = 1 : number_of_stat
                    stat_id      = stat_id_list(i_stat);
                    station_name = sched_data.stat(stat_id).name;
                    stat_label   = sched_data.stat(stat_id).label;
                    try
                        number_of_ref_frequ = eval(['length(vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_if);']);
                    catch error
                        error_code = 1;
                        error_msg = '"ref_if" parameter not available';  
                        return;
                    end
                    for i_index = 1 : number_of_ref_frequ
                        temp_str = eval(['vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_if(', num2str(i_index), ').str;']);
                        % fprintf(fid_vex, '    ref $IF = %s : %s;\n', temp_str, stat_label);
                        i_temp = i_temp + 1;
                        temp_str_array{i_temp,1} = temp_str;
                        temp_str_array{i_temp,2} = stat_label;
                    end % for i_index = 1 : number_of_ref_if
                end % for i_stat = 1 : number_of_stat
                print_references(fid_vex, temp_str_array, 'IF');


                % --- ref_bbc: ---
                temp_str_array = {};
                i_temp = 0;
                % ### Loop over all stations ###
                for i_stat = 1 : number_of_stat
                    stat_id      = stat_id_list(i_stat);
                    station_name = sched_data.stat(stat_id).name;
                    stat_label   = sched_data.stat(stat_id).label;
                    try
                        number_of_ref_bbc = eval(['length(vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_bbc);']);
                    catch error
                        error_code = 1;
                        error_msg = '"ref_bbc" parameter not available';  
                        return;
                    end
                    for i_index = 1 : number_of_ref_bbc
                        temp_str = eval(['vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_bbc(', num2str(i_index), ').str;']);
                        % fprintf(fid_vex, '    ref $BBC = %s : %s;\n', temp_str, stat_label);
                        i_temp = i_temp + 1;
                        temp_str_array{i_temp,1} = temp_str;
                        temp_str_array{i_temp,2} = stat_label;
                    end % for i_index = 1 : number_of_ref_bbc
                end % for i_stat = 1 : number_of_stat
                print_references(fid_vex, temp_str_array, 'BBC');


                % --- ref_tracks: ---
                temp_str_array = {};
                i_temp = 0;
                % ### Loop over all stations ###
                for i_stat = 1 : number_of_stat
                    stat_id      = stat_id_list(i_stat);
                    station_name = sched_data.stat(stat_id).name;
                    stat_label   = sched_data.stat(stat_id).label;
                    try
                        number_of_ref_tracks = eval(['length(vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_tracks);']);
                    catch error
                        error_code = 1;
                        error_msg = '"ref_tracks" parameter not available';  
                        return;
                    end
                    for i_index = 1 : number_of_ref_tracks
                        temp_str = eval(['vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_tracks(', num2str(i_index), ').str;']);
                        % fprintf(fid_vex, '    ref $TRACKS = %s : %s;\n', temp_str, stat_label);
                        i_temp = i_temp + 1;
                        temp_str_array{i_temp,1} = temp_str;
                        temp_str_array{i_temp,2} = stat_label;
                    end % for i_index = 1 : number_of_ref_tracks
                end % for i_stat = 1 : number_of_stat
                print_references(fid_vex, temp_str_array, 'TRACKS');


                % --- ref_roll: ---
                temp_str_array = {};
                i_temp = 0;
                % ### Loop over all stations ###
                for i_stat = 1 : number_of_stat
                    stat_id      = stat_id_list(i_stat);
                    station_name = sched_data.stat(stat_id).name;
                    stat_label   = sched_data.stat(stat_id).label;
                    try
                        number_of_ref_roll = eval(['length(vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_roll);']);
                    catch error
                        error_code = 1;
                        error_msg = '"ref_roll" parameter not available';  
                        return;
                    end
                    for i_index = 1 : number_of_ref_roll
                        temp_str = eval(['vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_roll(', num2str(i_index), ').str;']);
                        % fprintf(fid_vex, '    ref $ROLL = %s : %s;\n', temp_str, stat_label);
                        i_temp = i_temp + 1;
                        temp_str_array{i_temp,1} = temp_str;
                        temp_str_array{i_temp,2} = stat_label;
                    end % for i_index = 1 : number_of_ref_roll
                end % for i_stat = 1 : number_of_stat
                print_references(fid_vex, temp_str_array, 'ROLL');


                % --- ref_phase_cal_detect: ---
                temp_str_array = {};
                i_temp = 0;
                % ### Loop over all stations ###
                for i_stat = 1 : number_of_stat
                    stat_id      = stat_id_list(i_stat);
                    station_name = sched_data.stat(stat_id).name;
                    stat_label   = sched_data.stat(stat_id).label;
                    try
                        number_of_ref_phase_cal_detect = eval(['length(vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_phase_cal_detect);']);
                    catch error
                        error_code = 1;
                        error_msg = '"ref_phase_cal_detect" parameter not available';  
                        return;
                    end
                    for i_index = 1 : number_of_ref_phase_cal_detect
                        temp_str = eval(['vex_para.', station_name , '.mode(', num2str(i_mode), ').ref_phase_cal_detect(', num2str(i_index), ').str;']);
                        % fprintf(fid_vex, '    ref $PHASE_CAL_DETECT = %s : %s;\n', temp_str, stat_label);
                        i_temp = i_temp + 1;
                        temp_str_array{i_temp,1} = temp_str;
                        temp_str_array{i_temp,2} = stat_label;
                    end % for i_index = 1 : number_of_ref_phase_cal_detect
                end % for i_stat = 1 : number_of_stat
                print_references(fid_vex, temp_str_array, 'PHASE_CAL_DETECT');


                fprintf(fid_vex, 'enddef;\n');
                fprintf(fid_vex, '*\n');
                
            end % if flag_write_def

        end % if number_of_quasars > 0



        % ##### Additional notes, if available #####
        try
            fprintf(fid_vex, '*    %s\n', vex_para.mode_note);
        catch error

        end

        fprintf(fid_vex, '* ---------------------------------------------------\n');
        fprintf(fid_vex, '*\n');
        
        
    end % if isempty(stat_id_list) 
    
    return;
end


% ############### Sub routines ####################

function print_references(fid_vex, temp_str_array, block_name_str)
    % Get unique labels:
    unique_labels = unique(temp_str_array(:,1));
    % Loop over all unique labels:
    for i_label = 1 : length(unique_labels)
        % Get stations where this unique label is used:
        station_labels = temp_str_array(find(strcmp(unique_labels(i_label), temp_str_array(:,1))),2);
        temp_str = unique_labels{i_label}; 
        % Append station labels:
        for i_stat_lab = 1 : length(station_labels)
            temp_str = [temp_str, ' : ', station_labels{i_stat_lab}];
        end
        % Print the ref. line:
        fprintf(fid_vex, '    ref $%s = %s;\n', block_name_str, temp_str);
    end
end

