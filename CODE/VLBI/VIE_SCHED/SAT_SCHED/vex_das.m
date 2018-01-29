% -------------------------------------------------------------------------
%
%                              vex_das
%
%   Writes the $DAS Block to a VEX File.
%       - 
%
%   Author: 
%       2013-11-28 : Andreas Hellerschmied (heller182@gmx.at)
%   
%   changes       :
%   - 2014-01-13 : Andreas Hellerschmied: Parameter "tape_motion" added.
%   - 2015-06-18, A. Hellerschmied: Updated for VieVS 2.3
%   - 2015-10-20, A. Hellerschmied: Possibility added to write combined VEX
%               files: DAS blocks are written for all stations defined by the input
%               argument "stat_id_list". Multiple def. labels are skipped.
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

function [error_code, error_msg] = vex_das(fid_vex, vex_para, sched_data, stat_id_list)

    % ##### Init #####
    error_code = 0;
    error_msg = '';
    all_labels_str_array = {};
    
    % ##### Wite block #####
    if isempty(stat_id_list) % Error
        error_code = 2;
        error_msg = 'Station ID list is empty';  
        return;
    else
        
        fprintf(fid_vex, '$DAS;\n');
        fprintf(fid_vex, '*\n');
        
        
        % ##### loop over all stations in ID list ##### 
        for i_stat = 1 : length(stat_id_list)
            
            stat_id      = stat_id_list(i_stat);
            station_name = sched_data.stat(stat_id).name; 

            number_of_das = eval(['length(vex_para.', station_name , '.das);']);
    
            % #### Loop over all DAS definitions of the station ####
            for i_das = 1 : number_of_das

                % init:
                temp_str = '';
                temp_str_1 = '';

                % Definition:
                try
                    temp_str = eval(['vex_para.', station_name , '.das(', num2str(i_das), ').label;']);
                    
                    % Do not write the current def. to the VEX file, if the label was already written to the VEX file!
                    if logical(sum(ismember(all_labels_str_array, temp_str)))
                        continue;
                    end
                    all_labels_str_array = union(all_labels_str_array, temp_str);
                    
                    fprintf(fid_vex, 'def %s;\n', temp_str);
                catch error
                    error_code = 1;
                    error_msg = 'DAS Label not available';  
                    return;
                end


                % das_note:
                try
                    number_of_notes = eval(['length(vex_para.', station_name , '.das(', num2str(i_das), ').das_note)']);

                    for i_note = 1 : number_of_notes
                        temp_str = eval(['vex_para.', station_name , '.das(', num2str(i_das), ').das_note(', num2str(i_note), ').str']);
                        fprintf(fid_vex, '* %s\n', temp_str);
                    end

                catch error

                end

                % record_transport_type:
                try
                    temp_str = eval(['vex_para.', station_name , '.das(', num2str(i_das), ').record_transport_type.str;']);
                    fprintf(fid_vex, '    record_transport_type = %s;\n', temp_str);
                catch error
                    error_code = 1;
                    error_msg = '"record_transport_type" parameter not available';  
                    return;
                end

                % electronics_rack_type:
                try
                    temp_str = eval(['vex_para.', station_name , '.das(', num2str(i_das), ').electronics_rack_type.str;']);
                    fprintf(fid_vex, '    electronics_rack_type = %s;\n', temp_str);
                catch error
                    error_code = 1;
                    error_msg = '"electronics_rack_type" parameter not available';  
                    return;
                end

                % recording_system_ID:
                try
                    temp_str = eval(['vex_para.', station_name , '.das(', num2str(i_das), ').recording_system_ID.str;']);
                    fprintf(fid_vex, '    recording_system_ID = %s;\n', temp_str);
                catch error
                    % error_code = 1;
                    % error_msg = '"number_drives" parameter not available';  
                    % return;
                end

                % number_drives:
                try
                    temp_str = eval(['vex_para.', station_name , '.das(', num2str(i_das), ').number_drives.str;']);
                    fprintf(fid_vex, '    number_drives = %s;\n', temp_str);
                catch error
                    % error_code = 1;
                    % error_msg = '"number_drives" parameter not available';  
                    % return;
                end


                % headstack:

                if eval(['isfield(vex_para.',station_name  , '.das(1),', '''headstack_number'')']) % If a headstack is defined in the file "vex_das.cat"
                    number_of_headstack = eval(['length(vex_para.', station_name , '.das(', num2str(i_das), ').headstack_number);']);

                    for i_headstack = 1 : number_of_headstack

                        temp_str = '';
                        temp_str_1 = '';

                        % headstack_number (necessary parameter):
                        try
                            temp_str = eval(['vex_para.', station_name , '.das(', num2str(i_das), ').headstack_number(', num2str(i_headstack), ').str;']);
                        catch error
                            error_code = 1;
                            error_msg = '"headstack_number" parameter not available';  
                            return;
                        end

                        % headstack_function (optional parameter):
                        try
                            temp_str_1 = eval(['vex_para.', station_name , '.das(', num2str(i_das), ').headstack_function(', num2str(i_headstack), ').str;']);
                            temp_str = [temp_str, ' : ', temp_str_1];
                        catch error
                            temp_str = [temp_str, ' :  '];
                        end

                        % drive_number_offset (necessary parameter):
                        try
                            temp_str_1 = eval(['vex_para.', station_name , '.das(', num2str(i_das), ').drive_number_offset(', num2str(i_headstack), ').str;']);
                            temp_str = [temp_str, ' : ', temp_str_1];
                        catch error
                            error_code = 1;
                            error_msg = '"drive_number_offset" parameter not available';  
                            return;
                        end


                        % --- Write "headstack": ---

                        fprintf(fid_vex, '    headstack = %s;\n', temp_str);

                    end % for i_headstack = 1 : number_of_headstack
                end

                % --- tape_motion: ----

                temp_str = '';
                temp_str_1 = '';
                flag_tape_motion_parameter_complete = 1;

                % tape_controll:
                try
                    temp_str = eval(['vex_para.', station_name , '.das(', num2str(i_das), ').tape_controll.str;']);
                catch error
        %             error_code = 1;
        %             error_msg = '"tape_controll" parameter not available';  
        %             return;
                    flag_tape_motion_parameter_complete = 0;
                end

                % tape_early_start_time:
                try
                    temp_str_1 = eval(['vex_para.', station_name , '.das(', num2str(i_das), ').tape_early_start_time.str;']);
                    temp_str = [temp_str, ' : ', temp_str_1];
                catch error
                    flag_tape_motion_parameter_complete = 0;
        %             error_code = 1;
        %             error_msg = '"tape_early_start_time" parameter not available';  
        %             return;
                end

                % tape_late_stop_time:
                try
                    temp_str_1 = eval(['vex_para.', station_name , '.das(', num2str(i_das), ').tape_late_stop_time.str;']);
                    temp_str = [temp_str, ' : ', temp_str_1];
                catch error
                    flag_tape_motion_parameter_complete = 0;
        %             error_code = 1;
        %             error_msg = '"tape_late_stop_time" parameter not available';  
        %             return;
                end

                % min_gap_for_stopping_tape:
                try
                    temp_str_1 = eval(['vex_para.', station_name , '.das(', num2str(i_das), ').min_gap_for_stopping_tape.str;']);
                    temp_str = [temp_str, ' : ', temp_str_1];
                catch error
                    flag_tape_motion_parameter_complete = 0;
        %             error_code = 1;
        %             error_msg = '"min_gap_for_stopping_tape" parameter not available';  
        %             return;
                end

                % --- Write "tape_motion": ---
                if flag_tape_motion_parameter_complete
                    fprintf(fid_vex, '    tape_motion = %s;\n', temp_str);
                else
                    % Place error msg. here if required...
                end

                fprintf(fid_vex, 'enddef;\n');
                fprintf(fid_vex, '*\n');

            end % for i_das = 1 : das
            
        end % for i_stat = 1 : length(stat_id_list)


        % ##### Additional notes, if available: #####
        try
            fprintf(fid_vex, '*    %s\n', vex_para.das_note);
        catch error

        end

        fprintf(fid_vex, '* ---------------------------------------------------\n');
        fprintf(fid_vex, '*\n');
        
    end % if isempty(stat_id_list)
    
return;

