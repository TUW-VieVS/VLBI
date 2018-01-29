% -------------------------------------------------------------------------
%
%                              vex_procedures
%
%   Writes the $PROCEDURES Block to a VEX File.
%
%   Author: 
%       2013-11-28 : Andreas Hellerschmied (heller182@gmx.at)
%   
%   changes       :
%   - 2015-06-18, A. Hellerschmied: Updated for VieVS 2.3
%   - 2015-10-20, A. Hellerschmied: Possibility added to write combined VEX files (write multiple proc. definitions without repetitions of the same def. label)
%           
%
%   inputs        :
%   - fid_vex           : File ID of VEX File
%   - vex_para          : Data from VEX paramer file, structure
%   - sched_data        : scheduling data structure
%   - i_stat            : station count
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

function [error_code, error_msg] = vex_procedures(fid_vex, vex_para, sched_data, stat_id_list)


    % ##### Init #####
    error_code = 0;
    error_msg = '';
    all_labels_str_array = {};
    flag_first_proc_def = 1;
    

    % ##### Wite block #####
    if isempty(stat_id_list) % Error
        error_code = 2;
        error_msg = 'Station ID list is empty';  
        return;
    else
        
        % ##### loop over all stations in ID list ##### 
        for i_stat = 1 : length(stat_id_list)
            
            stat_id      = stat_id_list(i_stat);
            station_name = sched_data.stat(stat_id).name; 
            
            % If the label == "none" => do not write a procedure definiton for this station!
            try
                proc_1_label_str = eval(['vex_para.', station_name , '.procedures(1).label;']);
            catch error
                error_code = 1;
                error_msg = 'PROCEDURES Label not available';  
                return;
            end
            if strncmpi(proc_1_label_str, 'none', 4) % No $PROCEDURES Block should be written!
                continue;
            end
            
            if flag_first_proc_def % write the Block headline only once!
                fprintf(fid_vex, '$PROCEDURES;\n');
                fprintf(fid_vex, '*\n');
                flag_first_proc_def = 0;
            end
            
            number_of_proc = eval(['length(vex_para.', station_name , '.procedures);']);
            
            % #### Loop over all procedure definitions ####
            for i_proc = 1 : number_of_proc
                
                % init:
                temp_str = '';
                temp_str_1 = '';

                % Definition:
                try
                    temp_str = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').label;']);
                    % Do not write the current def. to the VEX file, if the label was already written to the VEX file!
                    if logical(sum(ismember(all_labels_str_array, temp_str)))
                        continue;
                    end
                    all_labels_str_array = union(all_labels_str_array, temp_str);
                    fprintf(fid_vex, 'def %s;\n', temp_str);
                catch error
                    error_code = 1;
                    error_msg = 'PROCEDURES Label not available';  
                    return;
                end

                % proceures_note:
                try
                    number_of_notes = eval(['length(vex_para.', station_name , '.procedures(', num2str(i_proc), ').procedures_note)']);

                    for i_note = 1 : number_of_notes
                        temp_str = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').procedures_note(', num2str(i_note), ').str']);
                        fprintf(fid_vex, '* %s\n', temp_str);
                    end

                catch error

                end


                % ------ Optional parameters: ------

                % tape_change:
                try
                    temp_str = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').tape_change.str;']);
                    fprintf(fid_vex, '    tape_change = %s;\n', temp_str);
                catch error
                    error_code = 1;
                    error_msg = '"tape_change" parameter not available';  
                    return;
                end

                % headstack_motion:
                try
                    temp_str = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').headstack_motion.str;']);
                    fprintf(fid_vex, '    headstack_motion = %s;\n', temp_str);
                catch error
                    error_code = 1;
                    error_msg = '"headstack_motion" parameter not available';  
                    return;
                end

                % new_source_command:
                try
                    temp_str = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').new_source_command.str;']);
                    fprintf(fid_vex, '    new_source_command = %s;\n', temp_str);
                catch error
                    error_code = 1;
                    error_msg = '"new_source_command" parameter not available';  
                    return;
                end

                % new_tape_setup:
                try
                    temp_str = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').new_tape_setup.str;']);
                    fprintf(fid_vex, '    new_tape_setup = %s;\n', temp_str);
                catch error
                    error_code = 1;
                    error_msg = '"new_tape_setup" parameter not available';  
                    return;
                end



                % ------ Optional parameters: ------

                % procedure_name_prefix:
                try
                    temp_str = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').procedure_name_prefix.str;']);
                    fprintf(fid_vex, '    procedure_name_prefix = %s;\n', temp_str);
                catch error
                    error_code = 1;
                    error_msg = '"procedure_name_prefix" parameter not available';  
                    return;
                end

                % setup_always:

                temp_str = '';
                temp_str_1 = '';
                flag_write_line = 1;

                % setup_system_for_each_obs:
                try
                    temp_str = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').setup_system_for_each_obs.str;']);
                catch error
                    flag_write_line = 0;
                end

                % time_to_setup_system:
                try
                    temp_str_1 = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').time_to_setup_system.str;']);
                    temp_str = [temp_str, ' : ', temp_str_1];
                catch error
                    flag_write_line = 0;
                end

                if (flag_write_line)
                    fprintf(fid_vex, '    setup_always = %s;\n', temp_str);
                end


                % parity_check:

                temp_str = '';
                temp_str_1 = '';
                flag_write_line = 1;

                % do_parity_check:
                try
                    temp_str = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').do_parity_check.str;']);
                catch error
                    flag_write_line = 0;
                end

                % time_to_do_parity_check:
                try
                    temp_str_1 = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').time_to_do_parity_check.str;']);
                    temp_str = [temp_str, ' : ', temp_str_1];
                catch error
                    flag_write_line = 0;
                end

                if (flag_write_line)
                    fprintf(fid_vex, '    parity_check = %s;\n', temp_str);
                end


                % tape_prepass:

                temp_str = '';
                temp_str_1 = '';
                flag_write_line = 1;

                % do_tape_prepass:
                try
                    temp_str = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').do_tape_prepass.str;']);
                catch error
                    flag_write_line = 0;
                end

                % time_to_do_tape_prepass:
                try
                    temp_str_1 = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').time_to_do_tape_prepass.str;']);
                    temp_str = [temp_str, ' : ', temp_str_1];
                catch error
                    flag_write_line = 0;
                end

                if (flag_write_line)
                    fprintf(fid_vex, '    tape_prepass = %s;\n', temp_str);
                end


                % preob_cal:

                temp_str = '';
                temp_str_1 = '';
                flag_write_line = 1;

                % do_preob_cal:
                try
                    temp_str = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').do_preob_cal.str;']);
                catch error
                    flag_write_line = 0;
                end

                % time_needed_for_preob_cal:
                try
                    temp_str_1 = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').time_needed_for_preob_cal.str;']);
                    temp_str = [temp_str, ' : ', temp_str_1];
                catch error
                    flag_write_line = 0;
                end

                % preob_prodecure_name:
                try
                    temp_str_1 = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').preob_prodecure_name.str;']);
                    temp_str = [temp_str, ' : ', temp_str_1];
                catch error
                    flag_write_line = 0;
                end

                if (flag_write_line)
                    fprintf(fid_vex, '    preob_cal = %s;\n', temp_str);
                end


                % midob_cal:

                temp_str = '';
                temp_str_1 = '';
                flag_write_line = 1;

                % do_midob_cal:
                try
                    temp_str = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').do_midob_cal.str;']);
                catch error
                    flag_write_line = 0;
                end

                % time_needed_for_midob_cal:
                try
                    temp_str_1 = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').time_needed_for_midob_cal.str;']);
                    temp_str = [temp_str, ' : ', temp_str_1];
                catch error
                    flag_write_line = 0;
                end

                % midob_prodecure_name:
                try
                    temp_str_1 = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').midob_prodecure_name.str;']);
                    temp_str = [temp_str, ' : ', temp_str_1];
                catch error
                    flag_write_line = 0;
                end

                if (flag_write_line)
                    fprintf(fid_vex, '    midob_cal = %s;\n', temp_str);
                end


                % postob_cal:

                temp_str = '';
                temp_str_1 = '';
                flag_write_line = 1;

                % do_postob_cal:
                try
                    temp_str = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').do_postob_cal.str;']);
                catch error
                    flag_write_line = 0;
                end

                % time_needed_for_postob_cal:
                try
                    temp_str_1 = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').time_needed_for_postob_cal.str;']);
                    temp_str = [temp_str, ' : ', temp_str_1];
                catch error
                    flag_write_line = 0;
                end

                % postob_prodecure_name:
                try
                    temp_str_1 = eval(['vex_para.', station_name , '.procedures(', num2str(i_proc), ').postob_prodecure_name.str;']);
                    temp_str = [temp_str, ' : ', temp_str_1];
                catch error
                    flag_write_line = 0;
                end

                if (flag_write_line)
                    fprintf(fid_vex, '    postob_cal = %s;\n', temp_str);
                end

                fprintf(fid_vex, 'enddef;\n');
                fprintf(fid_vex, '*\n');

            end % for i_proc = 1 : number_of_proc
              
        end % for i_stat = 1 : length(stat_id_list)

        
        % ##### Additional notes, if available: #####
        try
            fprintf(fid_vex, '*    %s\n', vex_para.procedures_note);
        catch error

        end
        
        fprintf(fid_vex, '* ---------------------------------------------------\n');
        fprintf(fid_vex, '*\n');
    
    end % if isempty(stat_id_list)
    
return;

