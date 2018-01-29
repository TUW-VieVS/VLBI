% -------------------------------------------------------------------------
%
%                              vex_if
%
%   Writes the $IF Block to a VEX File.
%       - 
%
%   Author: 
%       2013-11-19 : Andreas Hellerschmied (heller182@gmx.at)
%   
%   changes       :
%   - 2015-06-18, A. Hellerschmied: Updated for VieVS 2.3
%   - 2015-10-20, A. Hellerschmied: Possibility added to write combined VEX
%               files: IF blocks are written for all stations defined by the input
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

function [error_code, error_msg] = vex_if(fid_vex, vex_para, sched_data, stat_id_list)

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
        
        fprintf(fid_vex, '$IF;\n');
        fprintf(fid_vex, '*\n');
        
        
        % ##### loop over all stations in ID list ##### 
        for i_stat = 1 : length(stat_id_list)
            
            stat_id      = stat_id_list(i_stat);
            station_name = sched_data.stat(stat_id).name; 
    
            number_of_if = eval(['length(vex_para.', station_name , '.if);']);


            % #### Loop over all IF definitions of the station ####
            for i_if = 1 : number_of_if

                % init:
                temp_str = '';
                temp_str_1 = '';

                % Definition:
                try
                    temp_str = eval(['vex_para.', station_name , '.if(', num2str(i_if), ').label;']);
                    
                    % Do not write the current def. to the VEX file, if the label was already written to the VEX file!
                    if logical(sum(ismember(all_labels_str_array, temp_str)))
                        continue;
                    end
                    all_labels_str_array = union(all_labels_str_array, temp_str);
                    
                    fprintf(fid_vex, 'def %s;\n', temp_str);
                catch error
                    error_code = 1;
                    error_msg = 'if Label not available';  
                    return;
                end


                % if_note:
                try
                    number_of_notes = eval(['length(vex_para.', station_name , '.if(', num2str(i_if), ').if_note)']);
                    for i_note = 1 : number_of_notes
                        temp_str = eval(['vex_para.', station_name , '.if(', num2str(i_if), ').if_note(', num2str(i_note), ').str']);
                        fprintf(fid_vex, '* %s\n', temp_str);
                    end

                catch error

                end

                % if_def:
                number_of_if_def = eval(['length(vex_para.', station_name , '.if(', num2str(i_if), ').if_id);']);

                for i_if_def = 1 : number_of_if_def

                    % --- necessary parameters: ---

                    % if_id:
                    try
                        temp_str = eval(['vex_para.', station_name , '.if(', num2str(i_if), ').if_id(', num2str(i_if_def), ').str;']);
                    catch error
                        error_code = 1;
                        error_msg = '"if_id" parameter not available';  
                        return;
                    end

                    % physical_if_name:
                    try
                        temp_str_1 = eval(['vex_para.', station_name , '.if(', num2str(i_if), ').physical_if_name(', num2str(i_if_def), ').str;']);
                        temp_str = [temp_str, ' : ', temp_str_1];
                    catch error
                        error_code = 1;
                        error_msg = '"physical_if_name" parameter not available';  
                        return;
                    end

                    % polarization:
                    try
                        temp_str_1 = eval(['vex_para.', station_name , '.if(', num2str(i_if), ').polarization(', num2str(i_if_def), ').str;']);
                        temp_str = [temp_str, ' : ', temp_str_1];
                    catch error
                        error_code = 1;
                        error_msg = '"polarization" parameter not available';  
                        return;
                    end

                    % total_effective_lo_of_if:
                    try
                        temp_str_1 = eval(['vex_para.', station_name , '.if(', num2str(i_if), ').total_effective_lo_of_if(', num2str(i_if_def), ').str;']);
                        temp_str = [temp_str, ' : ', temp_str_1];
                    catch error
                        error_code = 1;
                        error_msg = '"total_effective_lo_of_if" parameter not available';  
                        return;
                    end

                    % net_sideband_of_if:
                    try
                        temp_str_1 = eval(['vex_para.', station_name , '.if(', num2str(i_if), ').net_sideband_of_if(', num2str(i_if_def), ').str;']);
                        temp_str = [temp_str, ' : ', temp_str_1];
                    catch error
                        error_code = 1;
                        error_msg = '"net_sideband_of_if" parameter not available';  
                        return;
                    end

                    % --- optional parameters: ---

                    % phase_cal_frequ_interva:
                    try
                        temp_str_1 = eval(['vex_para.', station_name , '.if(', num2str(i_if), ').phase_cal_frequ_interva(', num2str(i_if_def), ').str;']);
                        if (length(temp_str_1) > 0)
                            temp_str = [temp_str, ' : ', temp_str_1];
                        end
                    catch error

                    end

                    % phase_cal_base_frequency:
                    try
                        temp_str_1 = eval(['vex_para.', station_name , '.if(', num2str(i_if), ').phase_cal_base_frequency(', num2str(i_if_def), ').str;']);
                        if (length(temp_str_1) > 0)
                            temp_str = [temp_str, ' : ', temp_str_1];
                        end
                    catch error

                    end

                    % --- Write "if_def": ---

                    fprintf(fid_vex, '    if_def = %s;\n', temp_str);

                end % for i_if_def = 1 : number_of_if_def

                fprintf(fid_vex, 'enddef;\n');
                fprintf(fid_vex, '*\n');

            end % for i_if = 1 : number_of_if
            
        end % for i_stat = 1 : length(stat_id_list)


        % ##### Additional notes, if available: #####
        try
            fprintf(fid_vex, '*    %s\n', vex_para.if_note);
        catch error

        end

        fprintf(fid_vex, '* ---------------------------------------------------\n');
        fprintf(fid_vex, '*\n');
    
    end % if isempty(stat_id_list)
    
return;

