% -------------------------------------------------------------------------
%
%                              vex_bbc
%
%   Writes the $BBC Block to a VEX File.
%       - 
%
%   Author: 
%       2013-11-19 : Andreas Hellerschmied
%   
%   changes       :
%   - 2015-06-18, A. Hellerschmied: Updated for VieVS 2.3
%   - 2015-10-20, A. Hellerschmied: Possibility added to write combined VEX
%               files: BBC blocks are written for all stations defined by the input
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

function [error_code, error_msg] = vex_bbc(fid_vex, vex_para, sched_data, stat_id_list)

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
        
        fprintf(fid_vex, '$BBC;\n');
        fprintf(fid_vex, '*\n');
        
        
        % ##### loop over all stations in ID list ##### 
        for i_stat = 1 : length(stat_id_list)
            
            stat_id      = stat_id_list(i_stat);
            station_name = sched_data.stat(stat_id).name; 
    
            number_of_bbc = eval(['length(vex_para.', station_name , '.bbc);']);


            % #### Loop over all BBC definitions of the station ####
            for i_bbc = 1 : number_of_bbc

                % init:
                temp_str = '';
                temp_str_1 = '';

                % Definition:
                try
                    temp_str = eval(['vex_para.', station_name , '.bbc(', num2str(i_bbc), ').label;']);
                    
                    % Do not write the current def. to the VEX file, if the label was already written to the VEX file!
                    if logical(sum(ismember(all_labels_str_array, temp_str)))
                        continue;
                    end
                    all_labels_str_array = union(all_labels_str_array, temp_str);
                    
                    fprintf(fid_vex, 'def %s;\n', temp_str);
                catch error
                    error_code = 1;
                    error_msg = 'BBC Label not available';  
                    return;
                end


                % bbc_note:
                try
                    number_of_notes = eval(['length(vex_para.', station_name , '.bbc(', num2str(i_bbc), ').bbc_note)']);

                    for i_note = 1 : number_of_notes
                        temp_str = eval(['vex_para.', station_name , '.bbc(', num2str(i_bbc), ').bbc_note(', num2str(i_note), ').str']);
                        fprintf(fid_vex, '* %s\n', temp_str);
                    end

                catch error

                end

                % bbc_assign:
                number_of_bbc_assign = eval(['length(vex_para.', station_name , '.bbc(', num2str(i_bbc), ').logic_bbc_link);']);

                for i_bbc_assign = 1 : number_of_bbc_assign

                    % --- necessary parameters: ---

                    % logic_bbc_link:
                    try
                        temp_str = eval(['vex_para.', station_name , '.bbc(', num2str(i_bbc), ').logic_bbc_link(', num2str(i_bbc_assign), ').str;']);
                    catch error
                        error_code = 1;
                        error_msg = '"logic_bbc_link" parameter not available';  
                        return;
                    end

                    % physical_bbc_number:
                    try
                        temp_str_1 = eval(['vex_para.', station_name , '.bbc(', num2str(i_bbc), ').physical_bbc_number(', num2str(i_bbc_assign), ').str;']);
                        temp_str = [temp_str, ' : ', temp_str_1];
                    catch error
                        error_code = 1;
                        error_msg = '"physical_bbc_number" parameter not available';  
                        return;
                    end

                    % logic_if_link:
                    try
                        temp_str_1 = eval(['vex_para.', station_name , '.bbc(', num2str(i_bbc), ').logic_if_link(', num2str(i_bbc_assign), ').str;']);
                        temp_str = [temp_str, ' : ', temp_str_1];
                    catch error
                        error_code = 1;
                        error_msg = '"logic_if_link" parameter not available';  
                        return;
                    end


                    % --- Write "bbc_assign": ---

                    fprintf(fid_vex, '    BBC_assign = %s;\n', temp_str);

                end % for i_bbc_assign = 1 : number_of_bbc_assign

                fprintf(fid_vex, 'enddef;\n');
                fprintf(fid_vex, '*\n');

            end % for i_bbc = 1 : number_of_bbc
            
        end % for i_stat = 1 : length(stat_id_list)


        % ##### Additional notes, if available: #####
        try
            fprintf(fid_vex, '*    %s\n', vex_para.bbc_note);
        catch error

        end

        fprintf(fid_vex, '* ---------------------------------------------------\n');
        fprintf(fid_vex, '*\n');
        
    end % if isempty(stat_id_list)
    
return;

