% -------------------------------------------------------------------------
%
%                              vex_roll
%
%   Writes the $ROLL Block to a VEX File.
%       - 
%
%   Author: 
%       2013-11-20 : Andreas Hellerschmied (andreas.hellerschmied@geo.tuwien.ac.at)
%   
%   changes       :
%   - 2015-06-19, A. Hellerschmied: Updated for VieVS 2.3
%   - 2015-10-20, A. Hellerschmied: Possibility added to write combined VEX
%               files: ROLL blocks are written for all stations defined by the input
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

function [error_code, error_msg] = vex_roll(fid_vex, vex_para, sched_data, stat_id_list)

    % Init
    error_code = 0;
    error_msg = '';
    all_labels_str_array = {};
    
    % ##### Wite block #####
    if isempty(stat_id_list) % Error
        error_code = 2;
        error_msg = 'Station ID list is empty';  
        return;
    else
        
        fprintf(fid_vex, '$ROLL;\n');
        fprintf(fid_vex, '*\n');
        
        
        % ##### loop over all stations in ID list ##### 
        for i_stat = 1 : length(stat_id_list)
            
            stat_id      = stat_id_list(i_stat);
            station_name = sched_data.stat(stat_id).name; 
    
            number_of_roll = eval(['length(vex_para.', station_name , '.roll);']);


            % #### Loop over all ROLL definitions of the station ####
            for i_roll = 1 : number_of_roll

                % Definition:
                try
                    temp_str = eval(['vex_para.', station_name , '.roll(', num2str(i_roll), ').label;']);
                    
                    % Do not write the current def. to the VEX file, if the label was already written to the VEX file!
                    if logical(sum(ismember(all_labels_str_array, temp_str)))
                        continue;
                    end
                    all_labels_str_array = union(all_labels_str_array, temp_str);
                    
                    fprintf(fid_vex, 'def %s;\n', temp_str);
                catch error
                    error_code = 1;
                    error_msg = 'roll Label not available';  
                    return;
                end


                % roll_note:
                try
                    number_of_notes = eval(['length(vex_para.', station_name , '.roll(', num2str(i_roll), ').roll_note)']);
                    for i_note = 1 : number_of_notes
                        temp_str = eval(['vex_para.', station_name , '.roll(', num2str(i_roll), ').roll_note(', num2str(i_note), ').str']);
                        fprintf(fid_vex, '* %s\n', temp_str);
                    end

                catch error

                end

                % roll:
                try
                    temp_str = eval(['vex_para.', station_name , '.roll(', num2str(i_roll), ').roll.str;']);
                    fprintf(fid_vex, '    roll  = %s;\n', temp_str);
                catch error
                    error_code = 1;
                    error_msg = '"roll" parameters not available';  
                    return;
                end

                fprintf(fid_vex, 'enddef;\n');
                fprintf(fid_vex, '*\n');

            end % for i_roll = 1 : number_of_roll
            
        end % for i_stat = 1 : length(stat_id_list)


        % ##### Additional notes, if available: #####
        try
            fprintf(fid_vex, '*    %s\n', vex_para.roll_note);
        catch error

        end

        fprintf(fid_vex, '* ---------------------------------------------------\n');
        fprintf(fid_vex, '*\n');
        
    end % if isempty(stat_id_list)
            
return;

