% -------------------------------------------------------------------------
%
%                              vex_tracks
%
%   Writes the $TRACKS Block to a VEX File.
%       - 
%
%   Author: 
%       2013-11-19 : Andreas Hellerschmied
%   
%   changes       :
%   - 2015-06-18, A. Hellerschmied: Updated for VieVS 2.3
%   - 2015-10-20, A. Hellerschmied: Possibility added to write combined VEX
%               files: TRACKS blocks are written for all stations defined by the input
%               argument "stat_id_list". Multiple def. labels are skipped.
%   - 2016-05-17, A. Hellerschmied: fanout definitons are now optional
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

function [error_code, error_msg] = vex_tracks(fid_vex, vex_para, sched_data, stat_id_list)

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
        
        fprintf(fid_vex, '$TRACKS;\n');
        fprintf(fid_vex, '*\n');
        
        
        % ##### loop over all stations in ID list ##### 
        for i_stat = 1 : length(stat_id_list)
            
            stat_id      = stat_id_list(i_stat);
            station_name = sched_data.stat(stat_id).name; 
    
            number_of_tracks = eval(['length(vex_para.', station_name , '.tracks);']);


            % #### Loop over all TRACKS definitions of the station ####
            for i_tracks = 1 : number_of_tracks

                % init:
                temp_str = '';
                temp_str_1 = '';

                % --- Definition: ---
                try
                    temp_str = eval(['vex_para.', station_name , '.tracks(', num2str(i_tracks), ').label;']);
                    
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


                % --- tracks_note: ---
                try
                    number_of_notes = eval(['length(vex_para.', station_name , '.tracks(', num2str(i_tracks), ').tracks_note)']);

                    for i_note = 1 : number_of_notes
                        temp_str = eval(['vex_para.', station_name , '.tracks(', num2str(i_tracks), ').tracks_note(', num2str(i_note), ').str']);
                        fprintf(fid_vex, '* %s\n', temp_str);
                    end

                catch error

                end

                % --- track_frame_format: ---
                try
                    temp_str = eval(['vex_para.', station_name , '.tracks(', num2str(i_tracks), ').track_frame_format.str;']);
                    fprintf(fid_vex, '    track_frame_format = %s;\n', temp_str);
                catch error
                    error_code = 1;
                    error_msg = '"track_frame_format" parameter not available';  
                    return;
                end

                % --- data_modulation: ---
                try
                    temp_str = eval(['vex_para.', station_name , '.tracks(', num2str(i_tracks), ').data_modulation.str;']);
                    fprintf(fid_vex, '    data_modulation = %s;\n', temp_str);
                catch error
        %             error_code = 1;
        %             error_msg = '"data_modulation" parameter not available';  
        %             return;
                end


                % --- fanout_def: ---
                try
                    number_of_fanout_def = eval(['length(vex_para.', station_name , '.tracks(', num2str(i_tracks), ').can_id);']);
                    flag_tracks_available = true;
                catch
                    flag_tracks_available = false;
                end

                if flag_tracks_available
                    for i_fanout_def = 1 : number_of_fanout_def

                        temp_str = '';
                        temp_str_1 = '';

                        % --- optional parameters ---

                        % sub_pass_id:
                        try
                            temp_str = eval(['vex_para.', station_name , '.tracks(', num2str(i_tracks), ').sub_pass_id(', num2str(i_fanout_def), ').str;']);
                        catch error
                            temp_str = '  ';
                        end

                        % --- necessary parameters: ---

                        % can_id:
                        try
                            temp_str_1 = eval(['vex_para.', station_name , '.tracks(', num2str(i_tracks), ').can_id(', num2str(i_fanout_def), ').str;']);
                            temp_str = [temp_str, ' : ', temp_str_1];
                        catch error
                            error_code = 1;
                            error_msg = '"can_id" parameter not available';  
                            return;
                        end

                        % sign_mag:
                        try
                            temp_str_1 = eval(['vex_para.', station_name , '.tracks(', num2str(i_tracks), ').sign_mag(', num2str(i_fanout_def), ').str;']);
                            temp_str_1 = sprintf('%4s', temp_str_1);
                            temp_str = [temp_str, ' : ', temp_str_1];
                        catch error
                            error_code = 1;
                            error_msg = '"sign_mag" parameter not available';  
                            return;
                        end

                        % headstack_number:
                        try
                            temp_str_1 = eval(['vex_para.', station_name , '.tracks(', num2str(i_tracks), ').headstack_number(', num2str(i_fanout_def), ').str;']);
                            temp_str = [temp_str, ' : ', temp_str_1];
                        catch error
                            error_code = 1;
                            error_msg = '"headstack_number" parameter not available';  
                            return;
                        end

                        % multiplex_track
                        number_of_m_track = eval(['length(vex_para.', station_name , '.tracks(', num2str(i_tracks), ').multiplex_track(', num2str(i_fanout_def), ').track);']);

                        for i_m_track = 1 : number_of_m_track

                            try
                                temp_str_1 = eval(['vex_para.', station_name , '.tracks(', num2str(i_tracks), ').multiplex_track(', num2str(i_fanout_def), ').track(', num2str(i_m_track), ').str;']);
                                temp_str_1 = sprintf('%2s', temp_str_1);
                                temp_str = [temp_str, ' : ', temp_str_1];
                            catch error
                                error_code = 1;
                                error_msg = '"multiplex_track" parameter not available';  
                                return;
                            end

                        end % for i_m_track = 1 : number_of_m_track


                        % --- Write "fanout_def": ---

                        fprintf(fid_vex, '    fanout_def = %s;\n', temp_str);

                    end % for i_fanout_def = 1 : number_of_fanout_def
                    
                end % if flag_tracks_available

                fprintf(fid_vex, 'enddef;\n');
                fprintf(fid_vex, '*\n');

            end % for i_tracks = 1 : number_of_bbc

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

