% #########################################################################
%             TLEDataUpdate.m
% #########################################################################
%
%   This function updates the *.tle files, used for satellite scheduling 
%   by vie_sched from "raw" tle files from the internet (celestrak.org).
%
%   Author: 
%       Andreas Hellerschmied (heller182@gmx.at), 2013-10-13
%   
%   changes       :
%   - 2013-10-16 (Hellerschmied A.): Implementation of "update_log'
%                                      struct.
%   - 2015-06-11 (Hellerschmied A.): The Element Set Number check is ignored, because these numbers are not updated any more by "celestrak.com" for many satellites (constantly set to "999").
%   - 2017-02-22 (Hellerschmied A.): - Problems with Linux EOL symbols solved.
%                                    - "flag_check_set_number" added to activate/deactivate set number checks
%           
%
%   inputs        :
%   - config_filename   : Filename of configuration file.
%   - config_filepath   : Filepath of configurtaion file.
%   - path_tle          : Filepath of local TLE directory
%   - raw_tle_filepath  : Filepath of raw TLE directory
%     
%
%   outputs       :
%   - error_code    : Error Code.
%   - error_msg     : Error message.
%   - update_log    : Struct, containing detailled update information.
%    
%
%   locals        :
% 
%
%   coupling      :
%   
%
%   references    :
%
%-------------------------------------------------------------------------

function [error_code, error_msg, update_log] = TLEDataUpdate(config_filename, config_filepath, path_tle, raw_tle_filepath)

%% --- Options ---
flag_check_set_number = false ;   

%% --- Preallocation and Intialization ---
    
    % var
    error_code = 0;
    error_msg = '';
    state = 1;
    temp_line = '';
    n_lines = 0;
    i_begin = 0;
    length_tle_line = 0;
    total_num_raw_tle_files = 0;
    tle_file_line_counter = 0;
    raw_file_line_counter = 0;
    temp_line_tle = '';
    raw_tle_line_1 = '';
    raw_tle_line_2 = '';
    flag_found_dataset = 0;
    current_sat_name = '';
    tle_el_set_num = '';
    
    
    % struct update_log
    update_log = struct('tle_file', [], 'total_num_updated_sats', [], 'total_num_not_updated_sats', [], 'total_num_no_match_in_raw_tle', []);
    update_log.tle_file = struct('filename', [],'updated_sats', [],...
        'num_updated_sats', [], 'not_updated_sats', [], 'num_not_updated_sats', [], 'no_match_in_raw_tle', [], 'num_no_match_in_raw_tle', []);
    update_log.tle_file.updated_sats = struct('name', []);
    update_log.tle_file.not_updated_sats = struct('name', []);
    update_log.tle_file.no_match_in_raw_tle = struct('name', []);
    
    update_log.total_num_updated_sats = 0;
    update_log.total_num_not_updated_sats = 0;
    update_log.total_num_no_match_in_raw_tle = 0;

    % struct update_tle
    update_tle = struct('tle_filename', [], 'numb_raw_files', [], 'raw_tle_files', []);
    update_tle.raw_tle_files = struct('name', []);
    
    
    %% --- Read configuration file --- 
    
    % Open configuration file
    if ( ~isempty(config_filename) && ~isempty(config_filepath) )
        fid_config = fopen([config_filepath, config_filename], 'r');
        if (fid_config == -1)
            error_code = 1;
            [error_msg] = assign_error_msg(error_code);
            return;
        end;
    else
        error_code = 2;
        [error_msg] = assign_error_msg(error_code);
        return;
    end;
        
    % Read configuration file
     while(~feof(fid_config))
         temp_line = fgetl(fid_config);

         if (temp_line(1) == '#')   % Skip comment lines
             continue;
         end

         length_temp_line = length(temp_line);
         n_lines = n_lines + 1;
         state = 1;
         update_tle(n_lines).numb_raw_files = 0;

         for i = 1 : length_temp_line

             switch(state)  % Finite State Machine to read config file

                 case 1 % Read TLE filename

                     if (temp_line(i) == ' ')
                        update_tle(n_lines).tle_filename = temp_line(1 : (i-1));
                        state = 2;
                     end

                 case 2 % Skip delimiter

                    if ( (temp_line(i) ~= ' ') && (temp_line(i) ~= '-') )
                        i_begin = i;
                        state = 3;
                    end                          

                case 3 %Read filenames of tle raw data files

                     if (temp_line(i) == ' ')
                         update_tle(n_lines).numb_raw_files = update_tle(n_lines).numb_raw_files + 1;
                         update_tle(n_lines).raw_tle_files(update_tle(n_lines).numb_raw_files).name = temp_line(i_begin : (i-1));
                         state = 2;                             
                     end

                     if (i == length_temp_line)
                         update_tle(n_lines).numb_raw_files = update_tle(n_lines).numb_raw_files + 1;
                         update_tle(n_lines).raw_tle_files(update_tle(n_lines).numb_raw_files).name = temp_line(i_begin : (i));

                     end

             end    % switch(state)

         end

     end % while ...

     if ( (exist('fid_config', 'var') ) && (fid_config ~= -1) )
        fclose(fid_config);
    end

    numb_of_tle_files = length(update_tle);

    % Check if the config file contains enough info 
    if (numb_of_tle_files == 0)
        error_code = 3;
        [error_msg] = assign_error_msg(error_code);
        return;
    else
        for i = 1 :  numb_of_tle_files
           %total_num_raw_tle_files = total_num_raw_tle_files + update_tle(i).numb_raw_files;

            if (update_tle(i).numb_raw_files == 0)
                error_code = 3;
                [error_msg] = assign_error_msg(error_code);
                return;
            end
        end
    end


%% --- Update TLE datasets from raw TLE data files (*.txt) ---

    for i_tle_files = 1 : numb_of_tle_files

        % Open TLE file
        fid_tle = fopen([path_tle, update_tle(i_tle_files).tle_filename], 'r+'); % read and write permission

        if (fid_tle == -1)
            error_code = 4;
            [error_msg] = assign_error_msg(error_code);
            return;
        else

            tle_file_line_counter = 0;

            % Set initital log data
            update_log.tle_file(i_tle_files).filename = update_tle(i_tle_files).tle_filename;
            update_log.tle_file(i_tle_files).num_updated_sats = 0;
            update_log.tle_file(i_tle_files).num_not_updated_sats = 0;
            update_log.tle_file(i_tle_files).num_no_match_in_raw_tle = 0;

            % Read TLE file
            while(~feof(fid_tle))
                temp_line_tle = fgets(fid_tle);

                length_tle_line = length(temp_line_tle);
                numb_of_eol_symbols = length_tle_line - 69; % 69 = number of characters per TLE line
                

                tle_file_line_counter = tle_file_line_counter + 1;

                if  ( mod((tle_file_line_counter - 2), 3) ~= 0)   % Skip lines 1,3,4,6,7,9,.. 
                    if (temp_line_tle(1) ~= '2')
                       current_sat_name = temp_line_tle;
                    end
                    continue;
                elseif (temp_line_tle(1) == '1')    % Line 2 of each TLE dataset in file

                    tle_satellite_number = temp_line_tle(3:7);
                    tle_el_set_num = temp_line_tle(65:68);  % Get the Element Set Number (running count for all TLE: increased until it exceeds "999")
                    flag_found_dataset = 0;

                    for i_raw_data = 1 : update_tle(i_tle_files).numb_raw_files

                        % Open raw TLE file
                        fid_raw = fopen([raw_tle_filepath, update_tle(i_tle_files).raw_tle_files(i_raw_data).name], 'r');

                        if (fid_raw == -1)
                            error_code = 5;
                            [error_msg] = assign_error_msg(error_code);
                            return;
                        else

                        raw_file_line_counter = 0;

                            % Read raw TLE file
                            while(~feof(fid_raw))
                                temp_line_tle = fgetl(fid_raw);

                                raw_file_line_counter = raw_file_line_counter + 1;

                                if  ( mod((raw_file_line_counter - 2), 3) ~= 0)   % Skip lines 1,3,4,6,7,9,.. 
                                    continue;
                                elseif (temp_line_tle(1) == '1')    % Line 2 of each TLE dataset in file

                                    % Get NORAD satellite number of dataset in raw
                                    % TLE file:
                                    raw_satellite_number = temp_line_tle(3:7);

                                    % Check for matching NORAD satelite
                                    % numbers in both files:
                                    if (raw_satellite_number == tle_satellite_number)
                                        raw_tle_line_1 = temp_line_tle;
                                        raw_tle_line_2 = fgetl(fid_raw);
                                        flag_found_dataset = 1;
                                        break; % while loop
                                    end

                                else
                                    error_code = 7;
                                    [error_msg] = assign_error_msg(error_code);
                                    return;
                                end
                            end % while(~feof(fid_raw))

                            if ( (exist('fid_raw', 'var') ) && (fid_raw ~= -1) )
                                fclose(fid_raw);
                            end
                        end

                        if flag_found_dataset % Found matching datasets => leave for loop!
                            break; % for loop
                        end;

                    end % for i_raw_data = 1 : update_tle(i_tle_files).numb_raw_files

                else
                    error_code = 6;
                    [error_msg] = assign_error_msg(error_code);
                    return;
                end

                if flag_found_dataset % Found matching datasets? => yes => Copy content to *.tle file!

                    % Check TLE element set number => If nuber of
                    % datasets is not equal, the dataset changed.
                    if (str2double(tle_el_set_num) == str2double(raw_tle_line_1(65:68))) && flag_check_set_number

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% This check is not done any more, because the Element Set Number is consatntly set to "999" for many satellites in the data sets from "celestrak.com" (A. Hellerschmied, 2015-06-11)
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        update_log.tle_file(i_tle_files).num_not_updated_sats                                                           = update_log.tle_file(i_tle_files).num_not_updated_sats + 1;
                        update_log.tle_file(i_tle_files).not_updated_sats(update_log.tle_file(i_tle_files).num_not_updated_sats).name   = current_sat_name;
                        update_log.total_num_not_updated_sats                                                                           = update_log.total_num_not_updated_sats + 1;

                    else % no changes

                        status = fseek(fid_tle, -1*(length_tle_line), 0); 
                        fprintf(fid_tle, '%s' , raw_tle_line_1);
                        status = fseek(fid_tle, numb_of_eol_symbols, 0);
                        fprintf(fid_tle, '%s' , raw_tle_line_2);
                        status=fseek(fid_tle, numb_of_eol_symbols, 0);
                        tle_file_line_counter = tle_file_line_counter + 1;

                        update_log.total_num_updated_sats = update_log.total_num_updated_sats + 1;
                        update_log.tle_file(i_tle_files).num_updated_sats = update_log.tle_file(i_tle_files).num_updated_sats + 1;
                        update_log.tle_file(i_tle_files).updated_sats(update_log.tle_file(i_tle_files).num_updated_sats).name = current_sat_name;

                    end

                else % No matching satellite number found in raw tle files

                    update_log.tle_file(i_tle_files).num_no_match_in_raw_tle = update_log.tle_file(i_tle_files).num_no_match_in_raw_tle + 1;
                    update_log.tle_file(i_tle_files).no_match_in_raw_tle(update_log.tle_file(i_tle_files).num_no_match_in_raw_tle).name = current_sat_name;
                    update_log.total_num_no_match_in_raw_tle = update_log.total_num_no_match_in_raw_tle + 1;

                end;

                n_lines = n_lines + 1;

            end % while(~feof(fid_tle))

            if ( (exist('fid_tle', 'var') ) && (fid_tle ~= -1) )
                fclose(fid_tle);
            end
        end
    end % for i_tle_files = 1 : numb_of_tle_files    


    return;

end

%% --- Local Sub-Rountines ---

function [error_msg] = assign_error_msg(error_code)
    
    error_msg = '';
    
    if (error_code ~= 0) % In case of ERROR => Assignme Error Message!
        
        if (error_code == 1)
            error_msg = 'Failed to open config file.';

        elseif (error_code == 2)
            error_msg = 'Filenameof config file is missing.';

        elseif (error_code == 3)
            error_msg = 'Invalid config. file content.';

        elseif (error_code == 4)
            error_msg = 'Failed to open TLE file.';
            
        elseif (error_code == 5)
            error_msg = 'Failed to open raw TLE file.';
            
        elseif (error_code == 6)
            error_msg = 'Invalid TLE file content.';
            
        elseif (error_code == 7)
            error_msg = 'Invalid raw TLE file content.';
            
        else
            error_msg = 'Unknown Error.';
        end
    
    end
end

