% #########################################################################
%
%                              TLEDataDownload.m
%
% #########################################################################
%
%   This function updates the *.tle files, used for satellite scheduling 
%   by vie_sched from "raw" tle files from the internet (celestrak.org).
%
%   Author: 
%       Andreas Hellerschmied, 2013-10-16
%   
%   Changes       :
%   - 2017-08-31 - A. Hellerschmied: Minor bug fix.
%           
%
%   Inputs        :
%   - config_filename   : Filename of configuration file.
%   - config_filepath   : Filepath of configurtaion file.
%     
%
%   outputs       :
%   - error_code    : Error Code.
%   - error_msg     : Error message.
%   - download_log  : Struct, containing detailled update information.
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

function [error_code, error_msg, download_log] = TLEDataDownload(config_filename, config_filepath, path_raw_tle)

    %% --- Preallocation and Intialization ---
    
    % var
    state = 1;
    error_code = 0;
    error_msg = '';
    temp_line = '';
    n_lines = 0;
    i_begin = 0;  

    % struct update_tle
    download_tle = struct('filename_raw_tle', [], 'URL', []);
    
    download_log = struct('filepath_and_name', [], 'flag_file_written', []);
    
    
    %% --- Set file-paths ---
    % TLE folder
    % path_raw_tle = '../CRF/TLE/RAW_TLE/';
    
    
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

         for i = 1 : length_temp_line
             
             %disp(temp_line(i)); % TEST

             switch(state)  % Finite State Machine to parse and read config file

                 case 1 % Read URL

                     if (temp_line(i) == ' ')
                        download_tle(n_lines).URL = temp_line(1 : (i-1));
                        state = 2;
                     end

                 case 2 % Skip delimiter

                    if ( (temp_line(i) ~= ' ') && (temp_line(i) ~= '-') )
                        i_begin = i;
                        state = 3;
                    end                          

                case 3 %Read  

                     if (temp_line(i) == ' ') 
                         download_tle(n_lines).filename_raw_tle = temp_line(i_begin : (i-1));
                         state = 4;                             
                     elseif (i == length_temp_line) % End of Line
                         download_tle(n_lines).filename_raw_tle = temp_line(i_begin : (i));
                         state = 4;
                     end

                 case 4 % Do nothing and wait

             end    % switch(state)
             
         end % for i = 1 : length_temp_line
         
         if (state ~= 4) % ERROR: config file incompelete
             error_code = 3;
             [error_msg] = assign_error_msg(error_code);
             return;
         end

    end % while ...

    if ( (exist('fid_config', 'var') ) && (fid_config ~= -1) ) % close config file
        fclose(fid_config);
    end

    numb_of_tle_files = length(download_tle);



%% --- Download raw TLE data files (*.txt) ---

    for i_tle_files = 1 : numb_of_tle_files
         [download_log(i_tle_files).filepath_and_name, download_log(i_tle_files).flag_file_written] = urlwrite(download_tle(i_tle_files).URL, [path_raw_tle, download_tle(i_tle_files).filename_raw_tle]);
         
         if (download_log(i_tle_files).flag_file_written == 1) % urlwrite() successfull
            fid_raw_tle = fopen(download_log(i_tle_files).filepath_and_name, 'r');
            if (fid_raw_tle == -1)
                error_code = 5;
                [error_msg] = assign_error_msg(error_code);
                return;
            else
                temp_line = fgetl(fid_raw_tle);
                if (ischar(temp_line)) % Isn't the downloaded file empty?
                    temp_line = fgetl(fid_raw_tle);
                    if (temp_line(1) ~= '1') % Does the second TLE start with an "1"? 
                        error_code = 7;
                        [error_msg] = assign_error_msg(error_code);
                        return;
                    end
                else
                    error_code = 6;
                    [error_msg] = assign_error_msg(error_code);
                    return;
                end
                fclose(fid_raw_tle);
            end
             
         else
             error_code = 4;
             [error_msg] = assign_error_msg(error_code);
             return;
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
            error_msg = 'Invalid configuration file content.';

        elseif (error_code == 4)
            error_msg = 'Failed to download raw TLE file with "urlwrite()".';
            
        elseif (error_code == 5)
            error_msg = 'Failed to open raw TLE file.';
            
        elseif (error_code == 6)
            error_msg = 'Downloaded TLE file is empty.';
            
        elseif (error_code == 7)
            error_msg = 'Invalid raw TLE file content.';
            
        else
            error_msg = 'Unknown Error.';
        end
    
    end
end

