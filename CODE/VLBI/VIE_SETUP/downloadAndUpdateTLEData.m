% #########################################################################
% #     downloadAndUpdateTLEData
% #########################################################################
%
% DESCRITPION
% Function to download raw TLE data from defined web-resources and to update
% the local TLE files, based on the raw TLE data.
%
% AUTHOR 
%   Andreas Hellerschmied (2015.01.22)
%
% INPUT
%   handles      structure from the GUI (also containing data, eg residuals)
%   hObject      ...unused so far.
%
% OUTPUT
%
% CHANGES 
%  - 2016-12-22, A. Hellerschmied: Changed TLE directory from "/CATALOGS/TLE/" to "/ORBIT/TLE/" in the VieVS root dir.
% 

function downloadAndUpdateTLEData()

% #### Set paths and filenames: ####
% configuration files:
% - to define the update procedure (from raw TLE data to local TLE files):
filename_tle_update_config = 'tle_update_config.txt';
filepath_tle_update_config = '../ORBIT/TLE/CONFIG/';
% - to define the web resources for raw TLE files:
filename_tle_download_config = 'tle_download_config.txt';
filepath_tle_download_config = '../ORBIT/TLE/CONFIG/';
% - raw TLE folder:
path_raw_tle = '../ORBIT/TLE/RAW_TLE/';
% - local TLE folder:
path_local_tle = '../ORBIT/TLE/';

% #### Download raw TLE data ####
fprintf(1, '#### Download raw TLE data (config. file: %s) ####\n', filename_tle_download_config);
fprintf(1, '\n');

[error_code, error_msg, download_log] = TLEDataDownload(filename_tle_download_config, filepath_tle_download_config, path_raw_tle);
if (error_code ~= 0) % Error while downloading data
    fprintf(1, 'ERROR (TLE_data_download.m): %s\n', error_msg);
    fprintf(1, '\n');
else % No errors
    
    fprintf(1, ' => ...Done!\n');
    fprintf(1, '\n');
    fprintf(1, ' => Download log:\n');
    for i = 1 : length(download_log)
        fprintf(1, '   => %s\n', download_log(i).filepath_and_name);
    end

    fprintf(1, '\n');

    fprintf(1, '#### Update TLE dataset from raw tle data (config file: %s) ####\n', filename_tle_update_config);
    fprintf(1, '\n');
    
% #### Update TLE Datasets according to information in configurationfile: ####
    [error_code, error_msg, update_log] = TLEDataUpdate(filename_tle_update_config, filepath_tle_update_config, path_local_tle, path_raw_tle);

    % Command Window softcopy
    if (error_code == 0) % no errors

        fprintf(1, ' => ...Done!\n');
        fprintf(1, '\n');
        fprintf(1, ' - Configuration file:                      %s\n', [filepath_tle_update_config, filename_tle_update_config]);
        fprintf(1, ' - Number of TLE datasets:                  %d\n', length(update_log.tle_file));
        fprintf(1, ' - Total number of updated TLE datasets:     %d\n', update_log.total_num_updated_sats);
        fprintf(1, ' - Total number of non updated TLE datasets: %d\n', update_log.total_num_not_updated_sats);
        fprintf(1, ' - Total number of non matching datasets:   %d\n', update_log.total_num_no_match_in_raw_tle);
        fprintf(1, '\n');
        fprintf(1, ' => Update log:\n');
        fprintf(1, '\n');


        for i = 1 : length(update_log.tle_file)
            fprintf(1, '   => TLE file:    %s\n', update_log.tle_file(i).filename);
            fprintf(1, '       - Updated datasets (%d):\n', update_log.tle_file(i).num_updated_sats);
            for j = 1 : update_log.tle_file(i).num_updated_sats
                fprintf(1, '           - %s\n', update_log.tle_file(i).updated_sats(j).name);
            end
            fprintf(1, '       - Non updated datasets (%d):\n', update_log.tle_file(i).num_not_updated_sats);
            for j = 1 : update_log.tle_file(i).num_not_updated_sats
                fprintf(1, '           - %s\n', update_log.tle_file(i).not_updated_sats(j).name);
            end
            fprintf(1, '       - Non matching datasets (%d):\n', update_log.tle_file(i).num_no_match_in_raw_tle);
            for j = 1 : update_log.tle_file(i).num_no_match_in_raw_tle
                fprintf(1, '           - %s\n', update_log.tle_file(i).no_match_in_raw_tle(j).name);
            end  
        end

        fprintf(1, '\n');


    else % error case
        fprintf(1, 'ERROR (TLE_data_update.m): %s\n', error_msg);
        fprintf(1, '\n');
    end
end