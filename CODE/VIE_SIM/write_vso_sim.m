% #########################################################################
% #     write_vso_sim
% #########################################################################
%
% This function writes vso files with simulated observations for VIE_SIM.
% 
% VSO Format #6 (see function /VIE_INIT_VXX/read_vso.m):
% 6.) ##### Baseline delay + formal error, Ion. delay + formal error, cable corr, met. data (T, p, e): #####
% | Delay reference epoch            | Stat. 1 name  | Stat. 2 name  | source name | obs. type | baseline delay [nsec]/[sec] | delay formal error [nsec]/[sec] | Ion delay [nsec] | Ion formal error [nsec] | cable corr. stat. 1 [ns] | cable corr. stat. 2 [ns] | T stat. 1 [°C] | p stat. 1 [mbar] | Rel. humidity e stat. 1 |T stat. 2 [°C] | p stat. 2 [mbar] | Rel. humidity e stat. 2 |
% | YYYY MM DD hh mm ss.ssssssssssss | <max 8 char.> | <max 8 char.> | <string>    |   s/q     |      float                  |      float                      |      float       |      float              |  float                   |  float                   | float          | float            | float                   |float          | float            | float                   |    
%  
% Note 1:
% ### Units in the scan structure: ###
% delay                         : [sec]
% delay formal error            : [sec]
% ion. delay                    : [nsec]
% ion. delay formal error       : [nsec]
% cable delays                  : [nsec]
%
%
% AUTHOR 
%   A. Hellerschmied (2016-10-21)
%
% INPUT
%      session_name_str     -    string with session name
%      antenna              -    the VieVS antenna structure array
%      scan                 -    the VieVS scan structure array (with simulated observations in field "scan.obs.obs")
%      sources              -    the VieVS source structure array
%      year                 -    year (number)
%      n_days               -    number of days for which simulations were carried out = number of NGS files to be written
%      filename_index       -    first index for the running number of the NGS file
%      flag_zero_input      -    zero input? (1 = yes, 0 = no) set in vie_sim.m
%      sim_subdir_str       -    subdirectory in SIM
%
% OUTPUT
%
% CHANGES
%
%
function write_vso_sim(session_name_str, antenna, scan, sources, year, n_days, filename_index, flag_zero_input, sim_subdir_str)

% ##### Check if SIM sub-directory in the DATA directory exists and create it if it doesn't #####
if isempty(sim_subdir_str)
    out_dir_str = ['../DATA/SIM/', num2str(year), '/'];
    if ~isdir(out_dir_str)
        mkdir(out_dir_str);
    end
else
    out_dir_str = ['../DATA/SIM/', sim_subdir_str, '/',num2str(year), '/'];
    if ~isdir(out_dir_str)
        mkdir(out_dir_str);
    end
end

% ##### Loop over all days wich were simulated #####
for i_day = 1 : n_days
    
    % ##### Open output vso file #####
    
    % Name of output file:
    filename_index_str  = sprintf('_S%03d', filename_index + i_day -1);
    
    if strcmp(session_name_str(end-3:end), '.vso') % VSO input file
        filename_str        = [session_name_str(1:end-4), filename_index_str, '.vso'];
    else % other input file types
        filename_str        = [session_name_str, filename_index_str, '.vso']; 
    end
    
    fid_vso = fopen([out_dir_str ,filename_str], 'w');
    if fid_vso == -1
       error('Opening file %s failed!', [out_dir_str,'/',filename_str]);
    end
    fprintf('  - Writing file: %s', [out_dir_str ,filename_str])
    
    % ##### Write header to file #####
    fprintf(fid_vso, '# VSO format 6\n');
    fprintf(fid_vso, '# This file contains simulated observations for the session: %s\n', session_name_str);
    fprintf(fid_vso, '# Creation date: %s\n',date);
    fprintf(fid_vso, '# Simulated day: %d of %d\n', i_day, n_days);
    fprintf(fid_vso, '# ------------------------------+--------+---------+---------+----+----------------------------+------------------------------+----------------------+-----+-------+-----+-----+-------+-----+\n');
    fprintf(fid_vso, '# Observation epoch             | Station names    | Source  |type| Delay + formal error [ns]  | Ion. delay + formal error    | Cable corrections    | T   | p     | e   | T   | p     | e   |\n');
    fprintf(fid_vso, '#                               | Stat 1 | Stat 2  | Name    |q/s | Modeled delay              |  [ns]                        | [ns]                 |[°C] |[mbar] | [%%] |[°C] |[mbar] | [%%] |\n');
    fprintf(fid_vso, '#                               |        |         |         |    |                            |                              | Stat. 1   | Stat. 2  | Station 1         | Station 2         |\n');
    fprintf(fid_vso, '# ------------------------------+--------+---------+---------+----+----------------------------+------------------------------+----------------------+-------------------+-------------------+\n');
    
    % ##### Write observations to vso file: #####
    % Loop over all scans:
    for i_scan = 1 : length(scan)
        
        % Loop over all observations in scan
        for i_obs = 1 : length(scan(i_scan).obs)
            
            % If fields are empty => set them to zero:
            if isempty(scan(i_scan).obs(i_obs).delion)
                scan(i_scan).obs(i_obs).delion = 0;
            end
            if isempty(scan(i_scan).obs(i_obs).sgdion)
                scan(i_scan).obs(i_obs).sgdion = 0;
            end
            
             % remove delay-corrections applied in read_ngs
            cable_corr_ns = scan(i_scan).stat(scan(i_scan).obs(i_obs).i2).cab - scan(i_scan).stat(scan(i_scan).obs(i_obs).i1).cab; % cable cal of 2nd station - cable cal of 1st station, [ns]
            % zero input or simulated input?
            if flag_zero_input
                % if zero input -> o = c -> o-c = 0
                sim_delay_ns = 1.0d9 * scan(i_scan).obs(i_obs).com + scan(i_scan).obs(i_obs).delion - cable_corr_ns; % [ns]
            else
                sim_delay_ns = 1.0d9 * (scan(i_scan).obs(i_obs).com + scan(i_scan).obs(i_obs).obs(1,i_day)) + scan(i_scan).obs(i_obs).delion - cable_corr_ns; % [ns]
            end
            % the sigma of the simulated delay observable is set to the value of the simulated thermal noise (ionospheric formal error has to be taken into account, if available)
            sim_formal_delay_error_ns = sqrt((1.0d9 * scan(i_scan).obs(i_obs).sig)^2 - scan(i_scan).obs(i_obs).sgdion^2); % [ns]

            % 6.) ##### Baseline delay + formal error, Ion. delay + formal error, cable corr, met. data (T, p, e): #####
            % | Delay reference epoch            | Stat. 1 name  | Stat. 2 name  | source name | obs. type | baseline delay [nsec]/[sec] | delay formal error [nsec]/[sec] | Ion delay [nsec] | Ion formal error [nsec] | cable corr. stat. 1 [ns] | cable corr. stat. 2 [ns] | T stat. 1 [°C] | p stat. 1 [mbar] | Rel. humidity e stat. 1 |T stat. 2 [°C] | p stat. 2 [mbar] | Rel. humidity e stat. 2 |
            % | YYYY MM DD hh mm ss.ssssssssssss | <max 8 char.> | <max 8 char.> | <string>    |   s/q     |      float                  |      float                      |      float       |      float              |  float                   |  float                   | float          | float            | float                   |float          | float            | float                   |  
            fprintf(fid_vso,'%04d %02d %02d %02d %02d %015.12f %8s %8s %10s %1s %20.6f %10.6f %20.6f %10.6f %10.6f %10.6f %5.2f %7.2f %5.2f %5.2f %7.2f %5.2f\n', ...
                scan(i_scan).tim(1), scan(i_scan).tim(2), scan(i_scan).tim(3), scan(i_scan).tim(4), scan(i_scan).tim(5), scan(i_scan).tim(6),...
                antenna(scan(i_scan).obs(i_obs).i1).name, antenna(scan(i_scan).obs(i_obs).i2).name,...
                scan(i_scan).src_name, ...
                scan(i_scan).obs_type, ...
                sim_delay_ns, sim_formal_delay_error_ns,...                                                                 % BL delay + formal error [s]
                scan(i_scan).obs(i_obs).delion, scan(i_scan).obs(i_obs).sgdion,...                                          % Ion. delay + formal error[ns]
                scan(i_scan).stat(scan(i_scan).obs(i_obs).i1).cab, scan(i_scan).stat(scan(i_scan).obs(i_obs).i2).cab,...   % cable corrections [ns]
                scan(i_scan).stat(scan(i_scan).obs(i_obs).i1).temp, scan(i_scan).stat(scan(i_scan).obs(i_obs).i1).pres, scan(i_scan).stat(scan(i_scan).obs(i_obs).i1).e,...  % Station 1: T, p, e
                scan(i_scan).stat(scan(i_scan).obs(i_obs).i2).temp, scan(i_scan).stat(scan(i_scan).obs(i_obs).i2).pres, scan(i_scan).stat(scan(i_scan).obs(i_obs).i2).e);    % Station 2: T, p, e

        end % for i_obs = 1 : length(scan(i_scan).obs)
    end % for i_scan = 1 : length(scan)
    
    fprintf('  ...finished\n');

    % Close output file:
    fclose(fid_vso);

end


