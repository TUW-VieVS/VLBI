% #########################################################################
% #     vievs2vso
% #########################################################################
%
% DESCRITPION
% This function writes VSO files based on the VieVS data structures written by VIE_INIT.
%
%
% Possible VSO formats (currently only format 6 is supported!):
%
% 1.) ##### Baseline delay only: #####
% | Delay reference epoch            | Stat. 1 name  | Stat. 2 name  | source name | obs. type | baseline delay [nsec]/[sec] |
% | YYYY MM DD hh mm ss.ssssssssssss | <max 8 char.> | <max 8 char.> | <string>    |   s/q     |      float                  |
%
% 2.) ##### Baseline delay + formal error: #####
% | Delay reference epoch            | Stat. 1 name  | Stat. 2 name  | source name | obs. type | baseline delay [nsec]/[sec] | delay formal error [nsec]/[sec] |
% | YYYY MM DD hh mm ss.ssssssssssss | <max 8 char.> | <max 8 char.> | <string>    |   s/q     |      float                  |      float                      |
%
% 3.) ##### Baseline delay + formal error  and Ion. delay + formal error: #####
% | Delay reference epoch            | Stat. 1 name  | Stat. 2 name  | source name | obs. type | baseline delay [nsec]/[sec] | delay formal error [nsec]/[sec] | Ion delay [nsec] | Ion formal error [nsec] |
% | YYYY MM DD hh mm ss.ssssssssssss | <max 8 char.> | <max 8 char.> | <string>    |   s/q     |      float                  |      float                      |      float       |      float              |
%
% 4.) ##### No observationas (used e.g. to calculate a priori delays in vie_mod for the correlator input model) #####
% - If the input vso file has only 10 colums, this formati is selected automatically
% - If the second vso header line (starting with "#") contains the phrase "noDelay", all entries after colums 10 are skipped (= vso format 4 is selected)
% | Delay reference epoch            | Stat. 1 name  | Stat. 2 name  | source name | obs. type |
% | YYYY MM DD hh mm ss.ssssssssssss | <max 8 char.> | <max 8 char.> | <string>    |   s/q     |
%
% 5.) ##### Baseline delay + formal error, Ion. delay + formal error, cable corr: #####
% | Delay reference epoch            | Stat. 1 name  | Stat. 2 name  | source name | obs. type | baseline delay [nsec]/[sec] | delay formal error [nsec]/[sec] | Ion delay [nsec] | Ion formal error [nsec] | cable corr. stat. 1 [ns] | cable corr. stat. 2 [ns] |
% | YYYY MM DD hh mm ss.ssssssssssss | <max 8 char.> | <max 8 char.> | <string>    |   s/q     |      float                  |      float                      |      float       |      float              |  float                   |  float                   |
%
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
%   A. Hellerschmied, 2017-07-26
%
% INPUT
%  - session_name_str           - Name of the session which should be loaded (string)
%  - level0_path_str            - Path to LEVEL0 sub-dir in VieVS (string)
%  - vso_type                   - VSO output type (default = type 6) (int) (optional)
%
% OUTPUT
%  - 
%
% COUPLING
%  - 
%
% CHANGES
%  - 
%
function vievs2vso(session_name_str, level0_path_str, varargin)


% ##### Options #####
flag_add_clk_offset = false;  % Implementation not finished yet...
if flag_add_clk_offset
    add_clk_offset_antenna_name_str = {'HOBART12', 'YARRA12M', 'KATH12M'}; % Edit antenna list here
    add_clk_offset_sec_q = [0, -3.893e-6, -1.9924e-6]; % Clk offset for quasar scans [sec]
    add_clk_offset_sec_s = [0, -3.893e-6, -1.9924e-6]; % Clk offset for satellite scans [sec]
    add_clk_rate_sec = [4.3876e-13, 0, -3.2907e-13]; % [sec/sec]
    add_clk_ref_epoch_mjd = [modjuldat(2016, 11, 27, 5, 5, 0), 0, modjuldat(2016, 11, 27, 5, 5, 0)]; % Ref. epochs for clock rate [MJD]
end

% Output dir fro VSO files (relatice to VieVS root dir.):
vso_out_dir_str = 'OUT/VSO';

% VSO type:
switch(nargin)
    case 2
        vso_type = 6;
    case 3
        vso_type = varargin{1};
        if type ~= 6
            error('Sorry, currently only VSO format 6 is supported!');
        end
    otherwise
        error('Invalid number of input pamareters.');
end


% ##### Load VieVS structures #####

% Get vievs root dir:
vievs_root_dir = pwd;
vievs_root_dir = vievs_root_dir(1:strfind(vievs_root_dir, 'WORK') - 1);

% Clear old and load new structures:
clear antenna scan parameter sources
load([vievs_root_dir, 'DATA/LEVEL0/', level0_path_str, '/',session_name_str, '_antenna.mat']);
load([vievs_root_dir, 'DATA/LEVEL0/', level0_path_str, '/', session_name_str, '_scan.mat']);
load([vievs_root_dir, 'DATA/LEVEL0/', level0_path_str, '/', session_name_str, '_parameter.mat']);
load([vievs_root_dir, 'DATA/LEVEL0/', level0_path_str, '/', session_name_str, '_sources.mat']);


% Apply clock offsets and rates on observations:
if flag_add_clk_offset
    for i_stat = 1 : length(add_clk_offset_antenna_name_str)
        if strcmp(deblank(antenna(stat_1_id).name), add_clk_offset_antenna_name_str{i_stat}) % station 1 => subtract clk offset
            switch(scan(isc).obs_type)
                case 'q'
                    scan(isc).obs(iobs).obs = scan(isc).obs(iobs).obs - add_clk_offset_sec_q(i_stat);
                case 's'
                    scan(isc).obs(iobs).obs = scan(isc).obs(iobs).obs - add_clk_offset_sec_s(i_stat);
            end
            if add_clk_rate_sec(i_stat) ~= 0 % rate
                scan(isc).obs(iobs).obs = scan(isc).obs(iobs).obs - (mjd - add_clk_ref_epoch_mjd(i_stat))*86400 * add_clk_rate_sec(i_stat);
            end
        elseif strcmp(deblank(antenna(stat_2_id).name), add_clk_offset_antenna_name_str{i_stat})
            switch(scan(isc).obs_type)
                case 'q'
                    scan(isc).obs(iobs).obs = scan(isc).obs(iobs).obs + add_clk_offset_sec_q(i_stat);
                case 's'
                    scan(isc).obs(iobs).obs = scan(isc).obs(iobs).obs + add_clk_offset_sec_s(i_stat);
            end
            if add_clk_rate_sec(i_stat) ~= 0 % rate
                scan(isc).obs(iobs).obs = scan(isc).obs(iobs).obs + (mjd - add_clk_ref_epoch_mjd(i_stat))*86400 * add_clk_rate_sec(i_stat);
            end
        end
    end
end


% ##### Write VSO file #####

% remove ".vso" from session name:
if ~isempty(strfind(session_name_str, '.vso'))
    session_name_str = session_name_str(1 : strfind(session_name_str, '.vso') - 1);
end

% Add info about applied Clk reductions to filename:
if flag_add_clk_offset
    session_name_str_2 = [session_name_str_2, '_redClk'];
else
    session_name_str_2 = session_name_str;
end


% Open VSO file:
vso_absolute_out_dir_str = [vievs_root_dir, vso_out_dir_str, '/', level0_path_str];
if ~exist(vso_absolute_out_dir_str, 'dir')
    mkdir(vso_absolute_out_dir_str);
end

fid_vso = fopen([vso_absolute_out_dir_str, '/', session_name_str_2, '.vso'] , 'w');
if fid_vso < 0
    error('Cannot open output VSO file!');
end

station_list_str = '';
for i_stat = 1 : length(antenna)
    station_list_str = [station_list_str, deblank(antenna(i_stat).name), ' '];
end

% Header
fprintf(fid_vso, '################################################################################################\n');
fprintf(fid_vso, '# VSO file for experiment:  %s\n', session_name_str);
fprintf(fid_vso, '# VSO type:                 %d\n', vso_type);
fprintf(fid_vso, '# Created in VieVS:         %s\n', datestr(clock));
fprintf(fid_vso, '# Start:                    %s\n', sprintf('%04d-%02d-%02d %02d:%02d:%05.2f', scan(1).tim(1), scan(1).tim(2), scan(1).tim(3), scan(1).tim(4), scan(1).tim(5), scan(1).tim(6)));
fprintf(fid_vso, '# End:                      %s\n', sprintf('%04d-%02d-%02d %02d:%02d:%05.2f', scan(end).tim(1), scan(end).tim(2), scan(end).tim(3), scan(end).tim(4), scan(end).tim(5), scan(end).tim(6)));
fprintf(fid_vso, '# Stations:                 %s\n', station_list_str);
fprintf(fid_vso, '# Number of scans:          %d\n', length(scan));
if flag_add_clk_offset
fprintf(fid_vso, '# Observations reduced by station clock offsets:\n');
fprintf(fid_vso, '# - %8s: %+f nsec\n', stat_names_cell{1}, add_clk_offset_sec(strcmp(deblank(stat_names_cell{1}),  add_clk_offset_antenna_name_str))*1e9);
fprintf(fid_vso, '# - %8s: %+f nsec\n', stat_names_cell{2}, add_clk_offset_sec(strcmp(deblank(stat_names_cell{2}),  add_clk_offset_antenna_name_str))*1e9);
end
fprintf(fid_vso, '################################################################################################\n');

switch(vso_type)
%     case 1
%         fprintf(fid_vso, '# Reference epoch                | Baseline        | Source name               | Obs. type |   Delay [nsec] \n');
%         % Write lines (one line per obs.):
%         for i_obs = 1 : length(year)
%             fprintf(fid_vso, '%4d %02d %02d %02d %02d %015.12f  %8s %8s  %25s  %s           %f   %f  \n', year(i_obs), month(i_obs), day(i_obs), hour(i_obs), minu(i_obs), sec(i_obs), stat_1{i_obs}, stat_2{i_obs}, src_name{i_obs}, obs_type{i_obs}, delay_nsec(i_obs), delay_formal_error_nsec(i_obs));
%         end

    case 6 % VSO type 6
        fprintf(fid_vso, '# Reference epoch                | Baseline        | Source name               | Obs. type |   Delay [nsec] | del. sig [nsec] | Ion-del [nsec] | Ion-sig [nsec] | cable corr. stat. 1 [ns] | cable corr. stat. 2 [ns] | T stat. 1 [C] | p stat. 1 [mbar] | e stat.1 [] |T stat. 2 [C] | p stat. 2 [mbar] | e stat. 2 [] | \n');
        % Write lines (one line per obs.):
        for i_scan = 1 : length(scan)
            for i_obs = 1 : length(scan(i_scan).obs)
                stat1_id = scan(i_scan).obs(i_obs).i1;
                stat2_id = scan(i_scan).obs(i_obs).i2;
                sig_delay_ns    = sqrt((scan(i_scan).obs(i_obs).sig*1e9)^2 - scan(i_scan).obs(i_obs).sgdion^2);
                delay_ns        = scan(i_scan).obs(i_obs).obs*1e9 - scan(i_scan).stat(stat2_id).cab + scan(i_scan).stat(stat1_id).cab + scan(i_scan).obs(i_obs).delion;
                hum_1           = scan(i_scan).stat(stat1_id).e * 100/6.1078 * exp(-((17.1 * scan(i_scan).stat(stat1_id).temp)/(235 + scan(i_scan).stat(stat1_id).temp)));
                hum_2           = scan(i_scan).stat(stat2_id).e * 100/6.1078 * exp(-((17.1 * scan(i_scan).stat(stat2_id).temp)/(235 + scan(i_scan).stat(stat2_id).temp)));
                fprintf(fid_vso, '%4d %02d %02d %02d %02d %015.12f  %8s %8s  %25s  %7s %18.6f %15.6f %16.6f %16.6f %26.6f %26.6f %15.2f %18.2f %13.2f %15.2f %18.2f %13.2f\n', scan(i_scan).tim(1), scan(i_scan).tim(2), scan(i_scan).tim(3), scan(i_scan).tim(4), scan(i_scan).tim(5), scan(i_scan).tim(6), antenna(stat1_id).name, antenna(stat2_id).name, scan(i_scan).src_name, scan(i_scan).obs_type, delay_ns, sig_delay_ns, scan(i_scan).obs(i_obs).delion, scan(i_scan).obs(i_obs).sgdion, scan(i_scan).stat(stat1_id).cab, scan(i_scan).stat(stat2_id).cab, scan(i_scan).stat(stat1_id).temp, scan(i_scan).stat(stat1_id).pres, hum_1, scan(i_scan).stat(stat2_id).temp, scan(i_scan).stat(stat2_id).pres, hum_2);
            end
        end
        
    otherwise
        error('Selected VSO type (%d) not supported yet.', vso_type);
end

fclose(fid_vso);

fprintf('Writing VSO file %s finished!\n', [vso_absolute_out_dir_str, '/', session_name_str_2, '.vso']);




