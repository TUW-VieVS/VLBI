% #########################################################################
% #     read_sat_ephem_trf
% #########################################################################
%
% DESCRITPION
% This function reads ephemids files (TRF psotions and velocities) and outputs a MATLAB structure containing the data.
% 
% NOTE:
% - Currently only ephem data for ONE satellite is allowed per input file!!
% - The epoch information is taken from the MJD entry!
%
% The ephem. file needs to have one of the to following format:
%
% => Format 1:
% =====================================================================================================================================
% | Epoch                              | TRF position            | TRF velocity                   | sat. name                         | 
% -------------------------------------------------------------------------------------------------------------------------------------
% | yyyy | mm | dd | hh | mm| ss | MJD | x [m] | y [m] | z [m/s] | vx [m/s] | vy [m/s] | vz [m/s] | sat name string (without blanks!) |
% =====================================================================================================================================
% Example: 2016  7 20 10 35 50 57589.4415513  -3046810.4310   3995517.7040  -4635606.2040   -1721.915 5081.260188000 5531.866180000  APOD
%
% => Format 2:
% =====================================================================================================================================
% |Epoch| TRF position            | TRF velocity                   | sat. name                         | 
% ------------------------------------------------------------------------------------------------------
% | MJD | x [m] | y [m] | z [m/s] | vx [m/s] | vy [m/s] | vz [m/s] | sat name string (without blanks!) |
% ======================================================================================================
% Example: 57695.4236111775  -1544187.951  500632.576  -6642154.807  -4581.229739  5984.915100  1528.428682  APOD
%
% => Format 3:
% ===============================================================================================================================
% | Epoch                              | TRF position            | TRF velocity                   | sat. name                   | 
% -------------------------------------------------------------------------------------------------------------------------------
% | yyyy | mm | dd | hh | mm| ss | x [m] | y [m] | z [m/s] | vx [m/s] | vy [m/s] | vz [m/s] | sat name string (without blanks!) |
% ===============================================================================================================================
% Example: 2016  7 20 10 35 50  -3046810.4310   3995517.7040  -4635606.2040   -1721.915 5081.260188000 5531.866180000  APOD
%
%
% AUTHOR 
%   A. Hellerschmied, 2016-09-21
%
% INPUT
%  - ephem_file_path          - Path of ephem. file
%  - ephem_file_name          - Name of ephem. file
%  
% OUTPUT
%  - orbit_data              - structure containing the orbit data
%
% CHANGES
%  - yyyy-mm-dd, <first + second name>: Description
%  - 2016-11-08, A. Hellerschmied: Format 2 added
%  - 2016-11-08, A. Hellerschmied: added fields "epoch_sec_of_day", "year", "month", "day", "hour", "minu", "sec" to struct "orbit_data"
%  - 2017-06-14, A. Hellerschmied: Format 3 added
%
% EXTERNAL CALLS:
%  - modjuldat.m
%   
% Test call:
% [orbit_data] = read_sat_ephem_trf('../CATALOGS/SAT_EPHEM_TRF/', 'apod_2016-11-03_101000-104000.dat')
%
function [orbit_data] = read_sat_ephem_trf(ephem_file_path, ephem_file_name)


% ##### Init. #####

% ##### Preallocation #####
orbit_data = struct('file_type', 'sat_ephem_trf', 'file_name', ephem_file_name, 'file_path', ephem_file_path, 'mjd_start', [], 'mjd_end', [], 'sat_name_list', [], 'sat', [], 'epoch_mjd', [], 'sec_of_day', [], 'year', [], 'month', [], 'day', [], 'hour', [], 'minu', [], 'sec', []);
orbit_data.sat = struct('name', '', 'x_trf', [], 'y_trf', [], 'z_trf', [], 'vx_trf', [], 'vy_trf', [], 'vz_trf', [], 'flag_v_trf', []);

% ##### Check input arguments #####


% ##### Options #####


% ##### Open ephem file and read the data #####
fid_ephem = fopen([ephem_file_path, ephem_file_name], 'r');
if fid_ephem == -1
    error(['Cannot open satellite ephem. file: ', ephem_file_path, ephem_file_name]);
else
    fprintf(1, 'Loading orbit data from satellite ephem file: %s\n', [ephem_file_path, ephem_file_name]);
end 

% Check how many columns the file has to determine the format:
in_str = fgetl(fid_ephem);
while strcmp(in_str(1), '#')
    in_str = fgetl(fid_ephem);
end
ephem_data = textscan(in_str, '%s');
num_of_col = size(ephem_data{1}, 1);
frewind(fid_ephem);
switch num_of_col
    case 14
        % format = 1;
        % Read all the input data:
        ephem_data = textscan(fid_ephem, '%f %f %f %f %f %f %f %f %f %f %f %f %f %s', 'commentStyle', '#');
        % Prepare input data
        year            = ephem_data{1}; 
        month           = ephem_data{2}; 
        day             = ephem_data{3}; 
        hour            = ephem_data{4}; 
        minu            = ephem_data{5}; 
        sec             = ephem_data{6}; 
        mjd             = ephem_data{7}; 
        x_trf           = ephem_data{8};
        y_trf           = ephem_data{9};
        z_trf           = ephem_data{10};
        vx_trf          = ephem_data{11};
        vy_trf          = ephem_data{12};
        vz_trf          = ephem_data{13};
        sat_name_list   = ephem_data{14};
        sec_of_day      = 3600*hour + 60*minu + sec;
        
    case 8
        % format = 2 (just MJD epochs);
        % Read all the input data:
        ephem_data = textscan(fid_ephem, '%f %f %f %f %f %f %f %s', 'commentStyle', '#');
        % Prepare input data
        mjd             = ephem_data{1}; 
        x_trf           = ephem_data{2};
        y_trf           = ephem_data{3};
        z_trf           = ephem_data{4};
        vx_trf          = ephem_data{5};
        vy_trf          = ephem_data{6};
        vz_trf          = ephem_data{7};
        sat_name_list   = ephem_data{8};
        [year, month, day, hour, minu, sec] = mjd2date(mjd); % round to integer seconds!
        sec_of_day      = 3600*hour + 60*minu + sec;
        
    case 13
        % format = 3;
        % Read all the input data:
        ephem_data = textscan(fid_ephem, '%f %f %f %f %f %f %f %f %f %f %f %f %s', 'commentStyle', '#');
        % Prepare input data
        year            = ephem_data{1}; 
        month           = ephem_data{2}; 
        day             = ephem_data{3}; 
        hour            = ephem_data{4}; 
        minu            = ephem_data{5}; 
        sec             = ephem_data{6}; 
        x_trf           = ephem_data{7};
        y_trf           = ephem_data{8};
        z_trf           = ephem_data{9};
        vx_trf          = ephem_data{10};
        vy_trf          = ephem_data{11};
        vz_trf          = ephem_data{12};
        sat_name_list   = ephem_data{13};
        sec_of_day      = 3600*hour + 60*minu + sec;
        
        % Calc. MJD epochs:
        mjd = modjuldat(year, month, day, hour, minu, sec);
        
    otherwise
        error('Unknown TRF orbit file format.');
end


% Close ephem. file:
fclose(fid_ephem);




% Check if there is data for more than one satellite => Error msg.!
sat_name_str = unique(sat_name_list);
if length(sat_name_str) > 1
    error(['Only ephem data for ONE satellite is allowed in the ephem. file: ', ephem_file_path, ephem_file_name]);
end
 
% Save ephem. data in orbit_data structure:

% Satellite name:
orbit_data.sat_name_list            = sat_name_str;
if ~iscell(orbit_data.sat_name_list)
   orbit_data.sat_name_list = cellstr(orbit_data.sat_name_list); 
end
orbit_data.sat.name                 = sat_name_str{1};

i_sat = 1;

% Positions:
orbit_data.sat(i_sat).x_trf         = x_trf; % [m]
orbit_data.sat(i_sat).y_trf         = y_trf; % [m]
orbit_data.sat(i_sat).z_trf         = z_trf; % [m]

% Velocity
orbit_data.sat(i_sat).vx_trf        = vx_trf; % [m/s]
orbit_data.sat(i_sat).vy_trf        = vy_trf; % [m/s]
orbit_data.sat(i_sat).vz_trf        = vz_trf; % [m/s]

orbit_data.sat(i_sat).flag_v_trf    = true; 

% Epochs [MJD]:
orbit_data.epoch_mjd                = mjd;      % [mjd]
orbit_data.year                     = year;     %
orbit_data.month                    = month;    %
orbit_data.day                      = day;      %
orbit_data.hour                     = hour;  % 
orbit_data.minu                     = minu;     %
orbit_data.sec                      = sec;      % 
orbit_data.sec_of_day               = sec_of_day; % 

% First/last epoch:
orbit_data.mjd_end                  = orbit_data.epoch_mjd(end);
orbit_data.mjd_start                = orbit_data.epoch_mjd(1);

fprintf(1, ' - First epoch:      %s\n', mjd2datestr(orbit_data.mjd_start));
fprintf(1, ' - Last epoch:       %s\n', mjd2datestr(orbit_data.mjd_end));
fprintf(1, ' - Number of epochs: %d\n', length(orbit_data.epoch_mjd));
fprintf(1, ' - Sat. name:        %s\n', orbit_data.sat.name);


fprintf(1, '...finished loading orbit data from TRF ephem. file.\n');
fprintf(1, '\n');

