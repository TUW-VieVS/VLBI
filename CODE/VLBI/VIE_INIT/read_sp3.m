% #########################################################################
% #     read_sp3
% #########################################################################
%
% DESCRITPION
% This function reads SP3 orbit files and outputs a MATLAB structure containing the data.
% Options:
%  - Read data of a certain time span only
%  - Define a sub-set of satellites to read in the orbit data for them
%
% This function only read position record (lines starting with "P")!
% There has to be one entry in each position record for each satellite in the input file!
%
% AUTHOR 
%   A. Hellerschmied, 2016-07-28
%
% INPUT
%  - sp3_file_path          - Path of SP3 file
%  - sp3_file_name          - Name of SP3 file
%  - t_mjd_start            - Start time to read in data (optional)
%  - t_mjd_end              - End time to read in data (optional)
%  - sat_list_str           - List of satellites (cell array of strings)(optional)
%  
% OUTPUT
%  - orbit_data              - structure containing the orbit data
%
% CHANGES
%  - yyyy-mm-dd, <first + second name>: Description
%  - 2016-09-01, A. Hellerschmied: x,y,z coordinates converted to [m] from [km] (as provided in SP3 files)
%  - 2016-09-21, A. Hellerschmied: Added field vx_trf, vy_trf, vz_trf and flag_v_trf to orbit_data structure.
%  - 2016-12-06, A. Hellerschmied: - Added fields 'sec_of_day', 'year', 'month', 'day', 'hour', 'minu', 'sec' to orbit_data struct
%                                  - use modjuldat.m instead of juliandate.m
%
% EXTERNAL CALLS:
%  - tai_utc.m
%  - modjuldat.m
%  - mjd2date.m
%   
% Test call:
% [orbit_data] = read_sp3('../SP3/', '131a_igr18962.sp3',5.751804166666651e+04,5.751809375000000e+04)
%
function [orbit_data] = read_sp3(sp3_file_path, sp3_file_name, varargin)


% ##### Init. #####
flag_read_all_sat       = true;
flag_read_all_epochs    = true;

% ##### Preallocation #####
orbit_data = struct('file_type', 'sp3', 'file_name', sp3_file_name, 'file_path', sp3_file_path, 'mjd_start', [], 'mjd_end', [], 'sat_name_list', [], 'sat', [], 'epoch_mjd', [], 'sec_of_day', [], 'year', [], 'month', [], 'day', [], 'hour', [], 'minu', [], 'sec', []);
orbit_data.sat = struct('name', '', 'x_trf', [], 'y_trf', [], 'z_trf', [], 'vx_trf', [], 'vy_trf', [], 'vz_trf', [], 'flag_v_trf', []);

% ##### Check input arguments #####
switch(nargin)
    case 2
        
    case 3 % List of satellites
        flag_read_all_sat = false;
        orbit_data.sat_name_list = varargin{1};
        if ~iscell(orbit_data.sat_name_list)
           orbit_data.sat_name_list = cellstr(orbit_data.sat_name_list); 
        end
            
    case 4 % Start and end time
        flag_read_all_epochs    = false;
        orbit_data.mjd_start    = varargin{1};
        orbit_data.mjd_end      = varargin{2};
        
    case 5 % List of satellites & Start and end time
        flag_read_all_epochs    = false;
        orbit_data.mjd_start    = varargin{1};
        orbit_data.mjd_end      = varargin{2};
        flag_read_all_sat = false;
        orbit_data.sat_name_list = varargin{3};
        if ~iscell(orbit_data.sat_name_list)
           orbit_data.sat_name_list = cellstr(orbit_data.sat_name_list); 
        end
        
    otherwise % Error
        error('invalid number of input arguments!');

end

% ##### Options #####


% ##### Open SP3 file and read the data #####
fid_sp3 = fopen([sp3_file_path, sp3_file_name], 'r');
if fid_sp3 == -1
    error(['Cannot open SP3 file: ', sp3_file_path, sp3_file_name]);
else
    fprintf(1, 'Loading orbit data from SP3 file: %s\n', [sp3_file_path, sp3_file_name]);
end 

% Read all the data:
sp3_data = textscan(fid_sp3, '%s','Delimiter','');

% Close SP3 file:
fclose(fid_sp3);

sp3_data = sp3_data{1}; 

% ##### Epochs: #####
epoch_rec_ind = strncmp(sp3_data, '*', 1);
% number_of_epochs = sum(epoch_rec_ind);

epoch_rec_cell = sp3_data(epoch_rec_ind);
epoch_rec_char = char(epoch_rec_cell);

year    = str2double(cellstr(epoch_rec_char(:, 4:7)));
mon     = str2double(cellstr(epoch_rec_char(:, 9:10)));
day     = str2double(cellstr(epoch_rec_char(:, 12:13)));
hour    = str2double(cellstr(epoch_rec_char(:, 15:16)));
min     = str2double(cellstr(epoch_rec_char(:, 18:19)));
sec     = str2double(cellstr(epoch_rec_char(:, 21:31)));

% mjd_tmp = juliandate([year, mon, day, hour, min, sec]) -  2400000.5; 
mjd_tmp    = modjuldat(year, mon, day, hour, min, sec);

% Convert GPS time to UTC
 % - Get time diff. between UTC and GPS time:
 % tgps = UTC + leap_sec_tai_utc - 19 sec
 % => UTC = tgps + 19 sec - leap_sec_tai_utc
 leap_sec_tai_utc = tai_utc(mjd_tmp);
 leap_sec_gps_utc = 19 - leap_sec_tai_utc;
 sec = sec + leap_sec_gps_utc;
 
 unique_leap_sec_gps_utc = unique(leap_sec_gps_utc);
 fprintf(1, '  => Current difference between GPS time (ref. time in SP3 files) and UTC:\n');
 for i = 1 : length(unique_leap_sec_gps_utc)
     fprintf(1, '    => %s sec\n', num2str(unique_leap_sec_gps_utc(i)));
 end

% Convert epochs to MJD:
orbit_data.epoch_mjd = modjuldat(year, mon, day, hour, min, sec);

% Convert MJD to date/time (taking into account diff. bewtwwen GPS time and UTC!)
[year, month, day, hour, minu, sec] = mjd2date(orbit_data.epoch_mjd); % round to integer seconds!
sec_of_day = 3600*hour + 60*minu + sec;

orbit_data.year         = year;
orbit_data.month        = month;
orbit_data.day          = day;
orbit_data.hour         = hour;
orbit_data.minu         = minu;
orbit_data.sec          = sec;
orbit_data.sec_of_day   = sec_of_day;




% ##### Exclude mjd epochs, if they are out of the wanted time range #####
if ~flag_read_all_epochs
    excl_epochs_ind = (orbit_data.epoch_mjd < orbit_data.mjd_start) | (orbit_data.epoch_mjd > orbit_data.mjd_end);
    
    if length(excl_epochs_ind) == sum(excl_epochs_ind)
        error('All epochs in the SP3 file excluded by the defined time window!');
    end
    % Exclude data:
    orbit_data.year         = orbit_data.year(~excl_epochs_ind);
    orbit_data.month        = orbit_data.month(~excl_epochs_ind);
    orbit_data.day          = orbit_data.day(~excl_epochs_ind);
    orbit_data.hour         = orbit_data.hour(~excl_epochs_ind);
    orbit_data.minu         = orbit_data.minu(~excl_epochs_ind);
    orbit_data.sec          = orbit_data.sec(~excl_epochs_ind);
    orbit_data.sec_of_day   = orbit_data.sec_of_day(~excl_epochs_ind);
    orbit_data.epoch_mjd    = orbit_data.epoch_mjd(~excl_epochs_ind);
    
    
end


% ##### Satellite positions and names #####

% Indices of position Records:
pos_rec_ind = strncmp(sp3_data, 'P', 1);
number_of_pos_rec = sum(pos_rec_ind);

sp3_pos_rec_data = sp3_data(pos_rec_ind);
pos_rec_char = char(sp3_pos_rec_data);

% list of satellite names:
sat_names_char = pos_rec_char(:, 2:4);
if flag_read_all_sat
    sat_name_list = unique(cellstr(sat_names_char));
    orbit_data.sat_name_list = sat_name_list;
else
    sat_name_list = orbit_data.sat_name_list;
    % Exclude data (from sp3_pos_rec_data):
    incl_pos_rec_ind = false(number_of_pos_rec, 1);
    for i_sat = 1 : length(sat_name_list)
        incl_pos_rec_ind = incl_pos_rec_ind | strcmp(cellstr(sat_names_char), sat_name_list{i_sat});  
    end    
    sp3_pos_rec_data    = sp3_pos_rec_data(incl_pos_rec_ind);
    pos_rec_char        = char(sp3_pos_rec_data);
    sat_names_char      = pos_rec_char(:, 2:4);
end

% Number of satellites
number_of_sat = length(sat_name_list);

% If a time window to read data is defined, exclude pos. records out of this time range:
if ~flag_read_all_epochs
    excl_pos_rec_ind = kron(excl_epochs_ind, true(number_of_sat, 1));
    sp3_pos_rec_data = sp3_pos_rec_data(~excl_pos_rec_ind);
    pos_rec_char = char(sp3_pos_rec_data);
    sat_names_char = pos_rec_char(:, 2:4);
end

% Satellite names:
[orbit_data.sat(1:number_of_sat).name] = deal(sat_name_list{:});

% Positions:
x_char = pos_rec_char(:, 5:18);
y_char = pos_rec_char(:, 19:32);
z_char = pos_rec_char(:, 33:46);

x = str2double(cellstr(x_char));
y = str2double(cellstr(y_char));
z = str2double(cellstr(z_char));

% Convert x,y,z coordinates from [km] to [m]
x = x * 1e3;
y = y * 1e3;
z = z * 1e3;

% Loop over all satellites
for i_sat = 1 : number_of_sat
    % Get indices of current sat:
    sat_ind = strcmp(cellstr(sat_names_char), orbit_data.sat(i_sat).name);
    % Get positions:
    orbit_data.sat(i_sat).x_trf         = x(sat_ind); % [m]
    orbit_data.sat(i_sat).y_trf         = y(sat_ind); % [m]
    orbit_data.sat(i_sat).z_trf         = z(sat_ind); % [m]
    orbit_data.sat(i_sat).vx_trf        = []; % [m/s]
    orbit_data.sat(i_sat).vy_trf        = []; % [m/s]
    orbit_data.sat(i_sat).vz_trf        = []; % [m/s]
    orbit_data.sat(i_sat).flag_v_trf    = false; 
end



if flag_read_all_epochs
    orbit_data.mjd_end      = orbit_data.epoch_mjd(end);
    orbit_data.mjd_start    = orbit_data.epoch_mjd(1);
end

fprintf(1, '...finished loading orbit data from SP3 file.\n');
fprintf(1, '\n');

