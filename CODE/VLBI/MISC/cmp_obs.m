% #########################################################################
% #     cmp_obs
% #########################################################################
%
% DESCRITPION
% This function allows to compare the observation data (contained in the VieVS 
% data structures as written by VIE_INIT) of sessions, or of different versions
% of the same session (e.g. from NGS vs. vgosDB file). 
% Corresponding observations are matched by scan reference time and observing 
% baseline.
%
%
% AUTHOR 
%   2018-03-21: A. Hellerschmied
%
% INPUT
%   str_path_1            - Path of dataset 1
%   str_sess_name_1       - Name of dataset 1 (output data structures of VIE_MOD)
%   str_path_2            - Path of dataset 2
%   str_sess_name_2       - Name of dataset 2 (output data structures of VIE_MOD)
%   plot_type (optional)  - Plot differences, or delays of dataset 1/2 directly
%   ant_code (optional)   - Two letter antenna code, to restrict the output to observations of this antenna only
%
% OUTPUT
%
% CHANGES
%   
% EXAMPLE:
% cmp_obs('./NGS', '18JAN22XA_N004', './vgosDB', '18JAN22XA', 'diff', 'Ht')
%  => To plot the differences between the 2 datasets, only considering observations of stations Ht. 
%

function cmp_obs(str_path_1, str_sess_name_1, str_path_2, str_sess_name_2, varargin)

% Options:
flag_save_plot = false;

% Check input arguments
switch(nargin)
    case 4
        flag_all_ant = true;
        plot_type = 'diff';
    case 5
        plot_type = varargin{1};
        flag_all_ant = true;
    case 6
        plot_type = varargin{1};
        ant_code = varargin{2};
        flag_all_ant = false;
    otherwise
        error('Invalid number of input arguments!')
end


% ##### Load VieVS structures #####
% From dataset 1:
load([str_path_1, '/', str_sess_name_1, '_antenna.mat'])
load([str_path_1, '/', str_sess_name_1, '_parameter.mat'])
load([str_path_1, '/', str_sess_name_1, '_scan.mat'])
load([str_path_1, '/', str_sess_name_1, '_sources.mat'])
antenna_1 = antenna;
parameter_1 = parameter;
scan_1 = scan;
sources_1 = sources;
clear antenna parameter scan sources

% From dataset 2:
load([str_path_2, '/', str_sess_name_2, '_antenna.mat'])
load([str_path_2, '/', str_sess_name_2, '_parameter.mat'])
load([str_path_2, '/', str_sess_name_2, '_scan.mat'])
load([str_path_2, '/', str_sess_name_2, '_sources.mat'])
antenna_2 = antenna;
parameter_2 = parameter;
scan_2 = scan;
sources_2 = sources;
clear antenna parameter scan sources

% Sources:
num_of_q_1 = length(sources_1.q);
num_of_q_2 = length(sources_2.q);
num_of_s_1 = length(sources_1.s);
num_of_s_2 = length(sources_2.s);

% Number of scans:
num_of_scans_1 = length(scan_1);
num_of_scans_2 = length(scan_2);

% Number of observations:
num_of_obs_1 = length([scan_1.obs]);
num_of_obs_2 = length([scan_2.obs]);

% Scan reference epochs:
scan_epochs_1   = [scan_1.mjd];
scan_epochs_2   = [scan_2.mjd];
all_scan_epochs = unique([scan_epochs_1, scan_epochs_2]);
all_scan_epochs = sort(all_scan_epochs);
num_of_epochs = length(all_scan_epochs);

% antenna 
ant_names_1             = {antenna_1.name};
ant_names_2             = {antenna_2.name};
ant_code_1             = {antenna_1.code};
ant_code_2             = {antenna_2.code};
[all_ant_names, ic, ib ] = unique([ant_names_1, ant_names_2]);
unique_ant_id_1         = ib(1 : length(ant_names_1));
unique_ant_id_2         = ib(length(ant_names_1)+1 : end);

num_of_ant_1    = length(ant_names_1);
num_of_ant_2    = length(ant_names_2);
num_of_all_ant  = length(all_ant_names);
all_ant_codes   = [ant_code_1, ant_code_2];
all_ant_codes   = all_ant_codes(ic);


% baseline indices:
num_of_bl       = (num_of_all_ant * (num_of_all_ant - 1)/2);
bl_ind_list = nan(num_of_bl, 2);
i_bl = 0;
for i_ant = 1 : num_of_all_ant
    for i = i_ant + 1 : num_of_all_ant 
        i_bl = i_bl + 1;
        bl_ind_list(i_bl,1) = i_ant; % antenna 1 index
        bl_ind_list(i_bl,2) = i;     % antenna 2 index
    end
end

bl_ant_name_list = all_ant_names(bl_ind_list);
bl_ant_code_list = all_ant_codes(bl_ind_list);


% Print statistics to CW: 
fprintf(1, 'Session 1: %s from %s (%s):\n', parameter_1.session_name, str_path_1, parameter_1.data_type);
fprintf(1, '  Number of obs.:       %d\n', num_of_obs_1);
fprintf(1, '  Number of scans.:     %d\n', num_of_scans_1);
fprintf(1, '  Number of antennas.:  %d\n', num_of_ant_1);
fprintf(1, '  Number of quasars.:   %d\n', num_of_q_1);
fprintf(1, '\n');
fprintf(1, 'Session 2: %s from %s (%s):\n', parameter_2.session_name, str_path_2, parameter_2.data_type);
fprintf(1, '  Number of obs.:       %d\n', num_of_obs_2);
fprintf(1, '  Number of scans.:     %d\n', num_of_scans_2);
fprintf(1, '  Number of antennas.:  %d\n', num_of_ant_2);
fprintf(1, '  Number of quasars.:   %d\n', num_of_q_2);


% Preallocate data arrays for plotting
obs_bl_1 = nan(num_of_epochs, num_of_bl);
obs_bl_2 = nan(num_of_epochs, num_of_bl);

% dataset 1:
for i_scan = 1 : num_of_scans_1
    % Get epoch index:
    ind_epoch = scan_1(i_scan).mjd == all_scan_epochs;
    % checks:
    if sum(ind_epoch) > 1
        fprintf('WARNING: Multiple scan epichs found in one dataset!'\n);
        keyboard
    elseif sum(ind_epoch) == 0
       % That's OK => no data for ths epoch in the current dataset!
    end
    % Loop over observations:
    for i_obs = 1 : scan_1(i_scan).nobs
        % Get bl ID:
        bl_id = find((sum((bl_ind_list == [unique_ant_id_1(scan_1(i_scan).obs(i_obs).i1) unique_ant_id_1(scan_1(i_scan).obs(i_obs).i2)])') == 2) | (sum((bl_ind_list == [unique_ant_id_1(scan_1(i_scan).obs(i_obs).i2) unique_ant_id_1(scan_1(i_scan).obs(i_obs).i1)])') == 2));
        obs_bl_1(ind_epoch, bl_id) = scan_1(i_scan).obs(i_obs).obs;
    end
end

% dataset 2:
for i_scan = 1 : num_of_scans_2
    % Get epoch index:
    ind_epoch = scan_2(i_scan).mjd == all_scan_epochs;
    % checks:
    if sum(ind_epoch) > 1
        fprintf('WARNING: Multiple scan epichs found in one dataset!'\n);
        keyboard
    elseif sum(ind_epoch) == 0
       % That's OK => no data for ths epoch in the current dataset!
    end
    % Loop over observations:
    for i_obs = 1 : scan_2(i_scan).nobs
        % Get bl ID:
        bl_id = find((sum((bl_ind_list == [unique_ant_id_2(scan_2(i_scan).obs(i_obs).i1) unique_ant_id_2(scan_2(i_scan).obs(i_obs).i2)])') == 2) | (sum((bl_ind_list == [unique_ant_id_2(scan_2(i_scan).obs(i_obs).i2) unique_ant_id_2(scan_2(i_scan).obs(i_obs).i1)])') == 2));
        obs_bl_2(ind_epoch, bl_id) = scan_2(i_scan).obs(i_obs).obs;
    end
end

% Compute differences:
obs_bl_diff = obs_bl_1 - obs_bl_2;

% Select data to plot:
switch(plot_type)
    case 'diff'
        plot_data = obs_bl_diff;
        str_title = [parameter_1.session_name, ' (', parameter_1.data_type, ') - ', parameter_2.session_name, ' (', parameter_2.data_type ')'];
        str_ylabel = 'Delay differences [sec]';
    case '1'
        plot_data = obs_bl_1;
        str_title = [parameter_1.session_name, ' (', parameter_1.data_type, ')'];
        str_ylabel = 'Delay [sec]';
    case '2'
        plot_data = obs_bl_2;
        str_title = [parameter_2.session_name, ' (', parameter_2.data_type, ')'];
        str_ylabel = 'Delay [sec]';
end

% Plot data:
h_fig = figure;
hold on
leg_cell = {};
t_ref = all_scan_epochs(1);         % Reference epoch for the plot
cmap = colormap(jet(num_of_bl));    % Color map
% Loop over baselines:
for i_bl = 1 : num_of_bl
    if ~flag_all_ant
       if (~strcmp(bl_ant_code_list{i_bl, 1}, [' ', ant_code])) && (~strcmp(bl_ant_code_list{i_bl, 2}, [' ', ant_code]))
           continue
       end
    end
    plot_ind = ~isnan(plot_data(:, i_bl));
    plot((all_scan_epochs(plot_ind) - t_ref)*24, plot_data(plot_ind, i_bl), '-*', 'color', cmap(i_bl, :));
    if sum(plot_ind) ~= 0
        leg_cell = [leg_cell, {[bl_ant_code_list{i_bl,1}, bl_ant_code_list{i_bl,2}]}];
    end
end

xlabel(sprintf('hours since %s', mjd2datestr(t_ref)));
ylabel(str_ylabel)
title(str_title)
h_leg = legend(leg_cell);

if flag_save_plot
    save_pdf_a4_landscape(h_fig, './', [str_title, '.pdf'])
    fprintf(1, 'Plot saved as:   %s\n', ['./', str_title, '.pdf']);
end


return
