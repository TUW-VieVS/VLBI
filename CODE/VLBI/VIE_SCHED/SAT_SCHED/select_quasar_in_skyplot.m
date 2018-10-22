% #########################################################################
% #     select_quasar_in_skyplot
% #########################################################################
%
% DESCRIPTION
%   This function 
%
% CREATED  
%   2015-07-02     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
%   - 
%
%
% INPUT
%   - sched_handles                 : station data structure
%   - station_id_list               : List of IDs of the quasar observing stations (referring to the required entry of "stat_data.stat" and "obs_data.stat")
%
%
%   - stat_data                     : station data structure

%   - obs_data                      : Observation data structure (preallocated in "find_observation_times.m")
%   - quasar_id                     : Flag list, to define observable quasars at "t_epoch_jd"
%   - t_epoch_jd                    : Calculation epoch
%   - PARA                          : Global scheduling parameter strucutre
%
%
% OUTPUT
%   - quasar_id                     : ID of the selected source (used to ref. to "source" structure)
%   - error_code                    : Error Code (0 = no erros occured)
%   - error_msg                     : Error Message (empty, if no errors occured)
%
% CHANGES:
%

function [error_code, error_msg, quasar_id] = select_quasar_in_skyplot(sched_handles, station_id)

    % ##### Init #####
    error_code = 0;
    error_msg = '';
    
    
    % ##### Select source in skyplot #####
    
    % Set figure active
    figure(sched_handles.sky_plots(station_id).h_fig)
    fprintf(1, '=> Select quasar in the %s \n', sched_handles.sky_plots(station_id).h_fig.Name);
    
    %  Mouse cursor input 
    [x_in,y_in]=ginput(1);
    
    
    % ##### Select source by minimum distance #####

    % Init.:
    quasar_id_list = find(sched_handles.sky_plots(station_id).quasars.flag_observable_list);
    dist_list = zeros(length(quasar_id_list),2); % | quasar_id | distance |
    
    %Loop over all sources in the current skyplot:
    for i_quasar = 1 : length(quasar_id_list)
        % Calculate distance between
        x_plot = sched_handles.sky_plots(station_id).quasars.xx(i_quasar);
        y_plot = sched_handles.sky_plots(station_id).quasars.yy(i_quasar);
        dist_list(i_quasar, 2) = sqrt((x_plot - x_in)^2 + (y_plot - y_in)^2);
        dist_list(i_quasar, 1) = quasar_id_list(i_quasar);
    end
    % sort list (rows) by distance:
    dist_list = sortrows(dist_list, 2);
    quasar_id = dist_list(1,1);
    
return;