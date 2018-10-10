% #########################################################################
% #     update_elevation_plot
% #########################################################################
%
% DESCRIPTION
%   This function is used to update and modify the content of the
%   "elevation plot".
%
% CREATED  
%   2015-03-11     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
%   - 
%
%
% INPUT
%  
%   - sched_handles             : Graphic handles structure for scheduling
%   - jd_current_epoch_marker   : Currently treated epoch [JD], for setting the epoch marker
%
%
% OUTPUT
%   - sched_handles       : Graphic handles structure for scheduling (updated)
%   - error_code          : Error Code (0 = no erros occured)
%   - error_msg           : Error Message (empty, if no errors occured)
%
% CHANGES:
%

function [sched_handles, error_code, error_msg] = update_elevation_plot(sched_handles, jd_current_epoch_marker)

    % ##### Init #####
    error_code = 0;
    error_msg = '';
    
    % Check, if there is an elevation plot:
    if isempty(sched_handles.el_plot)
%         error_code = 1;
%         error_msg = 'No elevation plot available.';
        return; 
    end
    
    
    %% ##### Draw current epoch markers #####
    
    
    
    % ### Treat "station sub-plots" ###
    for i_stat = 1 : size(sched_handles.el_plot.stat,2) % loop over all station sub-plots
        
        % Set current axes:
        axes(sched_handles.el_plot.stat(i_stat).h_subplot);
        
        % Delete a previous time marker line, if it exists:
        if isfield(sched_handles.el_plot.stat(i_stat), 'h_current_epoch_marker')
            if ~isempty(sched_handles.el_plot.stat(i_stat).h_current_epoch_marker)
               delete(sched_handles.el_plot.stat(i_stat).h_current_epoch_marker); % Remove existing marker
            end
        end
        
        % Draw vertical line 
        temp = get(gca, 'YLim');
        y_min = temp(1);
        y_max = temp(2); 
        sched_handles.el_plot.stat(i_stat).h_current_epoch_marker = line([jd_current_epoch_marker, jd_current_epoch_marker], [y_min, y_max], 'color', 'red');

    end
    
    
    
    % ### Treat "Available obs. periods sub-plot" ###
    
    % Set current axes:
    axes(sched_handles.el_plot.h_available_obs_time_subplot);
    
    % Delete a previous time marker line, if it exists:
    if isfield(sched_handles.el_plot, 'h_available_obs_times_current_epoch_marker')
        if ~isempty(sched_handles.el_plot.h_available_obs_times_current_epoch_marker)
           delete(sched_handles.el_plot.h_available_obs_times_current_epoch_marker); % Remove existing marker
        end
    end
    
    % Draw vertical line 
    temp = get(gca, 'YLim');
    y_min = temp(1);
    y_max = temp(2);
    
    sched_handles.el_plot.h_available_obs_times_current_epoch_marker = line([jd_current_epoch_marker, jd_current_epoch_marker], [y_min, y_max], 'color', 'red');

    
    
return