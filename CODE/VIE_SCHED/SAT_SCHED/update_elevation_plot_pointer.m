% #########################################################################
% #     update_elevation_plot_pointer
% #########################################################################
%
% DESCRIPTION
%   This function is used to update and modify the content of the
%   "elevation plot".
%
% CREATED  
%   2015-04-21     Andreas Hellerschmied
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
%   - 

%   - PARA                      : Vie_Sched parameter vector

%
%
% OUTPUT
%   - sched_handles       : Graphic handles structure for scheduling (updated)
%   - error_code          : Error Code (0 = no erros occured)
%   - error_msg           : Error Message (empty, if no errors occured)
%
% CHANGES:
%   - 2018-12-20: A. Corbin       : unnecessary entries in legend removed
%

function [sched_handles, error_code, error_msg] = update_elevation_plot_pointer(sched_handles, obs_data, plot_opt, t_epoch_jd, PARA, current_station_id_list)

    % ##### Init #####
    error_code = 0;
    error_msg = '';
    
    reset_time_str = '--:--:--';
    
    % Flag init.:
    flag_delete_start_pointers_and_text = 0;
    flag_delete_end_pointers_and_text = 0;
    flag_delete_end_of_last_obs_pointers_and_text = 0;
    flag_update_start_pointers_and_text = 0;
    flag_update_end_pointers_and_text = 0;
    flag_update_end_of_last_obs_pointers_and_text = 0;
    
    flag_delete_current_epoch_marker = 0;
    flag_update_current_epoch_marker = 0;
    
    
    % Check, if there is an elevation plot:
    if isempty(sched_handles.el_plot)
%         error_code = 1;
%         error_msg = 'No elevation plot available.';
        return; 
    end
    
    % Check, if the ID list of the curretn sattion netword is not empty:
    if isempty(current_station_id_list)
        error_code = 1;
        error_msg = 'Station ID list is empty!';
        return; 
    end
    
 
    % ##### Draw current epoch markers #####
    
    
    % ##### Treat "station sub-plots" #####
    for i_stat = 1 : size(sched_handles.el_plot.stat,2) % loop over all station sub-plots
        
        station_id = sched_handles.el_plot.station_id_list(i_stat);
        
        if ~isempty(find(current_station_id_list == station_id)) % If the station in the elevation plot is inclusded in the current station network ("current_station_id_list")...

            % Set current axes:
            axes(sched_handles.el_plot.stat(i_stat).h_subplot);


            % ##### Set options and get antenna elevation angles + epoch according to "plot_opt" #####
            switch(plot_opt)

                case 0 % - 0 :   Delete all pointers only
                    flag_delete_start_pointers_and_text = 1;
                    flag_delete_end_pointers_and_text = 1;
                    flag_delete_end_of_last_obs_pointers_and_text = 1;
                    flag_delete_current_epoch_marker = 1;

                case 1 %- 1 :   Set the end_of_last_obs. (eolo) pointers to the epoch of the end of last observation
                    flag_delete_end_of_last_obs_pointers_and_text = 1;
                    pointer_el_rad_eolo = obs_data.stat(station_id).end_of_last_obs.el;
                    antenna_pos_time_eolo = obs_data.stat(station_id).end_of_last_obs.jd;
                    flag_update_end_of_last_obs_pointers_and_text = 1;
                    flag_delete_current_epoch_marker = 1;

                case 2 % - 2 :   Set start pointers to he epoch of the begin of the new scan (temp)
                    flag_delete_start_pointers_and_text = 1;
                    pointer_el_rad_start = obs_data.stat(station_id).begin_of_new_obs_temp.el;
                    antenna_pos_time_start = obs_data.stat(station_id).begin_of_new_obs_temp.jd;
                    flag_update_start_pointers_and_text = 1;
                    flag_delete_current_epoch_marker = 1;

                case 3 % - 3 :   Set end pointers to the epoch of end of the new scan (temp)
                    flag_delete_end_pointers_and_text = 1;
                    pointer_el_rad_end = obs_data.stat(station_id).end_of_new_obs_temp.el;
                    antenna_pos_time_end = obs_data.stat(station_id).end_of_new_obs_temp.jd;
                    flag_update_end_pointers_and_text = 1;
                    flag_delete_current_epoch_marker = 1;

                case 4 % - 4 :   Set end_of_last_obs. (eolo) pointers to he epoch of end of last obs and delete start/end pointers
                    flag_delete_start_pointers_and_text = 1;
                    flag_delete_end_pointers_and_text = 1;
                    pointer_el_rad_eolo = obs_data.stat(station_id).end_of_last_obs.el;
                    antenna_pos_time_eolo = obs_data.stat(station_id).end_of_last_obs.jd;
                    flag_update_end_of_last_obs_pointers_and_text = 1;
                    flag_delete_current_epoch_marker = 1;
                    flag_delete_end_of_last_obs_pointers_and_text = 1;

                case 5 % - 5 :   Set start pointers to he epoch of the begin of the new scan (temp) and delete end/current-epoch pointers
                    flag_delete_start_pointers_and_text = 1;
                    flag_delete_end_pointers_and_text = 1;
                    pointer_el_rad_start = obs_data.stat(station_id).begin_of_new_obs_temp.el;
                    antenna_pos_time_start = obs_data.stat(station_id).begin_of_new_obs_temp.jd;
                    flag_update_start_pointers_and_text = 1;
                    flag_delete_current_epoch_marker = 1;

                case 6 % - 6 :   Set start pointers to he epoch of the begin of the new scan (temp) and end pointers to az/el at the end of the new scan (temp)
                    flag_delete_start_pointers_and_text = 1;
                    flag_delete_end_pointers_and_text = 1;
                    pointer_el_rad_start = obs_data.stat(station_id).begin_of_new_obs_temp.el;
                    antenna_pos_time_start = obs_data.stat(station_id).begin_of_new_obs_temp.jd;
                    pointer_el_rad_end = obs_data.stat(station_id).end_of_new_obs_temp.el;
                    antenna_pos_time_end = obs_data.stat(station_id).end_of_new_obs_temp.jd;
                    flag_update_start_pointers_and_text = 1;
                    flag_update_end_pointers_and_text = 1;
                    flag_delete_current_epoch_marker = 1;

                case 9 % - 9 :   Set current epoch pointer to "t_epoch_jd"
                    flag_delete_current_epoch_marker = 1;
                    flag_update_current_epoch_marker = 1;

                    if isempty(t_epoch_jd)
                        error_code = 2;
                        error_msg = 't_epoch_jd is empty!';
                       return; 
                    end

                otherwise % Error
                    error_code = 1;
                    error_msg = 'Unknown plot option ("plot_opt".)';
                    return;

            end


            % ##### Delete previous pointers and reset texts, if they exist #####

            if flag_delete_current_epoch_marker
                if isfield(sched_handles.el_plot.stat(i_stat), 'h_current_epoch_marker') 
                    if ~isempty(sched_handles.el_plot.stat(i_stat).h_current_epoch_marker)
                       delete(sched_handles.el_plot.stat(i_stat).h_current_epoch_marker); % Remove an existing pointer
                    end
                end
            end



            if flag_delete_start_pointers_and_text
                if isfield(sched_handles.el_plot.stat(i_stat), 'h_epoch_marker_start') 
                    if ~isempty(sched_handles.el_plot.stat(i_stat).h_epoch_marker_start)
                       delete(sched_handles.el_plot.stat(i_stat).h_epoch_marker_start); % Remove an existing pointer
                    end
                end
                % Text
    %                 set(sched_handles.el_plot.stat(i_stat).h_text_antenna_el_start, 'String', '');
    %                 set(sched_handles.el_plot.stat(i_stat).h_text_antenna_pos_time_start, 'String', '');
            end


            if flag_delete_end_pointers_and_text
                if isfield(sched_handles.el_plot.stat(i_stat), 'h_epoch_marker_end') 
                    if ~isempty(sched_handles.el_plot.stat(i_stat).h_epoch_marker_end)
                       delete(sched_handles.el_plot.stat(i_stat).h_epoch_marker_end); % Remove an existing pointer
                    end
                end
                % Text
    %                 set(sched_handles.el_plot.stat(i_stat).h_text_antenna_el_end, 'String', '');
    %                 set(sched_handles.el_plot.stat(i_stat).h_text_antenna_pos_time_end, 'String', '');
            end


            if flag_delete_end_of_last_obs_pointers_and_text
                % Az:
                if isfield(sched_handles.el_plot.stat(i_stat), 'h_epoch_marker_end_of_last_obs') 
                    if ~isempty(sched_handles.el_plot.stat(i_stat).h_epoch_marker_end_of_last_obs)
                       delete(sched_handles.el_plot.stat(i_stat).h_epoch_marker_end_of_last_obs); % Remove an existing pointer
                    end
                end
                % Text
    %                 set(sched_handles.el_plot.stat(i_stat).h_text_end_of_last_obs_el, 'String', '');
    %                 set(sched_handles.el_plot.stat(i_stat).h_text_end_of_last_obs_time, 'String', '');
            end





            % ##### Update plots #####


            if flag_update_end_of_last_obs_pointers_and_text
                if isempty(antenna_pos_time_eolo)
                    antenna_pos_time_eolo = PARA.startmjd + 2400000.5;
                end
                % Draw epoch marker (vertical line) 
                temp = get(gca, 'YLim');
                y_min = temp(1);
                y_max = temp(2); 
                sched_handles.el_plot.stat(i_stat).h_epoch_marker_end_of_last_obs = line([antenna_pos_time_eolo, antenna_pos_time_eolo], [y_min, y_max], 'color', sched_handles.color_markers_end_of_last_scan, 'DisplayName', 'end epoch last scan');
            end


            if flag_update_start_pointers_and_text
                % Draw epoch marker (vertical line) 
                temp = get(gca, 'YLim');
                y_min = temp(1);
                y_max = temp(2); 
                sched_handles.el_plot.stat(i_stat).h_epoch_marker_start = line([antenna_pos_time_start, antenna_pos_time_start], [y_min, y_max], 'color', sched_handles.color_markers_scan_start, 'DisplayName', 'start epoch new scan'); % error: Vectors must be the same length.
            end



            if flag_update_end_pointers_and_text
                % Draw epoch marker (vertical line) 
                temp = get(gca, 'YLim');
                y_min = temp(1);
                y_max = temp(2); 
                sched_handles.el_plot.stat(i_stat).h_epoch_marker_end = line([antenna_pos_time_end, antenna_pos_time_end], [y_min, y_max], 'color', sched_handles.color_markers_scan_end, 'DisplayName', 'end epoch new scan');
            end

            if flag_update_current_epoch_marker
                % Draw epoch marker (vertical line) 
                temp = get(gca, 'YLim');
                y_min = temp(1);
                y_max = temp(2); 
                sched_handles.el_plot.stat(i_stat).h_current_epoch_marker = line([t_epoch_jd, t_epoch_jd], [y_min, y_max], 'color', sched_handles.color_current_epoch_marker, 'DisplayName','current epoch');
            end

        
        end % if ~isempty(find(current_station_id_list == station_id)) 
        
    end % for i_stat = 1 : size(sched_handles.el_plot.stat,2) =>  loop over all station sub-plots
    
    
    
    
    %% ##### Treat "Available obs. periods sub-plot" #####
    
    % Set current axes:
    axes(sched_handles.el_plot.h_available_obs_time_subplot);
    
    
    % ##### Set options and get antenna elevation angles + epoch according to "plot_opt" #####
    station_id = current_station_id_list(1);
    switch(plot_opt)

        case 0 % - 0 :   Delete all pointers only
            flag_delete_start_pointers_and_text = 1;
            flag_delete_end_pointers_and_text = 1;
            flag_delete_end_of_last_obs_pointers_and_text = 1;
            flag_delete_current_epoch_marker = 1;

        case 1 %- 1 :   Set the end_of_last_obs. (eolo) pointers to the epoch of the end of last observation
            flag_delete_end_of_last_obs_pointers_and_text = 1;
            pointer_el_rad_eolo = obs_data.stat(station_id).end_of_last_obs.el;
            antenna_pos_time_eolo = obs_data.stat(station_id).end_of_last_obs.jd;
            flag_update_end_of_last_obs_pointers_and_text = 1;
            flag_delete_current_epoch_marker = 1;

        case 2 % - 2 :   Set start pointers to he epoch of the begin of the new scan (temp)
            flag_delete_start_pointers_and_text = 1;
            pointer_el_rad_start = obs_data.stat(station_id).begin_of_new_obs_temp.el;
            antenna_pos_time_start = obs_data.stat(station_id).begin_of_new_obs_temp.jd;
            flag_update_start_pointers_and_text = 1;
            flag_delete_current_epoch_marker = 1;

        case 3 % - 3 :   Set end pointers to the epoch of end of the new scan (temp)
            flag_delete_end_pointers_and_text = 1;
            pointer_el_rad_end = obs_data.stat(station_id).end_of_new_obs_temp.el;
            antenna_pos_time_end = obs_data.stat(station_id).end_of_new_obs_temp.jd;
            flag_update_end_pointers_and_text = 1;
            flag_delete_current_epoch_marker = 1;

        case 4 % - 4 :   Set end_of_last_obs. (eolo) pointers to he epoch of end of last obs and delete start/end pointers
            flag_delete_start_pointers_and_text = 1;
            flag_delete_end_pointers_and_text = 1;
            pointer_el_rad_eolo = obs_data.stat(station_id).end_of_last_obs.el;
            antenna_pos_time_eolo = obs_data.stat(station_id).end_of_last_obs.jd;
            flag_update_end_of_last_obs_pointers_and_text = 1;
            flag_delete_current_epoch_marker = 1;
            flag_delete_end_of_last_obs_pointers_and_text = 1;

        case 5 % - 5 :   Set start pointers to he epoch of the begin of the new scan (temp) and delete end/current-epoch pointers
            flag_delete_start_pointers_and_text = 1;
            flag_delete_end_pointers_and_text = 1;
            pointer_el_rad_start = obs_data.stat(station_id).begin_of_new_obs_temp.el;
            antenna_pos_time_start = obs_data.stat(station_id).begin_of_new_obs_temp.jd;
            flag_update_start_pointers_and_text = 1;
            flag_delete_current_epoch_marker = 1;

        case 6 % - 6 :   Set start pointers to he epoch of the begin of the new scan (temp) and end pointers to az/el at the end of the new scan (temp)
            flag_delete_start_pointers_and_text = 1;
            flag_delete_end_pointers_and_text = 1;
            pointer_el_rad_start = obs_data.stat(station_id).begin_of_new_obs_temp.el;
            antenna_pos_time_start = obs_data.stat(station_id).begin_of_new_obs_temp.jd;
            pointer_el_rad_end = obs_data.stat(station_id).end_of_new_obs_temp.el;
            antenna_pos_time_end = obs_data.stat(station_id).end_of_new_obs_temp.jd;
            flag_update_start_pointers_and_text = 1;
            flag_update_end_pointers_and_text = 1;
            flag_delete_current_epoch_marker = 1;

        case 9 % - 9 :   Set current epoch pointer to "t_epoch_jd"
            flag_delete_current_epoch_marker = 1;
            flag_update_current_epoch_marker = 1;

            if isempty(t_epoch_jd)
                error_code = 2;
                error_msg = 't_epoch_jd is empty!';
               return; 
            end

        otherwise % Error
            error_code = 1;
            error_msg = 'Unknown plot option ("plot_opt".)';
            return;

    end

    
    % ##### Delete previous pointers and reset texts, if they exist #####
        
    if flag_delete_current_epoch_marker
        if isfield(sched_handles.el_plot, 'h_av_obs_times_marker_current_epoch') 
            if ~isempty(sched_handles.el_plot.h_av_obs_times_marker_current_epoch)
               delete(sched_handles.el_plot.h_av_obs_times_marker_current_epoch); % Remove an existing pointer
            end
        end
        set(sched_handles.el_plot.h_text_av_obs_time_current_epoch, 'String', reset_time_str);
    end
    
     if flag_delete_start_pointers_and_text
        if isfield(sched_handles.el_plot, 'h_av_obs_times_marker_start') 
            if ~isempty(sched_handles.el_plot.h_av_obs_times_marker_start)
               delete(sched_handles.el_plot.h_av_obs_times_marker_start); % Remove an existing pointer
            end
        end
    end


    if flag_delete_end_pointers_and_text
        if isfield(sched_handles.el_plot, 'h_av_obs_times_marker_end') 
            if ~isempty(sched_handles.el_plot.h_av_obs_times_marker_end)
               delete(sched_handles.el_plot.h_av_obs_times_marker_end); % Remove an existing pointer
            end
        end
    end


    if flag_delete_end_of_last_obs_pointers_and_text
        if isfield(sched_handles.el_plot, 'h_av_obs_times_marker_end_of_last_scan') 
            if ~isempty(sched_handles.el_plot.h_av_obs_times_marker_end_of_last_scan)
               delete(sched_handles.el_plot.h_av_obs_times_marker_end_of_last_scan); % Remove an existing pointer
            end
        end
        set(sched_handles.el_plot.h_text_av_obs_time_end_of_last_scan, 'String', reset_time_str);
    end

    
    
    
    % ##### Update plots #####
    
    if flag_update_current_epoch_marker
        % Draw epoch marker (vertical line) 
        temp = get(gca, 'YLim');
        y_min = temp(1);
        y_max = temp(2); 
        sched_handles.el_plot.h_av_obs_times_marker_current_epoch = line([t_epoch_jd, t_epoch_jd], [y_min, y_max], 'color', sched_handles.color_current_epoch_marker);
        % Text:
        [year, mon, day, hr, min, sec] = invjday(t_epoch_jd); % Time conversion
        temp_str = sprintf('%02.0f:%02.0f:%02.0f', hr, min, sec);
        set(sched_handles.el_plot.h_text_av_obs_time_current_epoch, 'String', temp_str);
    end
    
    if flag_update_start_pointers_and_text
        % Draw epoch marker (vertical line) 
        temp = get(gca, 'YLim');
        y_min = temp(1);
        y_max = temp(2); 
        if ~isempty(antenna_pos_time_start)
            sched_handles.el_plot.h_av_obs_times_marker_start = line([antenna_pos_time_start, antenna_pos_time_start], [y_min, y_max], 'color', sched_handles.color_markers_scan_start);
        end
    end
    
    if flag_update_end_pointers_and_text
        % Draw epoch marker (vertical line) 
        temp = get(gca, 'YLim');
        y_min = temp(1);
        y_max = temp(2); 
        if ~isempty(antenna_pos_time_end)
            sched_handles.el_plot.h_av_obs_times_marker_end = line([antenna_pos_time_end, antenna_pos_time_end], [y_min, y_max], 'color', sched_handles.color_markers_scan_end);
        end
    end
    
    
    if flag_update_end_of_last_obs_pointers_and_text
        if isempty(obs_data.end_of_last_scan_jd)
            antenna_pos_time_eolo = PARA.startmjd + 2400000.5;
        else
            antenna_pos_time_eolo = obs_data.end_of_last_scan_jd;
        end
        % Draw epoch marker (vertical line) 
        temp = get(gca, 'YLim');
        y_min = temp(1);
        y_max = temp(2); 
        sched_handles.el_plot.h_av_obs_times_marker_end_of_last_scan = line([antenna_pos_time_eolo, antenna_pos_time_eolo], [y_min, y_max], 'color', sched_handles.color_markers_end_of_last_scan);
        % Text:
        [year, mon, day, hr, min, sec] = invjday(antenna_pos_time_eolo); % Time conversion
        temp_str = sprintf('%02.0f:%02.0f:%02.0f', hr, min, sec);
        set(sched_handles.el_plot.h_text_av_obs_time_end_of_last_scan, 'String', temp_str);
    end
    

    
% Handles:
%     h_av_obs_times_marker_start
%     h_av_obs_times_marker_end
%     h_av_obs_times_marker_end_of_last_scan
%     h_av_obs_times_marker_current_epoch

%             % Update text:
%             
%             % time:
%             if isempty(antenna_pos_time_eolo)
%                 antenna_pos_time_str_end = '--:--:--';
%             else
%                 [year, mon, day, hr, min, sec] = invjday(antenna_pos_time_eolo); % Time conversion
%                 antenna_pos_time_eolo_str = sprintf('%02.0f:%02.0f:%02.0f', hr, min, sec);
%             end
%             set(sched_handles.el_plot.stat(i_stat).h_text_end_of_last_obs_time, 'String', antenna_pos_time_eolo_str);

    
    
return
