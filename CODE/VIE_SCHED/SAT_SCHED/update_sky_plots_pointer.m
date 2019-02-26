% #########################################################################
% #     update_sky_plots_pointer
% #########################################################################
%
% DESCRIPTION
%   This function updates the azimuth and elevation pointers of the sky plots included 
%   in the station list ("station_id_list").
%   Also the printed time is updated here.
%   
%   Plot options ( plot_opt = ...):
%       - 0 :   Delete pointers only
%       - 1 :   Set the end_of_last_obs. (eolo) pointers to az/el at the end of last observation
%       - 2 :   Set start pointers to az/el at the begin of the new scan (temp)
%       - 3 :   Set end pointers to az/el at the end of the new scan (temp)
%       - 4 :   Set end_of_last_obs. (eolo) pointers to az/el at the end of last obs and delete start/end pointers
%       - 5 :   Set start pointers to az/el at the begin of the new scan (temp) and delete end pointers
%       - 6 :   Set start pointers to az/el at the begin of the new scan (temp) and end pointers to az/el at the end of the new scan (temp)
%
% CREATED  
%   2015-04-09     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING       
%   - invjday
%
%
% INPUT
%
%   - sched_handles             : Graphic handles structure for scheduling
%   - jd_current_epoch_marker   : Currently treated epoch [JD], for setting the epoch marker
%   - station_id_list           : List of IDs of the observing stations (referring to the required entry of "stat_data.stat" and "obs_data.stat")
%   - obs_data                  : Observation data structure (preallocated in "find_observation_times.m")
%   - plot_opt                  : Plot options (see description above!)
%   - PARA                      : Vie_Sched parameter vector
%
% OUTPUT
%   - sched_handles       : Graphic handles structure for scheduling (updated)
%   - error_code          : Error Code (0 = no erros occured)
%   - error_msg           : Error Message (empty, if no errors occured)
%
% CHANGES:
%   - 2018-12-20: A. Corbin       : unnecessary entries in legend removed
%

function [sched_handles, error_code, error_msg] = update_sky_plots_pointer(sched_handles, station_id_list, obs_data, plot_opt, PARA)

    % ##### Init #####
    error_code = 0;
    error_msg = '';
    
    flag_delete_start_pointers_and_text = 0;
    flag_delete_end_pointers_and_text = 0;
    flag_delete_end_of_last_obs_pointers_and_text = 0;
    flag_update_start_pointers_and_text = 0;
    flag_update_end_pointers_and_text = 0;
    flag_update_end_of_last_obs_pointers_and_text = 0;
    
    reset_unaz_str = '---';
    reset_el_str = '---';
    reset_time_str = '--:--:--';
    
    

    % ##### Check, if there are sky plots #####
    if isempty(sched_handles.sky_plots)
        return; 
    end
    
 
    % Loop over all stations included in the station ID list:
    % => loop over all sky plots included in the current scan

    for i_stat = 1 : length(station_id_list)
        
        station_id = station_id_list(i_stat);
              
        % Check if current sky plot is available:
        if ~isempty(sched_handles.sky_plots(station_id).h_fig)
            

            % ##### Set current axes active #####
            axes(sched_handles.sky_plots(station_id).h_axis);


            % ##### Set options and get antenna azimuth and elevation angles + epoch #####
            switch(plot_opt)

                case 0 % - 0 :   Delete pointers only
                    flag_delete_start_pointers_and_text = 1;
                    flag_delete_end_pointers_and_text = 1;
                    flag_delete_end_of_last_obs_pointers_and_text = 1;

                case 1 %- 1 :   Set the end_of_last_obs. (eolo) pointers to az/el at the end of last observation
                    flag_delete_end_of_last_obs_pointers_and_text = 1;
                    pointer_az_rad_eolo = obs_data.stat(station_id).end_of_last_obs.un_az;
                    pointer_el_rad_eolo = obs_data.stat(station_id).end_of_last_obs.el;
                    antenna_pos_time_eolo = obs_data.stat(station_id).end_of_last_obs.jd;
                    flag_update_end_of_last_obs_pointers_and_text = 1;

                case 2 % - 2 :   Set start pointers to az/el at the begin of the new scan (temp)
                    flag_delete_start_pointers_and_text = 1;
                    pointer_az_rad_start = obs_data.stat(station_id).begin_of_new_obs_temp.un_az;
                    pointer_el_rad_start = obs_data.stat(station_id).begin_of_new_obs_temp.el;
                    antenna_pos_time_start = obs_data.stat(station_id).begin_of_new_obs_temp.jd;
                    flag_update_start_pointers_and_text = 1;

                case 3 % - 3 :   Set end pointers to az/el at the end of the new scan (temp)
                    flag_delete_end_pointers_and_text = 1;
                    pointer_az_rad_end = obs_data.stat(station_id).end_of_new_obs_temp.un_az;
                    pointer_el_rad_end = obs_data.stat(station_id).end_of_new_obs_temp.el;
                    antenna_pos_time_end = obs_data.stat(station_id).end_of_new_obs_temp.jd;
                    flag_update_end_pointers_and_text = 1;

                case 4 % - 4 :   Set end_of_last_obs. (eolo) pointers to az/el at the end of last obs and delete start/end pointers
                    flag_delete_start_pointers_and_text = 1;
                    flag_delete_end_pointers_and_text = 1;
                    pointer_az_rad_eolo = obs_data.stat(station_id).end_of_last_obs.un_az;
                    pointer_el_rad_eolo = obs_data.stat(station_id).end_of_last_obs.el;
                    antenna_pos_time_eolo = obs_data.stat(station_id).end_of_last_obs.jd;
                    flag_update_end_of_last_obs_pointers_and_text = 1;
                    flag_delete_end_of_last_obs_pointers_and_text = 1;

                case 5 % - 5 :   Set start pointers to az/el at the begin of the new scan (temp) and delete end pointers
                    flag_delete_start_pointers_and_text = 1;
                    flag_delete_end_pointers_and_text = 1;
                    pointer_az_rad_start = obs_data.stat(station_id).begin_of_new_obs_temp.un_az;
                    pointer_el_rad_start = obs_data.stat(station_id).begin_of_new_obs_temp.el;
                    antenna_pos_time_start = obs_data.stat(station_id).begin_of_new_obs_temp.jd;
                    flag_update_start_pointers_and_text = 1;
                    
                case 6 % - 6 :   Set start pointers to az/el at the begin of the new scan (temp) and end pointers to az/el at the end of the new scan (temp)
                    flag_delete_start_pointers_and_text = 1;
                    flag_delete_end_pointers_and_text = 1;
                    pointer_az_rad_start = obs_data.stat(station_id).begin_of_new_obs_temp.un_az;
                    pointer_el_rad_start = obs_data.stat(station_id).begin_of_new_obs_temp.el;
                    antenna_pos_time_start = obs_data.stat(station_id).begin_of_new_obs_temp.jd;
                    pointer_az_rad_end = obs_data.stat(station_id).end_of_new_obs_temp.un_az;
                    pointer_el_rad_end = obs_data.stat(station_id).end_of_new_obs_temp.el;
                    antenna_pos_time_end = obs_data.stat(station_id).end_of_new_obs_temp.jd;
                    flag_update_start_pointers_and_text = 1;
                    flag_update_end_pointers_and_text = 1;
                    
                otherwise % Error
                    error_code = 1;
                    error_msg = 'Unknown plot option ("plot_opt".)';

            end



            % ##### Delete previous pointers and reset texts, if they exist #####
            
            if flag_delete_start_pointers_and_text
            
                % Az:
                if isfield(sched_handles.sky_plots(station_id), 'h_azimuth_start_pointer') 
                    if ~isempty(sched_handles.sky_plots(station_id).h_azimuth_start_pointer)
                       delete(sched_handles.sky_plots(station_id).h_azimuth_start_pointer); % Remove an existing pointer
                    end
                end

                % El:
                if isfield(sched_handles.sky_plots(station_id), 'h_elevation_start_pointer') 
                    if ~isempty(sched_handles.sky_plots(station_id).h_elevation_start_pointer)
                       delete(sched_handles.sky_plots(station_id).h_elevation_start_pointer); % Remove an existing pointer
                    end
                end
                
                % Text
                set(sched_handles.sky_plots(station_id).h_text_antenna_un_az_start, 'String', reset_unaz_str);
                set(sched_handles.sky_plots(station_id).h_text_antenna_el_start, 'String', reset_el_str);
                set(sched_handles.sky_plots(station_id).h_text_antenna_pos_time_start, 'String', reset_time_str);
                
            end
            
            
            if flag_delete_end_pointers_and_text
            
                % Az:
                if isfield(sched_handles.sky_plots(station_id), 'h_azimuth_end_pointer') 
                    if ~isempty(sched_handles.sky_plots(station_id).h_azimuth_end_pointer)
                       delete(sched_handles.sky_plots(station_id).h_azimuth_end_pointer); % Remove an existing pointer
                    end
                end

                % El:
                if isfield(sched_handles.sky_plots(station_id), 'h_elevation_end_pointer') 
                    if ~isempty(sched_handles.sky_plots(station_id).h_elevation_end_pointer)
                       delete(sched_handles.sky_plots(station_id).h_elevation_end_pointer); % Remove an existing pointer
                    end
                end
                
                % Text
                set(sched_handles.sky_plots(station_id).h_text_antenna_un_az_end, 'String', reset_unaz_str);
                set(sched_handles.sky_plots(station_id).h_text_antenna_el_end, 'String', reset_el_str);
                set(sched_handles.sky_plots(station_id).h_text_antenna_pos_time_end, 'String', reset_time_str);
                
            end
            
            
            if flag_delete_end_of_last_obs_pointers_and_text
                
                % Az:
                if isfield(sched_handles.sky_plots(station_id), 'h_az_end_of_last_obs_pointer') 
                    if ~isempty(sched_handles.sky_plots(station_id).h_az_end_of_last_obs_pointer)
                       delete(sched_handles.sky_plots(station_id).h_az_end_of_last_obs_pointer); % Remove an existing pointer
                    end
                end

                % El:
                if isfield(sched_handles.sky_plots(station_id), 'h_el_end_of_last_obs_pointer') 
                    if ~isempty(sched_handles.sky_plots(station_id).h_el_end_of_last_obs_pointer)
                       delete(sched_handles.sky_plots(station_id).h_el_end_of_last_obs_pointer); % Remove an existing pointer
                    end
                end
                
                % Text
                set(sched_handles.sky_plots(station_id).h_text_end_of_last_obs_unaz, 'String', reset_unaz_str);
                set(sched_handles.sky_plots(station_id).h_text_end_of_last_obs_el, 'String', reset_el_str);
                set(sched_handles.sky_plots(station_id).h_text_end_of_last_obs_time, 'String', reset_time_str);
               
            end
            
            

            % ##### Update pointers #####
            
            if flag_update_end_of_last_obs_pointers_and_text
                
                % ##### Draw azimuth pointer #####
                pointer_el_rad_1 = 0;
                pointer_el_rad_2 = pi/2;

                % Transform elevation angle to a distance to the center of the plot:
                elSpherical_1 = 90*cos(pointer_el_rad_1); 
                elSpherical_2 = 90*cos(pointer_el_rad_2); 

                % Transform data to Cartesian coordinates:
                yy_pointer_1 = elSpherical_1 .* cos(pointer_az_rad_eolo); 
                xx_pointer_1 = elSpherical_1 .* sin(pointer_az_rad_eolo); 
                yy_pointer_2 = elSpherical_2 .* cos(pointer_az_rad_eolo); 
                xx_pointer_2 = elSpherical_2 .* sin(pointer_az_rad_eolo); 

                % Plot pointer
                sched_handles.sky_plots(station_id).h_az_end_of_last_obs_pointer = line([xx_pointer_1, xx_pointer_2], [yy_pointer_1, yy_pointer_2], ...
                                                                                        'color', sched_handles.color_markers_end_of_last_scan, ...
                                                                                        'DisplayName','end epoch last scan'); 

                % Update displayed azimuth value:
                text_antenna_un_az_str_start = sprintf('%5.1f', pointer_az_rad_eolo * 180/pi);
                set(sched_handles.sky_plots(station_id).h_text_end_of_last_obs_unaz, 'String', text_antenna_un_az_str_start);


                % ##### Draw elevation pointer and notification text #####
                if ~isempty(pointer_el_rad_eolo) 
                    % Transform elevation angle to a distance to the center of the plot:
                    elSpherical_3 = 90*cos(pointer_el_rad_eolo); 

                    % Transform data to Cartesian coordinates:
                    yy_pointer_3 = elSpherical_3 .* cos(pointer_az_rad_eolo); 
                    xx_pointer_3 = elSpherical_3 .* sin(pointer_az_rad_eolo);

                    % Plot pointer
                    sched_handles.sky_plots(station_id).h_el_end_of_last_obs_pointer = plot(xx_pointer_3, yy_pointer_3, 'x', 'color', sched_handles.color_markers_end_of_last_scan, 'HandleVisibility','off');

                    % Update displayed elevation value:
                    text_antenna_el_str_eolo = sprintf('%5.1f', pointer_el_rad_eolo * 180/pi);
                    set(sched_handles.sky_plots(station_id).h_text_end_of_last_obs_el, 'String', text_antenna_el_str_eolo)
                end


                % ##### Update antenna position time text #####
                if isempty(antenna_pos_time_eolo) % first observation
                    [year, mon, day, hr, min, sec] = invjday(PARA.startmjd + 2400000.5); % Time conversion
                    antenna_pos_time_str_start = sprintf('%02.0f:%02.0f:%02.0f', hr, min, sec);
                    % antenna_pos_time_str_start = '--:--:--';
                else
                    [year, mon, day, hr, min, sec] = invjday(antenna_pos_time_eolo); % Time conversion
                    antenna_pos_time_str_start = sprintf('%02.0f:%02.0f:%02.0f', hr, min, sec);
                end
                set(sched_handles.sky_plots(station_id).h_text_end_of_last_obs_time, 'String', antenna_pos_time_str_start);
                
            end
            
            
            if flag_update_start_pointers_and_text
                
                
                % ##### Draw azimuth pointer #####
                pointer_el_rad_1 = 0;
                pointer_el_rad_2 = pi/2;

                % Transform elevation angle to a distance to the center of the plot:
                elSpherical_1 = 90*cos(pointer_el_rad_1); 
                elSpherical_2 = 90*cos(pointer_el_rad_2); 

                % Transform data to Cartesian coordinates:
                yy_pointer_1 = elSpherical_1 .* cos(pointer_az_rad_start); 
                xx_pointer_1 = elSpherical_1 .* sin(pointer_az_rad_start); 
                yy_pointer_2 = elSpherical_2 .* cos(pointer_az_rad_start); 
                xx_pointer_2 = elSpherical_2 .* sin(pointer_az_rad_start); 

                % Plot pointer
                sched_handles.sky_plots(station_id).h_azimuth_start_pointer = line([xx_pointer_1, xx_pointer_2], [yy_pointer_1, yy_pointer_2], ...
                                                                                    'color', sched_handles.color_markers_scan_start, ...
                                                                                    'DisplayName', 'start epoch new scan'); 

                % Update displayed azimuth value:
                text_antenna_un_az_str_start = sprintf('%5.1f', pointer_az_rad_start * 180/pi);
                set(sched_handles.sky_plots(station_id).h_text_antenna_un_az_start, 'String', text_antenna_un_az_str_start);


                % ##### Draw elevation pointer and notification text #####
                if ~isempty(pointer_el_rad_start) 
                    % Transform elevation angle to a distance to the center of the plot:
                    elSpherical_3 = 90*cos(pointer_el_rad_start); 

                    % Transform data to Cartesian coordinates:
                    yy_pointer_3 = elSpherical_3 .* cos(pointer_az_rad_start); 
                    xx_pointer_3 = elSpherical_3 .* sin(pointer_az_rad_start);

                    % Plot pointer
                    sched_handles.sky_plots(station_id).h_elevation_start_pointer = plot(xx_pointer_3, yy_pointer_3, 'x', 'color', sched_handles.color_markers_scan_start, 'HandleVisibility','off');

                    % Update displayed elevation value:
                    text_antenna_el_str_start = sprintf('%5.1f', pointer_el_rad_start * 180/pi);
                    set(sched_handles.sky_plots(station_id).h_text_antenna_el_start, 'String', text_antenna_el_str_start)
                end


                % ##### Update antenna position time text #####
                if isempty(antenna_pos_time_start)
                    antenna_pos_time_str_start = '--:--:--';
                else
                    [year, mon, day, hr, min, sec] = invjday(antenna_pos_time_start); % Time conversion
                    antenna_pos_time_str_start = sprintf('%02.0f:%02.0f:%02.0f', hr, min, sec);
                end
                set(sched_handles.sky_plots(station_id).h_text_antenna_pos_time_start, 'String', antenna_pos_time_str_start);
                
            end
            
            
            
            if flag_update_end_pointers_and_text
                
                % ##### Draw azimuth pointer #####
                pointer_el_rad_1 = 0;
                pointer_el_rad_2 = pi/2;

                % Transform elevation angle to a distance to the center of the plot:
                elSpherical_1 = 90*cos(pointer_el_rad_1); 
                elSpherical_2 = 90*cos(pointer_el_rad_2); 

                % Transform data to Cartesian coordinates:
                yy_pointer_1 = elSpherical_1 .* cos(pointer_az_rad_end); 
                xx_pointer_1 = elSpherical_1 .* sin(pointer_az_rad_end); 
                yy_pointer_2 = elSpherical_2 .* cos(pointer_az_rad_end); 
                xx_pointer_2 = elSpherical_2 .* sin(pointer_az_rad_end); 

                % Plot pointer
                sched_handles.sky_plots(station_id).h_azimuth_end_pointer = line([xx_pointer_1, xx_pointer_2], [yy_pointer_1, yy_pointer_2], ...
                                                                                  'color', sched_handles.color_markers_scan_end, ...
                                                                                  'DisplayName', 'end epoch new scan'); 

                % Update displayed azimuth value:
                text_antenna_un_az_str_end = sprintf('%5.1f', pointer_az_rad_end * 180/pi);
                set(sched_handles.sky_plots(station_id).h_text_antenna_un_az_end, 'String', text_antenna_un_az_str_end);


                % ##### Draw elevation pointer and notification text #####
                if ~isempty(pointer_el_rad_end) 
                    % Transform elevation angle to a distance to the center of the plot:
                    elSpherical_3 = 90*cos(pointer_el_rad_end); 

                    % Transform data to Cartesian coordinates:
                    yy_pointer_3 = elSpherical_3 .* cos(pointer_az_rad_end); 
                    xx_pointer_3 = elSpherical_3 .* sin(pointer_az_rad_end);

                    % Plot pointer
                    sched_handles.sky_plots(station_id).h_elevation_end_pointer = plot(xx_pointer_3, yy_pointer_3, 'x', 'color', sched_handles.color_markers_scan_end, 'HandleVisibility','off');

                    % Update displayed elevation value:
                    text_antenna_el_str_end = sprintf('%5.1f', pointer_el_rad_end * 180/pi);
                    set(sched_handles.sky_plots(station_id).h_text_antenna_el_end, 'String', text_antenna_el_str_end)
                end


                % ##### Update antenna position time text #####
                if isempty(antenna_pos_time_end)
                    antenna_pos_time_str_end = '--:--:--';
                else
                    [year, mon, day, hr, min, sec] = invjday(antenna_pos_time_end); % Time conversion
                    antenna_pos_time_str_end = sprintf('%02.0f:%02.0f:%02.0f', hr, min, sec);
                end
                set(sched_handles.sky_plots(station_id).h_text_antenna_pos_time_end, 'String', antenna_pos_time_str_end);
                
                
                
            end
            

        end
        
    end % station loop
    

return;
