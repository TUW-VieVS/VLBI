% -------------------------------------------------------------------------
%
%                              create_skyplot.m
%
% Function to create skyplots for all stations, visualizing satellite tracks.  
%
%   Author: 
%       Andreas Hellerschmied, 25.02.2013
%   
%   changes       :
%   - 2014-08.27: A. Hellerschmied: Print Legends in Plots
%   - 2014-08.28: A. Hellerschmied: plot markers as "timescale" in skyplot. 
%   - 2015-07-29: A. Hellerschmied: Horizintal mask is plotted additionally to the cut-off elevation
%   - 2016-12-22, A. Hellerschmied: - PARA.INIT_PROP_INTERVAL in [sec] instead of [min]
%   - 2018-12-20: A. Corbin       : unnecessary entries in legend removed
%           
%
%   inputs        :
%   - stat_data         : station data structure
%   - obs_data          : observation data structure
%   - PARA              : VieVS parameter vector 
%   - sched_handles     : Graphic handles structure for scheduling
%     
%
%   outputs       :
%   - Screen Output (Command Window)
%   - error_code        : Error Code (0 = no erros occured)
%   - error_msg         : Error Message (empty, if no errors occured)
%   - sched_handles     : Graphic handles structure for scheduling
%
%
%   locals        :
% 
%
%   coupling      :
%   - skyplot.m
%   
%
%   references    :
%
%-------------------------------------------------------------------------

function [error_code, error_msg, sched_handles] = create_skyplot(stat_data, obs_data, PARA, sched_handles)

    % Init
    error_code = 0;
    error_msg = '';
    sat_names = '';


    % From stat_data struct:
    number_of_sat = stat_data.number_of_sat;
    number_of_stat = stat_data.number_of_stations;
    number_of_epochs = stat_data.number_of_epochs;
    
    % Get propagation interval [min]
    prop_int_min = PARA.INIT_PROP_INTERVAL/60;
    
    % ### Set Skyplot properties ###
    line_style = '-';
    % Set marker interval [min]
    marker_int_min = PARA.SKYPLOT_MARKER_INT;
    number_of_az_intervals_cut_off_el = 100;  % Number of azimuth intervals fot plotting the cut-off elevation marker (line plot)
    flag_plot_h_mask = 1;
    
    % Start pointer/status text color:
    color_start =  sched_handles.color_markers_scan_start;
    
    % End pointer/status text color:
    color_end =  sched_handles.color_markers_scan_end;

    color_last_scan_end = sched_handles.color_markers_end_of_last_scan;
   
    
    
    % Print notification, if the chosen marker interval is not valid!
    if mod(marker_int_min / prop_int_min, 1) ~= 0 % plot markers
        fprintf(1, '\n');
        fprintf(1, '--------------------------------------------------------------------------------------------\n');
        fprintf(1, '  Note: Skyplot marker interval = %3.1f min (= equal to initial propagation interval),\n', prop_int_min);
        fprintf(1, '        because the marker interval is not a multiple of the propagation interval!\n');
        fprintf(1, '--------------------------------------------------------------------------------------------\n');
    end
    
    
    
    for i_stat = 1 : number_of_stat
        
        % preallocation
        az_data = zeros(number_of_sat, number_of_epochs);
        el_data = zeros(number_of_sat, number_of_epochs);
        prn = zeros(1, number_of_sat);
       
        for i_sat = 1 : number_of_sat
            
            prn(1, i_sat) = i_sat;
            
            if (i_stat == 1)
                sat_names{i_sat} = [num2str(i_sat), ' - ', stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_header_line];
            end
            
            for i_epoch = 1 : number_of_epochs
                
                if stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).above_min_elevation
                    az = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).az;
                    el = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).el;
                else
                    az = NaN;
                    el = NaN;
                end
                
                az_data(i_sat, i_epoch) = az;
                el_data(i_sat, i_epoch) = el;
                
            end % i_epoch = 1 : number_of_sat

        end % i_sat = 1 : number_of_sat
        
        
        % ##### Create Skyplot #####
        
        % Create empty figure handles field for the current sky plot (preallocation):
        sched_handles.sky_plots(i_stat).h_fig = [];
        
        sched_handles.sky_plots(i_stat).color_start = color_start;
        sched_handles.sky_plots(i_stat).color_end =  color_end;
        
        
        
        % Check, if there is data to plot (at least one satellite is
        % visible from the station):
        if ~(sum(sum(isnan(el_data))) == numel(el_data)) % data is available
        
            
            sched_handles.sky_plots(i_stat).h_fig = figure;
            set(sched_handles.sky_plots(i_stat).h_fig, 'Name', ['Skyplot for Station: ', stat_data.stat(i_stat).name])

            % Set Figure size and position
            pos_fig = get(sched_handles.sky_plots(i_stat).h_fig, 'Position');
            x_fig = pos_fig(1);
            y_fig = pos_fig(2);
            set(sched_handles.sky_plots(i_stat).h_fig, 'Position', [x_fig (y_fig - 200) 600 400]);

            % Timescale markers:
            if mod(marker_int_min / prop_int_min, 1) == 0 % plot markers
                marker_int = marker_int_min / prop_int_min;
                [h_plot, h_axis, plotted_sats] = skyplot(az_data, el_data, prn, marker_int, line_style);
            else % do not plot time markers
                [h_plot, h_axis, plotted_sats] = skyplot(az_data, el_data, prn);
            end
            
            % Save handles to structure:
            sched_handles.sky_plots(i_stat).h_axis = h_axis;
            for i_sat = 1 : number_of_sat
                sched_handles.sky_plots(i_stat).sat(i_sat).h_plot = h_plot(i_sat);
                sched_handles.sky_plots(i_stat).sat(i_sat).flag_plotted = logical(plotted_sats(i_sat));
            end
            

            
            
            % ##### Plot horizontal mask and/or the cut-off elevation #####
            
            % Get cut-off elevation of the station:
            cut_off_el_rad = stat_data.stat(i_stat).min_elevation * pi / 180;
            az_interval = (2*pi) / number_of_az_intervals_cut_off_el;
            az_mask_rad = (0 : az_interval : (2*pi));
            
            if flag_plot_h_mask && (stat_data.stat(i_stat).horizontal_mask_num > 0)

                % Get h-mask data:
                hmasknum = stat_data.stat(i_stat).horizontal_mask_num;
                hmask    = stat_data.stat(i_stat).horizontal_mask;

                % Check the h-mask type:
                if (mod(hmasknum,2) ~= 0)
                    hmasktype = 1;   % step functions
                else
                    hmasktype = 2;   % line segments
                end

                el_mask_rad = zeros(1, length(az_mask_rad));

                for i_az = 1 : length(az_mask_rad)
                    az_rad = az_mask_rad(i_az);

                    % Loop over h-mask segments and get the elevation angle, defined for the given azimuth:
                    for i = 1 : floor(hmasknum/2)   
                        ib = i * 2;
                        if ((az_rad >= hmask(ib-1)) & (az_rad <= hmask(ib+1)))
                            if (hmasktype == 1)
                                el_mask_rad(i_az) = hmask(ib);
                            elseif (hmasktype == 2)
                                el_mask_rad(i_az) = ((hmask(ib+2) - hmask(ib)) / (hmask(ib+1) - hmask(ib-1))) * (az_rad - hmask(ib-1)) + hmask(ib);
                            end

%                             fprintf('el = %6.3f\n', el_mask_rad(i_az)*180/pi);
                            
                            if el_mask_rad(i_az) < cut_off_el_rad
                                el_mask_rad(i_az) = cut_off_el_rad;
                            end

                            break;
                        end
                    end
                end % for i_az = 1 : length(az_mask_rad)

                % Transform elevation angle to a distance to the center of the plot:
                elSpherical_1 = 90*cos(el_mask_rad); 

                % Transform data to Cartesian coordinates:
                yy_mask = elSpherical_1 .* cos(az_mask_rad); 
                xx_mask = elSpherical_1 .* sin(az_mask_rad); 

                
            else % Draw cut-off elevation angle only
                
                % Calc. cut-off el coordinates:
                elSpherical_cut_off_el = 90*cos(cut_off_el_rad); 
                yy_mask = elSpherical_cut_off_el .* cos(az_mask_rad); 
                xx_mask = elSpherical_cut_off_el .* sin(az_mask_rad);
                
            end % if flag_plot_h_mask
            
            % Plot:
            sched_handles.sky_plots(i_stat).h_cut_off_el_marker = line(xx_mask, yy_mask, 'Linestyle', '-.', 'Color', 'red', 'LineWidth', 0.1,  'HandleVisibility','off');
            
            
            % ##### Get observation type of the antenna #####
            switch(stat_data.stat(i_stat).obs_type)
               
                case {'sat_only'}
                    obs_type_str = 'satellite';
                    
                case {'quasar_only'}
                    obs_type_str = 'quasar';
                    
                case {'sat_and_quasar'}
                    obs_type_str = 'satellite & quasar';
                    
                case {'twin_sat'}
                    obs_type_str = ['satellite; twin with: ', stat_data.stat(i_stat).twin_partner];
                    
                case {'twin_quasar'}
                    obs_type_str = ['quasar; twin with: ', stat_data.stat(i_stat).twin_partner];
            end

 
            
            % Title:
            sched_handles.sky_plots(i_stat).h_title = title([stat_data.stat(i_stat).name, ' (obs. type: ', obs_type_str, ')']);
        
            % Legend:
            sched_handles.sky_plots(i_stat).h_legend = legend(sat_names);
            
            % Get positions:
            pos_leg = get(sched_handles.sky_plots(i_stat).h_legend, 'Position');      
            pos_axis = get(h_axis, 'Position');
            
            % Set new positions:
            set(h_axis, 'Position', [0.00 pos_axis(2) pos_axis(3) pos_axis(4)])
            set(sched_handles.sky_plots(i_stat).h_legend, 'Position', [1-pos_leg(3)-0.05   pos_leg(2)    pos_leg(3)    pos_leg(4)])
            

            
            % ##### Print antenna azimuth and elevation + epoch #####
            
            x_text_right_side = 170;
            x_text_left_side = 120;
            
            h_text_8  = text(x_text_left_side, 40, 'End of last obs.:', 'color', color_last_scan_end);
            h_text_9  = text(x_text_left_side, 30, 'unaz [deg] =');
            h_text_10 = text(x_text_left_side, 20, 'el [deg] =');
            h_text_11 = text(x_text_left_side, 10, 'time:');
            
            h_text_0 = text(x_text_left_side, -5, 'Scan start:', 'color', color_start);
            h_text_1 = text(x_text_left_side, -15, 'unaz [deg] =');
            h_text_2 = text(x_text_left_side, -25, 'el [deg] =');
            h_text_3 = text(x_text_left_side, -35, 'time:');
            
            h_text_4 = text(x_text_left_side, -50, 'Scan end:', 'color', color_end);
            h_text_5 = text(x_text_left_side, -60, 'unaz [deg] =');
            h_text_6 = text(x_text_left_side, -70, 'el [deg] =');
            h_text_7 = text(x_text_left_side, -80, 'time:');
            
            h_text_12 = text(x_text_left_side, -95, 'Curr. epoch:', 'color', sched_handles.color_current_epoch_marker);
            
            % End of last scan
            if isempty(obs_data.end_of_last_scan_jd)
                % text_end_of_last_obs_time_str = '--:--:--';
                [year, mon, day, hr, min, sec] = invjday(PARA.startmjd + 2400000.5); % Time conversion
                text_end_of_last_obs_time_str = sprintf('%02.0f:%02.0f:%02.0f', hr, min, sec);
            else
                [year, mon, day, hr, min, sec] = invjday(obs_data.end_of_last_scan_jd); % Time conversion
                text_end_of_last_obs_time_str = sprintf('%02.0f:%02.0f:%02.0f', hr, min, sec);
            end
            text_end_of_last_obs_unaz_str = sprintf('%5.1f', obs_data.stat(i_stat).end_of_last_obs.un_az * 180/pi);
            text_end_of_last_obs_el_str =    sprintf('%5.1f', obs_data.stat(i_stat).end_of_last_obs.el * 180/pi);
            
            sched_handles.sky_plots(i_stat).h_text_end_of_last_obs_unaz = text(x_text_right_side, 30, text_end_of_last_obs_unaz_str);
            sched_handles.sky_plots(i_stat).h_text_end_of_last_obs_el = text(x_text_right_side, 20, text_end_of_last_obs_el_str);
            sched_handles.sky_plots(i_stat).h_text_end_of_last_obs_time = text(x_text_right_side, 10, text_end_of_last_obs_time_str);
            
            % Scan start
            text_antenna_pos_time_start = '--:--:--';
            text_antenna_un_az_str_start = sprintf('---');
            text_antenna_el_str_start =    sprintf('---');
            sched_handles.sky_plots(i_stat).h_text_antenna_un_az_start = text(x_text_right_side,-15,text_antenna_un_az_str_start);
            sched_handles.sky_plots(i_stat).h_text_antenna_el_start = text(x_text_right_side,-25,text_antenna_el_str_start);
            sched_handles.sky_plots(i_stat).h_text_antenna_pos_time_start = text(x_text_right_side,-35,text_antenna_pos_time_start);
            
            % Scan end
            text_antenna_pos_time_end = '--:--:--';
            text_antenna_un_az_str_end = sprintf('---');
            text_antenna_el_str_end =    sprintf('---');
            sched_handles.sky_plots(i_stat).h_text_antenna_un_az_end = text(x_text_right_side,-60,text_antenna_un_az_str_end);
            sched_handles.sky_plots(i_stat).h_text_antenna_el_end = text(x_text_right_side,-70,text_antenna_el_str_end);
            sched_handles.sky_plots(i_stat).h_text_antenna_pos_time_end = text(x_text_right_side,-80,text_antenna_pos_time_end);
            
            % Current epoch:
            text_current_epoch =    sprintf('--:--:--');
            sched_handles.sky_plots(i_stat).h_current_epoch = text(x_text_right_side, -95, text_current_epoch);
            
            
            

            % ##### Draw azimuth pointer #####
            pointer_az_rad = obs_data.stat(i_stat).end_of_last_obs.un_az;

            pointer_el_rad_1 = 0;
            pointer_el_rad_2 = pi/2;

            % Transform elevation angle to a distance to the center of the plot:
            elSpherical_1 = 90*cos(pointer_el_rad_1); 
            elSpherical_2 = 90*cos(pointer_el_rad_2); 

            % Transform data to Cartesian coordinates:
            yy_pointer_1 = elSpherical_1 .* cos(pointer_az_rad); 
            xx_pointer_1 = elSpherical_1 .* sin(pointer_az_rad); 
            yy_pointer_2 = elSpherical_2 .* cos(pointer_az_rad); 
            xx_pointer_2 = elSpherical_2 .* sin(pointer_az_rad);  

            % Plot pointer
            sched_handles.sky_plots(i_stat).h_az_end_of_last_obs_pointer = line([xx_pointer_1, xx_pointer_2], [yy_pointer_1, yy_pointer_2], 'color', color_last_scan_end,  'HandleVisibility','off');
            
            
            

            % ##### Draw elevation pointer #####
            if ~isempty(obs_data.stat(i_stat).end_of_last_obs.el)

                pointer_el_rad = obs_data.stat(i_stat).end_of_last_obs.el;

                % Transform elevation angle to a distance to the center of the plot:
                elSpherical_3 = 90*cos(pointer_el_rad); 

                % Transform data to Cartesian coordinates:
                yy_pointer_3 = elSpherical_3 .* cos(pointer_az_rad); 
                xx_pointer_3 = elSpherical_3 .* sin(pointer_az_rad);

                % Plot pointer
                sched_handles.sky_plots(i_stat).h_el_end_of_last_obs_pointer = plot(xx_pointer_3, yy_pointer_3, 'x', 'color', color_last_scan_end,  'HandleVisibility','off');

            end

        end 
        
    end % i_stat = 1 : number_of_stat
   

return;



