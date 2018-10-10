% -------------------------------------------------------------------------
%
%                              create_elevation_plot.m
%
%   This function creates a plot window showing all necessary overpass data
%   of an arbitrary number of stations and satellites needed for the further
%   manual satellite scheduling approach.
%
%   Author: 
%       Andreas Hellerschmied (heller182@gmx.at), 27.10.2013
%   
%   changes       :
%   - 2014-01.19: A. Hellerschmied: Support of different antenna axis types
%   - 2014-01.19: A. Hellerschmied: Error treatment added.
%            Remark: Output of coupled procedures is currently not checked!
%   - 2015-02-16: A. Hellerschmied: Remanded (before: "TLE_overpass_plot.m")
%   - 2015-02-16: A. Hellerschmied: Remanded (before: "overpass_plot.m")
%   - 2015-08-18: A. Hellerschmied: Violations of axis limits are plotted
%   - 2016-11-08: A. Hellerschmied: Changed markers for observation condition violations
%   - 2016-11-25: A. Hellerschmied: Added session start time to the xlabel 
%   - 2016-12-02: A. Hellerschmied: - Adaptive x-tick labels
%                                   - optionial input arguments to set x-tick units and according interval
%   - 2017-01-23: A. Hellerschmied: Fixed small bug in defining x-tick separation
%           
%
%   inputs        :
%   - stat_data         : station data structure
%   - obs_data          : observation data structure
%   - PARA              : VieVS parameter vector 
%   - sched_handles     : Graphic handles structure for scheduling
%   - x_tick_unit_str   : Unit for x-axis tick labels ('minutes' or 'hours'), string (optional)
%   - int_unit          : Interval for x-axis tick labels (according unit defined in "x_tick_unit_str"), float (optional)
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
%   - mjd2datestr.m
%   
%
%   references    :
%
%-------------------------------------------------------------------------


function [sched_handles, error_code, error_msg] = create_elevation_plot(stat_data, obs_data, PARA, sched_handles, varargin)

    % Init
    error_code = 0;
    error_msg = '';
    
    
    % Set units for x-tick labels
    switch(nargin)
        case 4
            if PARA.duration <= 0.25
                x_tick_unit_str = 'minutes';
                int_unit = 1;
            elseif PARA.duration <= 1
                x_tick_unit_str = 'minutes';
                int_unit = 5;
            elseif PARA.duration <= 2
                x_tick_unit_str = 'minutes';
                int_unit = 15;
            elseif PARA.duration <= 12
                x_tick_unit_str = 'hours';
                int_unit = 0.5;
            else % > 12hours
                x_tick_unit_str = 'hours';
                int_unit = 1;
            end
        case 6
            x_tick_unit_str = varargin{1};
            int_unit        = varargin{2};
        otherwise
            error_code = 1;
            error_msg = sprintf('Invalid input: nargis = %d', nargin);
    end
    
    % preallocation
    t = zeros(1, stat_data.number_of_epochs);       % Time
    el = zeros(1, stat_data.number_of_epochs);      % Elevation
    sat_names = '';
    current_color = [0 0 0];
    
    station_id_list = [];
    station_counter = 0;
    
    reset_time_str = '--:--:--';
    
    
    obs_periods_el = zeros(stat_data.number_of_epochs, stat_data.number_of_sat, stat_data.number_of_stations);
    
    % Get start time string for the current session/plot:
    sess_start_time_str = mjd2datestr(PARA.startmjd);
    
                   
    %% --- Plot Sat.-Elevation over time for each station --- 
    % figure
    sched_handles.el_plot.h_fig = figure('Name','Elevation Timeseries');
    
    % Set window width
    %win_pos = get(gcf, 'Position');
    screen_size = get(0,'ScreenSize');  % ScreenSize = [left,bottom,width,height], [pixels]
    set(gcf,'Position', [screen_size(3)/10  screen_size(4)/8 screen_size(3)*8/10 screen_size(4)*6/8]);
    
    % Set Default Axis position
    shift_left = 0.08;
    pos = get(sched_handles.el_plot.h_fig, 'DefaultAxesPosition'); 
    % set(figureHandle, 'DefaultAxesPosition'. [0.08, 0.11, 0.70, 0.815]); => standard values
    set(sched_handles.el_plot.h_fig, 'DefaultAxesPosition', [pos(1) - shift_left, pos(2), pos(3), pos(4)]);

          
    
    for i_stat = 1 : stat_data.number_of_stations    % Stations
        
        if ~strcmp(stat_data.stat(i_stat).obs_type, 'quasar_only') && ~strcmp(stat_data.stat(i_stat).obs_type, 'twin_quasar')
            
            station_counter = station_counter + 1;
            station_id_list(station_counter) = i_stat;
        
            % Preset Axis dimensions
            y_min = stat_data.stat(i_stat).min_elevation;
            y_max = 90;
            x_min = stat_data.stat(1).sat(1).epoch(1).jd;
            x_max = stat_data.stat(1).sat(1).epoch(stat_data.number_of_epochs).jd;

            % Initialize Subplot
            sched_handles.el_plot.stat(i_stat).h_subplot = subplot(stat_data.number_of_stations + 1 - stat_data.number_of_stations_quasar, 1, i_stat); 

            exceed_max_axis1_rate = zeros(stat_data.number_of_sat, stat_data.number_of_epochs);
            exceed_max_axis2_rate = zeros(stat_data.number_of_sat, stat_data.number_of_epochs);
            exceed_min_sun_dist = zeros(stat_data.number_of_sat, stat_data.number_of_epochs);
            exceed_axis_limits = zeros(stat_data.number_of_sat, stat_data.number_of_epochs);

            for i_sat = 1 : stat_data.number_of_sat   % Satellites

                el = zeros(1, stat_data.number_of_epochs);

                % Get Satellite names for the legend
                if (i_stat == 1)
                   sat_names{i_sat} = [num2str(i_sat), ' - ', stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_header_line];
                end

                [number_of_obs_periods, col] = size(obs_data.sat(i_sat).obs_times);

                % Get additional data for plotting:
                for i_epoch = 1 : stat_data.number_of_epochs  % Propagation epochs

                    t(i_epoch) = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).jd;

                    % Elevation
                    if (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).above_min_elevation)
                        el(i_epoch) = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).el;
                    end;  

                    % Axis 2 Rate
                    if (PARA.PLOT_EL_RATE)
                        if (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_max_axis2_rate)
                            exceed_max_axis2_rate(i_sat, i_epoch) = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).el;
                        else
                            exceed_max_axis2_rate(i_sat, i_epoch) = NaN;
                        end
                    end

                    % Axis 1 Rate
                    if (PARA.PLOT_AZ_RATE)
                        if (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_max_axis1_rate)
                            exceed_max_axis1_rate(i_sat, i_epoch) = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).el;
                        else
                            exceed_max_axis1_rate(i_sat, i_epoch) = NaN;
                        end
                    end

                    % Sun dist.
                    if (PARA.PLOT_SUN_DIST)
                        if (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_min_sun_dist)
                            exceed_min_sun_dist(i_sat, i_epoch) = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).el;
                        else
                            exceed_min_sun_dist(i_sat, i_epoch) = NaN;
                        end
                    end
                    
                    % Axis limits
                    if (PARA.PLOT_AXIS_LIMITS)
                        if (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).exceed_axis_limits)
                            exceed_axis_limits(i_sat, i_epoch) = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).el;
                        else
                            exceed_axis_limits(i_sat, i_epoch) = NaN;
                        end
                    end
                    
                    

                    % Possible observation periods
                    for i_obs = 1 : number_of_obs_periods
                        t_start = obs_data.sat(i_sat).obs_times(i_obs, 1);
                        t_stop = obs_data.sat(i_sat).obs_times(i_obs, 2);

                        if (    (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).jd >= t_start)    && ...
                                (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).jd <= t_stop)             )
                            obs_periods_el(i_epoch, i_sat, i_stat) = stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).el;
                            break;
                        else
                            obs_periods_el(i_epoch, i_sat, i_stat) = NaN;
                        end
                    end


                end;    % i_epoch = 1 : number_of_epochs

                % Plot Elevation data
                sched_handles.el_plot.stat(i_stat).sat(i_sat).h_sat_el = plot(t, el);

                hold all 

                % Set dimensions
                set(gca, 'XLim', [x_min x_max]);
                set(gca, 'YLim', [y_min y_max]);

                switch(x_tick_unit_str)
                    case 'hours'
                        i_xtick = 1;
                        x(i_xtick) = x_min;
                        x_label_new{i_xtick} =  num2str((x(i_xtick) - x(1)) * 24);

                        while(x(i_xtick) <= (x_max - (1/24 * int_unit)))
                            i_xtick = i_xtick + 1;
                            x(i_xtick) = x(i_xtick - 1) + 1/24 * int_unit;
                            x_label_new{i_xtick} = num2str((x(i_xtick) - x(1)) * 24);
                        end

                    case 'minutes'
                        i_xtick = 1;
                        x(i_xtick) = x_min;
                        x_label_new{i_xtick} =  num2str((x(i_xtick) - x(1)) * 24);

                        while (x(i_xtick) <= (x_max - (1/(24*60) * int_unit)))
                            i_xtick = i_xtick + 1;
                            x(i_xtick) = x(i_xtick - 1) + (1/(24*60 )* int_unit);
                            x_label_new{i_xtick} = num2str((x(i_xtick) - x(1)) * 24*60);
                        end
                end

                set(gca,'XTick', x)
                set(gca,'XTickLabel', {})

                hold all;

                %Get current color
                current_color(i_sat, :) = get(sched_handles.el_plot.stat(i_stat).sat(i_sat).h_sat_el, 'Color');


                hold all;

            end;    % for i_sat = 1 : number_of_sat

            % Add Title to Plot
            title(['Station: ', stat_data.stat(i_stat).name]);
            ylabel('El [°]');
            % xlabel('Duration [hours]');



            % Add Legend
            if (i_stat == 1)
                leg_shift_right = 0.01;
                sched_handles.el_plot.h_legend = legend(sat_names);
                pos_leg = get(sched_handles.el_plot.h_legend, 'Position');
                pos_sub = get(sched_handles.el_plot.stat(i_stat).h_subplot, 'Position');
                set(sched_handles.el_plot.h_legend, 'Position', [pos_sub(1) + pos_sub(3) + leg_shift_right, pos_leg(2), pos_leg(3), pos_leg(4)]); 
            end

            % Plot additional data:
            for i_sat = 1 : stat_data.number_of_sat
                % El. Rate
                if (PARA.PLOT_EL_RATE)
                    sched_handles.el_plot.stat(i_stat).sat(i_sat).h_exceed_max_axis2_rate = plot(t, exceed_max_axis2_rate(i_sat, :), 'Color', current_color(i_sat,:));
                    set(sched_handles.el_plot.stat(i_stat).sat(i_sat).h_exceed_max_axis2_rate, 'Marker', 'd');
                end

                % Az. Rate
                if (PARA.PLOT_AZ_RATE)
                    sched_handles.el_plot.stat(i_stat).sat(i_sat).h_exceed_max_axis1_rate = plot(t, exceed_max_axis1_rate(i_sat, :), 'Color', current_color(i_sat,:));
                    set(sched_handles.el_plot.stat(i_stat).sat(i_sat).h_exceed_max_axis1_rate, 'Marker', 's');
                end

                % Sun dist.
                if (PARA.PLOT_SUN_DIST)
                    sched_handles.el_plot.stat(i_stat).sat(i_sat).h_exceed_min_sun_dist = plot(t, exceed_min_sun_dist(i_sat, :), 'Color', current_color(i_sat,:));
                    set(sched_handles.el_plot.stat(i_stat).sat(i_sat).h_exceed_min_sun_dist, 'Marker', '*');
                end
                
                % Axis limits
                if (PARA.PLOT_AXIS_LIMITS)
                    sched_handles.el_plot.stat(i_stat).sat(i_sat).h_exceed_axis_limits = plot(t, exceed_axis_limits(i_sat, :), 'Color', current_color(i_sat,:));
                    set(sched_handles.el_plot.stat(i_stat).sat(i_sat).h_exceed_axis_limits, 'Marker', 'o');
                end
                
            end;    % for i_sat = 1 : number_of_sat

        end

    end;    % for i_stat = 1 : stat_data.number_of_stations
    
    
    % Save "station id list":
    sched_handles.el_plot.station_id_list = station_id_list;
   
    
    
    %% --- Plot Rise / Peak / Set ---
    
    for i_stat = 1 : stat_data.number_of_stations    % Stations
        
        if ~strcmp(stat_data.stat(i_stat).obs_type, 'quasar_only') && ~strcmp(stat_data.stat(i_stat).obs_type, 'twin_quasar')
        
            % Set Current axes
            axes(sched_handles.el_plot.stat(i_stat).h_subplot);

            for i_sat = 1 : stat_data.number_of_sat    % Satellites

                sched_handles.el_plot.stat(i_stat).sat(i_sat).number_of_overpasses = stat_data.stat(i_stat).sat(i_sat).number_of_overpasses;

                % Plot Overpass data
                for i_overp = 1 : stat_data.stat(i_stat).sat(i_sat).number_of_overpasses

                    % Rise
                    if ~(isempty(stat_data.stat(i_stat).sat(i_sat).overpass(i_overp).rise.year))
                        sched_handles.el_plot.stat(i_stat).sat(i_sat).overpass(i_overp).h_rise = scatter(stat_data.stat(i_stat).sat(i_sat).overpass(i_overp).rise.jd, stat_data.stat(i_stat).sat(i_sat).overpass(i_overp).rise.el, 'CData', current_color(i_sat,:));
                    end

                    % Peak
                    for i_peak = 1 : stat_data.stat(i_stat).sat(i_sat).overpass(i_overp).number_of_peaks
                        sched_handles.el_plot.stat(i_stat).sat(i_sat).overpass(i_overp).h_peak = scatter(stat_data.stat(i_stat).sat(i_sat).overpass(i_overp).peak(i_peak).jd, stat_data.stat(i_stat).sat(i_sat).overpass(i_overp).peak(i_peak).el, 'CData', current_color(i_sat,:));
                    end

                    % Set
                    if ~(isempty(stat_data.stat(i_stat).sat(i_sat).overpass(i_overp).set.year))
                        sched_handles.el_plot.stat(i_stat).sat(i_sat).overpass(i_overp).h_set = scatter(stat_data.stat(i_stat).sat(i_sat).overpass(i_overp).set.jd, stat_data.stat(i_stat).sat(i_sat).overpass(i_overp).set.el, 'CData', current_color(i_sat,:));
                    end

                end

            end;    % for i_sat = 1 : number_of_sat
            
        end

    end;    % for i_stat = 1 : stat_data.number_of_stations
    
    
    
    
    
    %% --- Plot possible observation periods ---
    
    % Initialize Subplot
    %sched_handles.el_plot.stat(stat_data.number_of_stations + 1).h_subplot = subplot(stat_data.number_of_stations + 1, 1, stat_data.number_of_stations + 1); 
    sched_handles.el_plot.h_available_obs_time_subplot = subplot(stat_data.number_of_stations + 1 - stat_data.number_of_stations_quasar, 1, stat_data.number_of_stations + 1 - stat_data.number_of_stations_quasar); 
    
    % Plot Elevation of observation periods
    for i_sat = 1 : stat_data.number_of_sat   % Satellites
        for i_stat = 1 : (stat_data.number_of_stations - stat_data.number_of_stations_quasar) % Stations 
            sched_handles.el_plot.stat(i_stat).sat(i_sat).h_sat_el_obs_times = plot(t, obs_periods_el(:,i_sat, i_stat),  'Color', current_color(i_sat,:));
            hold all;
        end
    end
    
    % Set dimensions
    set(gca, 'YLim', [y_min y_max]);
    set(gca, 'XLim', [x_min x_max]);
    
    % Set x-Tick and x-Tick Labels [hours]
    switch(x_tick_unit_str)
        case 'hours'
            i_xtick = 1;
            x(i_xtick) = x_min;
            x_label_new{i_xtick} =  num2str((x(i_xtick) - x(1)) * 24);

            while(x(i_xtick) <= (x_max - (1/24 * int_unit)))
                i_xtick = i_xtick + 1;
                x(i_xtick) = x(i_xtick - 1) + (1/24) *  int_unit;
                x_label_new{i_xtick} = num2str((x(i_xtick) - x(1)) * 24);
            end
            xlabel(sprintf('[hours] since %s UTC (total duration: %2.2f hours)', sess_start_time_str,  PARA.duration));
            
        case 'minutes'
            i_xtick = 1;
            x(i_xtick) = x_min;
            x_label_new{i_xtick} =  num2str((x(i_xtick) - x(1)) * 24);

            while (x(i_xtick) <= (x_max - (1/(24*60) * int_unit)))
                i_xtick = i_xtick + 1;
                x(i_xtick) = x(i_xtick - 1) + (1/(24*60)*int_unit);
                x_label_new{i_xtick} = num2str((x(i_xtick) - x(1)) * 24*60);
            end
            unit_str = 'minutes';
            xlabel(sprintf('[minutes] since %s UTC (total duration: %2.2f min)', sess_start_time_str,  PARA.duration*60));
    end

    set(gca,'XTick', x)
    set(gca,'XTickLabel', x_label_new)
    
    % Add Title an Axes Labels to Plot
    title('Potential observation periods from station network');
    ylabel('El [°]');
    
    % Plot vertical lines and print time tags:
    
    y_text = (y_max - y_min) / 2;
    
    for i_sat = 1 : stat_data.number_of_sat   % Satellites
         
        [number_of_obs_periods, col] = size(obs_data.sat(i_sat).obs_times);
        
        sched_handles.el_plot.sat(i_sat).number_of_obs_periods = number_of_obs_periods;
        
        for i_obs = 1 : number_of_obs_periods
            
            % Vertical Lines:
            sched_handles.el_plot.sat(i_sat).obs_period(i_obs).h_line_start = line([obs_data.sat(i_sat).obs_times(i_obs, 1), obs_data.sat(i_sat).obs_times(i_obs, 1)], [y_min, y_max], 'color', current_color(i_sat,:), 'LineStyle', ':');
            sched_handles.el_plot.sat(i_sat).obs_period(i_obs).h_line_end = line([obs_data.sat(i_sat).obs_times(i_obs, 2), obs_data.sat(i_sat).obs_times(i_obs, 2)], [y_min, y_max], 'color', current_color(i_sat,:), 'LineStyle', ':');
            
            
            % Time tags:
            x_text = obs_data.sat(i_sat).obs_times(i_obs, 1);
            sched_handles.el_plot.sat(i_sat).obs_period(i_obs).h_text_start = text(x_text, y_text, [' t', num2str(obs_data.sat(i_sat).obs_times(i_obs, 3))], 'color', current_color(i_sat,:));
            x_text = obs_data.sat(i_sat).obs_times(i_obs, 2);
            sched_handles.el_plot.sat(i_sat).obs_period(i_obs).h_text_end = text(x_text, y_text, [' t', num2str(obs_data.sat(i_sat).obs_times(i_obs, 4))], 'color', current_color(i_sat,:));
            
        end
    end
    
    
    % ##### Text: End of last scan & Current epoch (marker) #####
    x_lim = get(gca, 'XLim');
    y_lim = get(gca, 'YLim');
    dx = x_lim(2) - x_lim(1);
    % dy = y_lim(2) - y_lim(1);
    h_text_1 = text(x_lim(2) + dx/35, y_lim(2), 'end of last scan:', 'color', sched_handles.color_markers_end_of_last_scan);
    sched_handles.el_plot.h_text_av_obs_time_end_of_last_scan = text(x_lim(2) + dx/30, y_lim(2) - 15, reset_time_str);
    
    h_text_2 = text(x_lim(2) + dx/35, y_lim(2) - 40, 'current epoch:', 'color', sched_handles.color_current_epoch_marker);
    sched_handles.el_plot.h_text_av_obs_time_current_epoch = text(x_lim(2) + dx/30, y_lim(2) - 55, reset_time_str);


return;



