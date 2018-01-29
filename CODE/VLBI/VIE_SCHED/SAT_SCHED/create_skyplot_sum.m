% -------------------------------------------------------------------------
%
%                              create_skyplot_sum.m
%
%  Description:
%   This function creates sky plots showing all schedules scans (quasars and satellites) for a defined time window (optional).
%   If no time window is defined, the whole schdule defined in the "sched_data" struct is visualized.
%   The observation time is color-coded (for sat. and quasar scans).
%   Plotted items (is scheduled for the defined time window):
%   - Satellite tracks (optional)
%   - Positions of quasar scans (at scans start)
%   - position of satellite scans (scan start to scan end)
%   - Sky-coverage segments (optional) => TO DO!
%   - Cut-off elevation and horizon mask (if available)
%   
%  
%   Author: 
%       Andreas Hellerschmied, 2016-04-28
%           
%
%   inputs        : (sched_data, PARA, flag_print_sattracks, flag_print_skycov, flag_save_fig, flag_save_pdf, t_min_jd, t_max_jd);
%   - sched_data            : schdeduling structure
%   - PARA                  : VieVS parameter vector 
%   - flag_print_sattracks  : Plot satellite tracks (optional, , default = true)
%   - flag_print_skycov     : Plot sky coverage segments (optional, default = false)
%   - flag_save_fig         : Save figure file to current /DATA/SCHED/ subdir.
%   - flag_save_pdf         : Save figure as PDF to current /DATA/SCHED/ subdir.
%   - t_min_jd              : Lower limit of time window (optional)
%   - t_max_jd              : Upper limit of time window (optional)
%     
%
%   outputs       :
%   - Screen Output (Command Window)
%   - error_code        : Error Code (0 = no erros occured)
%   - error_msg         : Error Message (empty, if no errors occured)
%   - sched_handles     : Graphic handles structure for scheduling
%
%
%   coupling      :
%   - skyplot
%   - get_station_statistics
%   - tle_2_topo
%   - save_pdf_a4_landscape
%   
%
%   references    :
%
%   changes       :
%   - 2016-05-09, A. Hellerschmied: Save sky_plot_sum figures as .fig and/or .pdf files
%   - 2016-05-11, A. Hellerschmied: Bug-fix
%   - 2016-05-12, A. Hellerschmied: Bug-fix
%   - 2016-07-13, A. Hellerschmied: Print PDF command adopted for MATLAB 2016a
%   - 2016-11-27, A. Hellerschmied: save_pdf_a4_landscape.m used to print PDF files
%   - 2016-12-22, A. Hellerschmied: PARA.INIT_PROP_INTERVAL in [sec] instead of [min]
%
%-------------------------------------------------------------------------

function [error_code, error_msg, skyplot_handles] = create_skyplot_sum(sched_data, PARA, varargin)

    % ##### Init #####
    error_code = 0;
    error_msg = '';
    skyplot_handles = [];
    flag_print_skycov       = false;
    flag_print_sattracks    = true;
    flag_check_time_window  = false;
    flag_save_fig           = false;
    flag_save_pdf           = false;
    sat_names               = '';
    t_min_jd                = sched_data.t_nominal_start_jd;
    t_max_jd                = sched_data.t_nominal_end_jd;
    
    
    % ##### Options #####
    flag_plot_h_mask                    = 1;
    number_of_az_intervals_cut_off_el   = 100;                          % Number of azimuth intervals fot plotting the cut-off elevation marker (line plot)
%     t_prop_int_min                      = 1;                          % Propagation interval for plotting satellite tracks
    t_prop_int_min                      = PARA.INIT_PROP_INTERVAL/60;   % Propagation interval for plotting satellite tracks
    line_style                          = '-';                          % Line type for satellite tracks
    marker_int_min                      = PARA.SKYPLOT_MARKER_INT;      % Marker interval [min]
    flag_plot_markers                   = true;
    flag_plot_horizon_mask              = true;                         % Plot cut-off elevation and horizon mask (if available)
    
    number_of_stat                      = length(sched_data.stat);
    
    
     % ##### Treat optinal input arguments #####
    switch(nargin)
        case 2
        case 3
            flag_print_sattracks    = varargin{1};
        case 4
            flag_print_sattracks    = varargin{1};
            flag_print_skycov       = varargin{2};
            
        case 5
            flag_print_sattracks    = varargin{1};
            flag_print_skycov       = varargin{2};
            flag_save_fig           = varargin{3};
        case 6
            flag_print_sattracks    = varargin{1};
            flag_print_skycov       = varargin{2};
            flag_save_fig           = varargin{3};
            flag_save_pdf           = varargin{4};
        case 8
            flag_check_time_window  = true;
            flag_print_sattracks    = varargin{1};
            flag_print_skycov       = varargin{2};
            flag_save_fig           = varargin{3};
            flag_save_pdf           = varargin{4};
            t_min_jd                = varargin{5};
            t_max_jd                = varargin{6};
        otherwise
            error_code = 1;
            error_msg = 'Invalid number of input arguments';
            return;
    end
    
    % ##### Get station statistics #####
    if flag_check_time_window
        [sched_data, error_code, error_msg] = get_station_statistics(sched_data, t_min_jd, t_max_jd);
        if error_code > 0
            error_msg = ['get_station_statistics: ', error_msg];
            return; 
        end
    else
        [sched_data, error_code, error_msg] = get_station_statistics(sched_data);
        if error_code > 0
            error_msg = ['get_station_statistics: ', error_msg];
            return; 
        end
    end
    
    
    % ##### Prepare "stat_data_tmp" #####
    % Input for: tle_2_topo.m
    for i_stat = 1 : number_of_stat
        stat_data_tmp.stat(i_stat).sat                          = sched_data.sat;
        stat_data_tmp.stat(i_stat).location.ellipsoid.long      = sched_data.stat(i_stat).long;
        stat_data_tmp.stat(i_stat).location.ellipsoid.lat       = sched_data.stat(i_stat).lat;
        stat_data_tmp.stat(i_stat).location.ellipsoid.altitude  = sched_data.stat(i_stat).altitude;
    end
    

    % ##### Create one sky-plot per station #####
    % loop over all stations
    
    if flag_print_sattracks
        % #### Get number of epochs: ####
        exp_duration_min = round((t_max_jd - t_min_jd) * (24*60)); % [min]
        number_of_epochs = exp_duration_min / t_prop_int_min;
        % Check "number_of_epochs":
        if mod(number_of_epochs, 1) ~= 0
            error_code = 1;
            error_msg = 'Number of epochs is not an integer number!';
            return;
        elseif number_of_epochs == 0;
            error_code = 2;
            error_msg = 'Number of time epochs for propagation = 0!';
            return;
        end
        % Print Warning, if the chosen marker interval is not valid!
        if mod(marker_int_min / t_prop_int_min, 1) ~= 0 % plot markers
            fprintf(1, '\n');
            fprintf(1, '--------------------------------------------------------------------------------------------\n');
            fprintf(1, '  Note: Time markers in the skyplot will not be printed,\n');
            fprintf(1, '        because the marker interval is not a multiple of the propagation interval!\n');
            fprintf(1, '--------------------------------------------------------------------------------------------\n');
            flag_plot_markers = false;
        end
        
    end
    
    % #### Loop over all stations ####
    for i_stat = 1 : number_of_stat
        
        % Create figure for each station:
        skyplot_handles.sky_plots_sum(i_stat).h_fig = figure;
        
        % ##### Create sky plot #####
        if flag_print_sattracks % Plot satellite tracks
            % #### Prepare Az/El values of satellite positions for plotting the satellite tracks ####
            number_of_sat = length(sched_data.stat(i_stat).statistics.sat_name_list);
            % ### preallocations ###
            az_data = zeros(number_of_sat, number_of_epochs);
            el_data = zeros(number_of_sat, number_of_epochs);
            prn = zeros(1, number_of_sat);
            for i_sat = 1 : number_of_sat
                prn(1, i_sat) = i_sat;
                sat_names{i_sat} = [num2str(i_sat), ' - ', deblank(sched_data.stat(i_stat).statistics.sat_name_list{i_sat})];
                sat_id = sched_data.stat(i_stat).statistics.sat_id_vector(i_sat);
                for i_epoch = 1 : (number_of_epochs + 1)
                    t_epoch_jd = t_min_jd + (i_epoch - 1) * t_prop_int_min / (24*60);
                    % ##### Check satellite visibility for the current epoch #####
                    number_of_obs_periods = size(sched_data.stat(i_stat).sat(sat_id).obs_times, 1);
                    flag_observable = false;
                    for i_obs_per = 1 : number_of_obs_periods
                        if (t_epoch_jd >=  sched_data.stat(i_stat).sat(sat_id).obs_times(i_obs_per,1)) && (t_epoch_jd <=  sched_data.stat(i_stat).sat(sat_id).obs_times(i_obs_per,2))
                           flag_observable = true;
                           break;
                        end
                    end % for i_obs_per = 1 : number_of_obs_periods
                    if flag_observable % Is observable
                        % ### Calc. Az and El for t_epoch_jd: ###
                        [az, el, range, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate, error_code, error_msg] = tle_2_topo(stat_data_tmp, i_stat, PARA, sat_id, t_epoch_jd);
                        if error_code > 0
                            error_msg = ['tle_2_topo: ', error_msg];
                            return; 
                        end
                    else % Not observable
                        az = NaN;
                        el = NaN;
                    end
                    az_data(i_sat, i_epoch) = az;
                    el_data(i_sat, i_epoch) = el;
                end % i_epoch = 1 : number_of_sat
            end % i_sat = 1 : number_of_sat
            
            % ##### Create sky-plot #####
            if flag_plot_markers
                marker_int = marker_int_min / t_prop_int_min;
                [h_plot, h_axis, plotted_sats, h_marker, h_sat_num_txt, h_sat_marker] = skyplot(az_data, el_data, prn, marker_int, line_style);    
            else % do not plot time markers
                [h_plot, h_axis, plotted_sats, h_marker, h_sat_num_txt, h_sat_marker] = skyplot(az_data, el_data, prn); 
            end
        else % Do not plot satellite tracks
            prn     = NaN;
            az_data = NaN;
            el_data = NaN;
            % ##### Create sky-plot (empty) #####
            [h_plot, h_axis, plotted_sats, h_marker, h_sat_num_txt, h_sat_marker] = skyplot(az_data, el_data, prn); 
        end
 
        % ##### Save handles to structure: #####
        skyplot_handles.sky_plots_sum(i_stat).h_axis = h_axis;
        for i_sat = 1 : number_of_sat
            skyplot_handles.sky_plots_sum(i_stat).sat(i_sat).h_plot = h_plot(i_sat);
            skyplot_handles.sky_plots_sum(i_stat).sat(i_sat).flag_plotted = logical(plotted_sats(i_sat));
        end
        
        % ##### Plot scans #####
        
        % Get number of scans
        number_of_scans = sched_data.stat(i_stat).statistics.num_of_all_scans;
        
        % Init:
        sat_count = 0;
        q_count = 0;
        
        % Preallocate AzEl values
        az_sat_rad      = zeros(sched_data.stat(i_stat).statistics.num_of_sat_scans, 1);
        el_sat_rad      = zeros(sched_data.stat(i_stat).statistics.num_of_sat_scans, 1);
        scan_num_sat    = zeros(sched_data.stat(i_stat).statistics.num_of_sat_scans, 1);
        az_q_rad        = zeros(sched_data.stat(i_stat).statistics.num_of_quasar_scans, 1);
        el_q_rad        = zeros(sched_data.stat(i_stat).statistics.num_of_quasar_scans, 1);
        scan_num_q      = zeros(sched_data.stat(i_stat).statistics.num_of_quasar_scans, 1);
        sat_mask        = zeros(number_of_scans, 1);
        q_mask          = zeros(number_of_scans, 1);

        % Loop over scans
        for i_scan = 1 : sched_data.number_of_scans
            
            % ##### Check, if the scan is within t_min_jd and t_max_jd #####
            if ((sched_data.scan(i_scan).t_start_jd < t_min_jd) || (sched_data.scan(i_scan).t_end_jd > t_max_jd) && flag_check_time_window)
                continue; 
            end
            
            % #### Check if station (i_stat) takes part in this scan
            if logical(sum([sched_data.scan(i_scan).stat.stat_id] == i_stat))
                
                % Scan epoch (start time)
                t_epoch_jd = sched_data.scan(i_scan).t_start_jd;
                
                % Distinguish between quasar and satellite scans
                switch(sched_data.scan(i_scan).obs_type)
                    case 'quasar'
                        q_count = q_count + 1;
                        q_mask(i_scan) = true;
                        % Station position (conversion deg => rad):
                        lon = sched_data.stat(i_stat).long * pi / 180;
                        lat = sched_data.stat(i_stat).lat * pi / 180;

                        [az, el, ha, dc] = zazel_s(t_epoch_jd - 2400000.5, lon, lat, sched_data.quasar(sched_data.scan(i_scan).quasar_id).ra_rad, sched_data.quasar(sched_data.scan(i_scan).quasar_id).de_rad);
                        az_q_rad(q_count)   = az;
                        el_q_rad(q_count)   = el;
                        scan_num_q(q_count) = i_scan;
                    case 'sat'
                        sat_count = sat_count + 1;
                        sat_mask(i_scan) = true;
                        [az, el, range, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate, error_code, error_msg] = tle_2_topo(stat_data_tmp, i_stat, PARA, sched_data.scan(i_scan).sat_id, t_epoch_jd);
                        if error_code > 0
                            error_msg = ['tle_2_topo: ', error_msg];
                            return; 
                        end
                        az_sat_rad(sat_count)   = az * pi/180;
                        el_sat_rad(sat_count)   = el * pi/180;
                        scan_num_sat(sat_count) = i_scan;
                end
            
            end % if logical(sum([sched_data.scan(i_scan).stat.stat_id] == i_stat))
            
        end % for i_scan = 1 : number_of_scans
        

        % ### Quasar scans ###
        if sched_data.stat(i_stat).statistics.num_of_quasar_scans > 0
            % Transform elevation angle to a distance to the center of the plot:
            elSpherical_1 = 90*cos(el_q_rad); 
            % Transform data to Cartesian coordinates:
            yy_q = elSpherical_1 .* cos(az_q_rad); 
            xx_q = elSpherical_1 .* sin(az_q_rad);
            % Plot
            skyplot_handles.sky_plots_sum(i_stat).h_quasar_scans = plot(xx_q, yy_q, 'x', 'Color', 'red');
            skyplot_handles.sky_plots_sum(i_stat).h_quasar_scans.MarkerSize     = 10;
            skyplot_handles.sky_plots_sum(i_stat).h_quasar_scans.LineWidth      = 1.5;

            % Scan numbers:
            xx_offset = 0;
            for i_scan_num = 1 : length(scan_num_q)
                temp_str = sprintf('%3d', scan_num_q(i_scan_num));
                skyplot_handles.sky_plots_sum(i_stat).h_scan_num_g(i_scan_num)  = text(xx_q(i_scan_num) + xx_offset, yy_q(i_scan_num), temp_str, 'color', 'red');
            end
        end
        
        
        % ### Satellite scans ###
        if sched_data.stat(i_stat).statistics.num_of_sat_scans > 0
            % Transform elevation angle to a distance to the center of the plot:
            elSpherical_1 = 90*cos(el_sat_rad); 
            % Transform data to Cartesian coordinates:
            yy_sat = elSpherical_1 .* cos(az_sat_rad); 
            xx_sat = elSpherical_1 .* sin(az_sat_rad);
            % Plot
            skyplot_handles.sky_plots_sum(i_stat).h_sat_scans = plot(xx_sat, yy_sat, 'o', 'Color', 'blue');
            skyplot_handles.sky_plots_sum(i_stat).h_sat_scans.MarkerSize    = 10;
            skyplot_handles.sky_plots_sum(i_stat).h_sat_scans.LineWidth     = 1.5;

            % Scan numbers:
            xx_offset = 0;
            for i_scan_num = 1 : length(scan_num_sat)
                temp_str = sprintf('%3d', scan_num_sat(i_scan_num));
                skyplot_handles.sky_plots_sum(i_stat).h_scan_num_sat(i_scan_num)  = text(xx_sat(i_scan_num) + xx_offset, yy_sat(i_scan_num), temp_str, 'color', 'blue');
            end
        end
        
%         % ##### Get observation type of the antenna #####
%         switch(sched_data.stat(i_stat).obs_type)
%             case {'sat_only'}
%                 obs_type_str = 'satellite';
%             case {'quasar_only'}
%                 obs_type_str = 'quasar';
%             case {'sat_and_quasar'}
%                 obs_type_str = 'satellite & quasar';
%             case {'twin_sat'}
%                 obs_type_str = ['satellite; twin with: ', stat_data.stat(i_stat).twin_partner];
%             case {'twin_quasar'}
%                 obs_type_str = ['quasar; twin with: ', stat_data.stat(i_stat).twin_partner];
%         end

 
        % ##### Title #####
        % skyplot_handles.sky_plots_sum(i_stat).h_title = title([sched_data.stat(i_stat).name, ' (obs. type: ', obs_type_str,', ', jd2datestr(t_min_jd), ' - ', jd2datestr(t_max_jd), ')']);
        skyplot_handles.sky_plots_sum(i_stat).h_title = title([sched_data.exper_name, ': ', sched_data.stat(i_stat).name, ' (', jd2datestr(t_min_jd), ' - ', jd2datestr(t_max_jd), ')']);
        skyplot_handles.sky_plots_sum(i_stat).h_fig.Name = ['Skyplot; ', sched_data.exper_name, ': ', sched_data.stat(i_stat).name, ' (', jd2datestr(t_min_jd), ' - ', jd2datestr(t_max_jd), ')'];

        % ##### Legend #####
        if flag_print_sattracks
            skyplot_handles.sky_plots_sum(i_stat).h_legend = legend(sat_names);

            % Get positions of legend
            pos_leg = get(skyplot_handles.sky_plots_sum(i_stat).h_legend, 'Position');      
            pos_axis = get(h_axis, 'Position');

            % Set new positions:
            set(h_axis, 'Position', [0.00 pos_axis(2) pos_axis(3) pos_axis(4)])
            set(skyplot_handles.sky_plots_sum(i_stat).h_legend, 'Position', [1-pos_leg(3)-0.05   pos_leg(2)    pos_leg(3)    pos_leg(4)])
        end  
        
        % ##### Plot horizontal mask and/or the cut-off elevation #####
        if flag_plot_horizon_mask
            % Get cut-off elevation of the station:
            cut_off_el_rad = PARA.MIN_CUTEL;
            az_interval = (2*pi) / number_of_az_intervals_cut_off_el;
            az_mask_rad = (0 : az_interval : (2*pi));

            if flag_plot_h_mask && (sched_data.stat(i_stat).horizontal_mask_num > 0)

                % Get h-mask data:
                hmasknum = sched_data.stat(i_stat).horizontal_mask_num;
                hmask    = sched_data.stat(i_stat).horizontal_mask;

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
            skyplot_handles.sky_plots_sum(i_stat).h_cut_off_el_marker = line(xx_mask, yy_mask, 'Linestyle', '-.', 'Color', 'red', 'LineWidth', 0.1);
            
        end % if flag_plot_horizon_mask

        
        % ##### Save figure #####
        if (flag_save_pdf || flag_save_fig)
            
            start_time_str = jd2datestr(t_min_jd);
            start_time_str(strfind(start_time_str,':'))='';
            start_time_str = start_time_str(12:15);
            end_time_str = jd2datestr(t_max_jd);
            end_time_str(strfind(end_time_str,':'))='';
            end_time_str = end_time_str(12:15);
            
            % Matlab figure
            if flag_save_fig
                savefig(skyplot_handles.sky_plots_sum(i_stat).h_fig, ['../DATA/SCHED/', PARA.OUTPUT_SUBDIR, '/sky_', sched_data.exper_name, '_', sched_data.stat(i_stat).label, '_', start_time_str, '-', end_time_str]);
                fprintf(1,'   ...Skyplot saved as: %s\n ', ['../DATA/SCHED/', PARA.OUTPUT_SUBDIR, '/sky_', sched_data.exper_name, '_', sched_data.stat(i_stat).label, '_', start_time_str, '-', end_time_str, '.fig']);
            end
            % PDF
            if flag_save_pdf
                save_pdf_a4_landscape(skyplot_handles.sky_plots_sum(i_stat).h_fig, ['../DATA/SCHED/', PARA.OUTPUT_SUBDIR, '/'], ['sky_', sched_data.exper_name, '_', sched_data.stat(i_stat).label, '_', start_time_str, '-', end_time_str])
            end

        end % if (flag_save_pdf || flag_save_fig)
        
    end % for i_stat = 1 : number_of_stat
        
return;



