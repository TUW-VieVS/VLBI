% #########################################################################
% #     update_sky_plots
% #########################################################################
%
% DESCRIPTION
%   This function is used to update and modify the content of the
%   "sky plots".
%
% CREATED  
%   2015-03-11     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
%   - tle_2_topo
%   - zazel_s                       : Conversion Ra/Dec => Az/El (Simplified approach for scheduling)
%   - preselect_observable_quasars  : Pre-selection of observable quasars
%   - invjday
%
% INPUT
%  
%   - sched_handles             : Graphic handles structure for scheduling
%   - jd_current_epoch_marker   : Currently treated epoch [JD], for setting the epoch marker
%   - stat_data                 : station data structure
%   - station_id_list           : List of IDs of the observing stations (referring to the required entry of "stat_data.stat" and "obs_data.stat")
%   - source                    : Source structure
%   - flag_plot_sources         : Flag, to choose, if sources should be plotted (1=yes, 0 =no)
%   - PARA                      : Global scheduling parameter structure
%   - obs_data                  : Observation data structure (preallocated in "find_observation_times.m")
%
%
% OUTPUT
%   - sched_handles       : Graphic handles structure for scheduling (updated)
%   - error_code          : Error Code (0 = no erros occured)
%   - error_msg           : Error Message (empty, if no errors occured)
%
% CHANGES:
%

function [sched_handles, error_code, error_msg] = update_sky_plots(sched_handles, stat_data, station_id_list, source, PARA, obs_data, flag_plot_sources, jd_current_epoch_marker, flag_observable_quasars_list)

    % ##### Const. #####
    rad2deg = 180/pi;
    deg2rad = pi/180;

    % ##### Init #####
    error_code = 0;
    error_msg = '';
%     pointer_az_rad = 0;
%     
%     % ##### Options #####
%     flag_draw_az_pointer = true;
    

    % ##### Check, if there are sky plots #####
    if isempty(sched_handles.sky_plots)
        return; 
    end
    
    
    % ##### Pre-select quasars, which are observable from the current network #####
    if flag_plot_sources
        if isempty(flag_observable_quasars_list)
            [source, flag_observable_quasars_list, error_code, error_msg] = preselect_observable_quasars(stat_data, station_id_list, source, jd_current_epoch_marker);
            if error_code > 0
                error_msg = ['preselect_observable_quasars: ', error_msg];
                return; 
            end
        else
            source = source(logical(flag_observable_quasars_list));
        end
    end
    
    
    
    
    % Loop over all stations included in the station ID list:
    % => loop over all sky plots included in the current scan

    for i_stat = 1 : length(station_id_list)
        
        station_id = station_id_list(i_stat);
              
        % Check if current sky plot is available:
        if ~isempty(sched_handles.sky_plots(station_id).h_fig)
            
            
            % ##### Print current epoch (Text) #####
            [year, mon, day, hr, min, sec] = invjday(jd_current_epoch_marker); % Time conversion
            current_epoch_str = sprintf('%02.0f:%02.0f:%02.0f', hr, min, sec);
            set(sched_handles.sky_plots(station_id).h_current_epoch, 'String', current_epoch_str);
            
            
            % ##### Set current axes #####
            axes(sched_handles.sky_plots(station_id).h_axis);

            
            % ##### Draw current epoch markers #####
            
            for i_sat = 1 : size(sched_handles.sky_plots(station_id).sat,2)
                
                if sched_handles.sky_plots(station_id).sat(i_sat).flag_plotted % If the track of the current satellite was plotted...
                    
                    % Set current plot color:
                    temp_color = get(sched_handles.sky_plots(station_id).sat(i_sat).h_plot, 'color');

                    % Delete a previous time marker line, if it exists:
                    if isfield(sched_handles.sky_plots(station_id).sat(i_sat), 'h_current_epoch_marker')
                        if ~isempty(sched_handles.sky_plots(station_id).sat(i_sat).h_current_epoch_marker)
                           delete(sched_handles.sky_plots(station_id).sat(i_sat).h_current_epoch_marker); % Remove existing marker
                        end
                    end
                    
                    % ### calculate antenna pointing data (az, el) for the given epoch (jd_current_epoch_marker) [JD] ###
                    
                    [az, el, range, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate, error_code, error_msg] = tle_2_topo(stat_data, station_id, PARA, i_sat, jd_current_epoch_marker);
                    if error_code > 0
                        error_msg = ['tle_2_topo: ', error_msg];
                        return; 
                    end
                       
                    % Check, if the satellite is above the cut-off elevation:
                    if el > stat_data.stat(station_id).min_elevation
                    
                        % Transform elevation angle to a distance to the center of the plot:
                        elSpherical = 90*cos(el * pi/180); 

                        % Transform data to Cartesian coordinates:
                        yy = elSpherical .* cos(az * pi/180); 
                        xx = elSpherical .* sin(az * pi/180); 

                        % ### Plot marker ###
                        sched_handles.sky_plots(station_id).sat(i_sat).h_current_epoch_marker = plot(xx', yy', '*', 'color', temp_color, 'MarkerSize', 20); 
                    
                    else % Satellite is below the cut-off elevation
                        
                                            

                    end
                    
                end
                
            end      
        
            
            
            % ##### Delete a previous source markers, if they exist #####
            if isfield(sched_handles.sky_plots(station_id), 'quasars')
                if isfield(sched_handles.sky_plots(station_id).quasars, 'h_quasar_markers')
                    if ~isempty(sched_handles.sky_plots(station_id).quasars.h_quasar_markers)
                       delete(sched_handles.sky_plots(station_id).quasars.h_quasar_markers); % Remove existing markers
                    end
                end
            end
        

            
            
            % ##### Plot quasars #####
            if flag_plot_sources
                
                
                % Delete a previous source markers, if they exist:
                if isfield(sched_handles.sky_plots(station_id), 'quasars')
                    if isfield(sched_handles.sky_plots(station_id).quasars, 'h_quasar_markers')
                        if ~isempty(sched_handles.sky_plots(station_id).quasars.h_quasar_markers)
                           delete(sched_handles.sky_plots(station_id).quasars.h_quasar_markers); % Remove existing markers
                        end
                    end
                end
                
                % Get station data for calculation of view angles (az, el)
                stat_lon = stat_data.stat(station_id).location.ellipsoid.long * pi / 180;
                stat_lat = stat_data.stat(station_id).location.ellipsoid.lat * pi / 180;

                % Preallocation:
                src_az = zeros(size(source,2), 1);
                src_el = zeros(size(source,2), 1);

                % ### Convert source coordinates (Ra/Dec) to Az/El ###
                for i_src = 1 : size(source,2)
                    [src_az(i_src), src_el(i_src), ha, dc] = zazel_s(jd_current_epoch_marker - 2400000.5, stat_lon, stat_lat, source(i_src).ra, source(i_src).de);
                end

                % Transform elevation angle to a distance to the center of the plot:
                elSpherical = 90*cos(src_el); 

                % Transform data to Cartesian coordinates:
                yy = elSpherical .* cos(src_az); 
                xx = elSpherical .* sin(src_az); 

                % ##### Plot quasar markers #####
                sched_handles.sky_plots(station_id).quasars.h_quasar_markers = plot(xx', yy', '*', 'color', sched_handles.color_quasars, 'MarkerSize', 3); 
                
                % Save data da sched handels structure:
                sched_handles.sky_plots(station_id).quasars.xx = xx;
                sched_handles.sky_plots(station_id).quasars.yy = yy;
                sched_handles.sky_plots(station_id).quasars.flag_observable_list = flag_observable_quasars_list;

            end
        
        end
        
    end % station loop
    

    
    
    
return;