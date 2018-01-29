% #########################################################################
% #     update_sky_plots_highlight_quasars
% #########################################################################
%
% DESCRIPTION
%   This function handles the markers to highlight quasars, defined in the "quasar_id_list".
%
%       Plot option (plot_opt):
%           1 - Highlighted sources
%           2 - Delete highlighted sources (markers)
%
% CREATED  
%   2015-04-22     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
%   - zazel_s                       : Conversion Ra/Dec => Az/El (Simplified approach for scheduling)
%
%
% INPUT
%  
%   - sched_handles             : Graphic handles structure for scheduling
%   - t_epoch_jd                : Currently treated epoch [JD], for setting the marker
%   - stat_data                 : station data structure
%   - station_id_list           : List of IDs of the observing stations (referring to the required entry of "stat_data.stat" and "obs_data.stat")
%   - source                    : Source structure
%   - plot_opt                  : Plot options
%
%
% OUTPUT
%   - sched_handles       : Graphic handles structure for scheduling (updated)
%   - error_code          : Error Code (0 = no erros occured)
%   - error_msg           : Error Message (empty, if no errors occured)
%
% CHANGES:
%

function [sched_handles, error_code, error_msg] = update_sky_plots_highlight_quasars(sched_handles, stat_data, station_id_list, source, t_epoch_jd, quasar_id_list, plot_opt)

    % ##### Const. #####
    rad2deg = 180/pi;
    deg2rad = pi/180;

    % ##### Init #####
    error_code = 0;
    error_msg = '';
    
    flag_delete_highlighted_quasars = 0;
    flag_plot_highlighted_quasars = 0;
    
    number_of_quasars = length(quasar_id_list);
    

    % ##### Check, if there are sky plots #####
    if isempty(sched_handles.sky_plots)
        return; 
    end
    
    
     % ##### Set plot options #####
    switch(plot_opt)
        
        case 1 % - 1 :   Highlighted sources
            flag_delete_highlighted_quasars = 1;
            flag_plot_highlighted_quasars = 1;
            
            if (number_of_quasars < 1)
                error_code = 1;
                error_msg = 'Quasar ID list does not contain entries!';
                return;
            end
            
            
        case 2 % - 2 :   Delete highlighted sources (markers)
            flag_delete_highlighted_quasars = 1;
    
    end
    
    

    
    
    
    
    % Loop over all stations included in the station ID list:
    % => loop over all sky plots included in the current scan

    for i_stat = 1 : length(station_id_list)
        
        station_id = station_id_list(i_stat);
              
        % Check if current sky plot is available:
        if ~isempty(sched_handles.sky_plots(station_id).h_fig)
            
            
            % ##### Set current axes #####
            axes(sched_handles.sky_plots(station_id).h_axis);
            
            
            
            % ##### Delete highlighted quasars #####
            if flag_delete_highlighted_quasars
                % Delete a previous source markers, if they exist:
                if isfield(sched_handles.sky_plots(station_id), 'highlighted_quasars')
                    if isfield(sched_handles.sky_plots(station_id).highlighted_quasars, 'h_quasar_markers')
                        if ~isempty(sched_handles.sky_plots(station_id).highlighted_quasars.h_quasar_markers)
                           delete(sched_handles.sky_plots(station_id).highlighted_quasars.h_quasar_markers); % Remove existing markers
                        end
                    end
                end
                
            end
           

            % ##### Plot highlighted quasars #####
            if flag_plot_highlighted_quasars
                
                % Get station data for calculation of view angles (az, el)
                stat_lon = stat_data.stat(station_id).location.ellipsoid.long * pi / 180;
                stat_lat = stat_data.stat(station_id).location.ellipsoid.lat * pi / 180;
                
                % Preallocation:
                src_az = zeros(number_of_quasars, 1);
                src_el = zeros(number_of_quasars, 1);

                % ### Convert source coordinates (Ra/Dec) to Az/El ###
                for i_src = 1 : number_of_quasars
                    quasar_id = quasar_id_list(i_src);
                    [src_az(i_src), src_el(i_src), ha, dc] = zazel_s(t_epoch_jd - 2400000.5, stat_lon, stat_lat, source(quasar_id).ra, source(quasar_id).de);
                end
                

                % Transform elevation angle to a distance to the center of the plot:
                elSpherical = 90*cos(src_el); 

                % Transform data to Cartesian coordinates:
                yy = elSpherical .* cos(src_az); 
                xx = elSpherical .* sin(src_az); 

                % ##### Plot quasar markers #####
                sched_handles.sky_plots(station_id).highlighted_quasars.h_quasar_markers = plot(xx', yy', 'o', 'color', sched_handles.color_highlighted_quasars, 'MarkerSize', 6); 
                                
                % Save data da sched handels structure:
                sched_handles.sky_plots(station_id).highlighted_quasars.xx = xx;
                sched_handles.sky_plots(station_id).highlighted_quasars.yy = yy;
                sched_handles.sky_plots(station_id).highlighted_quasars.quasar_id_list = quasar_id_list;

            end
        
        end
        
    end % station loop

    
return;