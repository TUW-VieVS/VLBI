% -------------------------------------------------------------------------
%
%                              print_available_sat_obs_periods.m
%
%   This function prints the available observation periods of all scheduled
%   satellites to the MATLAB Command Window.
%
%   Author: 
%       Andreas Hellerschmied, 19.03.2015
%   
%   changes       :
%   - 2016-11-25, A. Hellerschmied: Option added to print available observation persiods for each individual station
%           
%
%   inputs        :
%   - obs_data              : observation data structure
%   - flag_from_network     : flag, print availabe observations times for each satellite with common visibility from defined station network 
%   - flag_from_stations    : flag, print availabe observations times for each satellite for each station 
%     
%
%   outputs       :
%   - Screen Output (Command Window)
%   - error_code        : Error Code (0 = no erros occured)
%   - error_msg         : Error Message (empty, if no errors occured)
%    
%
%   locals        :
% 
%
%   coupling      :
%   - invjday.m
%   - jd2datestr.m
%   
%
%   references    :
%
%-------------------------------------------------------------------------


function [error_code, error_msg] = print_available_sat_obs_periods(obs_data, flag_from_network, flag_from_stations)

    % Init
    error_code = 0;
    error_msg = '';
    
    if flag_from_network
        fprintf(1, '\n\n###### Available observation time with common visibility from all stations ######\n');
        for i_sat = 1 : obs_data.number_of_sat
            [number_of_obs_periods, col] = size(obs_data.sat(i_sat).obs_times);
            fprintf(1, '\n  %d - %s', i_sat, obs_data.sat(i_sat).name);
            for i_obs = 1 :  number_of_obs_periods
                [year, mon, day, h, min, sec] = invjday(obs_data.sat(i_sat).obs_times(i_obs, 1));
                fprintf(1, '    t%d - %4.0d %2.0d %2.0d  %02.0f:%02.0f:%05.2f (start)\n', obs_data.sat(i_sat).obs_times(i_obs, 3), ...
                                            year, ...
                                            mon, ...
                                            day, ...
                                            h, ...
                                            min, ...
                                            sec     );
                [year, mon, day, h, min, sec] = invjday(obs_data.sat(i_sat).obs_times(i_obs, 2));
                fprintf(1, '    t%d - %4.0d %2.0d %2.0d  %02.0f:%02.0f:%05.2f (end)\n', obs_data.sat(i_sat).obs_times(i_obs, 4), ...
                                            year, ...
                                            mon, ...
                                            day, ...
                                            h, ...
                                            min, ...
                                            sec     );
                fprintf(1, '\n');
            end
        end
        fprintf(1, '\n');
    end
    
    if flag_from_stations
        fprintf(1, '\n\n###### Available observation time for individual stations ######\n');
        % Loop over all selected satellites
        for i_sat = 1 : obs_data.number_of_sat
            fprintf(1, '\n  %d - %s', i_sat, obs_data.sat(i_sat).name);
            % Loop over all stations:
            for i_stat = 1 : length(obs_data.stat)
                fprintf(1, '    Station: %s:\n',obs_data.stat(i_stat).name);
                % Loop over all available obs. periods:
                for i_obs = 1 : size(obs_data.sat(i_sat).stat(i_stat).obs_times, 1)
                    fprintf('     - %s => %s\n', jd2datestr(obs_data.sat(i_sat).stat(i_stat).obs_times(i_obs, 1)), jd2datestr(obs_data.sat(i_sat).stat(i_stat).obs_times(i_obs, 2)))
                end % for i_obs = 1 : size(obs_data.sat(i_sat).stat(i_stat).obs_times, 1)
                fprintf(1, '\n');
            end % for i_stat = 1 : length(obs_data.stat)
        end % for i_sat = 1 : obs_data.number_of_sat
    end % if flag_from_stations
        
        


return;