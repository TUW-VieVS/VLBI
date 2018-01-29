% -------------------------------------------------------------------------
%
%                              overpass_softcopy.m
%
%   This function produces screen output of the calculate antenna pointing
%   data/overpasses in the Matlab Command Window. 
%
%   Author: 
%       Andreas Hellerschmied, 22.9.2013
%   
%   changes       :
%       - (2013.10.01) Andreas Hellerschmied
%         - Fixed bug (error if there was no rise/set element in the 
%           overpass data).
%       - 2013.11.1, A. Hellerschmied: Output of topocentric Ra. and Dec.
%           added.
%   - 2014-01.19: A. Hellerschmied: Support of different antenna axis types
%   - 2014-01.19: A. Hellerschmied: Error treatment added.
%            Remark: Output of coupled procedures is currently not checked!
%   - 2015-02-16: A. Hellerschmied: Remanded (before: "TLE_overpass_softcopy.m")
%           
%
%   inputs        :
%       - stat_data: station data structure
%       - flag_print_whole_overpass: 
%           = 1: Print whole overpass time series.
%           = 0: Print only rise/peak/set of an overpass.
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
%   
%
%   references    :
%
%-------------------------------------------------------------------------


function [error_code, error_msg] = overpass_softcopy (stat_data, flag_print_whole_overpass)

    % Init
    error_code = 0;
    error_msg = '';

    if (flag_print_whole_overpass)
        
        fprintf(1,'\n###### Overpass Time Series ######\n');
        
        for i_stat = 1 : stat_data.number_of_stations    % Stations

            fprintf(1,'\n=> Station: %s\n', stat_data.stat(i_stat).name);

            for i_sat = 1 : stat_data.number_of_sat   % Satellites

                fprintf(1,'\n    => Satellite: %5.0d - %s', stat_data.stat(i_stat).sat(i_sat).TLE_data.sat_number,...
                    stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_header_line);

                next_line_flag = 0;

                for i_epoch = 1 : stat_data.number_of_epochs  % Propagation epochs

                    if (stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).above_min_elevation)

                        if next_line_flag
                            fprintf(1,'\n');
                            next_line_flag = 0;
                        end;

                            fprintf(1,'        => %4.0d %2.0d %2.0d  %02.0f:%02.0f:%05.2f - Az: %+8.3f %+6.4f  El: %6.3f %+6.4f  R: %9.3f  Ra: %+8.3f %+6.4f  Dec: %+7.3f %+6.4f  LHA: %+8.4f %+8.5f\n',...
                                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).year,...
                                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).mon,...
                                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).day,...
                                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).h,...
                                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).min,...
                                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).sec,...
                                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).az,...
                                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).az_rate,...
                                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).el,...
                                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).el_rate,...
                                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).range/1000,...
                                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).tra * 15,...
                                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).tra_rate * 15,...
                                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).tdec,...
                                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).tdec_rate,...
                                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).local_hour_angle * 15,...
                                stat_data.stat(i_stat).sat(i_sat).epoch(i_epoch).local_hour_angle_rate * 15);
                    else 
                        next_line_flag = 1;
                    end;  

                end;    % i_epoch = 1 : number_of_epochs

            end;    % for i_sat = 1 : number_of_sat

        end;    % for i_stat = 1 : stat_data.number_of_stations
        
    end;    % if (flag_print_whole_overpass)
    
    fprintf(1,'\n###### Overpass Rise / Peak / Set ######\n');
    
        for i_stat = 1 : stat_data.number_of_stations    % Stations

            fprintf(1,'\n=> Station: %s\n', stat_data.stat(i_stat).name);

            for i_sat = 1 : stat_data.number_of_sat   % Satellites

                fprintf(1,'\n    => Satellite: %5.0d - %s', stat_data.stat(i_stat).sat(i_sat).TLE_data.sat_number,...
                    stat_data.stat(i_stat).sat(i_sat).TLE_data.TLE_header_line);

                for i_overpass = 1 : stat_data.stat(i_stat).sat(i_sat).number_of_overpasses
                    
                    % Print "Rise" Data, if available:
                    if ~isempty(stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise.year)
                                                
                        fprintf(1,'     Rise: %4.0d %2.0d %2.0d  %02.0f:%02.0f:%05.2f - Az: %+8.3f  El: %9.6f  R: %9.3f\n',...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise.year,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise.mon,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise.day,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise.h,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise.min,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise.sec,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise.az,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise.el,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).rise.range / 1000);

                    end;
                    
                    % Print "Peak" Data, if available:
                    for i_peak = 1 : stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).number_of_peaks 
                        
                        fprintf(1,'     Peak: %4.0d %2.0d %2.0d  %02.0f:%02.0f:%05.2f - Az: %+8.3f  El: %9.6f  R: %9.3f\n',...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak(i_peak).year,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak(i_peak).mon,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak(i_peak).day,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak(i_peak).h,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak(i_peak).min,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak(i_peak).sec,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak(i_peak).az,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak(i_peak).el,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).peak(i_peak).range / 1000);

                    end;    %  for i_peak = 1 : stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).number_of_peaks
                    
                    % Print "Set" Data, if available:
                    if ~isempty(stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set.year)
                        
                        fprintf(1,'     Set:  %4.0d %2.0d %2.0d  %02.0f:%02.0f:%05.2f - Az: %+8.3f  El: %9.6f  R: %9.3f\n',...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set.year,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set.mon,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set.day,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set.h,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set.min,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set.sec,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set.az,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set.el,...
                        stat_data.stat(i_stat).sat(i_sat).overpass(i_overpass).set.range / 1000);
                        
                        fprintf(1,'\n');
                        
                    end;

                end;    % i_overpass = 1 : stat_data.stat(i_stat).sat(i_sat).number_of_overpasses
                
            end;    % for i_sat = 1 : number_of_sat

        end;    % for i_stat = 1 : stat_data.number_of_stations

return;



