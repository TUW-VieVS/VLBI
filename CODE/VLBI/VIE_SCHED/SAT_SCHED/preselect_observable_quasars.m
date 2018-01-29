% #########################################################################
% #     preselect_observable_quasars
% #########################################################################
%
% DESCRIPTION
%   This function does a pre-selection of observable quasars for the defined epoch.
%
% CREATED  
%   2015-04-13     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
%   - tle_2_topo
%   - zazel_s                  : Conversion Ra/Dec => Az/El (Simplified approach for scheduling)
%   - check_quasar_visibility_from_network
%
%
% INPUT
%  
%   - stat_data                 : station data structure
%   - station_id_list           : List of IDs of the observing stations (referring to the required entry of "stat_data.stat" and "obs_data.stat")
%   - source                    : Source structure
%   - t_epoch_jd                : Time epoch for calculations
%
%
% OUTPUT
%   - source_new                    : Observable sources (structure)
%   - flag_observable_quasars_list  : List of observable sources (logical)
%   - error_code                    : Error Code (0 = no erros occured)
%   - error_msg                     : Error Message (empty, if no errors occured)
%
% CHANGES:
%

function [source_new, flag_observable_quasars_list, error_code, error_msg] = preselect_observable_quasars(stat_data, station_id_list, source, t_epoch_jd)

    % ##### Init #####
    error_code = 0;
    error_msg = '';
    
    
    % ##### Pre-select quasars, which are observable from the current network #####

    % Loop init.:
    flag_observable_quasars_list = zeros(size(source,2), 1);

    % Loop over all provided sources:
    for i_src = 1 : size(source,2)

        % Check observability:
        [flag_observable_quasars_list(i_src), error_code, error_msg] = check_quasar_visibility_from_network(stat_data, station_id_list, source(i_src), t_epoch_jd);
        if error_code > 0
            error_msg = ['check_quasar_visibility_from_network: ', error_msg];
            return; 
        end
    end

    flag_observable_quasars_list = logical(flag_observable_quasars_list);
    
    
    % Preselection => Exclude all sources, which are currently not observable:
    source_new = source(flag_observable_quasars_list);
    
return;