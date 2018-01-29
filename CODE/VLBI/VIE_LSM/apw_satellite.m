% ************************************************************************
%   Description:
%   function to form the design matrix for satellite coordinates as
%   piecewise linear offset functions
%
%   Reference:
%
%   Input:
%       'per_satellite'     structure array  
%       'n_observ'          (1,1)               total number of observation in the session
%       'n_unk'             structure array     number of estimation intervals
%       'T_'                structure array     estimation epochs
%
%   Output:
%       'Apw_pos1'      (nobserv,number of pwlo unknowns)       design matrix for 1st satellite coord.
%       'Apw_pos2'      (nobserv,number of pwlo unknowns)       design matrix for 2nd satellite coord.
%       'Apw_pos3'      (nobserv,number of pwlo unknowns)       design matrix for 3rd satellite coord.
%
%   External calls:
%
%   Coded for VieVS:
%   04 May 2017 by A. Hellerschmied
%
%   Revision:
%
% ************************************************************************
function [Apw_pos1, Apw_pos2, Apw_pos3] = apw_satellite(per_satellite, n_observ, n_unk, T_)

% Init. an preallocate:
Apw_pos1 = zeros(n_observ, n_unk.sat + 1);
Apw_pos2 = zeros(n_observ, n_unk.sat + 1);
Apw_pos3 = zeros(n_observ, n_unk.sat + 1);

flag_within_interval = false;

minute              = per_satellite.minute;                % Estimation eouchs in minutes since mjd0
est_int_id_sat_pos  = per_satellite.est_int_id_sat_pos;    % Index of estimation interval that contains the current observation
nob                 = per_satellite.nob;                   % Onservation index  in session(absolut)

n_obs_in_int_list  = zeros(n_unk.sat, 1);


% Loop over all estimation intervals:
for i_inter = 1 : n_unk.sat
    n_obs_in_int = 0; % number of observations to the curretn target within the estimation interval with index "i_inter"
    
    % Loop over all observations to this satellite:
    for i_obs = 1 : per_satellite.total
        if i_inter == n_unk.sat
            flag_within_interval = (minute(i_obs) >= T_.sat(i_inter)) && (minute(i_obs) <= T_.sat(i_inter+1));
        else
            flag_within_interval = (minute(i_obs) >= T_.sat(i_inter)) && (minute(i_obs) < T_.sat(i_inter+1));
        end
        if flag_within_interval
            n_obs_in_int = n_obs_in_int + 1;
            n_obs_in_int_list(i_inter) = n_obs_in_int;
        end
    end
end

k = 0;

% Loop over all estimation intervals:
for i_inter = 1 : n_unk.sat
    
    % Loop over all observations to the target within the estimation interval with index "i_inter"
    for i_obs = 1 : n_obs_in_int_list(i_inter) % number of obs. of the source in an estimation interval
        k = k + 1;
        
        % pos 1
        try
        Apw_pos1(nob(k), i_inter)    = (1 - (minute(k) - T_.sat(est_int_id_sat_pos(k))) / (T_.sat(est_int_id_sat_pos(k) + 1) - T_.sat(est_int_id_sat_pos(k)))) * per_satellite.pd_pos1(k);
        Apw_pos1(nob(k) ,i_inter+1)  = ((minute(k) - T_.sat(est_int_id_sat_pos(k))) / (T_.sat(est_int_id_sat_pos(k) + 1) - T_.sat(est_int_id_sat_pos(k)))) * per_satellite.pd_pos1(k);
        catch
            keyboard
        end
        
        % pos 2
        Apw_pos2(nob(k), i_inter)    = (1 - (minute(k) - T_.sat(est_int_id_sat_pos(k))) / (T_.sat(est_int_id_sat_pos(k) + 1)-T_.sat(est_int_id_sat_pos(k)))) * per_satellite.pd_pos2(k);
        Apw_pos2(nob(k), i_inter+1)  = ((minute(k) - T_.sat(est_int_id_sat_pos(k))) / (T_.sat(est_int_id_sat_pos(k) + 1) - T_.sat(est_int_id_sat_pos(k)))) * per_satellite.pd_pos2(k);
        
        % pos 3:
        Apw_pos3(nob(k), i_inter)    = (1 - (minute(k) - T_.sat(est_int_id_sat_pos(k))) / (T_.sat(est_int_id_sat_pos(k) + 1) - T_.sat(est_int_id_sat_pos(k)))) * per_satellite.pd_pos3(k);
        Apw_pos3(nob(k), i_inter+1)  = ((minute(k) - T_.sat(est_int_id_sat_pos(k))) / (T_.sat(est_int_id_sat_pos(k) + 1) - T_.sat(est_int_id_sat_pos(k)))) * per_satellite.pd_pos3(k);
        
    end % for i_obs = 1 : n_obs_in_int_list(i_inter)
end % for i_inter = 1 : n_unk.sat






