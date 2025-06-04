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
%       'A'      wiht .pos1 / .pos2 / .pos3 -> design matrix for 3 satellite coord.
%
%   External calls:
%
%   Coded for VieVS:
%   04 May 2017 by A. Hellerschmied
%
%   Revision:
%	17 Dec 2024 by H. Wolf: per_satellite information is created in separate function
% ************************************************************************
function [A] = apw_satellite(per_satellite, n_observ)

    Apw_pos1 = zeros(n_observ, per_satellite.n_unk_sat_pos+ 1);
    Apw_pos2 = zeros(n_observ, per_satellite.n_unk_sat_pos + 1);
    Apw_pos3 = zeros(n_observ, per_satellite.n_unk_sat_pos + 1);
    
    flag_within_interval = false;
    
    minute              = per_satellite.minute_sat_pos;        % Estimation eouchs in minutes since mjd0
    est_int_id_sat_pos  = per_satellite.est_int_id_sat_pos;    % Index of estimation interval that contains the current observation
    nob                 = per_satellite.nob;                   % Onservation index  in session(absolut)
    T_                  = per_satellite.T_sat_pos;
    n_unk               = per_satellite.n_unk_sat_pos;

    n_obs_in_int_list  = zeros(n_unk, 1);
    
    for i_inter = 1 : n_unk
        n_obs_in_int = 0; % number of observations to the current target within the estimation interval with index "i_inter"
        
        for i_obs = 1 : per_satellite.total
            if i_inter == n_unk
                flag_within_interval = (minute(i_obs) >= T_(i_inter)) && (minute(i_obs) <= T_(i_inter+1));
            else
                flag_within_interval = (minute(i_obs) >= T_(i_inter)) && (minute(i_obs) < T_(i_inter+1));
            end
            if flag_within_interval
                n_obs_in_int = n_obs_in_int + 1;
                n_obs_in_int_list(i_inter) = n_obs_in_int;
            end
        end
    end
    
    if n_unk > 0
        k = 0;  
        % Loop over all estimation intervals:
        for i_inter = 1 : n_unk
            for i_obs = 1 : n_obs_in_int_list(i_inter) % number of obs. of the source in an estimation interval
                k = k + 1;
                
                Apw_pos1(nob(k), i_inter)    = (1 - (minute(k) - T_(est_int_id_sat_pos(k))) / (T_(est_int_id_sat_pos(k) + 1) - T_(est_int_id_sat_pos(k)))) * per_satellite.pd_pos1(k);
                Apw_pos1(nob(k) ,i_inter+1)  = ((minute(k) - T_(est_int_id_sat_pos(k))) / (T_(est_int_id_sat_pos(k) + 1) - T_(est_int_id_sat_pos(k)))) * per_satellite.pd_pos1(k);
    
                Apw_pos2(nob(k), i_inter)    = (1 - (minute(k) - T_(est_int_id_sat_pos(k))) / (T_(est_int_id_sat_pos(k) + 1)-T_(est_int_id_sat_pos(k)))) * per_satellite.pd_pos2(k);
                Apw_pos2(nob(k), i_inter+1)  = ((minute(k) - T_(est_int_id_sat_pos(k))) / (T_(est_int_id_sat_pos(k) + 1) - T_(est_int_id_sat_pos(k)))) * per_satellite.pd_pos2(k);
                
                Apw_pos3(nob(k), i_inter)    = (1 - (minute(k) - T_(est_int_id_sat_pos(k))) / (T_(est_int_id_sat_pos(k) + 1) - T_(est_int_id_sat_pos(k)))) * per_satellite.pd_pos3(k);
                Apw_pos3(nob(k), i_inter+1)  = ((minute(k) - T_(est_int_id_sat_pos(k))) / (T_(est_int_id_sat_pos(k) + 1) - T_(est_int_id_sat_pos(k)))) * per_satellite.pd_pos3(k);
            end 
        end 
    else
        for iObs = 1 : length(nob) % number of obs. of the source in an estimation interval
            Apw_pos1(nob(iObs))    = per_satellite.pd_pos1(iObs);
            Apw_pos2(nob(iObs))    = per_satellite.pd_pos2(iObs);
            Apw_pos3(nob(iObs))    = per_satellite.pd_pos3(iObs);
        end
    end

    A.pos1 = Apw_pos1;
    A.pos2 = Apw_pos2;
    A.pos3 = Apw_pos3;

    if any(sum(Apw_pos1,1)==0)
        error('*** The interval for estimating the satellite position is too short – longer intervals exist that contain no satellite observations to the satellite %s.', per_satellite.name)
    end
end