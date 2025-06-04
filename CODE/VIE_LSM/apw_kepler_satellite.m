% ************************************************************************
%   Description:
%   function to form the design matrix for keplerian elements as
%   piecewise linear offset functions
%
%   Input:										
%		per_satellite		struct with information per satellite
%   	n_observ			number of observations
%   	numKepEle			id of Keplerian Element
% 
%   Output:
%      'Apw_KepEle'      (nobserv,number of pwlo unknowns)       design matrix for keplerian elements
% 
%   External calls: 	
%       
%       
%   Coded for VieVS: 
%   12 July 2022 Helene Wolf
%
%   Revision: 
%
% ************************************************************************
function [Apw_KepEle] = apw_kepler_satellite(per_satellite, n_observ, numKepEle)
  
    minute              = per_satellite.('minute_KepEle' + string(numKepEle));                % Estimation epochs in minutes since mjd0
    estIntId            = per_satellite.('est_int_id_KepEle' + string(numKepEle));    % Index of estimation interval that contains the current observation
    nob                 = per_satellite.nob;                   % Onservation index  in session(absolut)
    T_                  = per_satellite.('T_KepEle' + string(numKepEle));
    n_unk               = per_satellite.('n_unk_KepEle' + string(numKepEle));
    n_obs_in_int_list   = zeros(n_unk, 1);

    Apw_KepEle = zeros(n_observ, n_unk + 1);
    flag_within_interval = false;

    % Loop over all estimation intervals:
    for iEstInt = 1 : n_unk
        n_obs_in_int = 0; % number of observations to the current target within the estimation interval with index "iEstInt"
        for iObs = 1 : per_satellite.total
            if iEstInt == n_unk
                flag_within_interval = (minute(iObs) >= T_(iEstInt)) && (minute(iObs) <= T_(iEstInt+1));
            else
                flag_within_interval = (minute(iObs) >= T_(iEstInt)) && (minute(iObs) < T_(iEstInt+1));
            end
            if flag_within_interval
                n_obs_in_int = n_obs_in_int + 1;
                n_obs_in_int_list(iEstInt) = n_obs_in_int;
            end
        end
    end
        
    if n_unk > 0
        k = 0;  
        for iEstInt = 1 : n_unk
            for iObs = 1 : n_obs_in_int_list(iEstInt) % number of obs. of the source in an estimation interval
                k= k + 1;
                Apw_KepEle(nob(k), iEstInt)    = (1 - (minute(k) - T_(estIntId(k))) / (T_(estIntId(k) + 1) - T_(estIntId(k)))) * per_satellite.('pd_KepEle' + string(numKepEle))(k);
                Apw_KepEle(nob(k) ,iEstInt+1)  = ((minute(k) - T_(estIntId(k))) / (T_(estIntId(k) + 1) - T_(estIntId(k)))) * per_satellite.('pd_KepEle' + string(numKepEle))(k);
            end 
        end
    else
        for iObs = 1 : length(nob) % number of obs. of the source in an estimation interval
            Apw_KepEle(nob(iObs))    = per_satellite.('pd_KepEle' + string(numKepEle))(iObs);
        end
    end

    if any(sum(Apw_KepEle,1)==0)
        error('*** The interval for estimating the orbital elements is too short – longer intervals exist that contain no satellite observations to the satellite %s.', per_satellite.name)
    end
end