% ************************************************************************
%   Description:
%   function to form the design matrix for orbital/solar radiation pressure parameters or satellite position as
%   piecewise linear offset functions
%
%   Input:										
%		per_satellite		struct with information per satellite
%   	n_observ			number of observations
%   	num		        	id of parameter
%       type                type (orb/srp/sat_pos),
%       num2                in case of satellite position: sat_pos1 / sat_pos2 / sat_pos3
% 
%   Output:
%      'Apw'      (nobserv,number of pwlo unknowns)       design matrix for
%      orbital/solar radiation /satellite_position
% 
%   External calls: 	
%       
%       
%   Coded for VieVS: 
%   12 July 2022 Helene Wolf
%
%   Revision: 
%   10 June 2026: edited function to have one function for all three kind
%   of satellite parameters (orbital/ solar radiation pressure parameters) 
% ************************************************************************

function [Apw] = apw_orb(per_satellite, n_observ, num, type, num2) 
    
    if nargin <5
        num2 = num; 
    end

    minute              = per_satellite.minute;                % Estimation epochs in minutes since mjd0
    estIntId            = per_satellite.('est_int_id_' + string(type) + string(num));    % Index of estimation interval that contains the current observation
    nob                 = per_satellite.nob;                   % Onservation index  in session(absolut)
    T                   = per_satellite.('T_' + string(type) + string(num));
    n_unk               = per_satellite.('n_unk_' + string(type) + string(num));

    Apw = zeros(n_observ, n_unk);
    pd = per_satellite.('pd_' + string(type) + string(num2));

    if n_unk == 1
        Apw(nob(1:numel(pd)), 1) = pd(:);
    else
         for i=1:length(T)
            valid = (estIntId == i);         
            if any(valid)
                obs_idx = nob(valid);
                int_idx = estIntId(valid);
                
                interval_len = T(int_idx + 1) - T(int_idx);
                interval_len(interval_len == 0) = 1;
                t = (minute(valid) - T(int_idx)) ./ interval_len;
                
                pd_val = pd(valid);
                int_idx = unique(int_idx);
                Apw(obs_idx, int_idx)    = (1 - t) .* pd_val;
                Apw(obs_idx, int_idx + 1) = t .* pd_val;
            end
        end
    end
    
    if any(sum(Apw, 1) == 0)
        error('*** The interval for estimating the orbital elements is too short – longer intervals exist that contain no satellite observations to the satellite %s.', per_satellite.name);
    end
end