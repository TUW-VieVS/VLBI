% ************************************************************************
%   Description:
%   function to determine the satellite position in the CRF
%
%   Reference:
%
%   Input:
%       'sources' structure array    sources structure 
%       'T2C_s'   (3,3,n)            terrrestrial to celestial matrices   
%
%   Output:
%       'sources' structure array    sources structure including CRF Position (and Velocity) 
%
%   External calls:
%
%   Coded for VieVS:
%   22 November 2021 by H. Wolf - created as external function of vie_mod
%
%   Revision:
%
% ************************************************************************
function [sources] = getCRF_PosVelSatellites(sources, T2C_s)
    numberOfOrbitEpochs = length([sources.s(1).mjd]);
    
    % loop over space crafts:
    for iSc = 1 : length(sources.s)
        xyzCRFtmp     = zeros(numberOfOrbitEpochs, 3);
        if  sources.s(iSc).flag_v_trf
            v_xyzCRFtmp   = zeros(numberOfOrbitEpochs, 3);
        end
        % loop over all orbit pos. epochs:
        for iOrbitEpoch = 1 : numberOfOrbitEpochs
            % Position:
            xyzCRFtmp(iOrbitEpoch, :) = (T2C_s(:, :, iOrbitEpoch) * [sources.s(iSc).x_trf(iOrbitEpoch); sources.s(iSc).y_trf(iOrbitEpoch); sources.s(iSc).z_trf(iOrbitEpoch)])';
            % Velocity:
            if  sources.s(iSc).flag_v_trf
                v_xyzCRFtmp(iOrbitEpoch, :) = (T2C_s(:, :, iOrbitEpoch) * ([sources.s(iSc).vx_trf(iOrbitEpoch); sources.s(iSc).vy_trf(iOrbitEpoch); sources.s(iSc).vz_trf(iOrbitEpoch)] + cross([0; 0; omega], [sources.s(iSc).x_trf(iOrbitEpoch); sources.s(iSc).y_trf(iOrbitEpoch); sources.s(iSc).z_trf(iOrbitEpoch)])))';
            end
        end
        % store results in sources structure
        sources.s(iSc).x_crf = xyzCRFtmp(:, 1);
        sources.s(iSc).y_crf = xyzCRFtmp(:, 2);
        sources.s(iSc).z_crf = xyzCRFtmp(:, 3);
        if  sources.s(iSc).flag_v_trf
            sources.s(iSc).vx_crf = v_xyzCRFtmp(:, 1);
            sources.s(iSc).vy_crf = v_xyzCRFtmp(:, 2);
            sources.s(iSc).vz_crf = v_xyzCRFtmp(:, 3);
        end
    end
end

