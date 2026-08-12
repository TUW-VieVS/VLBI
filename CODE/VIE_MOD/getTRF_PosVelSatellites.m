% ************************************************************************
%   Description:
%   function to determine the satellite position and velocity in the TRF
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

function [sources] = getTRF_PosVelSatellites(sources, T2C_s)
    global omega

    for iSc = 1 : length(sources.s)
        numberOfOrbitEpochs = length([sources.s(iSc).mjd]);
        src = sources.s(iSc);
        xyzTRFtmp     = zeros(numberOfOrbitEpochs, 3);
        if  src.flag_v_crf
            v_xyzTRFtmp   = xyzTRFtmp;
        end
        % loop over all orbit pos. epochs:
        for iOrbitEpoch = 1 : numberOfOrbitEpochs
            R = T2C_s(:, :, iOrbitEpoch)';
            xyzTRFtmp(iOrbitEpoch, :) = (R * [src.x_crf(iOrbitEpoch); src.y_crf(iOrbitEpoch); src.z_crf(iOrbitEpoch)])';
            if  sources.s(iSc).flag_v_crf
                v_xyzTRFtmp(iOrbitEpoch, :) = R * [src.vx_crf(iOrbitEpoch); src.vy_crf(iOrbitEpoch); src.vz_crf(iOrbitEpoch)] - cross([0; 0; omega], R*[src.x_crf(iOrbitEpoch); src.y_crf(iOrbitEpoch); src.z_crf(iOrbitEpoch)]);
            end
        end
        sources.s(iSc).x_trf = xyzTRFtmp(:, 1);
        sources.s(iSc).y_trf = xyzTRFtmp(:, 2);
        sources.s(iSc).z_trf = xyzTRFtmp(:, 3);
        if  src.flag_v_crf
            sources.s(iSc).vx_trf = v_xyzTRFtmp(:, 1);
            sources.s(iSc).vy_trf = v_xyzTRFtmp(:, 2);
            sources.s(iSc).vz_trf = v_xyzTRFtmp(:, 3);
            sources.s(iSc).flag_v_trf = 1;
        end
    end
end