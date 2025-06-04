% ************************************************************************
%   Description:
%	Computes the position of the satellite at the estimation time and adds 
%   the information to the sources struct.
%
%   Input:										
%     sources               sources struct
%     estInt				estimation intervals
%     i_sat					id of sat
%     mjd0					beginning of the session
%     scan					scan struct
%     type					type (position or KepEle)
%     numKepEle				number of Keplerian Element
%
% 
%   Output:
%     sources              sources struct
% 
%   External calls: 	
%       
%   Coded for VieVS: 
%   17 Dec 2024 by Helene Wolf
%
%   Revision: 
%
%
% ************************************************************************
function [sources] = addSatellitePositionAtEstimationInterval(sources, estInt, i_sat, mjd0, scan, type, numKepEle)

    mjdsat = sources.s(i_sat).mjd;
    if strcmp(type, 'position')
            sources.s(i_sat).posEstIntXtrf = zeros(length(estInt),3);
            sources.s(i_sat).posEstIntXcrf = zeros(length(estInt),3);
    elseif strcmp(type, 'KepEle')
            textInt = 'EstIntXtrfKepEle' + string(numKepEle);
            textIntcrf = 'EstIntXcrfKepEle' + string(numKepEle);
            textMjd = 'EstIntMjdKepEle' + string(numKepEle);
            sources.s(i_sat).(textInt) = zeros(length(estInt),3);
    end

    for i=1:length(estInt)
        mjd_est = mjd0 + estInt(i)/(24*60);
        [~, idx] = min(abs(mjdsat - mjd_est));
        [~, ~, ~, hour, minute, sec] = mjd2date(mjd_est);
        secOfDay = hour*3600 + minute *60 + sec;

        refidx           = idx(1)-5 : 1 : idx(1)+5; % suitable for interpolation with "lagint9.m" and "dt" < 1 sec
        tRefSecInterpol  = sources.s(i_sat).sec_of_day(refidx);
        tIntegerMjd      = floor(sources.s(i_sat).mjd(refidx));
        tRefMjd          = tIntegerMjd(1);
        offsetSec        = (tIntegerMjd - tRefMjd) * 86400; % full days since first interpolation epoch in [sec]
        tRefSecInterpol  = tRefSecInterpol + offsetSec; % add since first interpolation epoch in [sec]
        tRefSec          = tRefSecInterpol(1);
        tRefSecInterpol  = tRefSecInterpol - tRefSec;

        tIntegerMjdObs   = floor(mjd_est);
        tRefOffsetObs    = (tIntegerMjdObs - tRefMjd) * 86400;
        tRefSecObs       = secOfDay + tRefOffsetObs;
        tRefSecObs       = tRefSecObs -  tRefSec;

        trfScPosX = lagint9(tRefSecInterpol, sources.s(i_sat).x_trf(refidx), tRefSecObs);
        trfScPosY = lagint9(tRefSecInterpol, sources.s(i_sat).y_trf(refidx), tRefSecObs);
        trfScPosZ = lagint9(tRefSecInterpol, sources.s(i_sat).z_trf(refidx), tRefSecObs);
        crfScPosX = lagint9(tRefSecInterpol, sources.s(i_sat).x_crf(refidx), tRefSecObs);
        crfScPosY = lagint9(tRefSecInterpol, sources.s(i_sat).y_crf(refidx), tRefSecObs);
        crfScPosZ = lagint9(tRefSecInterpol, sources.s(i_sat).z_crf(refidx), tRefSecObs);
        
        if strcmp(type, 'position')
            sources.s(i_sat).posEstIntXtrf(i, :) = [trfScPosX trfScPosY trfScPosZ];
            sources.s(i_sat).posEstIntXcrf(i, :) = [crfScPosX crfScPosY crfScPosZ];
            sources.s(i_sat).posEstIntMjd(i,1) = mjd_est;
        elseif strcmp(type, 'KepEle') || strcmp(type, 'SatPara')
            sources.s(i_sat).(textInt)(i, :) = [trfScPosX trfScPosY trfScPosZ];
            sources.s(i_sat).(textIntcrf)(i, :) = [crfScPosX crfScPosY crfScPosZ];
            sources.s(i_sat).(textMjd)(i,1) = mjd_est;
        end
    end
    
    if ~isfield(sources.s(i_sat),'obs')
        obs = zeros(1,length([sources.s(i_sat).mjd]));
        
        for iScan = 1:length(scan)
            if scan(iScan).obs_type == 's'
                mjdScan = scan(iScan).mjd;
                for i=1:length([sources.s(i_sat).mjd])
                    if mjdScan > sources.s(i_sat).mjd(i) - 1/(24*60)*4 && mjdScan < sources.s(i_sat).mjd(i) + 1/(24*60)*4
                        obs(i) = 1;
                    end
                end
            else
                continue
            end
        end
        sources.s(i_sat).obs = obs; 
    end
end