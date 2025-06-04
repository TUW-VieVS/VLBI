% ************************************************************************
%   Description:
%   This function adds the orbit information of the "orbit_data" structure to the 
%   VieVS "sources" structure
%
%
% INPUT
%   orbit_data          - structure containing the orbit data
%   sources             - VieVS sources structure
%  
% OUTPUT
%  - sources             - VieVS sources structure (updates satellite positions)
%
%
%   Coded for VieVS: 
%   16 Aug 2016 by A. Hellerschmied
%
%   Revision: 
%  - yyyy-mm-dd, <first + second name>: Description
%  - 2016-09-21, A. Hellerschmied: Added field vx_trf, vy_trf, vz_trf and flag_v_trf to orbit_data and source.s structures 
%  - 2016-12-05, A. Hellerschmied: Added fields 'sec_of_day', 'year', 'month', 'day', 'hour', 'minu', 'sec' to orbit_data and source.s structures
%  - 2024-12-16, H. Wolf: Added the velocity computation in TRF 
% ************************************************************************

function [sources] = orbit_data2sources(orbit_data, sources)

for i_sat = 1 : length(sources.s)
    orbit_data_ind = strcmp({orbit_data.sat.name}, sources.s(i_sat).name);
    sources.s(i_sat).x_trf      = orbit_data.sat(orbit_data_ind).x_trf;
    sources.s(i_sat).y_trf      = orbit_data.sat(orbit_data_ind).y_trf;
    sources.s(i_sat).z_trf      = orbit_data.sat(orbit_data_ind).z_trf;
    sources.s(i_sat).mjd        = orbit_data.epoch_mjd;
    sources.s(i_sat).year       = orbit_data.year;
    sources.s(i_sat).month      = orbit_data.month;
    sources.s(i_sat).day        = orbit_data.day;
    sources.s(i_sat).hour       = orbit_data.hour;
    sources.s(i_sat).minu       = orbit_data.minu;
    sources.s(i_sat).sec        = orbit_data.sec;
    sources.s(i_sat).sec_of_day = orbit_data.sec_of_day;
    sources.s(i_sat).tosc = fix(orbit_data.tosc);
    sources.s(i_sat).flag_v_trf = false;
   
    % compute trf velocity based on interpolation and numerical derivative 
    sources.s(i_sat).vx_trf =zeros(length(sources.s(i_sat).x_trf),1);
    sources.s(i_sat).vy_trf =zeros(length(sources.s(i_sat).y_trf),1);
    sources.s(i_sat).vz_trf =zeros(length(sources.s(i_sat).z_trf),1);
    for idx= 100:700      
        dSec = 1; 
        refidx           = idx-7 : 1 : idx+6; % suitable for interpolation with "lagint9.m" and "dt" < 1 sec
        tRefSecInterpol = sources.s(i_sat).sec_of_day(refidx);
        tIntegerMjd   = floor(sources.s(i_sat).mjd(refidx));
        tRefMjd          = tIntegerMjd(1);
        offsetSec        = (tIntegerMjd - tRefMjd) * 86400; % full days since first interpolation epoch in [sec]
        tRefSecInterpol  = tRefSecInterpol + offsetSec; % add since first interpolation epoch in [sec]
        tRefSec          = tRefSecInterpol(1);
        tRefSecInterpol  = tRefSecInterpol - tRefSec;

        tIntegerMjdObs   = floor(sources.s(i_sat).mjd(idx));
        tRefOffsetObs    = (tIntegerMjdObs - tRefMjd) * 86400;
        tRefSecObs       = sources.s(i_sat).sec_of_day(idx) + tRefOffsetObs;
        tRefSecObs       = tRefSecObs -  tRefSec;

        tRefSecm = tRefSecObs - dSec;
        tRefSecp = tRefSecObs + dSec;

        trfScPosXm = lagint9(tRefSecInterpol, sources.s(i_sat).x_trf(refidx), tRefSecm);
        trfScPosYm = lagint9(tRefSecInterpol, sources.s(i_sat).y_trf(refidx), tRefSecm);
        trfScPosZm = lagint9(tRefSecInterpol, sources.s(i_sat).z_trf(refidx), tRefSecm);

        trfScPosXp = lagint9(tRefSecInterpol, sources.s(i_sat).x_trf(refidx), tRefSecp);
        trfScPosYp = lagint9(tRefSecInterpol, sources.s(i_sat).y_trf(refidx), tRefSecp);
        trfScPosZp = lagint9(tRefSecInterpol, sources.s(i_sat).z_trf(refidx), tRefSecp);

        trfScVelX =  (trfScPosXp - trfScPosXm)/(2*dSec);
        trfScVelY =  (trfScPosYp - trfScPosYm)/(2*dSec); 
        trfScVelZ =  (trfScPosZp - trfScPosZm)/(2*dSec); 

        sources.s(i_sat).vx_trf(idx) = trfScVelX;
        sources.s(i_sat).vy_trf(idx) = trfScVelY;
        sources.s(i_sat).vz_trf(idx) = trfScVelZ;
    end
    sources.s(i_sat).flag_v_trf = true;
end