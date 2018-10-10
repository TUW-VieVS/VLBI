% ************************************************************************
%   Description:
%   Computes an angular factor (psi) needed for calculation of the antenna 
%   axis offset altitude correction
%   Sovers, Fanselow and Jacobs, Eq. 3.203 - 3.206
% 
%   Reference: 
%   Sovers, Fanselow and Jacobs, Eq. 3.203 - 3.206
%
%   Input:										
%      ant                antenna coordinates   [m]
%      LHAe               source local hour angle east of the meridian [rad] 
%                
%   Output:
%      psifac             angular factor [-]
% 
%   External calls: 	
%   cart2phigd;
%
%   Coded for VieVS: 
%   04 Oct 2016 by Hana Krasna
%
%   Revision: 
%   02 Mar 2017 by Hana Krasna: bug fixed for Richmond
%   
% *************************************************************************
                  


function [psifac] = ao_altcorr(axtyp,ant,LHAe,azim,zd,corz)

phigd = cart2phigd(ant);


CosAZ=cos(azim);
SinAZ=sin(azim);
SinE = sin (pi/2-(zd-corz));        %   sin(elev)=sin(pi/2-(zenith-correction to zenith distance due to tropospheric refraction)))
CosE = cos (pi/2-(zd-corz));        %   cos(elev)=cos(pi/2-(zenith-correction to zenith distance due to tropospheric refraction)))


switch axtyp
    %  TREATING THE CASE OF AN EQUATORIAL ANTENNA
    case 'EQUA'
        psifac0 = cos(phigd)*cos(LHAe);
    %  TREATING THE CASE OF AN X-Y ANTENNA WITH FIXED AXIS NORTH-SOUTH
    case 'X-Y1'
        psifac0 = SinE/sqrt(1-CosAZ^2*CosE^2);
    case 'X-YN'
        psifac0 = SinE/sqrt(1-CosAZ^2*CosE^2);
    case 'XYNO'
        psifac0 = SinE/sqrt(1-CosAZ^2*CosE^2);
    %  TREATING THE CASE OF AN X-Y ANTENNA WITH FIXED AXIS EAST-WEST
    case 'X-Y2'
        psifac0 = SinE/sqrt(1-SinAZ^2*CosE^2);
    case 'X-YE'
        psifac0 = SinE/sqrt(1-SinAZ^2*CosE^2);
    case 'XYEA'
        psifac0 = SinE/sqrt(1-SinAZ^2*CosE^2);
    %  TREATING Richmond
    case 'RICHMOND'
        E=0.12/180*pi; % E=0.12deg W of N - azimuth misalignment (W of N XXX E of N - is the sign correct???)
        phiW=39.06/180*pi; % phiW = 39.06deg - latitude of Washington
        LHA_Rich = atan2(CosE*sin(azim-E),cos(phiW)*SinE-sin(phiW)*CosE*cos(azim+E));
        psifac0 = cos(phiW)*cos(LHA_Rich);
    case 'RICH'
        E=0.12/180*pi; % E=0.12deg W of N - azimuth misalignment (W of N XXX E of N - is the sign correct???)
        phiW=39.06/180*pi; % phiW = 39.06deg - latitude of Washington
        LHA_Rich = atan2(CosE*sin(azim-E),cos(phiW)*SinE-sin(phiW)*CosE*cos(azim+E));
        psifac0 = cos(phiW)*cos(LHA_Rich);

    %  TREATING other ANTENNA types
    otherwise
        psifac0 = 0;
end

psifac= psifac0;            

