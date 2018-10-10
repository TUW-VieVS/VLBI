% ************************************************************************
%   Description:
%   Corrections to the delay due to axis offset of the VLBI antenna
%   Corrections to the delay rate is not included! Based on occam
%   subroutine axis.
%
%   Input:
%      antPEF_ell  ellipsoial coordinates of antenna (lam,phi,hgt) [rad,rad,m]
%      azim        azimuth [rad]
%      zd          zenith distance [rad]
%      corz        correction to zenith distance due to trop. refract. [rad]
%      De          apparent declination of the source [rad]
%      axtyp       typ of the axis ['char']
%      offs        axis offset of the antenna [m]
%      aname       name of the antenna ['char']
% 
%   Output:
%      axkt        correction to the delay (and rate) due to axis offset
%                  of the antenna [sec]      
%
%   External calls: 	
%       ---					    											
%       
%   Coded for VieVS: 
%   01 Feb 2009 by Hana Spicakova
%
%   Revision:
%   25 May 2010 by Lucia Plank: Error message
%   23 Aug 2010 by Lucia Plank: antenna type 'RICH'
%   04 Oct 2013 by Hana Krasna: partials for antenna axis offset added
%
% ************************************************************************

function [axkt daxkt] = axis_stat (phi,azim,zd,corz,De,axtyp,offs,aname)
 
global c;

CAZ=cos(azim);
SAZ=sin(azim);
SZ = sin (zd-corz);        %   sin(zenith)
CZ = cos (zd-corz);        %   cos(zenith)
%CD = cos (De);            %   cos(declination)
%CHA = cos (LHA);
%SHA = sin (LHA);

switch axtyp
    %  TREATING THE CASE OF AN AZIMUTH-ELEVATION MOUNTING
    case 'AZEL'          
        axkt0 = - (offs * SZ)/c;
        daxkt0 = -SZ; % [-]
    %  TREATING THE CASE OF AN EQUATORIAL ANTENNA
    case 'EQUA'
        %  COMPUTE ANGLE BETWEEN NORTH DIRECTION AND SOURCE DIRECTION
        %  CORRECTED BY REFRACTION
        csx = SZ*CAZ*cos(phi) + CZ*sin(phi);
        sn  = acos(csx);    
        snx = sin(sn);
        %A0 = -offs * CD / c;
        axkt0 = -offs * snx/c;
        daxkt0 = -snx;
    %  TREATING THE CASE OF AN X-Y ANTENNA WITH FIXED AXIS NORTH-SOUTH
    case 'X-Y1'
        axkt0 = -(offs * sqrt(1 - (SZ*CAZ)^2))/c;
        daxkt0 =-(sqrt(1 - (SZ*CAZ)^2));
    case 'X-YN'
        axkt0 = -(offs * sqrt(1 - (SZ*CAZ)^2))/c;
        daxkt0 = -(sqrt(1 - (SZ*CAZ)^2));
    case 'XYNO'
        axkt0 = -(offs * sqrt(1 - (SZ*CAZ)^2))/c; 
        daxkt0 =-(sqrt(1 - (SZ*CAZ)^2));
    %  TREATING THE CASE OF AN X-Y ANTENNA WITH FIXED AXIS EAST-WEST
    case 'X-Y2'
        axkt0 = -(offs * sqrt(1 - (SZ*SAZ)^2)) / c;
        daxkt0 =-(sqrt(1 - (SZ*SAZ)^2)) ;
    case 'X-YE'
        axkt0 = -(offs * sqrt(1 - (SZ*SAZ)^2)) / c;
        daxkt0 =-(sqrt(1 - (SZ*SAZ)^2)) ;
    case 'XYEA'
        axkt0 = -(offs * sqrt(1 - (SZ*SAZ)^2)) / c;
        daxkt0 =-(sqrt(1 - (SZ*SAZ)^2)) ;
    case 'RICHMOND'
        axkt0  = - (offs * sqrt(1 - (CZ*sin(0.6817256)+ ...
                SZ*cos(0.6817256)*(CAZ*cos(0.0020944)- ...
                SAZ*sin(0.0020944)))^2))/c;
        daxkt0 =-(sqrt(1 - (CZ*sin(0.6817256)+ ...
                SZ*cos(0.6817256)*(CAZ*cos(0.0020944)- ...
                SAZ*sin(0.0020944)))^2));
    case 'RICH'
        axkt0  = - (offs * sqrt(1 - (CZ*sin(0.6817256)+ ...
                SZ*cos(0.6817256)*(CAZ*cos(0.0020944)- ...
                SAZ*sin(0.0020944)))^2))/c;
       daxkt0 = - (sqrt(1 - (CZ*sin(0.6817256)+ ...
                SZ*cos(0.6817256)*(CAZ*cos(0.0020944)- ...
                SAZ*sin(0.0020944)))^2));     
    %  TREATING other ANTENNA types
    otherwise
        axkt0 = 0;
        daxkt0 =0;
end
       
axkt= axkt0;             % correction to the delay
%axkt(2) = 0;                % correction to the delay rate is not computed!!!!

daxkt = daxkt0; % partial derivative
