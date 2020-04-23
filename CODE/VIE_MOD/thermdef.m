% ************************************************************************
%   Description:
%   Corrections to the delay due to thermal deformation of the VLBI antenna
%   (corr. to the rate not included). Based on Occam subroutine thermdef.f.
% 
%   Reference: 
%   Merit Standards, Appendix 11
%
%   Input:										
%       temp               measured temperatur at the station      [�C]
%       thermal            data from superstation file (antenna-info.txt)
%       azim               local source azimuth                    [rad]
%       zd                 local source zenith distance            [rad]
%       corz               correction to zd due to trop. refract.  [rad]
%       De                 apparent declination of the source      [rad]
%       axtyp              type of the axis                        ['char']
%       offs               axis offset of the antenna              [m]
%       aname              name of the antenna                     ['char']
%                
%   Output:
%       therm_d            correction to the delay due to thermal
%                          deformation of the antenna [sec]
% 
%   External calls: 	
%      global c (light velocity) 
%
%   Coded for VieVS: 
%   12 Feb 2009 by Hana Spicakova
%
%   Revision:
%   25 May 2010 by Lucia Plank: Error message
%   24 Aug 2010 by Lucia Plank: antenna type 'RICH'
%   08 Nov 2012 by Hana Kr�sn�: changed to coefficients from superstation
%   file
%   26 Feb 2015 by Lucia Plank: bug fixed for antenna types with fixed axes
%   1-(SZ*SAZ)--> 1-(SZ*SAZ)^2
%   23 Feb 2016 by Matthias Madzak: Bug fix: The sign (-) of the axis
%   offset part is changed to (+) according to Nothnagels publication.
%  
% *************************************************************************
function [therm_d] = thermdef (temp,thermal,azim,zd,corz,...
                                                    De,axtyp,offs,aname)

global c;

% original from Antenna Information File (Nothnagel)
gf=thermal.found_texp; % Foundation thermal expansion coefficient (1/K)
ga=thermal.fixedaxis_texp; % Fixed axis thermal expansion coefficient (1/K)
hf=thermal.found_h; % Height of foundation (m)
hp=thermal.fixedaxis; % Length of the fixed axis (m)
hv=thermal.avertex; % Distance from the movable axis to the antenna vertex (m)
hs=thermal.subref_h; % Height of the sub-reflector above the vertex (m)
hd=thermal.found_d; % Depth of foundation (m)
focus=thermal.focus; % Focus type: FO_PRIM -- primary, FO_SECN -- secondary
ref_temp=thermal.reftemp; % Reference temperature (degree C)


% temp0: reference temperatur
if (ref_temp ~= 999)
    temp0=ref_temp;            % ref_temp: from THERMAL.DEF
else
    temp0=temp;                % temp: measured temperatur at the station
end
 
% if no temperature is measured => temp is ref_temp
% => correction is zero
if (temp == -999)
	temp=temp0;
else
    % ### Correction for station WESTFORD: ### 
    % For calculating the thermal deformation of the antenna the
    % temperature within the radome is required. Therfor the measured or 
    % modelled exterior temp. is corrected by a formular by Arthur Niell.
    % (ref.: header of antenna-info.txt file):
    if strcmp(aname, 'WESTFORD')
       temp =  20 + 0.6 * (temp - 20); 
    end
    
end
       
% COPE WITH DIFFERENT FOCUS OF ANTENNAS
focfac=1.8;
% if antenna is prime focus => focfac=0.9
if strcmp(focus,'FO_PRIM')
	focfac=0.9;
end
% if antenna is secondary focus => focfac=1.8
if strcmp(focus,'FO_SECN')
	focfac=1.8;
end
       
     
CAZ=cos(azim);
SAZ=sin(azim);
SZ = sin (zd-corz);        %!    sin(zenith)
CZ = cos (zd-corz);        %!    cos(zenith)
CD = cos (De);             %!    cos(declination)

switch axtyp
    %  TREATING THE CASE OF AN AZIMUTH-ELEVATION MOUNTING
    case 'AZEL'
       therm_d = (gf*(temp-temp0)*hf*CZ+ ...
               ga*(temp-temp0)*(hp*CZ+hv-focfac*hs))/c ...
             + ga*(temp-temp0)*offs*SZ/c;          % Elena Skurihina's term 
    %  TREATING THE CASE OF AN EQUATORIAL ANTENNA
    case 'EQUA'
       therm_d = (gf*(temp-temp0)*hf*CZ+ ...
               ga*(temp-temp0)*(hp*CZ+hv-focfac*hs+hd*CD))/c ...
             + ga*(temp-temp0)*offs*CD/c;          % Elena Skurihina's term 
    %  TREATING THE CASE OF AN X/Y-fixed N-S axis ANTENNA
    case 'X-Y1'
        therm_d = (gf*(temp-temp0)*hf*CZ + ...
               ga*(temp-temp0)*(hp*CZ+hv-focfac*hs))/c ...
               + ga*(temp-temp0)*offs*sqrt(1-(SZ*CAZ)^2)/c;  % Elena Skurihina's term
%              - ga*(temp-temp0)*offs*sqrt(1-SZ*CAZ)/c;  % Elena Skurihina's term
    case 'X-YN'
        therm_d = (gf*(temp-temp0)*hf*CZ + ...
               ga*(temp-temp0)*(hp*CZ+hv-focfac*hs))/c ...
               + ga*(temp-temp0)*offs*sqrt(1-(SZ*CAZ)^2)/c;  % Elena Skurihina's term
%              - ga*(temp-temp0)*offs*sqrt(1-SZ*CAZ)/c;  % Elena Skurihina's term
    case 'XYNO'
        therm_d = (gf*(temp-temp0)*hf*CZ + ...
               ga*(temp-temp0)*(hp*CZ+hv-focfac*hs))/c ...
               + ga*(temp-temp0)*offs*sqrt(1-(SZ*CAZ)^2)/c;  % Elena Skurihina's term
%              - ga*(temp-temp0)*offs*sqrt(1-SZ*CAZ)/c;  % Elena Skurihina's term
 
%  TREATING THE CASE OF AN X/Y-fixed E-W axis ANTENNA
    case 'X-Y2'
        therm_d = (gf*(temp-temp0)*hf*CZ + ...
               ga*(temp-temp0)*(hp*CZ+hv-focfac*hs))/c ...
               + ga*(temp-temp0)*offs*sqrt(1-(SZ*SAZ)^2)/c;  % Elena Skurihina's term
%              - ga*(temp-temp0)*offs*sqrt(1-SZ*SAZ)/c;  % Elena Skurihina's term
    case 'X-YE'
        therm_d = (gf*(temp-temp0)*hf*CZ + ...
               ga*(temp-temp0)*(hp*CZ+hv-focfac*hs))/c ...
               + ga*(temp-temp0)*offs*sqrt(1-(SZ*SAZ)^2)/c;  % Elena Skurihina's term
%              - ga*(temp-temp0)*offs*sqrt(1-SZ*SAZ)/c;  % Elena Skurihina's term
    case 'XYEA'
        therm_d = (gf*(temp-temp0)*hf*CZ + ...
               ga*(temp-temp0)*(hp*CZ+hv-focfac*hs))/c ...
               + ga*(temp-temp0)*offs*sqrt(1-(SZ*SAZ)^2)/c;  % Elena Skurihina's term    
%              - ga*(temp-temp0)*offs*sqrt(1-SZ*SAZ)/c;  % Elena Skurihina's term         
    %  AXIS MODEL FOR RICHMOND   SHOULD BE APPLIED FOR !!!!
    case 'RICHMOND'
        therm_d = (gf*(temp-temp0)*hf*CZ + ...
                    ga*(temp-temp0)*(hp*CZ+hv-focfac*hs))/c ...
                  + ga*(temp-temp0)*offs* ...           % Elena Skurihina's term
                    sqrt(1-(CZ*sin(0.6817256)+SZ*cos(0.6817256)*...
                    (CAZ*cos(0.0020944)-SAZ*sin(0.0020944)))^2)/c;
    case 'RICH'
        therm_d = (gf*(temp-temp0)*hf*CZ + ...
                    ga*(temp-temp0)*(hp*CZ+hv-focfac*hs))/c ...
                  + ga*(temp-temp0)*offs* ...           % Elena Skurihina's term
                    sqrt(1-(CZ*sin(0.6817256)+SZ*cos(0.6817256)*...
                    (CAZ*cos(0.0020944)-SAZ*sin(0.0020944)))^2)/c;
    %  TREATING other ANTENNA types
    otherwise
        therm_d=0;
        
end




