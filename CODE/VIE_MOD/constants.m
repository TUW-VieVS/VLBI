% ************************************************************************
%   Description:
%   Constants needed in Vie_mod. Variables are used as globals.
%  
%   Reference: 
%   Various, see below             					    											
%       
%   Coded for VieVS: 
%   23 Nov 2009 by Lucia Plank
%
%   Revision: 
%   27 Sep 2011 by Hana Spicakova  parameter K_opl added as global parameter
%   (needed for computation of station displacement due to ocean pole tide loading)
%  
%
% *************************************************************************
function constants

% Reference
%   (REF:1) ..... IERS numerical standards, IERS Conv. Chap. 1 

global c 
global gms gme gmm  
global nstations ntimmax nobserv rad2uas au omega;
global rearthm massrSE massrME a_tidefree f_tidefree a_grs80 f_grs80
global K_opl

% light velocity in m/s (1)
c    = 299792458;   %[m/s], (1)
cinv = .4990047838061;
% astronomical unit 
au   = cinv * c;
%factors for matthews.m --> NO REFERENCE!
rearthm=6378136.55; %[m]
massrSE=332945.943062; %mass ratio Sun/Earth
massrME=0.012300034;   %mass ratio Moon/Earth

% IERS 2003 numerical standards
% ellipsoid parameters for xyz2ellip.m
a_tidefree = 6378136.6; %m      Equatorial radius of the Earth
f_tidefree = 1/298.25642;     % Flattening factor of the Earth
ge_tidefree = 9.7803278; %[m/s^2]
% GRS 80
% (http://www.bkg.bund.de/nn_164850/geodIS/EVRS/EN/References/...
%  Definitions/Def__GRS80-pdf,templateId=raw,property=publication...
%  File.pdf/Def_GRS80-pdf.pdf)
a_grs80    = 6378137;
f_grs80    = 0.00335281068118;

% Constant of gravitation
g   = 6.67428e-11; %[m^3/kg s^2]
% gravitational constant for the sun(G*M_SUN)[m^3/s^2]
gms = 1.32712442076e20; % [m^3/s^2] (REF: 1)
% gravitational constant for the earth(G*M_EARTH)[m^3/s^2]
gme = 3.986004418e14; % [m^3/s^2] (REF: 1)
% gravitational constant for the moon (G*M_MOON)[m^3/s^2]
mm  = 7.349e22;
gmm = g*mm;

% nominal earth rotation velocity [rad/sec]
omega = 7.292115e-5;

row = 1025; %[kg/m^3]  density of sea water
% Parameter K_opl is needed for computation of station displacement caused by
% ocean pole tide loading (IERS2010, Ch7.1.5)
Hp = sqrt(8*pi/15)*omega^2*a_tidefree^4/gme; %[m]
K_opl = 4*pi*g*a_tidefree*row*Hp/3/ge_tidefree; %[m]

% convert radians to microarcseconds
rad2uas = 180/pi*3600*1000*1000;

% limits from OCCAM_N.FI
nstations = 20;     % maximum number of station
ntimmax   = 2000;   % maximum number of time epochs, ie number of scans
nobserv   = 30000;  % number of observations


