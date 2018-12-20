% Purpose  
%   Calculate the sidereal time for Greenwich (0--2*pi).
% History  
%   2010-04-06   Jing SUN   Created
%   2016-05-11 Matthias SCHARTNER added mod function to find solution
%   between 0 and 24 instead of loops
%   2016-06-14 Matthias Schartner: vectorised function


function [gmst] = tgmst(mjd)

T = (floor(mjd) - 51544.5) / 36525;
hhs = (mjd - floor(mjd)) * 24;
% time seconds
p1 = 6*3600 + 41*60 + 50.54841;
p2 =            8640184.812866;  
p3 =                  0.093104;
theta = (p1 + p2*T + p3*T.^2) ./ 3600;  % hours

% find solution between 0 and 24 hours
theta = mod(theta,24);

theta = theta * 15 * pi / 180;         % Sidereal time Greenwich in radians at 0 UT
st = hhs*366.2422/365.2422*15*pi/180;  % sidereal time since midnight in radians
theta1 = theta + st;                   % sidereal time Greenwich in radians at the epoch
gmst = mod(theta1, 2*pi);


