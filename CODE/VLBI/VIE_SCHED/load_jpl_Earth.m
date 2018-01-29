% ************************************************************************
%   Description:
%   This function is only a part of load_eph.m to get the velocity of the 
%   Earth.
%   
%   Use load_eph.m to use all functions!
% 	2016-07-11 M. Schartner: created
%
% *************************************************************************
function [ jpl ] = load_jpl_Earth( jplnum )
      
% unit conversion
posu = 1e3;         % km --> m
velu = 1e3/86400;   % km/day --> m/s

% calculate geocentric positions
geo = 1;

% load ephemerides file
load(strcat('../EPHEM/',jplnum,'.mat'))

jpl.gm = [jpl.gm([3 10])];
jpl.order = 'Earth-Moon_Barycenter Moon';
jpl.nset = [jpl.nset([3 10])];
jpl.ncoeff = [jpl.ncoeff([3 10])];
removeField = {'merc','venu','mars','jupi','satu','uran','nept','plut','sun'};
jpl.rec = rmfield(jpl.rec,removeField);

end

