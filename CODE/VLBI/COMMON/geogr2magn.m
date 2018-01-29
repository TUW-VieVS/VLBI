% ************************************************************************
%   Description:
%   Transformation from Cartesian geographic x,y,z to Cartesian geomagnetic
%   coordinates X,Y,Z
%
%   Input:										
%      pos = [x,y,z]                 [m,m,m]
%               can be a matrix: number of rows = number of stations
%
%   Output:
%      pos_magn = [X,Y,Z]            [m,m,m]
%
%   Coded for VieVS: 
%   3 Oct 2012 by Benedikt Soja
% *************************************************************************
function pos_magn = geogr2magn(pos)

%geocentric lat/lon of the dipole North geomagnetic pole
lam_mp = 288.22*pi/180;
phi_mp = 79.75*pi/180;

%Rotation matrix
R = [cos(pi/2-phi_mp) 0 -sin(pi/2-phi_mp); 0 1 0; sin(pi/2-phi_mp) 0 cos(pi/2-phi_mp)]*...
    [cos(lam_mp) sin(lam_mp) 0; -sin(lam_mp) cos(lam_mp) 0; 0 0 1 ];

pos_magn = (R*pos')';