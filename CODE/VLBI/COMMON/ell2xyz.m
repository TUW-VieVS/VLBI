% ************************************************************************
%   Description:
%   Transformation from ellipsoidal coordinates lam,phi,elh to Cartesian coordinates X,Y,Z
%
%   Input:										
%      coor_ell = [lat,lon,h]      [rad,rad,m]
%               can be a matrix: number of rows = number of stations
%
%   Output:
%      pos = [x,y,z]                 [m,m,m]
%
%   Coded for VieVS: 
%   3 Oct 2012 by Benedikt Soja
% *************************************************************************
function pos=ell2xyz(coor_ell)


a = 6378137;  %m      Equatorial radius of the Earth
f = 0.003352810681225;      % Flattening factor of the Earth

e2=2*f-f^2;
ei2=e2/(1-f)^2;

c=a/(1-f);
V=sqrt(1+ei2*cos(coor_ell(:,1)).^2);

x=(c./V+coor_ell(:,3)).*cos(coor_ell(:,1)).*cos(coor_ell(:,2));
y=(c./V+coor_ell(:,3)).*cos(coor_ell(:,1)).*sin(coor_ell(:,2));
z=((1-e2)*c./V+coor_ell(:,3)).*sin(coor_ell(:,1));

pos = [x y z];