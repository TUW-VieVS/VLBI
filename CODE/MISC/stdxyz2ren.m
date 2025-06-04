% ************************************************************************
%   Description:
%   Transformation of a covariance matrix of station coordinates at a station (lam,phi)
%   from a geocentric coordinate system XYZ into the local system REN
%   XYZ --> REN               REN is positiv towards East and North
% 
%   Reference: 
%   Merit Standards, Appendix 11
%
%   Input:										
%      Qxyz = [dx,dy,dz]   covariance matrix of displacement in the geocentric system XYZ
%      phi                 latitude of the station    [rad]
%      lam                 longitude of the station   [rad]
%                
%   Output:
%      Quen = [dr,de,dn]     displacement in the local system REN 
%
%   Coded for VieVS: 
%   28-02-2023 Helene Wolf
% *************************************************************************
function [mx_uen]=stdxyz2ren(mxyz,phi,lam)


l1=size(mxyz);
l2=size(phi);
l3=size(lam);

l1=l1(1);
l2=l2(1);
l3=l3(1);

eq= isequal(l1,l2,l3);
if eq==0
    fprintf('Number of rows in dxyz, phi and lam must be the same! \n')
end

mx_uen = zeros(l1,3);

for i=1:l1
    Qxyz = diag(mxyz(i,:));
    rot =  [cos(lam(i))*cos(phi(i)), sin(lam(i))*cos(phi(i)), sin(phi(i));
            -sin(lam(i)),  cos(lam(i)), 0;
            -cos(lam(i))*sin(phi(i)), -sin(lam(i))*sin(phi(i)),  cos(phi(i))];

    Quen=rot*Qxyz*rot';

    mx_uen(i,:) = diag(Quen)'; 
end
   
