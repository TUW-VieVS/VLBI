% ************************************************************************
%   Description:
%   Transformation of a diplacement vector at a station (lam,phi)
%   from a geocentric coordinate system XYZ into the local system REN
%   XYZ --> REN               REN is positiv towards East and North
% 
%   Reference: 
%   Merit Standards, Appendix 11
%
%   Input:										
%      dxyz = [dx,dy,dz]   displacement in the geocentric system XYZ
%      phi                 latitude of the station    [rad]
%      lam                 longitude of the station   [rad]
%           It is possible to do the computation for more station.
%           e.q. 2 stations:
%                       dxyz = [5 4 6;
%                               2 4 3];
%                       phi = [0;   1.5];
%                       lam = [3.14;  0];
%                
%   Output:
%      dr = [dr,de,dn]     displacement in the local system REN 
%
%   Coded for VieVS: 
%   23 Mar 2009 by Hana Spicakova
%
%   Revision: 
%   20 Jan 2010 by Hana Spicakova
% *************************************************************************
function [dren]=xyz2ren(dxyz,phi,lam)

l1=size(dxyz);
l2=size(phi);
l3=size(lam);

l1=l1(1);
l2=l2(1);
l3=l3(1);

eq= isequal(l1,l2,l3);
if eq==0
    fprintf('Number of rows in dxyz, phi and lam must be the same! \n')
end


for i=1:l1
    phi1=phi(i);
    lam1=lam(i);
    rot =  [cos(phi1)*cos(lam1), -sin(lam1), -sin(phi1)*cos(lam1);
            cos(phi1)*sin(lam1),  cos(lam1), -sin(phi1)*sin(lam1);
            sin(phi1)          ,          0,  cos(phi1)          ];
        
    dren(i,:)=(inv(rot)*dxyz(i,:)')';
end    
   
