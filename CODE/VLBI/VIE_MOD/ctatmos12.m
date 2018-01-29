% ************************************************************************
%   Description:
%   Determination of the tidal S1 and S2 atmosphere loading site
%   displacement. Based on Occam subroutine CS1S2.
% 
%   Reference: 
%   Merit Standards, Appendix 11
%
%   Input:										
%      mjd                 Modified Julian Date [d]
%      ant                 station cartesian coordinates from catalogue, 
%                          TRF [m]
%      S12                 cos,sin component of S1 and S2 tidal wave in
%                          Up,North,East system , from S12_ATM.DEF [mm]
%                
%   Output:
%      cta                 station displacement vector (x,y,z) [m]
% 
%   External calls: 	
%      cart2phigd.m, ren2xyz.m   
%
%   Coded for VieVS: 
%   10 Jan 2009 by Hana Spicakova
%
%   Revision: 
%   29 Sep 2009 by Lucia Plank: change phi to geodetic latitude
% *************************************************************************
function cta = ctatmos12(mjd,ant,S12)

UTJ=(mjd-fix(mjd))*2*pi;        % [rad]

dr=S12(1)*cos(UTJ)+S12(2) *sin(UTJ)+S12(3) *cos(UTJ*2)+S12(4) *sin(UTJ*2);
dn=S12(5)*cos(UTJ)+S12(6) *sin(UTJ)+S12(7) *cos(UTJ*2)+S12(8) *sin(UTJ*2);
de=S12(9)*cos(UTJ)+S12(10)*sin(UTJ)+S12(11)*cos(UTJ*2)+S12(12)*sin(UTJ*2);
      
dren=[dr,de,dn]/1000;           % [m]

% Station coordinates from catalogue
X=ant(1);
Y=ant(2);
%Z=ant(3);

% geodetic latitude
phi=cart2phigd(ant); 

lam=atan2(Y,X);

% Correction to station vector
[cta]=ren2xyz(dren,phi,lam);    % [m]
