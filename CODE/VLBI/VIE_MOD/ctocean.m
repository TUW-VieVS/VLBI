% ************************************************************************
%   Description:
%   Ocean loading site displacement. Computes the station height
%   displacement for ocean loading (9 tides). Based on Occam subroutine
%   COCLD. 
% 
%   Reference: 
%   Merit Standards, April 81, second draft, appendices 7,11 and
%   C. C. Goads, J. Geoph. Res., 85, p. 2679 - 2683, May 1980
%
%   Input:	
%      mjd              Modified Julian Date [d]
%      leap             difference between UTC and TAI in [s]
%                       leap seconds(TAI-UTC))
%      ant              station cartesian coordinates from catalogue,TRF[m]
%      cto_data         apriori tidal ocean loading amplitudes and phases,
%                       from 'LOADING.OCE' [m,deg]
%                
%   Output:
%      cto              station displacement vector (x,y,z) [m]
% 
%   External calls: 	
%      ang_schwi.m, cart2phigd.m, ren2xyz.m                 					    											
%       
%   Coded for VieVS: 
%   10 Jan 2009 by Hana Spicakova
%
%   Revision: 
%   29 Sep 2009 by Lucia Plank: change phi to geodetic latitude
% *************************************************************************
function [cto]=ctocean(mjd,leap,ant,cto_data)
deg2pi=pi/180;

% Amplitude of the ocean tidal wave in REN system
rA=cto_data(1,:);
eA=cto_data(2,:);
nA=cto_data(3,:);

% Phase of the ocean tidal wave in REN system
rphase=cto_data(4,:)*deg2pi;
ephase=cto_data(5,:)*deg2pi;
nphase=cto_data(6,:)*deg2pi;

[angle]=ang_schwi(mjd,leap);

% Height displacement
dr=rA'.*cos(angle-rphase');
de=eA'.*cos(angle-ephase');         % + west  ??
dn=nA'.*cos(angle-nphase');         % + south ??

% Change the sign of horizontal contribution
dren=[dr,-de,-dn];
dren=sum(dren);

X=ant(1);
Y=ant(2);
%Z=ant(3);

%geodetic latitude
%phi=atan(Z/(sqrt(X^2+Y^2)));
phi=cart2phigd(ant); 

lam=atan2(Y,X);

[cto]=ren2xyz(dren,phi,lam);

