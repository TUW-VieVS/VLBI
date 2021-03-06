% Corrections to the station coordinates caused by ocean pole tide loading.
% The ocean pole tide is the ocean response to the variation of both the solid
% Earth and the oceans to the centrifugal potential that is generated by small
% perturbations to the Earth's rotation axis.
% IERS Conventions 2010, Chapter 7.1.5
% The ocean pole load tide coefficients are taken from Desai (2002):
% ftp://tai.bipm.org/iers/conv2010/chapter7/opoleloadcoefcmcor.txt.gz
%
%   Input:										
%      tim      (7,1)      year,month,day,hour,min,sec,doy
%      ant                 antenna coordinates xyz [m]
%      xp                  X Wobble [rad]
%      yp                  Y Wobble [rad]
%      opl                 ocean pole load tide coefficients from
%                          Desai (2002): u_r^R, u_r^I, u_n^R, u_n^I, u_e^R, u_e^I 
% 
%   Output:
%     ctop            
%     flgm_ctp
%
%   External calls: 	
%      meanpole.m                 					    											
%       
%   Coded for VieVS: 
%   10 Nov 2011 by Hana Spicakova
%%

function [ctop,flgm_ctp]=ctoceanpole(tim,ant,xp,yp,opl,ctpm)

global K_opl

% ant_ell - ellipsoial coordinates of antenna (lam,phi,hgt)
[phi,lam] = xyz2ell(ant);

yr  = tim(1);
doy = tim(7);
h   = tim(4);
m   = tim(5);
s   = tim(6);

IFAC = 365;
if (mod(yr,4)==0)
    IFAC=366;
end
t = yr + (doy + h/23.93447 + m/1440 + s/86400) / ((IFAC) + 0.2422);

% xpm, ypm: mean pole coordinates [mas]
[xpm,ypm,flgm_ctp] = meanpole(t,ctpm);        

% m1,m2 : time dependent offset of the instantaneous rotation
%         pole from the mean
m1=  rad2as(xp) - xpm/1000;   % [as] arcsecond
m2=-(rad2as(yp) - ypm/1000);  % [as]

% values from IERS Conv. 2010, Ch 7.1.5
% (appropriate to k2 and h2 for pole tide)
ga2R = 0.6870; %[-]
ga2I = 0.0036; %[-]

m1m2ga2 = m1*ga2R + m2*ga2I; %[as]
m2m1ga2 = m2*ga2R - m1*ga2I; %[as]

m1m2ga2 = m1m2ga2/3600/180*pi; %[rad]
m2m1ga2 = m2m1ga2/3600/180*pi; %[rad]

oplR = opl([1 3 5]); %[-]
oplI = opl([2 4 6]); %[-]

uRNE=K_opl*(m1m2ga2.*oplR + m2m1ga2.*oplI); 
uREN=[uRNE(1) uRNE(3) uRNE(2)]; % [m]

% Transformation of the displacement vector into geocentric system XYZ
[ctop]=ren2xyz(uREN,phi,lam);          % [m]


