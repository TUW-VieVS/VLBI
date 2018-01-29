% ************************************************************************
%   Description:
%   This function computes the displacement of the VLBI-antennas caused by
%   pole tide (polar motion). Based on Occam Subroutine CPOLTD.
% 
%   Reference: 
%   IERS Conventions 2003, Ch. 7.1.4
%
%   Input:
%      tim      (7,1)      year,month,day,hour,min,sec,doy           
%      ant                 station cartesian coordinates from catalogue,
%                          TRF [m]
%      xp                  X Wobble [rad]
%      yp                  Y Wobble [rad]
%                
%   Output:
%      ctp                 station displacement vector (x,y,z) [m]
%      flgm_ctp            flagmessage if a linear model instead of a cubic 
%                          one for the mean pole had to be used
%      phpole,plpole       partials w.r.t. the pole tide Love and Shida number
% 
%   External calls: 	
%      xyz2ell.m, meanpole.m, ren2xyz.m   
%
%   Coded for VieVS: 
%   15 Feb 2009 by Hana Spicakova
%
%   Revision:
%   17 Dec 2009 by Lucia Plank: replace xyz2ellip.m with xyz2ell.m
%   30 Mar 2010 by Lucia Plank: use tim as input, jul2dat not needed
%   13 Oct 2010 by Hana Spicakova: bug in the transformation from local
%      system into geocentrical fixed
%   13 Oct 2010 by Hana Spicakova: model for conventional mean pole added
%      according to IERS Conv. 2010 (cubic model over the period 1976-2010)
%   15 Feb 2012 by Hana Spicakova: wrong sign in the south/north displacement
%   15 Feb 2012 by Hana Spicakova: changed latitude to colatitude
%   15 Feb 2012 by Hana Spicakova: partials for pole tide Love and Shida numbers
%      added
%
% *************************************************************************
function [ctp,flgm_ctp,phpole,plpole]=ctpole(tim,ant,xp,yp,ctpm)

% ant_ell - ellipsoial coordinates of antenna (lam,phi,hgt)
[phi,lam] = xyz2ell(ant);
clt = pi/2 - phi; % colatitude (measured from North Pole <0 : 180 deg>)

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
% call meanpole: approximation by a linear trend
[xpm,ypm,flgm_ctp] = meanpole(t,ctpm);        

% m1,m2 : time dependent offset of the instantaneous rotation
%         pole from the mean
m1=  rad2as(xp) - xpm/1000;   % [as] arcsecond
m2=-(rad2as(yp) - ypm/1000);  % [as]


% According to IERS 2010, Ch:7.1.4:
% Using Love number values appropriate to the frequency of the pole tide
% (h2=0.6207, l2=0.0836) and r=a=6378km, displacement vector  for dr,de,dn:
h2=0.6207;
l2=0.0836;

omega = 7.292115e-5; %rad/s
re=6.378e6; %m
g=9.7803278; %m/s^2

dR_m = h2/g*(-omega^2*re^2/2); %m
dR = dR_m*pi/180/3600 ; % m/as

dT_m = 2*l2/g*(-omega^2*re^2/2); %m
dT = dT_m*pi/180/3600 ; % m/as

% dR = -0.033; %IERS 2010
% dT = -0.009; %IERS 2010

dr =  dR*sin(2*clt)*(m1*cos(lam) + m2*sin(lam)); % [m]   
de = -dT*cos(clt)  *(m1*sin(lam) - m2*cos(lam)); % [m] 
dn = -dT*cos(2*clt)*(m1*cos(lam) + m2*sin(lam)); % [m] 

dren=[dr,de,dn];

% Transformation of the displacement vector into geocentric system XYZ
[ctp]=ren2xyz(dren,phi,lam);          % [m]


%% partial derivatives w.r.t. pole love and shida numbers

dhr = 1/g*(-omega^2*re^2/2)*pi/180/3600* sin(2*clt)*(m1*cos(lam) + m2*sin(lam)); %[m]
dle = 2/g*(-omega^2*re^2/2)*pi/180/3600* cos(clt)  *(m1*sin(lam) - m2*cos(lam))*(-1);
dln = 2/g*(-omega^2*re^2/2)*pi/180/3600* cos(2*clt)*(m1*cos(lam) + m2*sin(lam))*(-1);

phren = [dhr,0,0];
plren = [0,dle,dln];

[phpole]=ren2xyz(phren,phi,lam);          % [m]
[plpole]=ren2xyz(plren,phi,lam);          % [m]

phpole=phpole.*100; %[cm]
plpole=plpole.*100; %[cm]
