% ************************************************************************
%   Description:
%   Computes local azimuth, zenith distance; a_grad: preparation for 
%   partial derivatives of atmosphere gradient, and correction to zenith 
%   distance due to tropospheric refraction. Translated from occam function
%   station.f.
% 
%   Reference: 
%      ---
%
%   Input:										
%      lam, phi           antenna ellipsoidal coordinates   [rad]
%      rqu        (1,3)   source vector in catalogue system [rad] 
%      t2c        (3,3)   transformation matrix TRS --> CRS       
%                
%   Output:
%      az                 local azimuth [rad]
%      zd                 local zenith distance   [rad]
%      a_grad             for partial derivatives of atmosphere gradient[1]
%      corz               correction to zenith distance due to trop. 
%                         refract. [rad]
%      de                 local declination [rad]
% 
%   External calls: 	
%      ---   
%
%   Coded for VieVS: 
%   10 Oct 2009 by Hana Spicakova
%
%   Revision: 
%   4 Oct 2016 by H. Krasna: Local Hour Angle east of the meridian added
%   
% *************************************************************************
function [az,zd,a_grad,corz,de,LHAe]=locsource(lam,phi,rqu,t2c)
  
   % source in TRS (c2t = t2c')
   rq_trs = t2c'*rqu';
   rq     = rq_trs/norm(rq_trs);

   % source in local system 
   my  = [1,0,0;0,-1,0;0,0,1];  % mirror matrix along y-axis
   g2l = my*rotm((pi/2-phi),2)*rotm(lam,3);

   lq = g2l*rq;
   
   % zenith distance
   zd = acos(lq(3));
   
   % azimuth
   saz = atan2(lq(2),lq(1));        % south azimuth
     if (lq(2) < 0) 
         saz =  2*pi + saz;
     end
   az  = saz+pi;
   az  = mod(az,2*pi);              % north azimuth
  
   de     = atan2(rq(3),sqrt((rq(1))^2 + (rq(2))^2)); 
   
   % Local Hour Angle east!!! of the meridian = Greenwich HA east of the meridian - station longitude 
   LHAe = atan2(rq(2),rq(1)) - lam; % [rad]

% for partial derivatives of atmosphere gradient
    a_grad(1)=cos(az)*tan(zd);
    a_grad(2)=sin(az)*tan(zd);

%% Correction to zenith distance due to tropospheric refraction
%  in OCCAM as option from Tesmer
%  corz=refract. factor / tangent of elevation
%  Default: YESTRC = 'Y'

corz = 313e-6*tan(zd); % Modest Handbook, p. 42, (2.200)
if corz > pi/2
   corz=-corz;
end

%zdr = zd-corz;
%lq = [sin(zdr)*cos(saz);sin(zdr)*sin(saz);cos(zdr)];
%rqref = t2c*(g2l'*lq);
