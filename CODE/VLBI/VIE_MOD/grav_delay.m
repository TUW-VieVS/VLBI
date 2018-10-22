% ************************************************************************
%   Description:
%   Calculates the gravitational delay due to celestial bodies.
% 
%   Reference: 
%   IERS Conventions Chapt. 11
%
%   Input:										
%      X1,X2    (3,1)     barycentric station vector [m]
%      vearth             barycentric earth velocity
%      b                  baseline GCRS
%      rqu      (3,1)     unit source vector barycentrum-source [m]
%      ephem    (struct)  ephemerides 
%      isc                number of scan
%                
%   Output:
%      Tgrav              gravitaional delay
% 
%   External calls: 	 
%
%   Coded for VieVS: 
%   25 Aug 2009 by Lucia Plank
%
%   Revision:
%   05 Aug 2010 by Lucia Plank: time of closest approach corrected;
%       calculate Tgrav for all bodies (sun, moon, planets)
%   25 Sep 2012 by Hana Krásná: partials w.r.t. relativistic parameter
%          gamma added (only Sun's contribution)
% *************************************************************************
function [Tgrav, pGammaSun] = grav_delay(X1,X2,vearth,b,rqu,ephem,isc,opt)

global c 
global gms gmm 


Tgrav = 0;
for body = 1 : 10 
    switch body
        case 1
            Xbody_bar = ephem.sun(isc).xbar;
            Vbody_bar = ephem.sun(isc).vbar;
            gg = gms;
        case 2
            Xbody_bar = ephem.moon(isc).xbar;
            Vbody_bar = ephem.moon(isc).vbar;
            gg = gmm;
        case 3
            Xbody_bar = ephem.merc(isc).xbar;
            Vbody_bar = ephem.merc(isc).vbar;
            gg = ephem.gmmerc;
        case 4
            Xbody_bar = ephem.venu(isc).xbar;
            Vbody_bar = ephem.venu(isc).vbar;
            gg = ephem.gmvenu;
        case 5
            Xbody_bar = ephem.mars(isc).xbar;
            Vbody_bar = ephem.mars(isc).vbar;
            gg = ephem.gmmars;
        case 6
            Xbody_bar = ephem.jupi(isc).xbar;
            Vbody_bar = ephem.jupi(isc).vbar;
            gg = ephem.gmjupi;
        case 7
            Xbody_bar = ephem.satu(isc).xbar;
            Vbody_bar = ephem.satu(isc).vbar;
            gg = ephem.gmsatu;
        case 8
            Xbody_bar = ephem.uran(isc).xbar;
            Vbody_bar = ephem.uran(isc).vbar;
            gg = ephem.gmuran;
        case 9
            Xbody_bar = ephem.nept(isc).xbar;
            Vbody_bar = ephem.nept(isc).vbar;
            gg = ephem.gmnept;
        case 10
            Xbody_bar = ephem.plut(isc).xbar;
            Vbody_bar = ephem.plut(isc).vbar;
            gg = ephem.gmplut;
    end
            
 % time of closest approach
 dt = norm(Xbody_bar-X1)/c;
xbody = Xbody_bar-dt*Vbody_bar; 
 dt = norm(xbody-X1)/c;
xbody = Xbody_bar-dt*Vbody_bar; 
 dt = norm(xbody-X1)/c;
xbody = Xbody_bar-dt*Vbody_bar;

    R1b = X1 - xbody;                                    % (eq. 4)
    R2b = X2 - vearth/c *(rqu*b) - xbody;                % (eq. 5)
    gamma = 1;
    Tgrav = Tgrav + (1+gamma)*gg/(c^3)*...
            log((norm(R1b)+rqu*R1b)/(norm(R2b)+rqu*R2b));    % (eq. 1)
end



%% partial derivative w.r.t. gamma (Sun only)
pGammaSun=0;
if opt.est_gamma ==1
    Xbody_bar = ephem.sun(isc).xbar;
    Vbody_bar = ephem.sun(isc).vbar;
    gg = gms;

     % time of closest approach
     dt = norm(Xbody_bar-X1)/c;
    xbody = Xbody_bar-dt*Vbody_bar; 
     dt = norm(xbody-X1)/c;
    xbody = Xbody_bar-dt*Vbody_bar; 
     dt = norm(xbody-X1)/c;
    xbody = Xbody_bar-dt*Vbody_bar;

    R1b = X1 - xbody;                                    % (eq. 4)
    R2b = X2 - vearth/c *(rqu*b) - xbody;                % (eq. 5)


    pGammaSun = gg/(c^3)* log((norm(R1b)+rqu*R1b)/(norm(R2b)+rqu*R2b));
end

        