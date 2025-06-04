% ************************************************************************
%   Description:
%	Computes the partial derivative of the time delay tau w.r.t. the 
%   orbital elements.
%
%   Input:										
%     time_UTC              time of observation
%     sources               sources struct
%     r                     position of satellite
%     v                     velocity of satellite
%
% 
%   Output:
%     drdpar               partial derivative (6x3) 
% 								line 1: partial derivative w.r.t. semi-major axis (1x3)
%                               line 2: partial derivative w.r.t. eccentricity (1x3)
% 								line 3: partial derivative w.r.t. inclination (1x3)
%                               line 4: partial derivative w.r.t. Omega (right ascension of ascending node) (1x3)								
%								line 5: partial derivative w.r.t. omega (argument of perigee) (1x3)
%                               line 6: partial derivative w.r.t. argument of latitude (1x3)
% 
%   External calls: 	
%
%       
%   Coded for VieVS: 
%   17 Dec 2024 by Helene Wolf
%
%   Revision: 
%
%
% ************************************************************************

function [ drdpar ] = drdorb_num(GM, time_UTC, tosc, r , v)

    [a,e,i,Omega,omega,t0] = xyzele(GM,time_UTC,r,v);
    u0 = getu0(GM, a,e,omega,t0,tosc);

    dele=[1,0.000777,mas2rad(100),mas2rad(100),deg2rad(3),mas2rad(100)];

    drdpar = zeros(6,3);

    t00p = gett0(GM,a+dele(1)/2,e,omega,tosc,u0);
    t00m = gett0(GM,a-dele(1)/2,e,omega,tosc,u0);
    rp = getPositionFromKeplerianElements(GM,a+dele(1)/2,e,i,Omega,omega,time_UTC,t00p);
    rm = getPositionFromKeplerianElements(GM,a-dele(1)/2,e,i,Omega,omega,time_UTC,t00m);
    drdpar(1,:) = (rp-rm)./dele(1);

    t00p = gett0(GM,a,e+dele(2)/2,omega,tosc,u0);
    t00m = gett0(GM,a,e-dele(2)/2,omega,tosc,u0);
    rp = getPositionFromKeplerianElements(GM,a,e+dele(2)/2,i,Omega,omega,time_UTC,t00p);
    rm = getPositionFromKeplerianElements(GM,a,e-dele(2)/2,i,Omega,omega,time_UTC,t00m);
    drdpar(2,:) = (rp-rm)./dele(2);

    t00 = gett0(GM,a,e,omega,tosc,u0);
    rp = getPositionFromKeplerianElements(GM,a,e,i+dele(3)/2,Omega,omega,time_UTC,t00);
    rm = getPositionFromKeplerianElements(GM,a,e,i-dele(3)/2,Omega,omega,time_UTC,t00);
    drdpar(3,:) = (rp-rm)./dele(3);

    t00 = gett0(GM,a,e,omega,tosc,u0);
    rp = getPositionFromKeplerianElements(GM,a,e,i,Omega+dele(4)/2,omega,time_UTC,t00);
    rm = getPositionFromKeplerianElements(GM,a,e,i,Omega-dele(4)/2,omega,time_UTC,t00);
    drdpar(4,:) = (rp-rm)./dele(4);

    t00p = gett0(GM,a,e,omega+dele(5)/2,tosc,u0);
    t00m = gett0(GM,a,e,omega-dele(5)/2,tosc,u0);
    rp = getPositionFromKeplerianElements(GM,a,e,i,Omega,omega+dele(5)/2,time_UTC,t00p);
    rm = getPositionFromKeplerianElements(GM,a,e,i,Omega,omega-dele(5)/2,time_UTC,t00m);
    drdpar(5,:) = (rp-rm)./dele(5);
    
    t00p = gett0(GM,a,e,omega,tosc,u0+dele(6)/2);
    t00m = gett0(GM,a,e,omega,tosc,u0-dele(6)/2);
    rp = getPositionFromKeplerianElements(GM,a,e,i,Omega,omega,time_UTC,t00p);
    rm = getPositionFromKeplerianElements(GM,a,e,i,Omega,omega,time_UTC,t00m);
    drdpar(6,:) = (rp-rm)./dele(6);
end