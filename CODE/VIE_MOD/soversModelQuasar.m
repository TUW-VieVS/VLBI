% ************************************************************************
%   Description:
%   function to calculate the time delay following sovers model
%   influence of planets not included yet
%
%   Reference:
%
%   Input:
%       'isc'           (1,1)             index of current scan 
%       'crsStation1'   (3,1)             CRS coordiantes of station 1
%       'crsStation2'   (3,1)             CRS coordinates of station 2
%       'rqu'           (1,3)             unit source vector barycentrum-source
%       'ephem'         structure array   Ephermerides
%       'v2'            (3,1)             velocity of station 2

%   Output:
%       'tau'           (1,1)             time delay 
%       'pGammaSun'     (1,1)             partial derivative of gravitaional delay w.r.t. gamma (Sun only)
%       'k1a'           (1,3)             aberrated source vector
%       'k2a'           (1,3)             aberrated source vector
%       'fac1'          (1,1)             Correction to the scale factor [sec]
%
%   External calls:
%
%   Coded for VieVS:
%   22 November 2021 by H. Wolf - created as external function of vie_mod
%
%   Revision:
%
% ************************************************************************

function [tau, pGammaSun, k1a, k2a, fac1] = soversModelQuasar(isc, crsStation1, crsStation2, rqu, ephem, v2)
    
    global gme c gms

    crsB = crsStation2 - crsStation1;     % baseline vector CRS
    % Ephemerids:
    earth   = ephem.earth(isc).xbar;
    vearth  = ephem.earth(isc).vbar;

    % input
    bet2s = v2/c;        % velocity st2 GCRS
    bet   = vearth/c;    % SBB
    bs    = crsB;       % baseline GCRS
    k     = -rqu;        % source vector SBB (1,3)

    %  transform to SBB frame
    gam  = 1/sqrt(1-bet'*bet);
    bet2 = (bet2s + (gam-1)*(bet2s'*bet)*bet/norm(bet)^2+gam*bet)/...
        (gam*(1 + bet2s'*bet)); % v2 SBB
    crsB    = bs +(gam-1)*(bs'*bet)*bet/norm(bet)^2 -...
        gam*bet2*(bs'*bet); % baseline SBB

    % calculate in SBB frame
    tausbb = (k*crsB)/(1-k*bet2);
    bt1t2  = crsB + bet2*tausbb;
    
    % transform back to GCRS frame
    tau = gam*tausbb - gam*(bt1t2'*bet);   % [m]
    tau = tau/c;                           % [sec]
    
    % gravitational correction (Sun)
    Rr   = ephem.sun(isc).xbar; % position SBB at reception time
    Vr   = ephem.sun(isc).vbar;
    trta = norm(earth)/c;
    Rn   = Rr - trta*Vr/c;
    
    % iterative solution for position at time of closest approach
    Rn = Rr - norm(Rn)*Vr/c;
    Rn = Rr - norm(Rn)*Vr/c;
    Rn = Rr - norm(Rn)*Vr/c;
    Ra = Rn;

    xb1 = earth + crsStation1;
    xb2 = earth + crsStation2;
    r1 = xb1 - Ra;
    r2 = xb2 + bet2*tausbb - Ra;

    dGp  = 2*gms/c^3*log((norm(r1)+rqu*r1)/(norm(r2)+rqu*r2));
    U    = gms/(norm(ephem.sun(isc).xgeo)*c^2);
    dGps = dGp - 2*U*tau;

    % differential gravitational delay due to the earth
    st2t2s = crsStation2 + tau*bet2s;
    dGpe = 2*gme/c^3*...
        log((norm(crsStation1)+rqu*crsStation1)/(norm(st2t2s)+rqu*st2t2s));

    tau = tau + dGps + dGpe;

    % aberrated source vector (yearly aberration)
    s0d  = rqu;
    coas = 1/(gam*(1+bet'*s0d'));
    coab = coas*((gam-1)*(bet'*s0d')/(bet'*bet)+gam);
    k1a  = (coas*s0d'+coab*bet)';
    k2a  = k1a;
    
    pGammaSun = 0;
    fac1 = 0;
end

