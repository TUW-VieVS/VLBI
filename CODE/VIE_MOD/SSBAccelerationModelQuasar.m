% ************************************************************************
%   Description:
%   function to calculate the time delay following the  calculation of the time delay following the consensus model + Titov 2011
%
%   Reference:
%
%   Input:
%       'isc'           (1,1)             index of current scan 
%       'crsStation1'   (3,1)             CRS coordiantes of station 1
%       'crsStation2'   (3,1)             CRS coordinates of station 2
%       'rqu'           (1,3)             unit source vector barycentrum-source
%       'ephem'         structure array   Ephermerides
%       'v1'            (3,1)             velocity of station 1
%       'v2'            (3,1)             velocity of station 2
%       'opt'           structure array   (for info. /DOC/opt.doc)
%       'delt_accSSB'   (1,1)             time since reference epoch for SSB acceleration [sec] 
%
%
%   Output:
%       'tau'           (1,1)             time delay 
%       'pGammaSun'     (1,1)             partial derivative of gravitaional delay w.r.t. gamma (Sun only)
%       'k1a'           (1,3)             aberrated source vector
%       'k2a'           (1,3)             aberrated source vector
%       'fac1'          (1,1)             Correction to the scale factor [sec]
%
%   External calls:
%        grav_delay.m     
%
%   Coded for VieVS:
%   22 November 2021 by H. Wolf - created as external function
%
%   Revision:
%
% ************************************************************************


function [tau, pGammaSun, k1a, k2a, fac1] = SSBAccelerationModelQuasar(isc, crsStation1, crsStation2, rqu, ephem, v1, v2 , opt, delt_accSSB)

    global c gms gme

    % acceleration of SSB - Galactocentric aberration
    % GA value:
    GA=5.8;                          %[uas/y]
    RAGC = (17+45/60+40/3600)/12*pi; %[rad] RA of Galactic center: 17h45min40sec
    DeGC = (-29-28/3600)/180*pi;     %[rad]  De of Galactic center: -29deg 00min 28sec

    % num of sec in one year
    GA_T = 60*60*24*365.25;
    % conversion rad -> as
    GA_k=180/pi*60*60 *GA_T; % as*s
    accMag = GA* 1e-6 *c  / GA_k ; % m/s^2 magnitude of acceleration (~ 2.44e-10 m/s^2)

    acc=[accMag*cos(RAGC)*cos(DeGC)
         accMag*sin(RAGC)*cos(DeGC)
         accMag*sin(DeGC)]; % [m/sec^2]


    crsB = crsStation2 - crsStation1;     % baseline vector CRS
    % Ephemerids:
    earth   = ephem.earth(isc).xbar;
    vearth  = ephem.earth(isc).vbar;
    sun     = ephem.sun(isc).xgeo;

    %(1) barycentric station vector (eq. 6)
    xb1 = earth + crsStation1;
    xb2 = earth + crsStation2;

    %(2)&(3) differential gravitational delay due to celestial bodies
    [Tgrav, pGammaSun] = grav_delay(xb1,xb2,vearth,crsB,rqu,ephem,isc,opt);

    %(4) differential gravitational delay due to the earth
    Tgrave = 2*gme/c^3*...
        log((norm(crsStation1)+rqu*crsStation1)/(norm(crsStation2)+rqu*crsStation2)); % (eq. 1)

    %(5) total differential gravitational delay (eq. 2)
    Tgrav = Tgrave + Tgrav;

    %(6) vacuum delay
    U     = gms/norm(sun); % gravitational potential at the geocenter
    gamma = 1;
    fac1  = (rqu*crsB)/c;
    term1 = 1-(1+gamma)*U/c^2-(norm(vearth))^2/(2*c^2)-(vearth'*v2)/c^2;
    fac2  = ((vearth+acc*delt_accSSB)'*crsB)/c^2;
    term2 = 1+(rqu*vearth)/(2*c);
    quot  = 1+(rqu*(vearth+v2+acc*delt_accSSB))/c;
    tau   = (Tgrav - fac1*term1 - fac2*term2)/quot;          %(eq. 9)

    pGammaSun = pGammaSun/quot;

    %(7) aberrated source vector (eq. 15)
    k1a = rqu + (vearth+v1)'/c - rqu*((rqu*(vearth+v1))')/c;
    k2a = rqu + (vearth+v2)'/c - rqu*((rqu*(vearth+v2))')/c; 
end

