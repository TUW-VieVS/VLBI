% ************************************************************************
%   Description:
%   function to calculate the time delay following the consensus model
%
%   Reference:
%
%   Input:
%       'isc'           (1,1)             index of current scan 
%       'crsStation1'   (3,1)             CRS coordiantes of station 1
%       'crsStation2'   (3,1)             CRS coordinates of station 2
%       'idStation1'    (1,1)             index of Station 1
%       'idStation2'    (1,1)             index of Station 2
%       'rqu'           (1,3)             unit source vector barycentrum-source
%       'ephem'         structure array   Ephermerides
%       'opt'           structure array   (for info. /DOC/opt.doc)
%       'antenna'       structure array   antenna structure 
%       'v1'            (3,1)             velocity of station 1
%       'v2'            (3,1)             velocity of station 2
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
%   22 November 2021 by H. Wolf - created as external function of vie_mod
%
%   Revision:
%
% ************************************************************************

function [tau, pGammaSun, k1a, k2a, fac1] = consensusModelQuasar(isc, crsStation1, crsStation2, idStation1, idStation2, rqu, ephem, opt, antenna, v1, v2)
    
    global gme c gms    
    
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
    if strcmp(antenna(idStation1).name,'GEOCENTR') || strcmp(antenna(idStation2).name,'GEOCENTR')
        Tgrav=0;
    end

    %(6) vacuum delay
    U     = gms/norm(sun); % gravitational potential at the geocenter
    gamma = 1;
    fac1  = (rqu*crsB)/c;
    term1 = 1-(1+gamma)*U/c^2-(norm(vearth))^2/(2*c^2)-(vearth'*v2)/c^2;
    fac2  = (vearth'*crsB)/c^2;
    term2 = 1+(rqu*vearth)/(2*c);
    quot  = 1+(rqu*(vearth+v2))/c;
    tau   = (Tgrav - fac1*term1 - fac2*term2)/quot;          %(eq. 9)

    pGammaSun = pGammaSun/quot;

    %(7) aberrated source vector (eq. 15)
    k1a = rqu + (vearth+v1)'/c - rqu*((rqu*(vearth+v1))')/c;
    k2a = rqu + (vearth+v2)'/c - rqu*((rqu*(vearth+v2))')/c;
end

