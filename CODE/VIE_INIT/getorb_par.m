function [SatPos,SatVel] = getorb_par(prn,tvec,coeff,t0,hstep,tbound,satnum)
% ------------------------------------------------------------------------
% [pos,vel] = getorb(prn,t,coeff,t0,hstep,tbound,satnum);
%
% Purpose: Read Formatted STD File
%
% Input:   filname        - GNSS precise orbit file name
%          prn            - satellite number
%          t              - epoch of pos, vel (mjd)
%          coeff(i,j,k,m) - Polynomial coefficients
%                           i=1,...,nint - integration interval
%                           j=1,...,nq+1 - degree
%                           k=1,2,3,     - component
%                           m=1,...,nsat - satellite nr
%          t0(i)          - epoch expansion point (mjd)
%                           i=1,...,nint - integration interval
%          hstep(i)       - integraiton interval (sec)
%                           i=1,...,nint - integration interval
%          tbound(i)      - integration interval boundaries
%                           i=1,...,nint+1
%          satnum(i)      - satellite number
%                           i=1,...,nsat - satellite nr
%                            nn: GPS
%                           1nn: GLONASS
%                           2nn: Galileo
%                           3nn: SBAS
%                           4nn: BeiDou
%                           5nn: QZSS
% Output:  pos(k)         - Position of satellite prn in ECI (m)
%          vel(k)         - Velocity of satellite prn in ECI (m/s)
%
% Author:  Urs Hugentobler, FESG
% Date  :  27-11-2022
% Changes: 
% -------------------------------------------------------------------------
format compact;
SatPos = zeros(length(tvec),3);
SatVel = zeros(length(tvec),3);
% find satellite
isat = find(satnum==prn);
if (isempty(isat)); error(sprintf('*** Satellite %d not found',prn)); end

% find integration interval
eps = 2/1440;
for idx = 1:length(tvec)
    t = tvec(idx);
    if (t < tbound(1)-eps)
        error(sprintf('*** Epoch earlier than first epoch %f',tbound(1)))
    end
    if (t > tbound(end)+eps)
        error(sprintf('*** Epoch later than last epoch %f',tbound(end)))
    end
    nint = size(coeff,1);
    ix   = find(tbound-eps < t);
    int  = min(ix(end),nint);
    
    % polynomial argument
    ti = (t-t0(int))*86400/hstep(int);
    
    vel = 0;
    
    % evaluation of polynomial
    nq = size(coeff,2)-1;
    % position
    pos = squeeze(coeff(int,end,:,isat));
    for iq=nq:-1:1
        pos = pos.*ti+ squeeze(coeff(int,iq,:,isat));
    end
    % velocity
    vel = squeeze(coeff(int,end,:,isat))*nq;
    for iq=nq:-1:2
        vel = vel.*ti + squeeze(coeff(int,iq,:,isat))*(iq-1);
    end
    vel = vel/hstep(int);
    SatPos(idx,:) = pos';
    SatVel(idx,:) = vel';
end