% ************************************************************************
%   Description:
%   Gives the Matrix needed to transform from the terrestrial to the
%   celestial system.
% 
%   Reference: 
%    - IERS Conventions plus updates (16 June 2009)
%    - MODEST-Handbook from JPL by Sovers and Jacobs (1994)
%    - Notice for the Fortran procedure nro_transf. by A.-M. Gontier for
%      the parital derivatives
%
%   Input:										
%      mjd             Modified Julian Date, observation time UTC [d]
%                
%   Output:
%      t2c(3,3,n)      terrrestrial to celestial matrices             
%
%   External calls: 	
%      tai_utc.m, as2rad.m, xys2000a.m, xys2006a.m, rotm.m, drotm.m
%
%   Coded for VieVS: 
%   03 Jun 2009 by Lucia Plank
%
%   Revision: 
%   26 May 2010 by Lucia Plank: tt instead of mjd for s'
%   23 Feb 2011 by Lucia: modified for vie_sched
% *************************************************************************
function [t2c] = ctrs2crs_sched(mjd)

% global globluc globnee
n = length(mjd);
% constants
p2  = 2*pi;
t2c = zeros(3,3,n);

% prepare time argument
tmu = tai_utc(mjd);
tt  = mjd + (32.184 + tmu)./86400;
tjc = (tt-51544.5)./36525;  % time since J2000 in jul .centuries

% quantity s'
ss = as2rad(-47e-6*tjc);

% earth rotation angle
dut1 = 0;
ut   = mjd + dut1./86400;
tu   = ut - 51544.5;             % days since fundamental epoch
frac = ut - floor(ut) + 0.5;     % UT1 Julian day fraction
fac  = 0.00273781191135448;       
era  = p2 * (frac + 0.7790572732640 + fac * tu );
era  = mod (era,p2);              % [rad]

% prec./nut. X,Y,S
nutmod = 'IAU_2006/2000A';
switch nutmod
    case 'IAU_2000A'
        [X,Y,S] = xys2000a (tt);       %[rad]
    case 'IAU_2006/2000A'
        [X,Y,S] = xys2006a (tt);       %[rad]
end
% % apply nutation corrections of EOP series
%    X = X + (dX');   % [rad]
%    Y = Y + (dY');
   
% build rotational matrices
for i = 1 : n
%     % polar motion matrix:
%     W     = rotm(-ss(i),3)* rotm(xp(i),2)* rotm(yp(i),1);
%     dWdx  = rotm(-ss(i),3)*drotm(xp(i),2)* rotm(yp(i),1);
%     dWdy  = rotm(-ss(i),3)* rotm(xp(i),2)*drotm(yp(i),1);
    % rotation around pole axis
    R     = rotm(-era(i),3);
    dRdut = - drotm(-era(i),3)*(fac+1);
    % precession/nutation matrix:
    E     = atan(Y(i)/X(i));
    d     = asin(Y(i)/sin(E));
    PN    = rotm(-E,3)*rotm(-d,2)*rotm(E,3)*rotm(S(i),3);   
    % whole transformation matrix:
    t2c(:,:,i) = PN * R;
end


