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
%      mjd                 Modified Julian Date, observation time UTC [d]
%      xp, yp, dut1        Earth rotation parameters in [rad] resp. [s]
%      dX,dY               Nutation corrections in [rad]
%      nutmod              character: 'IAU_2000A' or IAU_2006A' 
%                
%   Output:
%      t2c    (3,3,n)      terrrestrial to celestial matrices             
%      dQdx   (3,3,n)      partial derivative of t2c w.r.t. pole x
%      dQdy   (3,3,n)      partial derivative of t2c w.r.t. pole y
%      dQdut  (3,3,n)      partial derivative of t2c w.r.t. dut1
%      dQdX   (3,3,n)      partial derivative of t2c w.r.t. celestial X
%      dQdY   (3,3,n)      partial derivative of t2c w.r.t. celestial Y
%      X      (n,1)        total value of celestial pole X [rad]
%                           (X = nutationmodel X + dX)
%      Y      (n,1)        total value of celestial pole Y [rad]
%                           (Y = nutationmodel Y + dY)
%      era    (n,1)        Earth rotation angle [rad]
%
%   External calls: 	
%      tai_utc.m, as2rad.m, xys2000a.m, xys2006a.m, rotm.m, drotm.m
%
%   Coded for VieVS: 
%   03 Jun 2009 by Lucia Plank
%
%   Revision: 
%   26 May 2010 by Lucia Plank: tt instead of mjd for s'
%   22 Jun 2012 by Lucia Plank: fatal sign error in dQdX corrected
%   26 Mar 2021 by Sigrid Boehm: in dMciodY the last element  
%   -Y(i)-X(i)^2*Y(i)/4 was changed to .../2 due to a comment by Axel Nothnagel
% *****************************************************************************
function [t2c,dQdx,dQdy,dQdut,dQdX,dQdY,X,Y,era] = ...
                                      trs2crs(mjd,xp,yp,dut1,dX,dY,nutmod)
% global globluc globnee
n = length(mjd);
% constants
p2    = 2*pi;
t2c   = zeros(3,3,n);
dQdx  = zeros(3,3,n);
dQdy  = zeros(3,3,n);
dQdut = zeros(3,3,n);
dQdX  = zeros(3,3,n);
dQdY  = zeros(3,3,n);


% prepare time argument

tmu = tai_utc(mjd);
tt  = mjd + (32.184 + tmu)./86400;

tjc =(tt-51544.5)./36525;  % time since J2000 in jul .centuries

% quantity s'
ss = as2rad(-47e-6*tjc);

% earth rotation angle

ut   = mjd + dut1./86400;
tu   = ut - 51544.5;             % days since fundamental epoch
frac = ut - floor(ut) + 0.5;     % UT1 Julian day fraction
fac  = 0.00273781191135448;
            
era = p2 * (frac + 0.7790572732640 + fac * tu );
era = mod (era,p2);              % [rad]

% prec./nut. X,Y,S
switch nutmod
    case 'IAU_2000A'
        disp('IAU_2000A');
        [X,Y,S] = xys2000a (tt);       %[rad]
    case 'IAU_2006/2000A'
        disp('IAU_2006A');
        [X,Y,S] = xys2006a (tt);       %[rad]
end
% apply nutation corrections of EOP series
   X = X + (dX');   % [rad]
   Y = Y + (dY');
   
% build rotational matrices
for i=1:n
    % polar motion matrix:
    W     = rotm(-ss(i),3)* rotm(xp(i),2)* rotm(yp(i),1);
    dWdx  = rotm(-ss(i),3)*drotm(xp(i),2)* rotm(yp(i),1);
    dWdy  = rotm(-ss(i),3)* rotm(xp(i),2)*drotm(yp(i),1);
    % rotation around pole axis
    R     = rotm(-era(i),3);
    dRdut = - drotm(-era(i),3)*(fac+1);
    % precession/nutation matrix:
    v     = sqrt(X(i)^2+Y(i)^2);
    E     = atan2(Y(i),X(i));
    z     = sqrt(1-(X(i)^2+Y(i)^2));
    d     = atan2(v,z);
    PN    = rotm(-E,3)*rotm(-d,2)*rotm(E,3)*rotm(S(i),3);
    dEdX  = -Y(i)/(X(i)^2+Y(i)^2);
    dEdY  =  X(i)/(X(i)^2+Y(i)^2);
       
    dddX  =  X(i)/(z*v); 
    dddY  =  Y(i)/(z*v); 
    
    dSdX  = -Y(i)/2;
    dSdY  = -X(i)/2;
    
    dPNdX = drotm(-E,3)* rotm(-d,2)* rotm(E,3)* rotm(S(i),3)* -dEdX + ...
             rotm(-E,3)*drotm(-d,2)* rotm(E,3)* rotm(S(i),3)* -dddX + ...
             rotm(-E,3)* rotm(-d,2)*drotm(E,3)* rotm(S(i),3)*  dEdX + ...
             rotm(-E,3)* rotm(-d,2)* rotm(E,3)*drotm(S(i),3)*  dSdX;
    dPNdY = drotm(-E,3)* rotm(-d,2)* rotm(E,3)* rotm(S(i),3)* -dEdY + ...
             rotm(-E,3)*drotm(-d,2)* rotm(E,3)* rotm(S(i),3)* -dddY + ...
             rotm(-E,3)* rotm(-d,2)*drotm(E,3)* rotm(S(i),3)*  dEdY + ...
             rotm(-E,3)* rotm(-d,2)* rotm(E,3)*drotm(S(i),3)*  dSdY;

        % 
%     a = 1/(1+cos(d));
%     Mcio = [1- a*X(i)^2 , -a*X(i)*Y(i), X(i);...
%             -a*X(i)*Y(i), 1-a*Y(i)^2  , Y(i);...
%             -X(i)       , -Y(i)       , 1-a*(X(i)^2+Y(i)^2)];
%     dMciodX = [-X(i)-X(i)^3/2         , -Y(i)/2-3*X(i)^2*Y(i)/8,   1;...
%                -Y(i)/2-3*X(i)^2*Y(i)/8, 0                      ,   0;...
%                -1                     , 0                      , -X(i)-X(i)^3/2];
%     dMciodY = [-X(i)^2*Y(i)/4   , -X(i)/2-X(i)^3/8   , 0;...
%                -X(i)/2-X(i)^3/8 , -Y(i)-X(i)^2*Y(i)/4, 1;...
%                0                , -1                 , -Y(i)-X(i)^2*Y(i)/2];

  
%     dPNdXneu= dMciodX * rotm(S(i),3) + Mcio*drotm(S(i),3)*dSdX
%     dPNdX
%     dPNdYneu= dMciodY * rotm(S(i),3) + Mcio*drotm(S(i),3)*dSdY
%   %
%     
%     dPNdY 
    
    % whole transformation matrix:
    t2c(:,:,i) = PN * R * W;
    
    % partial derivatives
    dQdx(:,:,i)  = PN    *     R * dWdx;
    dQdy(:,:,i)  = PN    *     R * dWdy;
    dQdut(:,:,i) = PN    * dRdut * W;
    dQdX(:,:,i)  = dPNdX *     R * W;
    dQdY(:,:,i)  = dPNdY *     R * W;
end
