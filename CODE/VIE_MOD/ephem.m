% ************************************************************************
%   Description:
%   COMPUTE POSITION X AND VELOCITY XP OF A CELESTIAL
%   BODY WHOSE OSCULATING ELEMENTS ARE GIVEN
%   (NO PARABOLAS, HYPERBOLAS)
%
%   Input:										
%              GM     : GRAVITY - CONSTANT                  R*8
%               A      : SEMIMAJOR AXIS                      R*8
%               E      : NUMERICAL EXCENTRICITY              R*8
%               I      : INCLINATION (RADIAN)                R*8
%               KN     : RIGHT ASCENSION OF ASCENDING NODE   R*8
%                        (RADIAN)
%               PER    : ARGUMENT OF PERICENTRE (RADIAN)     R*8
%               T0     : PERICENTRE-PASSING-TIME             R*8
%               T      : EPOCH OF EPHEMERIS COMPUTATION      R*8
%
% 
%   Output:
%               POS    : POSITION OF SATELLITE               R*8
%               VEL    : VELOCITY OF SATELLITE               R*8
%                        POS(K), VEL(K),K=1,2,3
% 
%   External calls: 	
%
%       
%   FROM BERNESE
%
%   Revision: 
%
%
% ************************************************************************


function [pos,vel] = ephem(GM,a,e,i,kn,per,t0,t)
 
% P=PARAMETER OF CONIC SECTION
        p = a*(1-e^2);

% N = MEAN MOTION, M = MEAN ANOMALY
        n = sqrt(GM/a^3);
        M = n*(t-t0)*86400;

% SOLVE KEPLER'S EQUATION
        EX1 = M;
        for ii = 1:10
          EX = EX1+(M+e*sin(EX1)-EX1)/(1-e*cos(EX1));
          if (abs(EX-EX1)<1e-12); break; end
          EX1 = EX;
        end

% V = TRUE ANOMALY
        v = 2*atan(sqrt((1+e)/(1-e))*tan(EX/2));
        r = a*(1-e*cos(EX));
        beta = sqrt(GM/p);
        x1 = r*cos(v);
        x2 = r*sin(v);
        v1 =-beta*sin(v);
        v2 = beta*(e+cos(v));

% SINES AND COSINES OF INCLINATION I, NODE K, PERIGEE O
        CK = cos(kn);
        SK = sin(kn);
        CI = cos(i);
        SI = sin(i);
        CO = cos(per);
        SO = sin(per);

% VECTORS P AND Q
        P1 = CK*CO-SK*CI*SO;
        P2 = SK*CO+CK*CI*SO;
        P3 = SI*SO;

        Q1 =-CK*SO-SK*CI*CO;
        Q2 =-SK*SO+CK*CI*CO;
        Q3 = SI*CO;

% COMPUTE POSITION AND VELOCITY
        pos(1,1) = P1*x1+Q1*x2;
        pos(2,1) = P2*x1+Q2*x2;
        pos(3,1) = P3*x1+Q3*x2;

        vel(1,1) = P1*v1+Q1*v2;
        vel(2,1) = P2*v1+Q2*v2;
        vel(3,1) = P3*v1+Q3*v2;
end