%----------------------------------------------------------------------
% PURPOSE    :  COMPUTE PARTIAL DERIVATIVES OF KEPLERIAN ORBIT WITH
%               RESPECT TO THE FOLLOWING ORBITAL ELEMENTS
%
%               NUMBER   ELEMENT   EXPLANATION
%
%                 1        A       SEMI MAJOR AXIS
%                 2        E       EXCENTRICITY
%                 3        I       INCLINATION
%                 4        KN      R.A. OF ASCENDING NODE
%                 5        PER     ARGUMENT OF PERIGEE
%                 6        U0      ARGUMENT OF LATITUDE AT TIME TOSC
%
%               MANY PARTS OF THE COMPUTATION HAVE TO BE
%               PERFORMED ONLY AT THE FIRST TIME A SPECIAL ARC IS
%               ACTUALLY USED. THEREFORE THESE ITEMS ARE ONLY
%               RECOMPUTED IF THE ARC IDENTIFICATION IS CHANGING.
%
% PARAMETERS :
%         IN :  GM     : GRAVITY CONSTANT                    R*8
%               T      : TIME                                R*8
%               TOSC   : OSCULATION EPOCH                    R*8
%               A      : SEMI MAJOR AXIS                     R*8
%               E      : NUMERICAL ECCENTRICITY              R*8
%               I      : INCLINATION                         R*8
%               KN     : R.A. OF ASCENDING NODE              R*8
%               PER    : PERIGEE                             R*8
%               T0     : PERIGEE PASSING TIME                R*8
%        OUT :  DRDELE : ARRAY CONTAINING PARTIALS WITH      R*8(3,6)
%                        RESPECT TO ELEMENTS A,E,I,KN,PER,U0
%----------------------------------------------------------------------

function [drdele] = rpartn(GM, t,tosc,a,e,i,kn,per,t0)

      % 1. COMPUTE ROTATION MATRICES
      CK = cos(kn);
      SK = sin(kn);
      CI = cos(i);
      SI = sin(i);
      CP = cos(per);
      SP = sin(per);

% R3(-KN)*R1(-I)*R3(-PER) (FIRST TWO COLUMNS ONLY)
      ABC(1,1) = CK*CP-SK*CI*SP;
      ABC(2,1) = SK*CP+CK*CI*SP;
      ABC(3,1) = SI*SP;
      ABC(1,2) =-CK*SP-SK*CI*CP;
      ABC(2,2) =-SK*SP+CK*CI*CP;
      ABC(3,2) = SI*CP;

% R3(-KN)*D(R1(-I))/DI*R3(-PER)
      D(1,1) =-SK*CP-CK*CI*SP;
      D(2,1) = CK*CP-SK*CI*SP;
      D(3,1) = 0;
      D(1,2) = SK*SP-CK*CI*CP;
      D(2,2) =-CK*SP-SK*CI*CP;
      D(3,2) = 0;

% D(R3(-KN))/DI*R1(-I)*R3(-PER)
      E1(1,1) = SK*SI*SP;
      E1(2,1) =-CK*SI*SP;
      E1(3,1) = CI*SP;
      E1(1,2) = SK*SI*CP;
      E1(2,2) =-CK*SI*CP;
      E1(3,2) = CI*CP;

% R3(-KN)*R1(-I)*D(R3(-PER))/DPER
      F(1,1) =-CK*SP-SK*CI*CP;
      F(2,1) =-SK*SP+CK*CI*CP;
      F(3,1) = SI*CP;
      F(1,2) =-CK*CP+SK*CI*SP;
      F(2,2) =-SK*CP-CK*CI*SP;
      F(3,2) =-SI*SP;

%C COMPUTE MEAN MOTION
      n = sqrt(GM/a^3);

% SQRT(1-E**2)
      ew = sqrt(1-e^2);

% COMPUTE E0, V0
      M  = n*(tosc-t0)*86400;
      EX = M+e*sin(M);
      for iter=1:3
        dE=(M-EX+e*sin(EX))/(1-e*cos(EX));
        EX=EX+dE;
      end
      CEX = cos(EX);
      SEX = sin(EX);
      v0  = 2*atan(sqrt((1+e)/(1-e))*tan(EX/2));

% COMPUTE PARTIAL OF T0 WITH RESPECT TO A
      DT0DA = -1.5*(EX-e*SEX)/(n*a);

% COMPUTE PARTIAL OF T0 WITH RESPECT TO E
      DE0DE = -2/(1-e^2)*tan(EX/2)/(1+tan(EX/2)^2);
      DT0DE = -(DE0DE*(1-e*CEX)-SEX)/n;

% COMPUTE PARTIAL OF T0 WITH RESPECT TO U0
      HE1 = 1/sqrt(1+e);
      HE2 = 1/sqrt(1-e);
      DE0DU0 = HE1/HE2*(1+(tan(v0/2))^2)/(1+(tan(EX/2))^2);
      DT0DU0 = -(1-e*CEX)/n*DE0DU0;

% COMPUTE PARTIALS
% ----------------
% 1. KEPLERIAN POSITION
      M = n*(t-t0)*86400;
      EX = M+e*sin(M);
      for iter=1:3
        dE = (M-EX+e*sin(EX))/(1-e*cos(EX));
        EX = EX+dE;
      end
      SEX = sin(EX);
      CEX = cos(EX);
      X1 = a*(CEX-e);
      X2 = a*ew*SEX;
      RR = sqrt(X1^2+X2^2);

% 2. PARTIALS WITH RESPECT TO THE THREE EULERIAN ANGLES
%    (FIRST PART FOR PERIGEE ONLY)
      for k=1:3
        drdele(3,k) = E1(k,1)*X1+E1(k,2)*X2;
        drdele(4,k) = D(k,1)*X1+D(k,2)*X2;
        drdele(5,k) = F(k,1)*X1+F(k,2)*X2;
      end

% 3. PARTIAL WITH RESPECT TO A
      DEDA = n*(-1.5/a*(t-t0)*86400-DT0DA)/RR*a;
      Y1 = X1/a-a*SEX*DEDA;
      Y2 = X2/a+a*ew*CEX*DEDA;
      for k=1:3
        drdele(1,k) = ABC(k,1)*Y1+ABC(k,2)*Y2;
      end

% 4. PARTIAL WITH RESPECT TO ECCENTRICITY
      DEDE = (SEX-n*DT0DE)/RR*a;
      Y1 = -a*(1+SEX*DEDE);
      Y2 =  a*(-e/ew*SEX+ew*CEX*DEDE);
      for k=1:3
        drdele(2,k) = ABC(k,1)*Y1+ABC(k,2)*Y2;
      end
        
% 5. SECOND PART OF PARTIAL WITH RESPECT TO PER, PARTIAL WITH RESPECT TO U0
      DEDU0 = -n*a/RR*DT0DU0;
      Y1 = -a*SEX*DEDU0;
      Y2 = a*ew*CEX*DEDU0;
      for k=1:3
        drdele(6,k) = ABC(k,1)*Y1+ABC(k,2)*Y2;
        drdele(5,k) = drdele(5,k)-(ABC(k,1)*Y1+ABC(k,2)*Y2);
      end      
end