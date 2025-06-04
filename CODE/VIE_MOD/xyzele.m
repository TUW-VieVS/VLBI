%----------------------------------------------------------------------
% PURPOSE    :  THIS SUBROUTINE CALCULATES THE (OSCULATING)
%               ORBITAL ELEMENTS OF A CELESTIAL BODY , WHOSE
%               CARTESIAN POSITION- AND VELOCITY COMPONENTS ARE
%               GIVEN IN THE ARRAYS X(K) , V(K) , K=1,2,3
%
% PARAMETERS :
%         IN :  GM     : GRAVITY-CONSTANT (GAUSSIAN CONSTANT R*8
%                        PLANETARY APPLICATIONS)
%               T      : EPOCH TO WHICH X,V (AND ELEMENTS)   R*8
%                        REFER (SEC)
%               POS    : POSITION OF SATELLITE AT TIME T     R*8(3)
%               VEL    : VELOCITY OF SATELLITE AT TIME T     R*8(3)
%                        BODY AT TIME T
%        OUT :  A      : SEMI MAJOR AXIS                     R*8
%               E      : NUMERICAL EXCENTRICITY              R*8
%               I      : INCLINATION OF ORBITAL PLANE WITH   R*8
%                        RESPECT TO FUNDAMENTAL PLANE OF
%                        COORDINATE SYSTEM
%               KN     : LONGITUDE OF ASCENDING NODE         R*8
%               PER    : ARGUMENT OF PERIGEE/HELION          R*8
%               T0     : TIME OF PERIGEE/HELION PASSING      R*8
%----------------------------------------------------------------------

function [a,e,i,kn,per,t0] = xyzele(GM,t,pos,vel)

      h(1) = pos(2)*vel(3)-pos(3)*vel(2);
      h(2) =-pos(3)*vel(1)+pos(1)*vel(3);
      h(3) = pos(1)*vel(2)-pos(2)*vel(1);

      kn = atan2(h(1),h(2));
      i  = atan2(sqrt(h(1)^2+h(2)^2),h(3));

      p = (h(1)^2+h(2)^2+h(3)^2)/GM;
      r = sqrt(pos(1)^2+pos(2)^2+pos(3)^2);
      
      ecv = p/r-1;
      esv = sqrt(p/GM)/r*(pos(1)*vel(1)+pos(2)*vel(2)+pos(3)*vel(3));
      
      v1 = atan2(esv,ecv);
      e  = sqrt(ecv^2+esv^2);

      CK = cos(kn);
      SK = sin(kn);
      CI = cos(i);
      SI = sin(i);
      
      x(1) =    CK*pos(1)   +SK*pos(2);
      x(2) =-CI*SK*pos(1)+CI*CK*pos(2)+SI*pos(3);
      x(3) = SI*SK*pos(1)-SI*CK*pos(2)+CI*pos(3);
      
      u = atan2(x(2),x(1));
      per = u-v1;
      
      EX = 2*atan(sqrt((1-e)/(1+e))*tan(v1/2));
      a  = p/(1-e^2);
      t0 = t-(EX-e*sin(EX))/sqrt(GM/a^3)/86400; 
end
      
