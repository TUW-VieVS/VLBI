% ************************************************************************
%   Description:
%	   a, e, omega, tosc, t0   --> u0
%
%   Input:										
%              GM      : GRAVITY PARAMETER                         R*8
%              a       : SEMIMAJOR AXIS                            R*8
%              e       : ECCENTRICITY                              R*8
%              omega   : ARGUMENT OF PERIGEE                       R*8
%              t0      : PERIGEE PASSING TIME                      R*8
%              tosc    : OSCULATION TIME (SEC)                     R*8
%
% 
%   Output:
%              u0      : ARGUMENT OF LATITUDE AT TIME TOSC         R*8
% 
%   External calls: 	
%
%       
%   FROM Bernese
%
%   Revision: 
%
%
% ************************************************************************


function [u0] = getu0(GM, a,e,omega,t0,tosc)

      n =sqrt(GM/a^3);
      
      M = n*(tosc-t0)*86400;
      E = M;
      for i=1:10
        EX = E+(M+e*sin(E)-E)/(1-e*cos(E));
        if (abs(EX-E)<1e-12)
            break
        end
        E = EX;
      end
      v0 = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
      u0 = omega+v0;
end