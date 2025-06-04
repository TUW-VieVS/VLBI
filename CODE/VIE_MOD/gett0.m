% ************************************************************************
%   Description:
%	   PURPOSE    : a, e, omega, tosc, u0   --> t0
%
%   Input:										
%              GM      : GRAVITY PARAMETER                         R*8
%              a       : SEMIMAJOR AXIS                            R*8
%              e       : ECCENTRICITY                              R*8
%              omega   : ARGUMENT OF PERIGEE                       R*8
%              tosc    : OSCULATION TIME (SEC)                     R*8
%              u0      : ARGUMENT OF LATITUDE AT TIME TOSC         R*8
% 
%   Output:
%              t0      : PERIGEE PASSING TIME                      R*8
% 
%   External calls: 	
%       
%   FROM BERNESE
%
%   Revision: 
%
%
% ************************************************************************

function t0 = gett0(GM,a,e,omega,tosc,u0)
        n  = sqrt(GM/a^3);
        v0 = u0-omega;
        E = 2*atan(sqrt((1-e)/(1+e))*tan(v0/2));
        M  = E-e*sin(E);
        while M<0 
            M = M+2*pi;
        end
        while M>2*pi 
            M =  M-2*pi;
        end
        t0 = tosc-M/n/86400; 
end