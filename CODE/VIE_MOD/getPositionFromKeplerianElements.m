% ************************************************************************
%   Description:
%	Computes the position and velocity for satellites based on the fso data
%   in a 5 minute interval. The position and velocity is then stored in the 
%   sources.s data for each satellite included in the NGS file.
%	The calculation is done using GPS Time - but in the sources.s struct the
%   the corresponding UTC time is stored.
%
%   Input:										
%      GM           GRAVITY PARAMETER     
%      a		    semi-major axis
%      e			eccentricity
%      i			inclination
%      Omega		Omega (right ascension of ascending node)
%      argp        argument of perigee
%      t			time
%      t0 			perigee passing time
% 
%   Output:
%     rrn           cartesian position of satellite
% 
%   External calls: 	
%
%       
%   Coded for VieVS: 
%   17 Dec 2024 by Helene Wolf
%
%   Revision: 
%
%
% ************************************************************************

function [rrn] = getPositionFromKeplerianElements(GM, a,e,i,Omega,argp,t,t0)

	% N = MEAN MOTION, M = MEAN ANOMALY
    n = sqrt(GM/a^3);
    M = n*(t-t0)*86400;

    % SOLVE KEPLER'S EQUATION
    E = M;
    for ii = 1:10
      EX = E+(M+e*sin(E)-E)/(1-e*cos(E));
      if (abs(EX-E)<1e-12)
		break
      end 
      E = EX;
    end

    %E = atan(tan(f/2)* (sqrt(1-e))/(sqrt(1+e)))*2;
    if E<0
        E=E+2*pi;
    end

    R1 = rotm(Omega,3)';
    R2 = rotm(i,1)';
    R3 = rotm(argp,3)';
    R = R1*R2*R3;
    q = [a*(cos(E)-e)
         a*sqrt(1-e^2)*sin(E)
         0];

    rrn = R*q;
end