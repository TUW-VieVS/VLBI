% ************************************************************************
%   Description:
%   Calculates positions and velocities per planet
%
%   Input:										
%		list_ib    	 specificator for interpolation
		%			       = 0 --> no interpolation for body ib
%   			    	   = 1 --> position only
%			       		   = 2 --> position and velocity
%		ib           number identificator of body
%		set			 data field
%		tc			 normalized chebyshev time
%		pc			 variable initialized in the load_eph.m
%		ncoeff		 number of coefficients per set
%		vc			 variable initialized in the load_eph.m
%		vfac		 variable initialized in the load_eph.m
%                
%   Output:
%      ephem              structure array with barycentric coordinates and
%                         velocities (if selected, also geocentric). Units
%                         are[m, m/s]. Stored at ../EPHEM/ephem.
%                         GM for the planets [m^3/s^2]
% 
%
%   Coded for VieVS: 
%   08 Aug 2016 by A. Girdiuk
%
%   Revision: 
%
% *************************************************************************
function [xbar,vbar]=eph_planet(list_ib,ib,set,tc,pc,ncoeff,vc,vfac)

twot = 0;
if list_ib~=0
	% calculate polynomial values
    if tc(ib) ~= pc(2)
    	np    = 2; nv = 3;
        pc(2) = tc(ib);
        twot  = tc(ib) + tc(ib);
	end
    if np < ncoeff(ib)
    	for i = np+1:ncoeff(ib)
        	pc(i) = twot*pc(i-1)- pc(i-2);
        end
        np=ncoeff;
	end
    % interpolate and get position for each component  
    xbar = (pc'* set)'*1e3;
    % Velocity
    if list_ib == 2   
    	vc(3)= twot + twot;    
		if nv < ncoeff(ib)
        	for i=nv+1:ncoeff(ib)     
        		vc(i) = twot*vc(i-1) + pc(i-1) + pc(i-1) - vc(i-2);
        	end
        	nv = ncoeff(ib);
      	end
	end 
    %  iterpolate and get velocity for each component
    vel = (vc'*set)';
    vbar = vfac * vel*1e3/86400;    
end
