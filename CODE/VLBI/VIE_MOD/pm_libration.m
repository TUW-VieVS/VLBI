% ************************************************************************
%   Description:
%   provides libration corrections for polar motion, 10 daily terms
%   according to table 5.1a
%   Reference: Conventions 5.5.1 (update 2 July 2010)
%
%    These corrections should be added to "average"
%    EOP values to get estimates of the instantaneous values.
%           
%               * vectorised *
%
%   Input:
%     rjd           epoch of interest given in mjd (n,1) - dynamical
%                   time or TT
%                
%   Output:
%     cor_x         tidal correction in x (radians of arc)
%     cor_y         tidal correction in y (radians of arc)
%    
% 
%   External calls: 	
%     fund_arg.m, as2rad.m
%
%
%   Coded for VieVS:
%   13 Aug 2010 by Lucia Plank
%
%   Revision: 
%
% *************************************************************************
function [cor_x,cor_y]=pm_libration(rjd)

dim=size(rjd);
if dim(2)>dim(1); rjd=rjd'; end;
      
      halfpi = 1.5707963267948966;

%  Oceanic tidal terms present in x (microas),y(microas),ut1(microseconds)       
%  NARG(j,6) : Multipliers of GMST+pi and Delaunay arguments. 

NARG=[...
1,-1, 0,-2, 0,-1
1,-1, 0,-2, 0,-2
1, 1, 0,-2,-2,-2
1, 0, 0,-2, 0,-1
1, 0, 0,-2, 0,-2
1,-1, 0, 0, 0, 0
1, 0, 0,-2, 2,-2
1, 0, 0, 0, 0, 0
1, 0, 0, 0, 0,-1
1, 1, 0, 0, 0, 0];

coeff= [...
  -.4,   .3,   -.3,  -.4
 -2.3,  1.3,  -1.3, -2.3
  -.4,   .3,   -.3,  -.4
 -2.1,  1.2,  -1.2, -2.1
-11.4,  6.5,  -6.5,-11.4
   .8,  -.5,    .5,   .8
 -4.8,  2.7,  -2.7, -4.8
 14.3, -8.2,   8.2, 14.3
  1.9, -1.1,   1.1,  1.9
   .8,  -.4,    .4,   .8];

XSIN =coeff(:,1);
XCOS =coeff(:,2);
YSIN =coeff(:,3);
YCOS =coeff(:,4);

      T = (rjd - 51544.5)/36525.0;  % julian century

      
% Arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega
      ARG = fund_arg(T,3); % [rad]

% Corrections
    
    agt=NARG*ARG';
    agt=mod(agt,4.*halfpi);
   
    cor_x  =  cos(agt')*XCOS    + sin(agt')*XSIN;
    cor_y  =  cos(agt')*YCOS    + sin(agt')*YSIN;
  
    cor_x   = as2rad(cor_x*1e-6);   %  radians
    cor_y   = as2rad(cor_y*1e-6);   %  radians
    
      
      
 
      