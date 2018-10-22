% ************************************************************************
%   Description:
%   provides libration corrections for UT1, 11 semi-diurnal terms
%   according to table 5.1b
%   Reference: Conventions 5.5.3 (update 2 July 2010)
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
%     cor_ut         tidal correction in ut (seconds)
% 
%   External calls: 	
%     fund_arg.m, as2rad.m
%
%
%   Coded for VieVS:
%   13 Aug 2010 by Lucia Plank
%
%   Revision: 
% *************************************************************************
function [cor_ut]=ut_libration(rjd)

dim=size(rjd);
if dim(2)>dim(1); rjd=rjd'; end;
      
      halfpi = 1.5707963267948966;

%  Oceanic tidal terms present in x (microas),y(microas),ut1(microseconds)       
%  NARG(j,6) : Multipliers of GMST+pi and Delaunay arguments. 

NARG=[...
2,-2, 0,-2, 0,-2
2, 0, 0,-2,-2,-2
2,-1, 0,-2, 0,-2
2, 1, 0,-2,-2,-2
2, 0, 0,-2, 0,-1
2, 0, 0,-2, 0,-2 
2, 1, 0,-2, 0,-2 
2, 0,-1,-2, 2,-2
2, 0, 0,-2, 2,-2
2, 0, 0, 0, 0, 0
2, 0, 0, 0, 0,-1];

coeff= [...
 0.05, -0.03
 0.06, -0.03
 0.35, -0.20
 0.07, -0.04
-0.07,  0.04
 1.75, -1.01
-0.05,  0.03
 0.05, -0.03
 0.76, -0.44
 0.21, -0.12
 0.06, -0.04];

XSIN =coeff(:,1);
XCOS =coeff(:,2);

      T = (rjd - 51544.5)/36525.0;  % julian century
      
% Arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega
      ARG = fund_arg(T,3); % [rad]

% Corrections
    
    agt=NARG*ARG';
    agt=mod(agt,4.*halfpi);
   
    cor_ut  =  cos(agt')*XCOS    + sin(agt')*XSIN; % [musec]
  
    cor_ut   = (cor_ut*1e-6);   %  sec

    
      
      
 
      