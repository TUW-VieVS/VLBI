% ************************************************************************
%   Function par_tiderp 
%   This function provides the fundamental arguments for the diurnal/subdiurnal
%   tidal terms to finally extract the partial derivatives for m
%   amplitudes of diurnal/subdiurnal tidal variations in ERP.
%
%   Input:
%     rjd            epoch of interest given in mjd (n,1) - dynamical
%                    time or TT
%                
%   Output:
%     cosag          cosine of tidal argument (m,n)
%     sinag          sine of tidal argument (m,n)
% 
%   External calls: 	
%     OUT/GLOB/tide.mat
%     fund_arg.m
%
%   Coded for VieVS:
%   03 Nov 2011 by Sigrid Boehm
%
% *************************************************************************


function [cosag1, sinag1, cosag2, sinag2] = par_tiderp(rjd)

 
%   NARG(j,6) : Multipliers of GMST+pi and Delaunay arguments for ocean
%   tidal terms as in IERS Conventions 2010: GMST+pi l l' F D OM
	
     fileID = fopen('../DATA/GLOB/tidalERPvar_list.txt');
     C = textscan(fileID,'%s %d %d %d %d %d %d %d','Delimiter',',','CommentStyle','%');
     fclose(fileID);
     
     tide.name = cell2mat(C{1,1});
     tide.num = double([C{1,8}]);
     tide.gmstpi = double([C{1,2}]);
     tide.l = double([C{1,3}]);
     tide.lp = double([C{1,4}]);
     tide.F = double([C{1,5}]);
     tide.D = double([C{1,6}]);
     tide.OM = double([C{1,7}]);
     
     save('../DATA/GLOB/tide.mat','tide');
     
     ind = (tide.gmstpi<2); 

     NARG1 = [tide.gmstpi(ind) tide.l(ind) tide.lp(ind) ...
              tide.F(ind) tide.D(ind) tide.OM(ind)];

     NARG2 = [tide.gmstpi(~ind) tide.l(~ind) tide.lp(~ind) ...
              tide.F(~ind) tide.D(~ind) tide.OM(~ind)];
 
         
% ------------------------
%  TIME Vector
% ------------------------
 
      T  = (rjd - 51544.5D0)/36525.0D0; % julian centuries
      
% ------------------------
% FUNDAMENTAL ARGUMENTS
% ------------------------
      
      ARG = fund_arg(T,3); % [rad]

% --------------------------
% CALCULATE SIN and COS
% --------------------------
    
    arg1  = NARG1*ARG';
    arg1  = mod(arg1,2*pi);
    
    arg2  = NARG2*ARG';
    arg2  = mod(arg2,2*pi);

    cosag1   = cos(arg1);
    sinag1   = sin(arg1);
    
    cosag2   = cos(arg2);
    sinag2   = sin(arg2);
  
      	
