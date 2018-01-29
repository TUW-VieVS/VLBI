% ************************************************************************
%   Description:
%   function to form the design matrix for the XPOL estimates as piecewise
%   linear ofsets
%
%   Reference: 
%
%   Input:										
%       'scan'        structure array           (for info. /DOC/scan.doc)
%       'iscan'       (1,1)                     Number of the scan that is being analysed
%       'opt'         structure array           (for info. /DOC/opt.doc)
%       'c'           (1,1)                     velocity of light in vacuum
%       rad2mas       (1,1)                     coefficient for conversion from radians to milli arc second 
%       't'           (1,num. of scans)         minutes from midnight to the time of the scans
%
%   Output:
%       'Apwxpol'     (nobserv/scan x number of xpol pwlo unknowns)  design matrix for xpol 
%       'T_xpol'      (1,2)                                          Estimation epochs for xpol  
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   25 July 2012 by Claudia Tierno Ros
%
%   Revision: 
%   
% ************************************************************************

function [Apwxpol,T_xpol] = ahp_xpol_scan_A(scan,iscan,opt,c,rad2mas,t)

% determining the xpol estimation intervals and puting observation times in to them 
 t10 = floor(t(iscan)/opt.xpol.int)*opt.xpol.int;
 t20 = ceil(t(iscan)/opt.xpol.int)*opt.xpol.int;

 T_xpol = t10:opt.xpol.int:t20; % The estimation intervals for xpol
 
 if length(T_xpol)==1
      T_xpol = [t10 opt.xpol.int+t20]; % The estimation intervals for xpol
 end
    
 for j=1:scan(iscan).nobs
      Apwxpol(j,1)= (1-(t(iscan)-T_xpol(1,1))/(T_xpol(1,2)-T_xpol(1,1)))...
             *scan(iscan).obs(j).ppol(1)*c*100*(1/rad2mas); % [sec/rad ----cm/mas]
         
      Apwxpol(j,2)= ((t(iscan)-T_xpol(1))/(T_xpol(2)-T_xpol(1)))...
             *scan(iscan).obs(j).ppol(1)*c*100*(1/rad2mas); % [sec/rad ----cm/mas]
 end



end






