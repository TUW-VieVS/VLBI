% ************************************************************************
%   Description:
%   function to form the design matrix for the YPOL estimates as piecewise
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
%       'Apwypol'     (nobserv/scan x number of ypol pwlo unknowns)  design matrix for ypol 
%       'T_ypol'      (1,2)                                          Estimation epochs for ypol  
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   25 July 2012 by Claudia Tierno Ros
%
%   Revision: 
%   
% ************************************************************************

function [Apwypol,T_ypol] = ahp_ypol_scan_A(scan,iscan,c,rad2mas,opt,t)


    % determining the xpol estimation intervals and puting observation times in to them 
      t10 = floor(t(iscan)/opt.ypol.int)*opt.ypol.int;
      t20 = ceil(t(iscan)/opt.ypol.int)*opt.ypol.int;

      T_ypol = t10:opt.ypol.int:t20; % The estimation intervals for xpol
      
      if length(T_ypol)==1
      T_ypol = [t10 opt.ypol.int+t20]; % The estimation intervals for xpol
      end
    
    
     for j=1:scan(iscan).nobs
         
         
         Apwypol(j,1)= (1-(t(iscan)-T_ypol(1,1))/(T_ypol(1,2)-T_ypol(1,1)))...
             *scan(iscan).obs(j).ppol(2)*c*100*(1/rad2mas); % [sec/rad ----cm/mas]
         
         Apwypol(j,2)= ((t(iscan)-T_ypol(1))/(T_ypol(2)-T_ypol(1)))...
             *scan(iscan).obs(j).ppol(2)*c*100*(1/rad2mas); % [sec/rad ----cm/mas]
        
            
     end
    



