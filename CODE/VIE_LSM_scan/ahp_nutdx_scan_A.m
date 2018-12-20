% ************************************************************************
%   Description:
%   function to form the design matrix for the nutdx estimates as piecewise
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
%       'Apwnutdx'     (nobserv/scan x number of nutdx pwlo unknowns)  design matrix for nutdx 
%       'T_nutdx'      (1,2)                                          Estimation epochs for nutdx  
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   25 July 2012 by Claudia Tierno Ros
%
%   Revision: 
%   
% ************************************************************************



function [Apwnutdx,T_nutdx] = ahp_nutdx_scan_A(scan,iscan,c,rad2mas,opt,t)


    % determining the xpol estimation intervals and puting observation times in to them 
      t10 = floor(t(iscan)/opt.nutdx.int)*opt.nutdx.int;
      t20 = ceil(t(iscan)/opt.nutdx.int)*opt.nutdx.int;

      T_nutdx = t10:opt.nutdx.int:t20; % The estimation intervals for xpol
      
      if length(T_nutdx)==1
      T_nutdx = [t10 opt.nutdx.int+t20]; % The estimation intervals for xpol
      end
    
    
     for j=1:scan(iscan).nobs
         
         
         Apwnutdx(j,1)= (1-(t(iscan)-T_nutdx(1,1))/(T_nutdx(1,2)-T_nutdx(1,1)))...
             *scan(iscan).obs(j).pnut(1)*c*100*(1/rad2mas); % [sec/rad ----cm/mas]
         
         Apwnutdx(j,2)= ((t(iscan)-T_nutdx(1))/(T_nutdx(2)-T_nutdx(1)))...
             *scan(iscan).obs(j).pnut(1)*c*100*(1/rad2mas); % [sec/rad ----cm/mas]
        
            
     end

