% ************************************************************************
%   Description:
%   function to form the design matrix for the nutdy estimates as piecewise
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
%       'Apwnutdy'     (nobserv/scan x number of nutdy pwlo unknowns)  design matrix for nutdy 
%       'T_nutdy'      (1,2)                                           Estimation epochs for nutdy  
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   25 July 2012 by Claudia Tierno Ros
%
%   Revision: 
%   
% ************************************************************************



function [Apwnutdy,T_nutdy] = ahp_nutdy_scan_A(scan,iscan,c,rad2mas,opt,t)

    % determining the xpol estimation intervals and puting observation times in to them 
      t10 = floor(t(iscan)/opt.nutdy.int)*opt.nutdy.int;
      t20 = ceil(t(iscan)/opt.nutdy.int)*opt.nutdy.int;

      T_nutdy = t10:opt.nutdy.int:t20; % The estimation intervals for xpol
      
      if length(T_nutdy)==1
      T_nutdy = [t10 opt.nutdy.int+t20]; % The estimation intervals for xpol
      end
    
    
     for j=1:scan(iscan).nobs
         
         
         Apwnutdy(j,1)= (1-(t(iscan)-T_nutdy(1))/(T_nutdy(2)-T_nutdy(1)))...
             *scan(iscan).obs(j).pnut(2)*c*100*(1/rad2mas); % [sec/rad ----cm/mas]
         
         Apwnutdy(j,2)= ((t(iscan)-T_nutdy(1))/(T_nutdy(2)-T_nutdy(1)))...
             *scan(iscan).obs(j).pnut(2)*c*100*(1/rad2mas); % [sec/rad ----cm/mas]
        
     end
 