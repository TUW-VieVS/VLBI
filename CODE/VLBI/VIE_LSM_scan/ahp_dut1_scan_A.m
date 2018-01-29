% ************************************************************************
%   Description:
%   function to form the design matrix for the DUT1 estimates as piecewise
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
%       'Apwdut1'     (nobserv/scan x number of dut1 pwlo unknowns)  design matrix for dut1 
%       'T_dut1'      (1,2)                                          Estimation epochs for dut1  
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   25 July 2012 by Claudia Tierno Ros
%
%   Revision: 
%   
% ************************************************************************


function [Apwdut1,T_dut1] = ahp_dut1_scan_A(scan,iscan,c,rad2mas,opt,t)


    % determining the xpol estimation intervals and puting observation times in to them 
      t10 = floor(t(iscan)/opt.dut1.int)*opt.dut1.int;
      t20 = ceil(t(iscan)/opt.dut1.int)*opt.dut1.int;

      T_dut1 = t10:opt.dut1.int:t20; % The estimation intervals for xpol
      
      if length(T_dut1)==1
      T_dut1 = [t10 opt.dut1.int+t20]; % The estimation intervals for xpol
      end
    
    
     for j=1:scan(iscan).nobs
         
         
         Apwdut1(j,1)= (1-(t(iscan)-T_dut1(1,1))/(T_dut1(1,2)-T_dut1(1,1)))...
             *scan(iscan).obs(j).ppol(3)*c*100*(1/rad2mas); % [sec/rad ----cm/mas]
         
         Apwdut1(j,2)= ((t(iscan)-T_dut1(1))/(T_dut1(2)-T_dut1(1)))...
             *scan(iscan).obs(j).ppol(3)*c*100*(1/rad2mas); % [sec/rad ----cm/mas]
        
            
     end


end