% ************************************************************************
%   Description:
%   function to form the design matrix for the zenith wet delays as  
%   piecewise linear ofset functions
%
%   Reference: 
%
%   Input:										
%       'per_stat'   structure array     (for info. /DOC/per_stat.doc)
%       'T_'         structure array      Estimation epochs for clocks, zwd, ngr, egr, xyz
%       'na'         (1,1)                Number of antennas
%       't'          (1,num. of scans)    minutes from midnight to the time of the scans
%       'iscan'      (1,1)                Number of the scan that is being analysed
%
%   Output:
%       'Azwd'       structure array      Design matrix for pwlo zwd estimates         
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   July 2012 by Claudia Tierno Ros
%
%   Revision: 
%   
% ************************************************************************ 

function [Azwd] = apw_zwd_scan(per_stat,T_,na,t,iscan)
                               
    for j=1:na
        
        for k=1:per_stat(j).total
  
           Azwd(j).sm(k,1)=per_stat(j).first(k)*...
             (1-(t(iscan)-T_(j).zwd(1))/(T_(j).zwd(2)...
              -T_(j).zwd(1)))*per_stat(j).mf(k);
            
          
             
            Azwd(j).sm(k,2)=per_stat(j).first(k)*...
               ((t(iscan)-T_(j).zwd(1))/(T_(j).zwd(2)...
                -T_(j).zwd(1)))*per_stat(j).mf(k);
             
              
        end
    end

