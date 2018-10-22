% ************************************************************************
%   Description:
%   function to form the design matrix for the troposphere east gradients 
%   as piecewise linear offset functions
%
%   Reference: 
%
%   Input:										
%       'per_stat'     structure array     (for info. /DOC/per_stat.doc)
%       'T_'           structure array      Estimation epochs for clocks, zwd, ngr, egr, xyz
%       'na'           (1,1)                Number of antennas
%
%   Output:
%       'Apwegr'       structure array      Design matrix for east gradients         
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   July 2012 by Claudia Tierno Ros
%
%   Revision: 
%   
% ************************************************************************ 


function [Apwegr] = apw_egr_scan(per_stat,T_,na)

    for j=1:na
        
        for k=1:per_stat(j).total

          el=(pi/2)-per_stat(j).zd(k);

          Apwegr(j).sm(k,1)= per_stat(j).first(k)*...
              (1-((per_stat(j).minute-T_(j).egr(1))/(T_(j).egr(2)...
              -T_(j).egr(1))))/(tan(el)*sin(el)+0.0032)*sin(per_stat(j).az); 

         
          Apwegr(j).sm(k,2)= per_stat(j).first(k)*...
              (((per_stat(j).minute-T_(j).egr(1))/(T_(j).egr(2)...
              -T_(j).egr(1))))/(tan(el)*sin(el)+0.0032)*sin(per_stat(j).az); 
        

        end
    end

