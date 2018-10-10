% ************************************************************************
%   Description:
%   function to form the design matrix for the troposphere north gradients 
%   as piecewise linear ofset functions
%
%   Reference: 
%
%   Input:										
%       'per_stat'     structure array     (for info. /DOC/per_stat.doc)
%       'T_'           structure array      Estimation epochs for clocks, zwd, ngr, egr, xyz
%       'na'           (1,1)                Number of antennas
%
%   Output:
%       'Apwngr'       structure array      Design matrix for north gradients         
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   July 2012 by Claudia Tierno Ros
%
%   Revision: 
%   
% ************************************************************************ 


function [Apwngr] = apw_ngr_scan(per_stat,T_,na)

    for j=1:na
   
        for k=1:per_stat(j).total

            el=(pi/2)-per_stat(j).zd(k);

            Apwngr(j).sm(k,1)= per_stat(j).first(k)*...
               (1-((per_stat(j).minute-T_(j).ngr(1))/(T_(j).ngr(2)...
               -T_(j).ngr(1))))/(tan(el)*sin(el)+0.0032)*cos(per_stat(j).az); 
           
            
         
            Apwngr(j).sm(k,2)= per_stat(j).first(k)*...
               (((per_stat(j).minute-T_(j).ngr(1))/(T_(j).ngr(2)...
                -T_(j).ngr(1))))/(tan(el)*sin(el)+0.0032)*cos(per_stat(j).az); 
       
            
         end
    end



