% ************************************************************************
%   Description:
%   function to form the design matrix for the TRF coordinates of antennas 
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
%       'Apwx'      structure array         Design matrix for antenna dx coordinates 
%       'Apwy'      structure array         Design matrix for antenna dy coordinates
%       'Apwz'      structure array         Design matrix for antenna dz coordinates
%  
%   External calls: 	
%   
%   Coded for VieVS: 
%   July 2012 by Claudia Tierno Ros
%
%   Revision: 
%   
% ************************************************************************


function [Apwx, Apwy, Apwz] = apw_xyz_scan(per_stat,T_,na)

    for j=1:na
        for k=1:per_stat(j).total

            Apwx(j).sm(k,1)=per_stat(j).first(k)*...
                (1-(per_stat(j).minute-T_(j).xyz(1))/(T_(j).xyz(2)...
                 -T_(j).xyz(1)))*per_stat(j).dx(k);
            
           Apwx(j).sm(k,2)=per_stat(j).first(k)*...
                ((per_stat(j).minute-T_(j).xyz(1))/(T_(j).xyz(2)...
                 -T_(j).xyz(1)))*per_stat(j).dx(k);  
             
             Apwy(j).sm(k,1)=per_stat(j).first(k)*...
                (1-(per_stat(j).minute-T_(j).xyz(1))/(T_(j).xyz(2)...
                 -T_(j).xyz(1)))*per_stat(j).dy(k);
            
           Apwy(j).sm(k,2)=per_stat(j).first(k)*...
                ((per_stat(j).minute-T_(j).xyz(1))/(T_(j).xyz(2)...
                 -T_(j).xyz(1)))*per_stat(j).dy(k); 
             
             Apwz(j).sm(k,1)=per_stat(j).first(k)*...
                (1-(per_stat(j).minute-T_(j).xyz(1))/(T_(j).xyz(2)...
                 -T_(j).xyz(1)))*per_stat(j).dz(k);
            
           Apwz(j).sm(k,2)=per_stat(j).first(k)*...
                ((per_stat(j).minute-T_(j).xyz(1))/(T_(j).xyz(2)...
                 -T_(j).xyz(1)))*per_stat(j).dz(k); 
             
        end
    end



