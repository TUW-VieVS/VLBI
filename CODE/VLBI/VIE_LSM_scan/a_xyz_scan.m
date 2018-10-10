% ************************************************************************
%   Description:
%   function to form the design matrix of station coordinates for single and global solution 
%   (one estimate per session)
%
%   Reference: 
%
%   Input:										
%       'per_stat'     structure array     (for info. /DOC/per_stat.doc)
%       'na'           (1,1)                Number of antennas
%
%   Output:
%      'Ax'    structure array              Design matrix of the station dx coordinate ofset estimate
%      'Ay'    structure array              Design matrix of the station dy coordinate offset estimate
%      'Az'    structure array              Design matrix of the station dz coordinate offset estimate
%      
%   External calls: 	
%   
%   Coded for VieVS: 
%   July 2012 by Claudia Tierno Ros 
%
%   Revision: 
%   
%
% ************************************************************************


function [Ax, Ay, Az] = a_xyz_scan(per_stat,na)

    for j=1:na
        for k=1:per_stat(j).total
            
            Ax(j).sm(k,1) = per_stat(j).first(k)*per_stat(j).dx(k);
            Ay(j).sm(k,1) = per_stat(j).first(k)*per_stat(j).dy(k);
            Az(j).sm(k,1) = per_stat(j).first(k)*per_stat(j).dz(k);
    
        end
    end



