% ************************************************************************
%   Description:
%   function to exclude models from zwd design matrix, weight matrix,
%   and o-c vector of constraints
%
%   Reference: 
%
%   Input:	
%       'HX'               structure array     design matrix of the pseudo-observation equations (constraints)
%       'PhX'              structure array     weight matrix of the pseudo-observation equations 
%       'ochX'             structure array     o-c vector of constraints (zero vector)
%       'opt'              structure array     (for info. /DOC/opt.doc)
%       'na'               (1,1)               number of antennas
%       'num_inter_zwd'    (1,na)              number of estimates of zwd
%
%   Output:
%       'HX'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'PhX'         structure array     weight matrix of the pseudo-observation equations 
%       'ochX'        structure array     o-c vector of constraints (zero vector)
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   27 Aug 2012 by Claudia Tierno Ros
%
%   Revision: 
%
% ************************************************************************




function [H3,Ph3,och3]=delparam_scan_const_zwd(opt,na,num_inter_zwd,H3,Ph3,och3)

for istat = na:-1:1 
    if opt.stat(istat).zwd_inc == 0
      
        H3(:,sum(num_inter_zwd(1:istat-1))+istat:sum(num_inter_zwd(1:istat))+istat) = [];
        H3(sum(num_inter_zwd(1:istat-1))+1:sum(num_inter_zwd(1:istat)),:) = [];
        
        Ph3(:,sum(num_inter_zwd(1:istat-1))+1:sum(num_inter_zwd(1:istat))) = []; 
        Ph3(sum(num_inter_zwd(1:istat-1))+1:sum(num_inter_zwd(1:istat)),:) = []; 
        
        if opt.constr_zwd == 1
            och3(sum(num_inter_zwd(1:istat-1))+1:sum(num_inter_zwd(1:istat))) = []; 
        end

    end
    
end
    
