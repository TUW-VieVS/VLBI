% ************************************************************************
%   Description:
%   function to exclude  models given below as stationwise 
%          - zenith wet delay
%          - tropospheric gradients
%          
%
%   Reference: 
%
%   Input:										
%       'na'             (1,1)                   number of antennas
%       'opt'             structure array       (for info. /DOC/opt.doc)
%       'num_inter_zwd'  (1, num. antennas)      number of estimation intervals of zwd
%       'num_inter_ngr'  (1, num. antennas)      number of estimation intervals of ngr
%       'num_inter_egr'  (1, num. antennas)      number of estimation intervals of egr
%       'A_session'       structure array        design matrix of the real-observation equations
%
%   Output:
%       'A_session'          structure array     design matrix of the real-observation equations
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   31 July 2012 by Claudia Tierno Ros
%
%   Revision: 
%  
% ************************************************************************  


function [A_session] = delparam_scan_A(na,opt,A_session,num_inter_zwd,num_inter_ngr,num_inter_egr)

for istat = na :-1: 1 
    
    if opt.stat(istat).zwd_inc == 0
               
        A_session(3).sm(:,sum(num_inter_zwd(1:istat-1))+istat:sum(num_inter_zwd(1:istat))+istat) = [];
       
    end
    
    if opt.stat(istat).ngr_inc == 0
              
        A_session(4).sm(:,sum(num_inter_ngr(1:istat-1))+istat:sum(num_inter_ngr(1:istat))+istat) = [];
        
    end
    
    if opt.stat(istat).egr_inc == 0
               
        A_session(5).sm(:,sum(num_inter_egr(1:istat-1))+istat:sum(num_inter_egr(1:istat))+istat) = [];
       
    end
    
    if opt.stc == 0
       opt.stat(istat).xyz_inc = 0;
    end
    if opt.stat(istat).xyz_inc == 0
       
        A_session(13).sm(:,istat) = [];
        A_session(14).sm(:,istat) = [];
        A_session(15).sm(:,istat) = [];
    end
end
    


end


