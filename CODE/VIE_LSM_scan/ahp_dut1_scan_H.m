% ************************************************************************
%   Description:
%   function to form the design matrix, weight matrix, and o-c vector of
%   dut1
%
%   Reference: 
%
%   Input:	
%       'num_inter_dut1'        (1,1)                   number of estimation intervals dut1
%       'opt'                  structure array         (for info. /DOC/opt.doc)
%
%   Output:
%       'HX'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'PhX'         structure array     weight matrix of the pseudo-observation equations 
%       'ochX'        structure array     o-c vector of constraints (zero vector)
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   25 July 2012 by Claudia Tierno Ros
%
%   Revision: 
%  
% ************************************************************************
function [H8,Ph8,och8] = ahp_dut1_scan_H(opt,num_inter_dut1)



 for inter = 1:num_inter_dut1 
     H8(inter, inter) = +1;  % design matrix for the xpol constraints as pseudo observations
     H8(inter, inter+1) = -1; % design matrix for the xpol constraints as pseudo observations 
     Ph8(inter, inter) =1./opt.dut1.coef^2; % [1/mas^2] weight matrix of the design matrix for the dut1 constraints
 end
     

och8=zeros(num_inter_dut1,1); % o-c vector for the xpol constraints

end






