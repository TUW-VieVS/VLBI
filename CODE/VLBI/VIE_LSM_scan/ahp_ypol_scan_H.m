% ************************************************************************
%   Description:
%   function to form the design matrix, weight matrix, and o-c vector of
%   ypol
%
%   Reference: 
%
%   Input:	
%       'num_inter_ypol'        (1,1)                   number of estimation intervals ypol
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
function [H7,Ph7,och7] = ahp_ypol_scan_H(opt,num_inter_ypol)



 for inter = 1:num_inter_ypol 
     H7(inter, inter) = +1;  % design matrix for the xpol constraints as pseudo observations
     H7(inter, inter+1) = -1; % design matrix for the xpol constraints as pseudo observations 
     Ph7(inter, inter) = 1./opt.ypol.coef^2; % [1/mas^2] weight matrix of the design matrix for the xpol constraints
 end
     

och7=zeros(num_inter_ypol,1); % o-c vector for the xpol constraints

end






