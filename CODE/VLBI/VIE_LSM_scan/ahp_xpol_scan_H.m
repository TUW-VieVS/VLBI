% ************************************************************************
%   Description:
%   function to form the design matrix, weight matrix, and o-c vector of
%   xpol
%
%   Reference: 
%
%   Input:	
%       'num_inter_xpol'        (1,1)                   number of estimation intervals xpol
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
function [H6,Ph6,och6] = ahp_xpol_scan_H(opt,num_inter_xpol)



 for inter = 1:num_inter_xpol 
     H6(inter, inter) = +1;  % design matrix for the xpol constraints as pseudo observations
     H6(inter, inter+1) = -1; % design matrix for the xpol constraints as pseudo observations 
     Ph6(inter, inter) = 1./opt.xpol.coef^2; % [1/mas^2] weight matrix of the design matrix for the xpol constraints
 end
     

och6=zeros(num_inter_xpol,1); % o-c vector for the xpol constraints

end






