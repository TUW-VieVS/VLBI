% ************************************************************************
%   Description:
%   function to form the design matrix, weight matrix, and o-c vector of
%   nutdx
%
%   Reference: 
%
%   Input:	
%       'num_inter_nutdx'        (1,1)                   number of estimation intervals nutdx
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

function [H9,Ph9,och9] = ahp_nutdx_scan_H(opt,num_inter_nutdx)



 for inter = 1:num_inter_nutdx 
     H9(inter, inter) = +1;  % design matrix for the xpol constraints as pseudo observations
     H9(inter, inter+1) = -1; % design matrix for the xpol constraints as pseudo observations 
     Ph9(inter, inter) = 1./opt.nutdx.coef^2; % [1/mas^2] weight matrix of the design matrix for the nutdx constraints
 end
     

och9=zeros(num_inter_nutdx,1); % o-c vector for the xpol constraints

end






