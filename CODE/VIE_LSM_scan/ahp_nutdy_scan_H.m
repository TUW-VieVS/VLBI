% ************************************************************************
%   Description:
%   function to form the design matrix, weight matrix, and o-c vector of
%   nutdy
%
%   Reference: 
%
%   Input:	
%       'num_inter_nutdy'        (1,1)                   number of estimation intervals nutdy
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
function [H10,Ph10,och10] = ahp_nutdy_scan_H(opt,num_inter_nutdy)



 for inter = 1:num_inter_nutdy 
     H10(inter, inter) = +1;  % design matrix for the xpol constraints as pseudo observations
     H10(inter, inter+1) = -1; % design matrix for the xpol constraints as pseudo observations 
     Ph10(inter, inter) =1./opt.nutdy.coef^2; % [1/mas^2] weight matrix of the design matrix for the nutdy constraints
 end
     

och10=zeros(num_inter_nutdy,1); % o-c vector for the xpol constraints

end






