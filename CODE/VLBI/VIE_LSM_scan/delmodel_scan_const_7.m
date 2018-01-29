% ************************************************************************
%   Description:
%   function to exclude models from ypol design matrix, weight matrix,
%   and o-c vector of constraints
%
%   Reference: 
%
%   Input:	
%       'HX'               structure array     design matrix of the pseudo-observation equations (constraints)
%       'PhX'              structure array     weight matrix of the pseudo-observation equations 
%       'ochX'             structure array     o-c vector of constraints (zero vector)
%       'opt'              structure array     (for info. /DOC/opt.doc)
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

function [H7,Ph7,och7] = delmodel_scan_const_7(opt,H7,Ph7,och7)

if opt.ypol.model ~= 1 
     H7= []; Ph7 = []; och7 = []; 
 else if opt.ypol.model == 1 && opt.ypol.constrain == 0    
         H7(1:size(H7,1),:) = 0; 
         Ph7(1:size(Ph7,1),:) = 0;
         och7 = [];
     end
end
