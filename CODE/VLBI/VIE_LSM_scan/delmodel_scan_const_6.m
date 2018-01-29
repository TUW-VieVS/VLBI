% ************************************************************************
%   Description:
%   function to exclude models from xpol design matrix, weight matrix,
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

function [H6,Ph6,och6]=delmodel_scan_const_6(opt,H6,Ph6,och6)

if opt.xpol.model ~= 1
     H6= []; Ph6= []; och6= []; 
 else if opt.xpol.model == 1 && opt.xpol.constrain == 0
         H6(1:size(H6,1),:) = 0; 
         Ph6(1:size(Ph6,1),:) = 0;
         och6 = [];
     end   
end