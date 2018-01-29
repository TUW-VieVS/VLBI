% ************************************************************************
%   Description:
%   function to exclude models from nutdy design matrix, weight matrix,
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


function [H10,Ph10,och10]=delmodel_scan_const_10(opt,H10,Ph10,och10)

if opt.nutdy.model ~= 1
    H10 = []; Ph10 = []; och10 = []; 
 else if opt.dut1.model == 1 && opt.nutdy.constrain == 0
         H10(1:size(H10,1),:) = 0; 
         Ph10(1:size(Ph10,1),:) = 0;
         och10 = [];
     end
end