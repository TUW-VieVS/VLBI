% ************************************************************************
%   Description:
%   function to exclude models from nutdx design matrix, weight matrix,
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


function [H9,Ph9,och9]=delmodel_scan_const_9(opt,H9,Ph9,och9)

if opt.nutdx.model ~= 1
     H9 = []; Ph9 = []; och9 = [];
 else if opt.dut1.model == 1 && opt.nutdx.constrain == 0
         H9(1:size(H9,1),:) = 0; 
         Ph9(1:size(Ph9,1),:) = 0;
         och9 = [];
     end
end