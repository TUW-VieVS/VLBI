% ************************************************************************
%   Description:
%   function to exclude models from dut1 design matrix, weight matrix,
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


function [H8,Ph8,och8]=delmodel_scan_const_8(opt,H8,Ph8,och8)

if opt.dut1.model ~= 1
     H8 = []; Ph8 = []; och8 = []; 
 else if opt.dut1.model == 1 && opt.dut1.constrain == 0 
         H8(1:size(H8,1),:) = 0; 
         Ph8(1:size(Ph8,1),:) = 0;  
         och8 = [];
     end
end