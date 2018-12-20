% ************************************************************************
%   Description:
%   function to exclude models from design matrix, weight matrix,
%   and o-c vector of constraints
%
%   Reference: 
%
%   Input:	
%       'HX'               structure array      design matrix of the pseudo-observation equations (constraints)
%       'PhX'              structure array      weight matrix of the pseudo-observation equations 
%       'ochX'             structure array      o-c vector of constraints (zero vector)
%       'opt'              structure array      (for info. /DOC/opt.doc)
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



function [H1,Ph1,och1,H2,Ph2] = delmodel_scan_const_1_2(opt,H1,Ph1,och1,H2,Ph2)

if opt.pw_clk == 1  % Deleting rate and quadratic terms from clock function - only pwl offsets (main solution)
   H2= []; Ph2 = []; 
    
end

if opt.pw_clk == 0  % Deleting pwl offsets + rate + quadratic terms from clock function - (main solution)
    H1 = []; Ph1= []; och1 = []; 
    H2 = []; Ph2 = []; 
    
end
