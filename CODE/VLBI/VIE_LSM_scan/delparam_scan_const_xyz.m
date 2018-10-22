% ************************************************************************
%   Description:
%   function to exclude models from antenna coord. design matrix, weight matrix,
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


function [H13,H14,H15,Ph13,Ph14,Ph15,och13,och14,och15] =delparam_scan_const_xyz(opt,H13,H14,H15,Ph13,Ph14,Ph15,och13,och14,och15,na)
         
for istat=na:-1:1   
if opt.stc == 0 || opt.stat(istat).xyz_inc==0
       
              
        H13(:,istat) = [];
               
        H14 = H13; 
        H15 = H13;


        if opt.constr_xyz == 1
            och13 = [];
            och14= och13; 
            och15 = och13;
        end

        
end
end
