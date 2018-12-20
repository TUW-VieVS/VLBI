% ************************************************************************
%   Description:
%   function to form the design matrix, weight matrix, and o-c vector of
%   zwd
%
%   Reference: 
%
%   Input:	
%       'num_inter_zwd'        (1,num. of antennas)     number of estimation intervals zwd
%       'opt'                  structure array         (for info. /DOC/opt.doc)
%       'na'                   (1,1)                    number of antennas
%       'c'                    (1,1)                    velocity of light in vacuum  
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

function [H,Ph,och] = hpoc_zwd_scan(num_inter_zwd,opt,na,c)

H = []; Ph = [];
t_zwd = sum(num_inter_zwd); % total estimate interval of WZDs in a session

for istat = 1:na  
    % FORMING THE CONSTRAINTS of WZDs
    coef_zwd = opt.stat(istat).coef_zwd;
    int_zwd = opt.stat(istat).int_zwd;
    H_zwd(t_zwd,num_inter_zwd(istat)+1) = 0;
    P_hzwd(t_zwd,num_inter_zwd(istat))= 0;
        sumzwd = 0;
        for i = 1:istat-1
            sumzwd = sumzwd + num_inter_zwd(i);
        end
        for inter = 1:num_inter_zwd(istat) 
            H_zwd(sumzwd + inter, inter) = +1;  % design matrix for the zwd pseudo observtaions as constarints
            H_zwd(sumzwd + inter, inter+1) = -1; % design matrix for the zwd pseudo observtaions as constarints
            P_hzwd(sumzwd + inter, inter) = 1./coef_zwd^2; % [1/cm^2] weight matrix of the design matrix for the zwd constraints
        end
    H = horzcat(H,H_zwd); % Concatenating
    Ph = horzcat(Ph,P_hzwd); % Concatenating
    clear H_zwd P_hzwd
end
% FORMING THE O-C VECTOR FOR THE CONSTRAINTS
och(size(H,1),1) = 0;  % O-C vector for the zwd constraints

if opt.constr_zwd == 0
    H(1:size(H,1),1:size(H,2)) = 0; 
    Ph(1:size(Ph,1),1:size(Ph,2)) = 0;
    och = [];
end

