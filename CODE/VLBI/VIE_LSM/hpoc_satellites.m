% ************************************************************************
%   Description:
%   function to form the design matrix, weight matrix, and o-c vector of
%   constraints of satellite coordinates
%
%   Reference: 
%
%   Input:	
%       'H'                     structure array     design matrix of the pseudo-observation equations (constraints)
%       'Ph'                    structure array     weight matrix of the pseudo-observation equations 
%       'och'                   structure array     o-c vector of constraints (zero vector)
%       'number_pwlo_per_sat'   (1, n_sat)          contains number of source coordinate offsets for each source
%       'n_sat'                 (1,1)               number of satellites
%       'opt'                   structure array     (for info. /DOC/opt.doc)
%
%   Output:
%       'H'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'Ph'         structure array     weight matrix of the pseudo-observation equations 
%       'och'        structure array     o-c vector of constraints (zero vector)
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   04 May 2017 by A. Hellerschmied
%
%   Revision: 
%   
% ************************************************************************
function [H, Ph, och] = hpoc_satellites(H, Ph, och, number_pwlo_per_sat, n_sat, opt)

% Init. and preallocate:
H_pos1      = []; 
H_pos2      = [];
H_pos3      = [];
Ph_pos1     = []; 
Ph_pos2     = []; 
Ph_pos3     = []; 
oc_pos1     = [];
oc_pos2     = [];
oc_pos3     = [];


% -------------------------------------------------------------------------
% FORMING THE CONSTRAINTS for SATELLITE COORDINATES
if opt.pw_sat == 1
    
    tot_numb_of_pwlo_estimates = sum(number_pwlo_per_sat);
    add_index = 0;
    
    % Loop over all satellites
    for i_sat = 1 : n_sat
        
        H_pos1_tmp(i_sat).h     = zeros(tot_numb_of_pwlo_estimates - n_sat, number_pwlo_per_sat(i_sat));
        Ph_pos1_tmp(i_sat).h   = zeros(tot_numb_of_pwlo_estimates - n_sat, number_pwlo_per_sat(i_sat) - 1);
        
        % Loop over all pwl estimation intervals:
        for i_inter = 1 : (number_pwlo_per_sat(i_sat) - 1)
            H_pos1_tmp(i_sat).h(add_index + i_inter, i_inter)         = +1;                                         % design matrix for the right ascension pseudo observtaions as constraints
            H_pos1_tmp(i_sat).h(add_index + i_inter, i_inter + 1)     = -1;                                         % design matrix for the declination pseudo observtaions as constraints
            Ph_pos1_tmp(i_sat).h(add_index + i_inter, i_inter)        = 1. / opt.satellite(i_sat).sat_pos_coef^2;   % weight matrix coefficients of the design matrix for the satellite coordinate constraints (H) [1/cm^2]
        end
        
        H_pos1 = horzcat(H_pos1, H_pos1_tmp(i_sat).h); % Concatenating
        Ph_pos1 = horzcat(Ph_pos1, Ph_pos1_tmp(i_sat).h); % Concatenating
        
        add_index = add_index + number_pwlo_per_sat(i_sat) - 1;
    end

    % -------------------------------------------------------------------------
    % FORMING THE O-C VECTOR FOR THE CONSTRAINTS
    
    if opt.constr_sat == 1
        oc_pos1 = zeros(size(H_pos1, 1), 1);
        oc_pos2 = oc_pos1;
        oc_pos3 = oc_pos1;
    else
        % Set H and P matrices to zero, if no constraints should be applied:
        H_pos1 = zeros(size(H_pos1, 1), size(H_pos1, 2));
        Ph_pos1 = zeros(size(Ph_pos1, 1), size(Ph_pos1, 2));
    end
    
    % In general, the same constraints apply for all three coordinates:
    H_pos2  = H_pos1;
    H_pos3  = H_pos1;
    Ph_pos2 = Ph_pos1;
    Ph_pos3 = Ph_pos1;
    
end

% Store results:
H(16).sm = H_pos1; 
H(17).sm = H_pos2; 
H(18).sm = H_pos3; 
Ph(16).sm = Ph_pos1; 
Ph(17).sm = Ph_pos2; 
Ph(18).sm = Ph_pos3; 
och(16).sv = oc_pos1; 
och(17).sv = oc_pos2;
och(18).sv = oc_pos3;

