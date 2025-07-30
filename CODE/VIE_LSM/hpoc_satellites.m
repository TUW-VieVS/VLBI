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
%   2025-02-02 by H.Wolf: added contraint matrix options for orbital elements
%
% ************************************************************************
function [H, Ph, och] = hpoc_satellites(H, Ph, och, nso, n_sat, opt)

    if opt.SatPos.pw_sat
        number_pwlo_per_sat = nso.sat_pos;

        H_pos1      = [];     H_pos2      = [];    H_pos3      = [];
        Ph_pos1     = [];     Ph_pos2     = [];    Ph_pos3     = []; 
        oc_pos1     = [];     oc_pos2     = [];    oc_pos3     = [];

        H_pos1_tmp.h    = zeros(number_pwlo_per_sat-1, number_pwlo_per_sat);
        Ph_pos1_tmp.h   = zeros(number_pwlo_per_sat-1, number_pwlo_per_sat - 1);
        
        % Loop over all pwl estimation intervals:
        for i_inter = 1 : (number_pwlo_per_sat - 1)
            H_pos1_tmp.h(i_inter, i_inter)         = +1;                                         % design matrix for the right ascension pseudo observtaions as constraints
            H_pos1_tmp.h(i_inter, i_inter + 1)     = -1;                                         % design matrix for the declination pseudo observtaions as constraints
            Ph_pos1_tmp.h(i_inter, i_inter)        = 1. / opt.satellite(1).SatPos.sat_pos_coef^2;   % weight matrix coefficients of the design matrix for the satellite coordinate constraints (H) [1/cm^2]
        end
        
        H_pos1 = horzcat(H_pos1, H_pos1_tmp.h); % Concatenating
        Ph_pos1 = horzcat(Ph_pos1, Ph_pos1_tmp.h); % Concatenating
    
        % FORMING THE O-C VECTOR FOR THE CONSTRAINTS
        if opt.SatPos.constr_sat == 1
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
    
        if opt.SatPos.constrRadialComponent == 1
            constr_RComp = opt.SatPos.constrRadialComponentValue;
            
            if strcmp(opt.SatPos.sat_pos_est_ref_frame, 'rsw') || strcmp(opt.SatPos.sat_pos_est_ref_frame, 'ntw') 
                H_pos1_fix = diag(ones(1,size(H_pos1,2)));
                Ph_pos1_fix = diag((1/constr_RComp)*ones(1,size(H_pos1,2)));
                H_pos1 = [H_pos1; H_pos1_fix];
                Ph_pos1 =[blkdiag(Ph_pos1,Ph_pos1_fix)];
                oc_pos1 = [oc_pos1;[zeros(size(H_pos1_fix,1),1)]];
            else
                disp('WARNING: Radial/Normal component cannot be fixed because satellite position is not estimated in the rsw-frame or ntw-frame!')
            end
        end

        for i=1:n_sat
            H(16).sm = blkdiag(H(16).sm, H_pos1);
            H(17).sm = blkdiag(H(17).sm, H_pos2);
            H(18).sm = blkdiag(H(18).sm, H_pos3);
            Ph(16).sm = blkdiag(Ph(16).sm, Ph_pos1);
            Ph(17).sm = blkdiag(Ph(17).sm, Ph_pos2);
            Ph(18).sm = blkdiag(Ph(18).sm, Ph_pos3);
            och(16).sv = vertcat(och(16).sv, oc_pos1);
            och(17).sv = vertcat(och(17).sv, oc_pos2);
            och(18).sv = vertcat(och(18).sv, oc_pos3);
        end
    end


    if opt.KepEle.estKepEle  
        for iKep=1:6
            for isat=1:n_sat
                if opt.KepEle.('estKepEle' + string(iKep))
                    n_pwlo = nso(isat).('KepEle' + string(iKep));
                    if  opt.KepEle.('relConstrKepEle' + string(iKep))
                        % to be implemented
                        %coef_p1 = opt.('relConstrValKepEle' + string(numKepEle));
                        %Ph_KepEle = diag(ones(1,n_unk)).*1./coef_p1^2;
                        %mat = diag(ones(1,n_unk+1)) - diag(ones(1,n_unk),1);
                        %H_KepEle = mat(1:n_unk,1:n_unk+1);
                        %oc_KepEle = zeros(size(H_p1,1),1);
                    else
                        H_KepEle = zeros(n_pwlo-1, n_pwlo);
                        Ph_KepEle= zeros(n_pwlo-1, n_pwlo-1);
                        oc_KepEle = [];
                    end
                        H(20+iKep).sm = blkdiag(H(20+iKep).sm, H_KepEle) ;
                        Ph(20+iKep).sm = blkdiag(Ph(20+iKep).sm, Ph_KepEle);
                        och(20+iKep).sv = vertcat(och(20+iKep).sv, oc_KepEle);
                end
             end
        end
    end
end