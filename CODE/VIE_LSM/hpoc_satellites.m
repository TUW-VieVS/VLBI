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
%   2026-06-01 by H.Wolf: optimized code, added constraints for SRP parameters
%
% ************************************************************************
function [H, Ph, och] = hpoc_satellites(H, Ph, och, nsat, opt)

    n_sat = length(nsat);

    if opt.SatPos.pw_sat
        n_pw = nsat(1).pos;
        n_int  = n_pw - 1;
        
        if opt.SatPos.constr_sat
            H_pos1 = -diff(eye(n_int+1),1,2)';
            Ph_pos1 = eye(n_int) / (opt.satellite(1).SatPos.sat_pos_coef^2);
            oc_pos = zeros(n_int, 1);
        else
            H_pos1 = zeros(n_int,n_int+1);
            Ph_pos1 = zeros(n_int,n_int);
            oc_pos = [];
        end

        if opt.SatPos.constrRadialComponent
            if strcmp(opt.SatPos.sat_pos_est_ref_frame, 'rsw') || strcmp(opt.SatPos.sat_pos_est_ref_frame, 'ntw')
                H_rad_  = eye(n_pw);
                Ph_rad_ = eye(n_pw) / opt.SatPos.constrRadialComponentValue;
                
                H_pos1 = [H_pos1; H_rad_];
                Ph_pos1 = blkdiag(Ph_pos1, Ph_rad_);
                oc_pos = [oc_pos; zeros(n_pw, 1)];
            else
                warning('Radiale Fixierung nur im rsw- oder ntw-Rahmen möglich!');
            end
        end

        if n_sat > 1
            H_full = kron(eye(n_sat), H_pos1);
            Ph_full = kron(eye(n_sat), Ph_pos1);
            oc_full = repmat(oc_pos, n_sat, 1);
        else
            H_full = H_pos1; Ph_full = Ph_pos1; oc_full = oc_pos;
        end

        H(16).sm = blkdiag(H(16).sm, H_full);
        H(17).sm = blkdiag(H(17).sm, H_full);
        H(18).sm = blkdiag(H(18).sm, H_full);
        Ph(16).sm = blkdiag(Ph(16).sm, Ph_full);
        Ph(17).sm = blkdiag(Ph(17).sm, Ph_full);
        Ph(18).sm = blkdiag(Ph(18).sm, Ph_full);
        och(16).sv = vertcat(och(16).sv, oc_full);
        och(17).sv = vertcat(och(17).sv, oc_full);
        och(18).sv = vertcat(och(18).sv, oc_full);            
    end

    if opt.ORB.estORB  
        est_iors = find([opt.ORB.params(:).estimate]);
        num_est = numel(est_iors);
        for k=1:num_est
            iorb = est_iors(k);
            n_off_vals = zeros(n_sat, 1);
            for i = 1:n_sat
                n_off_vals(i) = nsat(i).(sprintf('orb%d', iorb));
            end
            total_dim = sum(n_off_vals);
            
            H_orb = zeros(total_dim, total_dim);
            Ph_orb = zeros(total_dim, total_dim);
            oc_orb = [];  

            H(20+iorb).sm = blkdiag(H(20+iorb).sm, H_orb) ;
            Ph(20+iorb).sm = blkdiag(Ph(20+iorb).sm, Ph_orb);
            och(20+iorb).sv = vertcat(och(20+iorb).sv, oc_orb);
        end
    end

    if opt.SRP.estSRP
        est_iors = find([opt.SRP.params(:).estimate]);
        num_est = numel(est_iors);
        for k=1:num_est
            isrp = est_iors(k);
            n_off_vals = zeros(n_sat, 1);
            for i = 1:n_sat
                n_off_vals(i) = nsat(i).(sprintf('srp%d', isrp));
            end
            total_dim = sum(n_off_vals);

            H_srp = zeros(total_dim, total_dim);
            Ph_srp= zeros(total_dim, total_dim);
            oc_srp = [];

            H(32+isrp).sm = blkdiag(H(32+isrp).sm, H_srp) ;
            Ph(32+isrp).sm = blkdiag(Ph(32+isrp).sm, Ph_srp);
            och(32+isrp).sv = vertcat(och(32+isrp).sv, oc_srp);
        end
    end
end