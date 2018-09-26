% ************************************************************************
%   Description:
%   function to exclude  models given below as stationwise 
%          - zenith wet delay
%          - tropospheric gradients
%          - antenna TRF coordinates
%
%   Reference: 
%
%   Input:										
%       'na'        (1,1)                number of antennas
%       't'          structure array     estimation intervals of clk, zwd, ngr, egr, xyz
%       'opt'        structure array     (for info. /DOC/opt.doc)
%       'n_'         structure array     number of estimates (pwlo or one offset)
%       'sum_'       structure array     sum vector of number of estimates vector (n_)
%       'A'          structure array     design matrix of the real-observation equations
%       'H'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'Ph'         structure array     weight matrix of the pseudo-observation equations 
%       'och'        structure array     o-c vector of constraints (zero vector)
%
%   Output:
%       'n_'         structure array     number of estimates (pwlo or one offset)
%       't'          structure array     estimation intervals of clk, zwd, ngr, egr, xyz
%       'A'          structure array     design matrix of the real-observation equations
%       'H'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'Ph'         structure array     weight matrix of the pseudo-observation equations 
%       'och'        structure array     o-c vector of constraints (zero vector)
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   12 May 2009 by Kamil Teke
%
%   Revision: 
%   06 Dec 2009 by Kamil Teke: header added
%   06 Mar 2010 by Kamil Teke: absolute constraints for troposphere gradients corres. parts added
%   26 Sep 2018 by Daniel Landskron: bug corrected with excluding absolute constraints on gradients
% ************************************************************************ 
function [n_,t,A,H,Ph,och] = delparam(na,t,opt,n_,sum_,A,H,Ph,och)
         
inc_.zwd = 0; inc_.ngr = 0; inc_.egr = 0; inc_.xyz = 0; 
countzwd = 0; countngr = 0; countegr = 0; countxyz = 0;
del.zwd = []; del.ngr = []; del.egr = []; del.xyz = [];
vecm = 0:1:na;
for istat = 1 : na 
    if opt.stat(na+1-istat).zwd_inc == 0
        countzwd = countzwd + 1;
        del.zwd(countzwd) = na+1-istat;
        inc_.zwd = inc_.zwd + 1;
        summ_zwd = sum_.zwd-vecm;
        
        A(3).sm(:,sum_.zwd(del.zwd(countzwd))+1:sum_.zwd(del.zwd(countzwd)+1)) = [];
        H(3).sm(:,sum_.zwd(del.zwd(countzwd))+1:sum_.zwd(del.zwd(countzwd)+1)) = [];
        H(3).sm(summ_zwd(del.zwd(countzwd))+1:summ_zwd(del.zwd(countzwd)+1),:) = [];
        
        Ph(3).sm(:,summ_zwd(del.zwd(countzwd))+1:summ_zwd(del.zwd(countzwd)+1)) = []; 
        Ph(3).sm(summ_zwd(del.zwd(countzwd))+1:summ_zwd(del.zwd(countzwd)+1),:) = []; 
        
        if opt.constr_zwd == 1
            och(3).sv(summ_zwd(del.zwd(countzwd))+1:summ_zwd(del.zwd(countzwd)+1)) = []; 
        end
        
        n_(:,del.zwd(countzwd)).zwd = 0;
        t(del.zwd(countzwd)).zwd = 0;
    end
    
    if opt.stat(na+1-istat).ngr_inc == 0
        countngr = countngr + 1;
        del.ngr(countngr) = na+1-istat;
        inc_.ngr = inc_.ngr + 1;
        summ_ngr = sum_.ngr-vecm;
        
        A(4).sm(:,sum_.ngr(del.ngr(countngr))+1:sum_.ngr(del.ngr(countngr)+1)) = [];
        
        if opt.constr_rel_ngr == 1 & opt.constr_abs_ngr == 0
            H(4).rel_sm(:,sum_.ngr(del.ngr(countngr))+1:sum_.ngr(del.ngr(countngr)+1)) = [];
            H(4).rel_sm(summ_ngr(del.ngr(countngr))+1:summ_ngr(del.ngr(countngr)+1),:) = [];

            Ph(4).rel_sm(:,summ_ngr(del.ngr(countngr))+1:summ_ngr(del.ngr(countngr)+1)) = [];
            Ph(4).rel_sm(summ_ngr(del.ngr(countngr))+1:summ_ngr(del.ngr(countngr)+1),:) = [];

            och(4).rel_sv(summ_ngr(del.ngr(countngr))+1:summ_ngr(del.ngr(countngr)+1)) = [];
        end

        if opt.constr_rel_ngr == 0 & opt.constr_abs_ngr == 1
            H(4).abs_sm(:,sum_.ngr(del.ngr(countngr))+1:sum_.ngr(del.ngr(countngr)+1)) = [];
            H(4).abs_sm(summ_ngr(del.ngr(countngr))+1:summ_ngr(del.ngr(countngr)+1),:) = [];

            Ph(4).abs_sm(:,summ_ngr(del.ngr(countngr))+1:summ_ngr(del.ngr(countngr)+1)) = [];
            Ph(4).abs_sm(summ_ngr(del.ngr(countngr))+1:summ_ngr(del.ngr(countngr)+1),:) = [];

            och(4).abs_sv(summ_ngr(del.ngr(countngr))+1:summ_ngr(del.ngr(countngr)+1)) = [];
        end
        
        if opt.constr_rel_ngr == 1 & opt.constr_abs_ngr == 1            
            H(4).rel_sm(:,sum_.ngr(del.ngr(countngr))+1:sum_.ngr(del.ngr(countngr)+1)) = [];
            H(4).rel_sm(summ_ngr(del.ngr(countngr))+1:summ_ngr(del.ngr(countngr)+1),:) = [];

            Ph(4).rel_sm(:,summ_ngr(del.ngr(countngr))+1:summ_ngr(del.ngr(countngr)+1)) = [];
            Ph(4).rel_sm(summ_ngr(del.ngr(countngr))+1:summ_ngr(del.ngr(countngr)+1),:) = [];

            och(4).rel_sv(summ_ngr(del.ngr(countngr))+1:summ_ngr(del.ngr(countngr)+1)) = [];

            H(4).abs_sm(:,sum_.ngr(del.ngr(countngr))+1:sum_.ngr(del.ngr(countngr)+1)) = [];
            H(4).abs_sm(summ_ngr(del.ngr(countngr))+1:summ_ngr(del.ngr(countngr)+1),:) = [];

            Ph(4).abs_sm(:,summ_ngr(del.ngr(countngr))+1:summ_ngr(del.ngr(countngr)+1)) = [];
            Ph(4).abs_sm(summ_ngr(del.ngr(countngr))+1:summ_ngr(del.ngr(countngr)+1),:) = [];

            och(4).abs_sv(summ_ngr(del.ngr(countngr))+1:summ_ngr(del.ngr(countngr)+1)) = [];
        end
                
        n_(:,del.ngr(countngr)).ngr = 0;
        t(del.ngr(countngr)).ngr = 0;
    end
    
    if opt.stat(na+1-istat).egr_inc == 0
        countegr = countegr + 1;
        del.egr(countegr) = na+1-istat;
        inc_.egr = inc_.egr + 1;
        summ_egr = sum_.egr-vecm;
        
        A(5).sm(:,sum_.egr(del.egr(countegr))+1:sum_.egr(del.egr(countegr)+1)) = [];
        
        if opt.constr_rel_egr == 1 & opt.constr_abs_egr == 0
            H(5).rel_sm(:,sum_.egr(del.egr(countegr))+1:sum_.egr(del.egr(countegr)+1)) = [];
            H(5).rel_sm(summ_egr(del.egr(countegr))+1:summ_egr(del.egr(countegr)+1),:) = [];

            Ph(5).rel_sm(:,summ_egr(del.egr(countegr))+1:summ_egr(del.egr(countegr)+1)) = []; 
            Ph(5).rel_sm(summ_egr(del.egr(countegr))+1:summ_egr(del.egr(countegr)+1),:) = []; 

            och(5).rel_sv(summ_egr(del.egr(countegr))+1:summ_egr(del.egr(countegr)+1)) = [];
        end
        
        if opt.constr_rel_egr == 0 & opt.constr_abs_egr == 1
            H(5).abs_sm(:,sum_.egr(del.egr(countegr))+1:sum_.egr(del.egr(countegr)+1)) = [];
            H(5).abs_sm(summ_egr(del.egr(countegr))+1:summ_egr(del.egr(countegr)+1),:) = [];

            Ph(5).abs_sm(:,summ_egr(del.egr(countegr))+1:summ_egr(del.egr(countegr)+1)) = []; 
            Ph(5).abs_sm(summ_egr(del.egr(countegr))+1:summ_egr(del.egr(countegr)+1),:) = []; 

            och(5).abs_sv(summ_egr(del.egr(countegr))+1:summ_egr(del.egr(countegr)+1)) = [];
        end
        
        if opt.constr_rel_egr == 1 & opt.constr_abs_egr == 1   
            H(5).rel_sm(:,sum_.egr(del.egr(countegr))+1:sum_.egr(del.egr(countegr)+1)) = [];
            H(5).rel_sm(summ_egr(del.egr(countegr))+1:summ_egr(del.egr(countegr)+1),:) = [];

            Ph(5).rel_sm(:,summ_egr(del.egr(countegr))+1:summ_egr(del.egr(countegr)+1)) = []; 
            Ph(5).rel_sm(summ_egr(del.egr(countegr))+1:summ_egr(del.egr(countegr)+1),:) = []; 

            och(5).rel_sv(summ_egr(del.egr(countegr))+1:summ_egr(del.egr(countegr)+1)) = [];
            
            H(5).abs_sm(:,sum_.egr(del.egr(countegr))+1:sum_.egr(del.egr(countegr)+1)) = [];
            H(5).abs_sm(summ_egr(del.egr(countegr))+1:summ_egr(del.egr(countegr)+1),:) = [];

            Ph(5).abs_sm(:,summ_egr(del.egr(countegr))+1:summ_egr(del.egr(countegr)+1)) = []; 
            Ph(5).abs_sm(summ_egr(del.egr(countegr))+1:summ_egr(del.egr(countegr)+1),:) = []; 

            och(5).abs_sv(summ_egr(del.egr(countegr))+1:summ_egr(del.egr(countegr)+1)) = [];
        end
        
        n_(:,del.egr(countegr)).egr = 0;
        t(del.egr(countegr)).egr = 0;
    end
    
    if opt.stc == 0
       opt.stat(na+1-istat).xyz_inc = 0;
    end
    if opt.stat(na+1-istat).xyz_inc == 0
        
        countxyz = countxyz + 1;
        del.xyz(countxyz) = na+1-istat;
        
        inc_.xyz = inc_.xyz + 1;
        summ_xyz = sum_.xyz-vecm;

        A(13).sm(:,sum_.xyz(del.xyz(countxyz))+1:sum_.xyz(del.xyz(countxyz)+1)) = [];
        A(14).sm(:,sum_.xyz(del.xyz(countxyz))+1:sum_.xyz(del.xyz(countxyz)+1)) = [];
        A(15).sm(:,sum_.xyz(del.xyz(countxyz))+1:sum_.xyz(del.xyz(countxyz)+1)) = [];

        H(13).sm(:,sum_.xyz(del.xyz(countxyz))+1:sum_.xyz(del.xyz(countxyz)+1)) = [];
        H(13).sm(summ_xyz(del.xyz(countxyz))+1:summ_xyz(del.xyz(countxyz)+1),:) = [];

        H(14).sm = []; H(15).sm = [];
        H(14).sm = H(13).sm; 
        H(15).sm = H(13).sm;

        Ph(13).sm(:,summ_xyz(del.xyz(countxyz))+1:summ_xyz(del.xyz(countxyz)+1)) = []; 
        Ph(13).sm(summ_xyz(del.xyz(countxyz))+1:summ_xyz(del.xyz(countxyz)+1),:) = [];

        Ph(14).sm = []; Ph(15).sm = [];
        Ph(14).sm = Ph(13).sm; Ph(15).sm = Ph(13).sm;

        if opt.constr_xyz == 1
            och(14).sv = []; och(15).sv = []; 
            och(13).sv(summ_xyz(del.xyz(countxyz))+1:summ_xyz(del.xyz(countxyz)+1)) = [];
            och(14).sv = och(13).sv; 
            och(15).sv = och(13).sv;
        end

        n_(:,del.xyz(countxyz)).xyz = 0;
        t(del.xyz(countxyz)).xyz = 0;
    end
end
    
if opt.constr_rel_ngr == 1 & opt.constr_abs_ngr == 0
    H(4).sm = H(4).rel_sm; Ph(4).sm = Ph(4).rel_sm; och(4).sv = och(4).rel_sv;
end
if opt.constr_rel_ngr == 0 & opt.constr_abs_ngr == 1
    H(4).sm = H(4).abs_sm; Ph(4).sm = Ph(4).abs_sm; och(4).sv = och(4).abs_sv;
end
   
if opt.constr_rel_egr == 1 & opt.constr_abs_egr == 0
    H(5).sm = H(5).rel_sm; Ph(5).sm = Ph(5).rel_sm; och(5).sv = och(5).rel_sv;
end
if opt.constr_rel_egr == 0 & opt.constr_abs_egr == 1
    H(5).sm = H(5).abs_sm; Ph(5).sm = Ph(5).abs_sm; och(5).sv = och(5).abs_sv;
end

if opt.constr_rel_ngr == 1 & opt.constr_abs_ngr == 1
    H(4).sm = vertcat(H(4).rel_sm,H(4).abs_sm);
    Ph(4).sm = blkdiag(Ph(4).rel_sm,Ph(4).abs_sm);
    och(4).sv = vertcat(och(4).rel_sv,och(4).abs_sv);
end

if opt.constr_rel_egr == 1 & opt.constr_abs_egr == 1
    H(5).sm = vertcat(H(5).rel_sm,H(5).abs_sm);
    Ph(5).sm = blkdiag(Ph(5).rel_sm,Ph(5).abs_sm);
    och(5).sv = vertcat(och(5).rel_sv,och(5).abs_sv);
end

if opt.constr_rel_ngr == 0 & opt.constr_abs_ngr == 0
    if size(A(4).sm,2) == 0
        H(4).sm = [];
    else
        H(4).sm(1,size(A(4).sm,2)) = 0;
    end
    Ph(4).sm = [];
    och(4).sv = [];
end

if opt.constr_rel_egr == 0 & opt.constr_abs_egr == 0
    if size(A(5).sm,2) == 0
        H(5).sm = [];
    else
        H(5).sm(1,size(A(5).sm,2)) = 0;
    end   
    Ph(5).sm = [];
    och(5).sv = [];
end
