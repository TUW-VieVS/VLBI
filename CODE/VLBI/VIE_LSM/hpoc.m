% ************************************************************************
%   Description:
%   function to form the design matrix, weight matrix, and o-c vector of
%   constraints of clocks, zwd, ngr, egr
%
%   Reference: 
%		based on functions by Kamil Teke
%
%   Input:	
%       'H'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'Ph'         structure array     weight matrix of the pseudo-observation equations 
%       'och'        structure array     o-c vector of constraints (zero vector)
%       'n_'         structure array     number of estimates (pwlo or one offset)
%       'opt'        structure array     (for info. /DOC/opt.doc)
%       'na'         (1,1)               number of antennas
%       'size_A_2'   dimension           number of estimated parameters of the real-observation equations
%
%   Output:
%       'H'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'Ph'         structure array     weight matrix of the pseudo-observation equations 
%       'och'        structure array     o-c vector of constraints (zero vector)
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   10 Oct 2016 by A. Girdiuk
%
%   Revision: 
%   
% ************************************************************************
function [H,Ph,och] = hpoc(H,Ph,och,n_,opt,na,size_A_2)

Hclk  = []; Phclk  = [];  
Hzwd  = []; Phzwd  = [];
Hharm = []; Phharm =[];

Hrelngr = []; Phrelngr = [];
Habsngr = []; Phabsngr = [];

Hrelegr = []; Phrelegr = [];
Habsegr = []; Phabsegr = [];

for istat = 1:na
    
    coef_clk     = opt.stat(istat).coef_clk; 
    coef_zwd     = opt.stat(istat).coef_zwd;
    coef_rel_ngr = opt.stat(istat).coef_rel_ngr;
    coef_abs_ngr = opt.stat(istat).coef_abs_ngr;
    coef_rel_egr = opt.stat(istat).coef_rel_egr;
    coef_abs_egr = opt.stat(istat).coef_abs_egr;
    
    mat_clk = diag(ones(1,n_(istat).clk)) - diag(ones(1,n_(istat).clk-1),1);
    Hclk = blkdiag(Hclk,mat_clk(1:n_(istat).clk-1,1:n_(istat).clk));
    Phclk = blkdiag(Phclk,diag(ones(1,n_(istat).clk-1).*1./coef_clk^2));

    mat_zwd = diag(ones(1,n_(istat).zwd)) - diag(ones(1,n_(istat).zwd-1),1);
    Hzwd = blkdiag(Hzwd,mat_zwd(1:n_(istat).zwd-1,1:n_(istat).zwd));
    Phzwd = blkdiag(Phzwd,diag(ones(1,n_(istat).zwd-1).*1./coef_zwd^2));
    
	mat_ngr_rel = diag(ones(1,n_(istat).ngr)) - diag(ones(1,n_(istat).ngr-1),1);
    Hrelngr = blkdiag(Hrelngr,mat_ngr_rel(1:n_(istat).ngr-1,1:n_(istat).ngr));
    Phrelngr = blkdiag(Phrelngr,diag(ones(1,n_(istat).ngr-1).*1./coef_rel_ngr^2));
    
	mat_ngr_abs = diag(ones(1,n_(istat).ngr));
    Habsngr = blkdiag(Habsngr,mat_ngr_abs(1:n_(istat).ngr-1,1:n_(istat).ngr));
    Phabsngr = blkdiag(Phabsngr,diag(ones(1,n_(istat).ngr-1).*1./coef_abs_ngr^2));
    
	mat_egr_rel = diag(ones(1,n_(istat).egr)) - diag(ones(1,n_(istat).egr-1),1);
    Hrelegr = blkdiag(Hrelegr,mat_egr_rel(1:n_(istat).egr-1,1:n_(istat).egr));
    Phrelegr = blkdiag(Phrelegr,diag(ones(1,n_(istat).egr-1).*1./coef_rel_egr^2));
    
	mat_egr_abs = diag(ones(1,n_(istat).egr));
    Habsegr = blkdiag(Habsegr,mat_egr_abs(1:n_(istat).egr-1,1:n_(istat).egr));
    Phabsegr = blkdiag(Phabsegr,diag(ones(1,n_(istat).egr-1).*1./coef_abs_egr^2));
    
end

% FORMING THE O-C VECTOR FOR THE CONSTRAINTS
oc_hclk(size(Hclk,1),1) = 0;  % O-C vector for the clock constraints
oc_hzwd(size(Hzwd,1),1) = 0;  % O-C vector for the zwd constraints

if opt.constr_clk == 0
    Hclk(1:size(Hclk,1),1:size(Hclk,2)) = 0; 
    Phclk(1:size(Phclk,1),1:size(Phclk,2)) = 0;
    oc_hclk = [];
end

if opt.constr_zwd == 0
    Hzwd(1:size(Hzwd,1),1:size(Hzwd,2)) = 0; 
    Phzwd(1:size(Phzwd,1),1:size(Phzwd,2)) = 0;
    oc_hzwd = [];
end

H(1).sm = Hclk;
%H(2).sm(size(A(2).sm,2)-na+1,size(A(2).sm,2)) = 0;
H(2).sm(size_A_2-na+1,size_A_2) = 0;

Ph(1).sm = Phclk; 
Ph(2).sm(size_A_2-na+1,size_A_2-na+1) = 0;
och(1).sv = oc_hclk; 

H(3).sm = Hzwd;
Ph(3).sm = Phzwd; 
och(3).sv = oc_hzwd;

oc_rel_hngr(size(Hrelngr,1),1) = 0;  % O-C vector for the relative constraints of  North gradients
oc_abs_hngr(size(Habsngr,1),1) = 0;  % O-C vector for the absolute constraints of  North gradients

H(4).rel_sm = Hrelngr;
Ph(4).rel_sm = Phrelngr;
och(4).rel_sv = oc_rel_hngr;

H(4).abs_sm = Habsngr;
Ph(4).abs_sm = Phabsngr;
och(4).abs_sv = oc_abs_hngr;

oc_rel_hegr(size(Hrelegr,1),1) = 0;  % O-C vector for the relative constraints of  North gradients
oc_abs_hegr(size(Habsegr,1),1) = 0;  % O-C vector for the absolute constraints of  North gradients

H(5).rel_sm  = Hrelegr;
Ph(5).rel_sm  = Phrelegr;
och(5).rel_sv = oc_rel_hegr;

H(5).abs_sm = Habsegr;
Ph(5).abs_sm = Phabsegr;
och(5).abs_sv = oc_abs_hegr;

