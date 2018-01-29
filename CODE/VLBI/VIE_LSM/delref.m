% ************************************************************************
%   Description:
%   function to delete reference clock
%
%   Reference: 
%
%   Input:										
%       'na'        (1,1)                number of antennas
%       'opt'        structure array     (for info. /DOC/opt.doc)
%       'sum_'       structure array     sum vector of number of estimates vector (n_)
%       'A'          structure array     design matrix of the real-observation equations
%       'H'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'Ph'         structure array     weight matrix of the pseudo-observation equations 
%       'och'        structure array     o-c vector of constraints (zero vector)
%       'n_'         structure array     number of estimates (pwlo or one offset)
%       't'          structure array     estimation intervals of clk, zwd, ngr, egr, xyz
%
%   Output:
%       'A'          structure array     design matrix of the real-observation equations
%       'H'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'Ph'         structure array     weight matrix of the pseudo-observation equations 
%       'och'        structure array     o-c vector of constraints (zero vector)
%       'nistat'     (1,1)               the number of reference clock (in the order of "antenna")
%       'n_'         structure array     number of estimates (pwlo or one offset)
%       't'          structure array     estimation intervals of clk, zwd, ngr, egr, xyz
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   12 May 2009 by Kamil Teke
%
%   Revision: 
%   06 Dec 2009 by Kamil Teke: header added
% ************************************************************************
function [A,H,Ph,och,nistat,n_,t] = delref(na,opt,sum_,A,H,Ph,och,n_,t)

    vecm = 0:1:na;
    for istat = 1 : na
       if opt.stat(istat).ref == 1
           nistat = istat;
            break 
       end
    end
     
    % from piecewise linear functions
    A(1).sm(:,sum_.clk(nistat)+1:sum_.clk(nistat+1)) = [];

    summ_clk = sum_.clk-vecm;
    H(1).sm(:,sum_.clk(nistat)+1:sum_.clk(nistat+1)) = [];
    H(1).sm(summ_clk(nistat)+1:summ_clk(nistat+1),:) = [];

    Ph(1).sm(:,summ_clk(nistat)+1:summ_clk(nistat+1)) = [];
    Ph(1).sm(summ_clk(nistat)+1:summ_clk(nistat+1),:) = [];

    if opt.constr_clk == 1
        och(1).sv(summ_clk(nistat)+1:summ_clk(nistat+1)) = [];
    end

    % from quadratic clock functions
    A(2).sm(:,sum_.qclk(nistat)+1:sum_.qclk(nistat+1)) = [];

    summ_qclk = sum_.qclk-vecm;
    H(2).sm(:,sum_.qclk(nistat)+1:sum_.qclk(nistat+1)) = [];
    H(2).sm(summ_qclk(nistat)+1:summ_qclk(nistat+1),:) = [];

    Ph(2).sm(:,summ_qclk(nistat)+1:summ_qclk(nistat+1)) = [];
    Ph(2).sm(summ_qclk(nistat)+1:summ_qclk(nistat+1),:) = [];

    % -------------------------------------------------------------------------
    n_(nistat).clk = 0; n_(nistat).qclk = 0;
    t(nistat).clk = 0;
    
end


