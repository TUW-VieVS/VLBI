% ************************************************************************
%   Description:
%   function to form the design matrix, weight matrix, and o-c vector of
%   constraints of clock parameters (only clock rates) at the co-location sites
%
%   Reference: 
%
%   Input:	
%       'H'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'Ph'         structure array     weight matrix of the pseudo-observation equations 
%       'och'        structure array     o-c vector of constraints (zero vector)
%       'n_'         structure array     number of estimates (pwlo or one offset)
%       'opt'        structure array     (for info. /DOC/opt.doc)
%       'na'         (1,1)               number of antennas
%       'c'          (1,1)               velocity of light in vacuum  
%       'antenna'    structure array     antenna information of the session
%
%   Output:
%       'H'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'Ph'         structure array     weight matrix of the pseudo-observation equations 
%       'och'        structure array     o-c vector of constraints (zero vector)
%       't_'         structure array     total estimate intervals of a specific estimate in a session
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   28 Sep 2015 by Younghee Kwak
%
%   Revision: 
%   15 Feb 2016 by Younghee Kwak modified for the standard VieVS version
% ************************************************************************
function [H,Ph,och,t_] = hpoc_comclk(H,Ph,och,n_,opt,na,c,antenna)

%-----------------------------------------%
% should be a separate module or added to the superstation file%
% now only hobart antennas % There might be twin/sibling telescopes in Wettzell, Onsala, Hartebeesthoek, Yebes, Kashima ...
% please set 'one offset per clock' at the first solution (Run -> VieVS estimation settings)
% recommended to set 'pwl offset & one rate per clock' at the Clock estimation (Estimation -> Least squares -> Clock)
% so that only clock rates are set as common parameters without quadratic terms
for istat = 1:na
    if antenna(istat).name == 'HOBART26'
        ncombi(1,1)=istat;
    end
    if antenna(istat).name == 'HOBART12'
        ncombi(1,2)=istat;
    end
%    if antenna(istat).name == 'HARTRAO '
%        ncombi(2,1)=istat;
%    end
%    if antenna(istat).name == 'HART15M '
%        ncombi(2,2)=istat;
%    end
end 
%-----------------------------------------%

Hcomqclk = []; Phcomqclk = []; ochcomqclk = [];
t_.clk = sum([n_.clk]);   % total number of clk paramters in a session % pwl clock offsets
t_.qclk = sum([n_.qclk]); % total number of clk paramters in a session
ncolo=1;                  % number of co-location sites % currently only 1, hobart
for icom = 1:ncolo % number of co-location site
    sumqclkj = 0;
    sumqclkk = 0;

    j=ncombi(icom,1);  % istat for co-location site  !! the order at antenna.mat
    k=ncombi(icom,2);  % istat for co-location site  !! the order at antenna.mat
    H_comqclk(icom).h(n_(j).qclk,t_.qclk) = 0;     % initialization for constraint matrix        !! the matrix size: (1) row    = no. of co-location sites * number of clk per
                                                   %                                                                 (2) column = no. of clk parameters 
    Ph_comqclk(icom).h(n_(j).qclk,n_(j).qclk)= 0;  % initialization for constraint weight matrix !!diagonal square matrix with same number of row of H_clk               
    oc_hcomqclk(icom).h(n_(j).qclk,1) = 0;

    for i = 1:j-1 % for both antennas j & k
        sumqclkj = sumqclkj + n_(i).qclk;  % counting number of clk estimates right before the estimate to use
                                           % if it is the 6th site, count clk numbers till end of 5th site
    end
    for i = 1:k-1 
        sumqclkk = sumqclkk + n_(i).qclk;
    end

    for inter = 1:n_(j).qclk       % for convenience, assume the numbers of clk are same for both antennas n_(j).qclk = n_(k).qclk
        H_comqclk(icom).h(inter,sumqclkj+inter) = -1;  %  e.g. clk of 'HOBART26(Ho)' % 1*clk_Hb - 1*clk_Ho = 0   ! j-th is reference
        H_comqclk(icom).h(inter,sumqclkk+inter) = +1;  %  e.g. clk of 'HOBART12(Hb)'                                                
        Ph_comqclk(icom).h(inter,inter) = 1./...
            (10^2); % [1/cm^2] 10cm/day weight matrix of the design matrix for the clk tie constarints
        oc_hcomqclk(icom).h(inter,1) = 1e-20;  % O-C vector for the clk constraints 
                                               % non-zero value but almost 0 to avoid disappearance after sparse (0-element row/column squeezed out) command in vie_lsm 
    end
    % Concatenating
    Hcomqclk = vertcat(Hcomqclk,H_comqclk(icom).h);
    Phcomqclk = blkdiag(Phcomqclk,Ph_comqclk(icom).h);
    ochcomqclk = vertcat(ochcomqclk,oc_hcomqclk(icom).h);
end

%save h.mat H_comqclk   % test
%save ph.mat Ph_comqclk %test

% add to original constrain matrix
H(2).sm   = vertcat(H(2).sm,   Hcomqclk);   
Ph(2).sm  = blkdiag(Ph(2).sm,  Phcomqclk); 
och(2).sv = vertcat(och(2).sv, ochcomqclk); 
