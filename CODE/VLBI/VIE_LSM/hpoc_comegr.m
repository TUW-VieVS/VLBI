% ************************************************************************
%   Description:
%   function to form the design matrix, weight matrix, and o-c vector of
%   constraints of common troposphere east gradients at the co-location sites
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
%       't'          structure array     estimation intervals of clk, zwd, ngr, egr, xyz
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
%   based on Kamil's hpoc_* codes
%   08 Jun 2015 by Younghee Kwak
%   
%   Revision: 
%   04 Feb 2016 by Younghee Kwak modified for the standard VieVS version
%   13 Dec 2016 by Younghee Kwak added 3 more sites: Hobart + Hartebeesthoek, Wettzell, Yebes
%   31 Jan 2017 by Younghee Kwak check co-located sites
%   13 Feb 2017 by Younghee Kwak skip the non-matching intervals 
% ************************************************************************
function [H,Ph,och,t_] = hpoc_comegr(H,Ph,och,n_,opt,na,c,antenna,t)

%-----------------------------------------%
% should be a separate module or added to the superstation file%
% four co-located VLBI sites with twin/sibling telescopes, i.e. Hobart, Hartebeesthoek, Wettzell, Yebes
% check if there are co-located sites
icolo=0;
for istat = 1:na
    if antenna(istat).name == 'HOBART26'
        for jstat = 1:na
            if antenna(jstat).name == 'HOBART12'
                % if both antennas are included in the list, treat them as a co-located site
                icolo=icolo+1;
                ncombi(icolo,1)=istat;
                ncombi(icolo,2)=jstat;
            end
        end
    end
    if antenna(istat).name == 'HARTRAO '
        for jstat = 1:na
            if antenna(jstat).name == 'HART15M '
                icolo=icolo+1;
                ncombi(icolo,1)=istat;
                ncombi(icolo,2)=jstat;
            end
        end
    end
    if antenna(istat).name == 'WETTZELL'
        for jstat = 1:na
            if antenna(jstat).name == 'WETTZ13N'
                icolo=icolo+1;
                ncombi(icolo,1)=istat;
                ncombi(icolo,2)=jstat;
            end
        end
    end
    if antenna(istat).name == 'YEBES40M'
        for jstat = 1:na
            if antenna(jstat).name == 'RAEGYEB '
                icolo=icolo+1;
                ncombi(icolo,1)=istat;
                ncombi(icolo,2)=jstat;
            end
        end
    end
end 
%-----------------------------------------%

Hcomegr = []; Phcomegr = []; ochcomegr = [];
t_.egr = sum([n_.egr]); % total number of egr paramters in a session
% total number of egr parameters of co-located sites in a session should be same in this version!!! 
% need to be fixed for epoch-wise comparison
% ncolo=4; % number of co-location sites % currently only 4, Hobart, Hartebeesthoek, Wettzell, Yebes
ncolo=icolo; % number of co-location sites 
sigma=66.7e-12; % sigma in weight matrix in seconds e.g. 66ps (around 2cm)  

for icom = 1:ncolo
    sumegrj = 0;
    sumegrk = 0;

    j=ncombi(icom,1);  % istat for co-location site  !! the order at antenna.mat
    k=ncombi(icom,2);  % istat for co-location site  !! the order at antenna.mat
    H_comegr(icom).h(n_(j).egr,t_.egr) = 0;     % initialization for constraint matrix        !! the matrix size: (1) row    = no. of co-location sites * number of egr per  
                                                %                                                                 (2) column = no. of egr parameters 
%    Ph_comegr(icom).h(n_(j).egr,n_(j).egr)= 0;  % initialization for constraint weight matrix !! diagonal square matrix with same number of row of H_egr               
%    oc_hcomegr(icom).h(n_(j).egr,1) = [];    

    for i = 1:j-1 % for both antennas j & k
        sumegrj = sumegrj + n_(i).egr;  % counting number of egr estimates right before the estimate to use
                                        % if it is the 6th site, count egr numbers till end of 5th site
    end
    for i = 1:k-1 
        sumegrk = sumegrk + n_(i).egr;
    end

    inter=0; % inter: the counts of fictitious observations for constraints 
    for jnter = 1:n_(j).egr
        for knter = 1:n_(k).egr
            if (t(j).egr(jnter) == t(k).egr(knter))
                inter=inter+1;
                H_comegr(icom).h(inter,sumegrj+jnter) = -1;  %  e.g. egr of 'HOBART26(Ho)' % 1*egr_Hb - 1*egr_Ho = 0   ! j-th is reference 
                H_comegr(icom).h(inter,sumegrk+knter) = +1;  %  e.g. egr of 'HOBART12(Hb)'
                Ph_comegr(icom).h(inter,inter) = 1./...
                    ((sigma)^2*c^2*100^2); % [1/cm^2] weight matrix of the design matrix for the egr tie constarints
                oc_hcomegr(icom).h(inter,1) = 1e-20;  % O-C vector for the egr constraints 
                                                  % non-zero value but almost 0 to avoid disappearance after sparse (0-element row/column squeezed out) command in vie_lsm 
            end
        end
    end
    % Concatenating
    Hcomegr = vertcat(Hcomegr,H_comegr(icom).h);
    Phcomegr = blkdiag(Phcomegr,Ph_comegr(icom).h);
    ochcomegr = vertcat(ochcomegr,oc_hcomegr(icom).h);
end

%save h.mat H_comegr   % test
%save ph.mat Ph_comegr %test

% add to original constrain matrix
H(5).rel_sm   = vertcat(H(5).rel_sm,Hcomegr);     
Ph(5).rel_sm  = blkdiag(Ph(5).rel_sm,Phcomegr); 
och(5).rel_sv = vertcat(och(5).rel_sv,ochcomegr); 
