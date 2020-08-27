% ************************************************************************
%   Description:
%   function to form the design matrix, weight matrix, and o-c vector of
%   constraints of sources' coordinates
%
%   Reference: 
%
%   Input:	
%       'H'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'Ph'         structure array     weight matrix of the pseudo-observation equations 
%       'och'        structure array     o-c vector of constraints (zero vector)
%       'nso'        structure array     the number of source coordinate offsets for each source
%       'ns'         (1,1)               number of sources
%       'opt'        structure array     (for info. /DOC/opt.doc)
%       'ts'        (1,1)                total number of estimates for all sources AFTER eleminating
%                                        non-observed sources --- ts.sources = sum([nso.sources]);
%
%   Output:
%       'H'          structure array     design matrix of the pseudo-observation equations (constraints)
%       'Ph'         structure array     weight matrix of the pseudo-observation equations 
%       'och'        structure array     o-c vector of constraints (zero vector)
%       'ts'        (1,1)                total number of estimates for all sources AFTER eleminating
%                                        non-observed sources --- ts.sources = sum([nso.sources]);
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   26 Aug 2009 by Kamil Teke
%
%   Revision: 
%   06 Dec 2009 by Kamil Teke: header added
%   10 Aug 2012 by Kamil Teke: units of constraints are now in mas & cm
%   18 Jun 2014 by Hana Krasna: sources can be estimated in vie_lsm with
%   NNR condition
%   17 Jul 2014 by Hana Krasna: very loose constrains of 10 mas are applied
%   on all sources together with the NNR condition to prevent the matrix
%   being singular
%   06 Jul 2016 by David Mayer: Increased performance 
%   11 Jul 2016 by David Mayer: loose constraints can now be set in GUI
% ************************************************************************
function [H,Ph,och,ts] = hpoc_sources(H,Ph,och,nso,ns,opt,ts)

Hra = []; Phra = []; Hde = []; Phde = []; 

% FORMING THE CONSTRAINTS of SOURCES COORDINATES
if opt.pw_sou == 1 % piecewise
    for isou = 1 : ns
        
        coef_rade = opt.source(isou).coef_rade;
        int_rade = opt.source(isou).int_rade;
        H_ra(isou).h(ts.sources-ns,nso(isou).sources) = 0;
        P_hra(isou).h(ts.sources-ns,nso(isou).sources-1)= 0;
        sumrade = 0;
        for i = 1:isou-1
            sumrade = sumrade + nso(i).sources-1;
        end
        for inter = 1:nso(isou).sources-1
            H_ra(isou).h(sumrade + inter, inter) = +1;  % design matrix for the right ascension pseudo observtaions as constraints
            H_ra(isou).h(sumrade + inter, inter+1) = -1; % design matrix for the declination pseudo observtaions as constraints
            P_hra(isou).h(sumrade + inter, inter) = 1./coef_rade^2; % [1/mas^2] weight matrix of the design matrix for the sources coordinates constraints
        end
        Hra = horzcat(Hra,H_ra(isou).h); % Concatenating
        Phra = horzcat(Phra,P_hra(isou).h); % Concatenating
        Hde = Hra;
        Phde = Phra;
    end
end

if opt.est_sourceNNR==1
    if opt.UseSourceAbsConstrNNR    
        Hra=diag(ones(ns,1));
        Hde=diag(ones(ns,1));
        Phra=eye(ns);
        Phde=eye(ns);
        p=1/opt.sourceAbsConstrNNR^2; % absolute constraints
        Phra=Phra.*p;
        Phde=Phde.*p;
    else
        Hra = zeros(ns); Hde = zeros(ns); 
        Phra = zeros(ns); Phde = zeros(ns); 
    end
end

% -------------------------------------------------------------------------
% FORMING THE O-C VECTOR FOR THE CONSTRAINTS
oc_hra(size(Hra,1),1) = 0;  % O-C vector for the constraints of  right ascension coordinates
oc_hde(size(Hde,1),1) = 0;  % O-C vector for the constraints of  declination coordinates

if opt.pw_sou == 1 & opt.constr_sou == 0
    Hra(1:size(Hra,1),1:size(Hra,2)) = 0; Hde = Hra; 
    Phra(1:size(Phra,1),1:size(Phra,2)) = 0; Phde = Phra;
    oc_hra = []; oc_hde = [];
end

if opt.est_sourceNNR==1 && opt.UseSourceAbsConstrNNR==0
     oc_hra = []; oc_hde = [];
end

H(11).sm = Hra; H(12).sm = Hde; 
Ph(11).sm = Phra; Ph(12).sm = Phde; 
och(11).sv = oc_hra; och(12).sv = oc_hde;
