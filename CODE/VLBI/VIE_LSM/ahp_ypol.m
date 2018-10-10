% ************************************************************************
%   Description:
%   function to form the design matrix for the YPOL estimates as piecewise
%   linear ofsets
%
%   Reference: 
%
%   Input:										
%       'scan'        structure array           (for info. /DOC/scan.doc)
%       'mjd0'        (1,1)                     the midnight (the beginning of the day of the session), UTC
%       'ntim'        (1,1)                     number of scans in the session     
%       'opt'         structure array           (for info. /DOC/opt.doc)
%       'T'           structure array           estimation intervals for dUT1, xpol, ypol, nutdx, and nutdy
%       'c'           (1,1)                     velocity of light in vacuum
%       rad2mas       (1,1)                     coefficient for conversion from radians to milli arc second 
%
%   Output:
%       'Apwypol'      (nobserv x number of ypol pwlo unknowns)  design matrix for ypol 
%       'T'            structure array                           estimation intervals for dUT1, xpol, ypol, nutdx, and nutdy
%       'n_unk_ypol'   (1,1)                                     number of ypol estimation intervals in a session  
%       'Hypol'        (size(Apwypol,2)-1,size(Apwypol,2))       design matrix for the ypol constraints as pseudo observations  
%       'Phypol'       (size(Apwypol,2)-1,size(Apwypol,2)-1)     weight matrix of the ypol constraints
%       'oc_hypol'     (size(Hypol,1),1)                         o-c vector for the ypol constraints              
%      
%   External calls: 	
%   
%   Coded for VieVS: 
%   12 May 2009 by Kamil Teke
%
%   Revision: 
%   05 Dec 2009 by Kamil Teke: header added
%   10 Aug 2012 by Kamil Teke: units of constraints are now in mas & cm
%   10 Oct 2016 by A. Girdiuk: code-optimized: loop removed
% ************************************************************************
function [Apwypol,T,Hypol,Phypol,oc_hypol] = ahp_ypol(scan,mjd0,opt,T,c,rad2mas,obs_mjd)

int_ypol = opt.ypol.int;

mjd1 = min([scan.mjd]); % The time of the first scan in mjd [day]
mjd2 = max([scan.mjd]); % The time of the last scan in mjd [day]

t1 = (mjd1-mjd0)*24*60; % [minute] mjd0(midnight)--------mjd1(start)------(Session)------mjd2(end)-----
t2 = (mjd2-mjd0)*24*60; % [minute]
t =  ([scan.mjd] - mjd0)*24*60; % [minute]

% determining the ypol estimation intervals and puting observation times in to them 
t10 = floor(t1/int_ypol)*int_ypol;
t20 = ceil(t2/int_ypol)*int_ypol;

T.ypol = t10:int_ypol:t20; % The estimation intervals for ypol
n_unk_ypol = length(T.ypol)-1; % The number of estimation intervals for ypol

% number of scans per clock estimate interval [1]
stm_ypol(1,n_unk_ypol) = 0; 
for inter = 1:n_unk_ypol
    if inter == n_unk_ypol
        tst= (obs_mjd >= T.ypol(inter)).*(obs_mjd <= T.ypol(inter+1));
    else
        tst= (obs_mjd >= T.ypol(inter)).*(obs_mjd < T.ypol(inter+1));
    end
	stm_ypol(inter) = sum(tst); % number of scans per xpol estimation interval [1]
end

temp = [scan.obs];

k = 0; Apwypol(length(obs_mjd),n_unk_ypol+1)=0;
for inter = 1:n_unk_ypol % number of ypol estimate intervals in a session 1-24
    for iobs = 1:stm_ypol(inter) % number of observations in a ypol est. interval 241 229 
        k = k + 1;
        Apwypol(k,inter) = (1-(obs_mjd(k)-T.ypol(inter))/(T.ypol(inter+1)-T.ypol(inter)))...
            *temp(k).ppol(2)*c*100*(1/rad2mas); % [sec/rad ----cm/mas]
        Apwypol(k,inter+1) = ((obs_mjd(k)-T.ypol(inter))/(T.ypol(inter+1)-T.ypol(inter)))...
            *temp(k).ppol(2)*c*100*(1/rad2mas); % [sec/rad ----cm/mas]
    end
end

% FORMING THE CONSTRAINTS of YPOL
coef_ypol = opt.ypol.coef;

mat = diag(ones(1,n_unk_ypol+1)) - diag(ones(1,n_unk_ypol),1);
Hypol = mat(1:n_unk_ypol,1:n_unk_ypol+1);
Phypol = diag(ones(1,n_unk_ypol).*1./coef_ypol^2);

% FORMING THE O-C VECTOR FOR THE CONSTRAINTS
oc_hypol(size(Hypol,1),1) = 0; % o-c vector for the ypol constraints
