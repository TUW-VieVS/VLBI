% ************************************************************************
%   Description:
%   function to form the design matrix for the XPOL estimates as piecewise
%   linear ofsets
%
%   Reference: 
%
%   Input:										
%       'scan'        structure array           (for info. /DOC/scan.doc)
%       'mjd0'        (1,1)                     the midnight (the beginning of the day of the session), UTC
%       'ntim'        (1,1)                     number of scans in the session     
%       'opt'         structure array           (for info. /DOC/opt.doc)
%       'c'           (1,1)                     velocity of light in vacuum
%       rad2mas       (1,1)                     coefficient for conversion from radians to milli arc second 
%
%   Output:
%       'Apwxpol'      (nobserv x number of xpol pwlo unknowns)  design matrix for xpol 
%       'T'            structure array                           estimation intervals for dUT1, xpol, ypol, nutdx, and nutdy
%       'n_unk_xpol'   (1,1)                                     number of xpol estimation intervals in a session  
%       'Hxpol'        (size(Apwxpol,2)-1,size(Apwxpol,2))       design matrix for the xpol constraints as pseudo observations  
%       'Phxpol'       (size(Apwxpol,2)-1,size(Apwxpol,2)-1)     weight matrix of the xpol constraints
%       'oc_hxpol'     (size(Hxpol,1),1)                         o-c vector for the xpol constraints              
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
function [Apwxpol,T,Hxpol,Phxpol,oc_hxpol] = ahp_xpol(scan,mjd0,opt,c,rad2mas,obs_mjd)

int_xpol = opt.xpol.int;

mjd1 = min([scan.mjd]); % The time of the first scan in mjd [day]
mjd2 = max([scan.mjd]); % The time of the last scan in mjd [day]

t1 = (mjd1-mjd0)*24*60; % [minute] mjd0(midnight)--------mjd1(start)------(Session)------mjd2(end)-----
t2 = (mjd2-mjd0)*24*60; % [minute]
%t =  ([scan.mjd] - mjd0)*24*60; % [minute]

% determining the xpol estimation intervals and puting observation times in to them 
t10 = floor(t1/int_xpol)*int_xpol;
t20 = ceil(t2/int_xpol)*int_xpol;

[T.xpol] = t10:int_xpol:t20; % The estimation intervals for xpol
n_unk_xpol = length(T.xpol)-1; % The number of estimation intervals for xpol

% number of scans per clock estimate interval [1]
stm_xpol(1,n_unk_xpol) = 0; 
for inter = 1:n_unk_xpol
    if inter == n_unk_xpol
        tst= (obs_mjd >= T.xpol(inter)).*(obs_mjd <= T.xpol(inter+1));
    else
        tst= (obs_mjd >= T.xpol(inter)).*(obs_mjd < T.xpol(inter+1));
    end
	stm_xpol(inter) = sum(tst); % number of scans per xpol estimation interval [1]
end

temp = [scan.obs];

k = 0; Apwxpol(length(obs_mjd),n_unk_xpol+1)=0;
for inter = 1:n_unk_xpol % number of xpol estimate intervals in a session 1-24
    for iobs = 1:stm_xpol(inter) % number of observations in a xpol est. interval 241 229 
        k = k + 1;
        Apwxpol(k,inter) = (1-(obs_mjd(k)-T.xpol(inter))/(T.xpol(inter+1)-T.xpol(inter)))...
            *temp(k).ppol(1)*c*100*(1/rad2mas); % [sec/rad ----cm/mas]
        Apwxpol(k,inter+1) = ((obs_mjd(k)-T.xpol(inter))/(T.xpol(inter+1)-T.xpol(inter)))...
            *temp(k).ppol(1)*c*100*(1/rad2mas); % [sec/rad ----cm/mas]
    end
end

% FORMING THE CONSTRAINTS of XPOL
coef_xpol = opt.xpol.coef;

mat = diag(ones(1,n_unk_xpol+1)) - diag(ones(1,n_unk_xpol),1);
Hxpol = mat(1:n_unk_xpol,1:n_unk_xpol+1);
Phxpol = diag(ones(1,n_unk_xpol).*1./coef_xpol^2);

% FORMING THE O-C VECTOR FOR THE CONSTRAINTS
oc_hxpol(size(Hxpol,1),1) = 0; % o-c vector for the xpol constraints
