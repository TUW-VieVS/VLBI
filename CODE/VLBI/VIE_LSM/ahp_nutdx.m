% ************************************************************************
%   Description:
%   function to form the design matrix for the nutation dx estimates as
%   piecewise linear ofsets
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
%       'Apwnutdx'      (nobserv x number of nutation dx, pwlo unknowns)  design matrix for nutation 
%       'T'             structure array                                   estimation intervals for dUT1, xpol, ypol, nutdx, and nutdy
%       'n_unk_nutdx'   (1,1)                                             number of nutation dx estimation intervals in a session  
%       'Hnutdx'        (size(Apwnutdx,2)-1,size(Apwnutdx,2))             design matrix for the nutation dx constraints as pseudo observations  
%       'Phnutdx'       (size(Apwnutdx,2)-1,size(Apwnutdx,2)-1)           weight matrix of the nutation constraints
%       'oc_hnutdx'     (size(Hnutdx,1),1)                                o-c vector for the nutation dx constraints              
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
function [Apwnutdx,T,Hnutdx,Phnutdx,oc_hnutdx] = ahp_nutdx(scan,mjd0,opt,T,c,rad2mas,obs_mjd)

int_nutdx = opt.nutdx.int;

mjd1 = min([scan.mjd]); % The time of the first scan in mjd [day]
mjd2 = max([scan.mjd]); % The time of the last scan in mjd [day]

t1 = (mjd1-mjd0)*24*60; % [minute] mjd0(midnight)--------mjd1(start)------(Session)------mjd2(end)-----
t2 = (mjd2-mjd0)*24*60; % [minute]
t =  ([scan.mjd] - mjd0)*24*60; % [minute]

% determining the nutation dx estimation intervals and puting observation times in to them 
t10 = floor(t1/int_nutdx)*int_nutdx;
t20 = ceil(t2/int_nutdx)*int_nutdx;

T.nutdx = t10:int_nutdx:t20; % The estimation intervals for nutdx
n_unk_nutdx = length(T.nutdx)-1; % The number of estimation intervals for nutation dx

% number of scans per clock estimate interval [1]
stm_nutdx(1,n_unk_nutdx) = 0; 
for inter = 1:n_unk_nutdx
    if inter == n_unk_nutdx
        tst= (obs_mjd >= T.nutdx(inter)).*(obs_mjd <= T.nutdx(inter+1));
    else
        tst= (obs_mjd >= T.nutdx(inter)).*(obs_mjd < T.nutdx(inter+1));
    end
	stm_nutdx(inter) = sum(tst); % number of scans per xpol estimation interval [1]
end

temp = [scan.obs];

k = 0; Apwnutdx(length(obs_mjd),n_unk_nutdx+1)=0;
for inter = 1:n_unk_nutdx % number of nutation dx estimate intervals in a session 
    for iobs = 1:stm_nutdx(inter) % number of observations in a nutation dx est. interval 
        k = k + 1;
        Apwnutdx(k,inter) = (1-(obs_mjd(k)-T.nutdx(inter))/(T.nutdx(inter+1)-T.nutdx(inter)))...
            *temp(k).pnut(1)*c*100*(1/rad2mas); % [sec/rad ----cm/mas]
        Apwnutdx(k,inter+1) = ((obs_mjd(k)-T.nutdx(inter))/(T.nutdx(inter+1)-T.nutdx(inter)))...
            *temp(k).pnut(1)*c*100*(1/rad2mas); % [sec/rad ----cm/mas];
    end
end

% FORMING THE CONSTRAINTS of nutation dx
coef_nutdx = opt.nutdx.coef; 

mat = diag(ones(1,n_unk_nutdx+1)) - diag(ones(1,n_unk_nutdx),1);
Hnutdx = mat(1:n_unk_nutdx,1:n_unk_nutdx+1);
Phnutdx = diag(ones(1,n_unk_nutdx).*1./coef_nutdx^2);

% FORMING THE O-C VECTOR FOR THE CONSTRAINTS
oc_hnutdx(size(Hnutdx,1),1) = 0; % o-c vector for the nutation dx constraints
