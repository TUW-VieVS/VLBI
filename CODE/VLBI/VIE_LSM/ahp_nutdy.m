% ************************************************************************
%   Description:
%   function to form the design matrix for the nutation dy estimates as
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
%       'Apwnutdy'      (nobserv x number of nutation dx, pwlo unknowns)   design matrix for nutation dy
%       'T'             structure array                                    estimation intervals for dUT1, xpol, ypol, nutdx, and nutdy
%       'n_unk_nutdy'   (1,1)                                              number of nutation dy estimation intervals in a session  
%       'Hnutdy'        (size(Apwnutdy,2)-1,size(Apwnutdy,2))              design matrix for the nutation dy constraints as pseudo observations  
%       'Phnutdy'       (size(Apwnutdy,2)-1,size(Apwnutdy,2)-1)            weight matrix of the nutation dy constraints
%       'oc_hnutdy'     (size(Hnutdy,1),1)                                 o-c vector for the nutation dy constraints              
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
function [Apwnutdy,T,Hnutdy,Phnutdy,oc_hnutdy] = ahp_nutdy(scan,mjd0,opt,T,c,rad2mas,obs_mjd)

int_nutdy = opt.nutdy.int;

mjd1 = min([scan.mjd]); % The time of the first scan in mjd [day]
mjd2 = max([scan.mjd]); % The time of the last scan in mjd [day]

t1 = (mjd1-mjd0)*24*60; % [minute] mjd0(midnight)--------mjd1(start)------(Session)------mjd2(end)-----
t2 = (mjd2-mjd0)*24*60; % [minute]
t =  ([scan.mjd] - mjd0)*24*60; % [minute]

% determining the celestial pole offsets dy (~dPSI) estimation intervals and puting observation times in to them 
t10 = floor(t1/int_nutdy)*int_nutdy;
t20 = ceil(t2/int_nutdy)*int_nutdy;

T.nutdy = t10:int_nutdy:t20; % The estimation intervals for nutation dy
n_unk_nutdy = length(T.nutdy)-1; % The number of estimation intervals for nutation dy

% number of scans per clock estimate interval [1]
stm_nutdy(1,n_unk_nutdy) = 0; 
for inter = 1:n_unk_nutdy
    if inter == n_unk_nutdy
        tst= (obs_mjd >= T.nutdy(inter)).*(obs_mjd <= T.nutdy(inter+1));
    else
        tst= (obs_mjd >= T.nutdy(inter)).*(obs_mjd < T.nutdy(inter+1));
    end
	stm_nutdy(inter) = sum(tst); % number of scans per xpol estimation interval [1]
end

temp = [scan.obs];

k = 0; Apwnutdy(length(obs_mjd),n_unk_nutdy+1)=0;
for inter = 1:n_unk_nutdy % number of PSI estimate intervals in a session
    for iobs = 1:stm_nutdy(inter) % number of observations in a PSI est. interval
        k = k + 1;
        Apwnutdy(k,inter) = (1-(obs_mjd(k)-T.nutdy(inter))/(T.nutdy(inter+1)-T.nutdy(inter)))...
            *temp(k).pnut(2)*c*100*(1/rad2mas); % [sec/rad ----cm/mas]
        Apwnutdy(k,inter+1) = ((obs_mjd(k)-T.nutdy(inter))/(T.nutdy(inter+1)-T.nutdy(inter)))...
            *temp(k).pnut(2)*c*100*(1/rad2mas); % [sec/rad ----cm/mas]
    end
end

% FORMING THE CONSTRAINTS of PSI
coef_nutdy = opt.nutdy.coef;

mat = diag(ones(1,n_unk_nutdy+1)) - diag(ones(1,n_unk_nutdy),1);
Hnutdy = mat(1:n_unk_nutdy,1:n_unk_nutdy+1);
Phnutdy = diag(ones(1,n_unk_nutdy).*1./coef_nutdy^2);

% FORMING THE O-C VECTOR FOR THE CONSTRAINTS
oc_hnutdy(size(Hnutdy,1),1) = 0; % o-c vector for the nutation dy constraints
