% ************************************************************************
%   Description:
%   function to form the design matrix for the dUT1 estimates as piecewise
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
%       'Apwdut1'      (nobserv x number of dUT1 pwlo unknowns)  design matrix for dUT1 
%       'T'            structure array                           estimation intervals for dUT1, xpol, ypol, nutdx, and nutdy
%       'n_unk_dut1'   (1,1)                                     number of dUT1 estimation intervals in a session  
%       'Hdut1'        (size(Apwdut1,2)-1,size(Apwdut1,2))       design matrix for the dUT1 constraints as pseudo observations  
%       'Phdut1'       (size(Apwdut1,2)-1,size(Apwdut1,2)-1)     weight matrix of the dut1 constraints
%       'oc_hdut1'     (size(Hdut1,1),1)                         o-c vector for the dUT1 constraints              
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
function [Apwdut1,T,Hdut1,Phdut1,oc_hdut1] = ahp_dut1(scan,mjd0,opt,T,c,rad2mas,obs_mjd)

int_dut1 = opt.dut1.int;

mjd1 = min([scan.mjd]); % The time of the first scan in mjd [day]
mjd2 = max([scan.mjd]); % The time of the last scan in mjd [day]

t1 = (mjd1-mjd0)*24*60; % [minute] mjd0(midnight)--------mjd1(start)------(Session)------mjd2(end)-----
t2 = (mjd2-mjd0)*24*60; % [minute]
t =  ([scan.mjd] - mjd0)*24*60; % [minute]

% determining the dut1 estimation intervals and puting observation epochs in to them 
t10 = floor(t1/int_dut1)*int_dut1;
t20 = ceil(t2/int_dut1)*int_dut1;

T.dut1 = t10:int_dut1:t20; % The estimation intervals for dut1
n_unk_dut1 = length(T.dut1)-1; % The number of estimation intervals for dut1

% number of scans per clock estimate interval [1]
stm_dut1(1,n_unk_dut1) = 0; 
for inter = 1:n_unk_dut1
    if inter == n_unk_dut1
        tst= (obs_mjd >= T.dut1(inter)).*(obs_mjd <= T.dut1(inter+1));
    else
        tst= (obs_mjd >= T.dut1(inter)).*(obs_mjd < T.dut1(inter+1));
    end
	stm_dut1(inter) = sum(tst); % number of scans per xpol estimation interval [1]
end

temp = [scan.obs];

k = 0; Apwdut1(length(obs_mjd),n_unk_dut1+1)=0;
for inter = 1:n_unk_dut1 % number of dut1 estimate intervals in a session
    for iobs = 1:stm_dut1(inter) % number of observations in a dut1 est. interval 
        k = k + 1;
        Apwdut1(k,inter) = (1-(obs_mjd(k)-T.dut1(inter))/(T.dut1(inter+1)-T.dut1(inter)))...
            *temp(k).ppol(3)*c*100*(1/rad2mas); % [sec/rad ----cm/mas]
        Apwdut1(k,inter+1) = ((obs_mjd(k)-T.dut1(inter))/(T.dut1(inter+1)-T.dut1(inter)))...
            *temp(k).ppol(3)*c*100*(1/rad2mas); % [sec/rad ----cm/mas]
    end
end

% FORMING THE CONSTRAINTS of dUT1
coef_dut1 = opt.dut1.coef;

mat = diag(ones(1,n_unk_dut1+1)) - diag(ones(1,n_unk_dut1),1);
Hdut1 = mat(1:n_unk_dut1,1:n_unk_dut1+1);
Phdut1 = diag(ones(1,n_unk_dut1).*1./coef_dut1^2);

% FORMING THE O-C VECTOR FOR THE CONSTRAINTS
oc_hdut1(size(Hdut1,1),1) = 0; % o-c vector for the dUT1 constraints
