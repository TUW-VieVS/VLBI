% ************************************************************************
%   Description:
%   function to form the design matrix for the troposphere north gradients 
%   as piecewise linear ofset functions
%
%   Reference: 
%
%   Input:										
%       'obs_per_stat'        structure array     (for info. /DOC/obs_per_stat.doc)
%       'nobserv'             (1,1)               number of observation in the session
%	    'nob'			vector (1 x n)	   ordered number of observation per station
%       'n_unk'               structure array     number of estimation intervals for clocks, zwd, ngr, egr, xyz    
%       'T_'                  structure array     estimation epochs for clocks, zwd, ngr, egr, xyz
%
%   Output:
%       'Apwngr'      (nobserv,number of tropo. north gradient pwlo unknowns)  design matrix for north gradients           
%      
%   External calls: 	
%   
%   Coded for VieVS: 
%   12 May 2009 by Kamil Teke
%
%   Revision: 
%   05 Dec 2009 by Kamil Teke: header added
%   01 Sep 2010 by J. Boehm: gradient formulation changed from MacMillan (1995) to Chen and Herring (1997)
%   10 Oct 2016 by A. Girdiuk: code-optimized: loop removed, 
%							   function assignment is changed in accordance with stwisepar.m
% ************************************************************************ 
function [Apwngr] = apw_ngr(obs_per_stat,nobserv,nob,n_unk,T_)

Apwngr(nobserv,n_unk.ngr+1) = 0;

minute = obs_per_stat.minute;
% mf = obs_per_stat.mf;
first = obs_per_stat.first;
no_int_ngr = obs_per_stat.no_int_ngr;
%total = obs_per_stat.total;
%nob = obs_per_stat.nob;
az = obs_per_stat.az;
zd = obs_per_stat.zd;
el = pi/2-zd;

% number of scans per wet zenith delay estimate interval [1]

% for inter = 1:n_unk.ngr
%     if inter == n_unk.ngr
%         tst = (minute >= T_.ngr(inter)).*(minute <= T_.ngr(inter+1));
%     else
%         tst = (minute >= T_.ngr(inter)).*(minute < T_.ngr(inter+1));
%     end
%     stm_ngr(inter) = sum(tst); % number of observation to the corresponding source in the 'inter' estimation interval
% end

k = 0;
for inter = 1:n_unk.ngr % number of north gradient estimate intervals in a session
    for iobs = 1:n_unk.stm_ngr(inter) % number of observations of the station in a scan
        k = k + 1;
        Apwngr(nob(k),inter) = first(k)*(1-(minute(k)-T_.ngr(no_int_ngr(k)))/...
            (T_.ngr(no_int_ngr(k)+1)-T_.ngr(no_int_ngr(k))))/(tan(el(k))*sin(el(k))+0.0032)*cos(az(k));
        Apwngr(nob(k),inter+1) = first(k)*((minute(k)-T_.ngr(no_int_ngr(k)))/...
            (T_.ngr(no_int_ngr(k)+1)-T_.ngr(no_int_ngr(k))))/(tan(el(k))*sin(el(k))+0.0032)*cos(az(k));
    end
end
