% ************************************************************************
%   Description:
%   function to form the design matrix for the zenith wet delays as  
%   piecewise linear ofset functions
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
%       'Azwd'      (nobserv,number of zenith wet delay pwlo unknowns)  design matrix for pwlo zwd estimates           
%      
%   External calls: 	
%   
%   Coded for VieVS: 
%   12 May 2009 by Kamil Teke
%
%   Revision: 
%   05 Dec 2009 by Kamil Teke: header added
%   10 Oct 2016 by A. Girdiuk: code-optimized: loop removed, 
%							   function assignment is changed in accordance with stwisepar.m
% ************************************************************************ 
function [Azwd] = apw_zwd(obs_per_stat,nobserv,nob,n_unk,T_)

Azwd(nobserv,n_unk.zwd+1) = 0;

minute = obs_per_stat.minute;
mf = obs_per_stat.mf;
first = obs_per_stat.first;
no_int_zwd = obs_per_stat.no_int_zwd;
%total = obs_per_stat.total;

% number of scans per wet zenith delay estimate interval [1]

% for inter = 1:n_unk.zwd
%     if inter == n_unk.zwd
%         tst = (minute >= T_.zwd(inter)).*(minute <= T_.zwd(inter+1));
%     else
%         tst = (minute >= T_.zwd(inter)).*(minute < T_.zwd(inter+1));
%     end
%     stm_zwd(inter) = sum(tst); % number of observation to the corresponding source in the 'inter' estimation interval
% end

k = 0;
for inter = 1:n_unk.zwd % number of zwd estimate intervals in a session
    for iobs = 1:n_unk.stm_zwd(inter) % number of observations of the station in a scan
        k = k + 1;
        Azwd(nob(k),inter) = first(k)*(1-(minute(k)-T_.zwd(no_int_zwd(k)))/(T_.zwd(no_int_zwd(k)+1)-T_.zwd(no_int_zwd(k))))*mf(k);
        Azwd(nob(k),inter+1) = first(k)*((minute(k)-T_.zwd(no_int_zwd(k)))/(T_.zwd(no_int_zwd(k)+1)-T_.zwd(no_int_zwd(k))))*mf(k);
    end
end
