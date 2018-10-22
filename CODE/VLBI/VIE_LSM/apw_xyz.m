% ************************************************************************
%   Description:
%   function to form the design matrix for the TRF coordinates of antennas 
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
%       'Apwx'      (nobserv,number of TRF dx coordinate pwlo unknowns) design matrix for antenna dx coordinates 
%       'Apwy'      (nobserv,number of TRF dy coordinate pwlo unknowns) design matrix for antenna dy coordinates
%       'Apwz'      (nobserv,number of TRF dz coordinate pwlo unknowns) design matrix for antenna dz coordinates
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
function [Apwx,Apwy,Apwz] = apw_xyz(obs_per_stat,nobserv,nob,n_unk,T_)

Apwx(nobserv,n_unk.xyz+1) = 0;
Apwy(nobserv,n_unk.xyz+1) = 0;
Apwz(nobserv,n_unk.xyz+1) = 0;

minute = obs_per_stat.minute;
first = obs_per_stat.first;
no_int_xyz = obs_per_stat.no_int_xyz;
%total = obs_per_stat.total;
%nob = obs_per_stat.nob;
dx = obs_per_stat.dx;
dy = obs_per_stat.dy;
dz = obs_per_stat.dz;

% number of scans per coordinate estimate interval [1]
% for inter = 1:n_unk.xyz
%     if inter == n_unk.xyz
%         tst = (minute >= T_.xyz(inter)).*(minute <= T_.xyz(inter+1));
%     else
%         tst = (minute >= T_.xyz(inter)).*(minute < T_.xyz(inter+1));
%     end
%     stm_xyz(inter) = sum(tst); % number of observation to the corresponding source in the 'inter' estimation interval
% end

k = 0;
for inter = 1:n_unk.xyz % number of coordinate estimation intervals in a session
    for iobs = 1:n_unk.stm_xyz(inter) % number of observations of the station in a scan
        k = k + 1;
        Apwx(nob(k),inter) = first(k)*(1-(minute(k)-T_.xyz(no_int_xyz(k)))/...
            (T_.xyz(no_int_xyz(k)+1)-T_.xyz(no_int_xyz(k))))*dx(k);
        Apwx(nob(k),inter+1) = first(k)*((minute(k)-T_.xyz(no_int_xyz(k)))/...
            (T_.xyz(no_int_xyz(k)+1)-T_.xyz(no_int_xyz(k))))*dx(k);
        
        Apwy(nob(k),inter) = first(k)*(1-(minute(k)-T_.xyz(no_int_xyz(k)))/...
            (T_.xyz(no_int_xyz(k)+1)-T_.xyz(no_int_xyz(k))))*dy(k);
        Apwy(nob(k),inter+1) = first(k)*((minute(k)-T_.xyz(no_int_xyz(k)))/...
            (T_.xyz(no_int_xyz(k)+1)-T_.xyz(no_int_xyz(k))))*dy(k);
        
        Apwz(nob(k),inter) = first(k)*(1-(minute(k)-T_.xyz(no_int_xyz(k)))/...
            (T_.xyz(no_int_xyz(k)+1)-T_.xyz(no_int_xyz(k))))*dz(k);
        Apwz(nob(k),inter+1) = first(k)*((minute(k)-T_.xyz(no_int_xyz(k)))/...
            (T_.xyz(no_int_xyz(k)+1)-T_.xyz(no_int_xyz(k))))*dz(k);
    end
end

