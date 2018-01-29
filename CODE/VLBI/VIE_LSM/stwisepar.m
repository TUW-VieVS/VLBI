% ************************************************************************
%   Description:
%   function to rearrange data from scan based to station based (antenna)
%
%   Reference: 
%
%   Input:	
%       'per_stat'           structure array        variable prepared in loop over stattions in vie_lsm
%       'mjdstat'    (num. of obs. per stat,1)      epochs of the observations carried out by the antenna in mjd
%       'istat'             (1,1)                   ith station
%       'mjd0'              (1,1)                   mjd of the first day of the session at 0:00 UTC 
%       'int'               (1,optional)            estimation intervals of zwd, clk, egr, ngr, xyz(TRF)   
%
%   Output:
%       'obs_per_stat'     structure array     (for info. /DOC/obs_per_stat.doc)
%       'n_unk'            structure array     number of estimation intervals for clocks, zwd, ngr, egr, xyz, sou  
%       'T_'               structure array     estimation epochs for clocks, zwd, ngr, egr, xyz, sou
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   12 May 2009 by Kamil Teke
%
%   Revision: 
%   06 Dec 2009 by Kamil Teke: header added 
%   11 Oct 2011 by Hana Spicakova: amplitudes of seasonal variation in station positions
%   04 Oct 2013 by Hana Krasna: estimation of antenna axis offset added
%   05 Dec 2013 by Hana Krasna: APL regression coefficients added
%   10 Oct 2016 by A. Girdiuk: code-optimized:
%							   assignment of derivations is moved to loop over stations in vie_lsm:
%							   % ASSIGNING THE INPUT PARAMETERS TO EACH STATION
% ************************************************************************

% function to rearrange the data based on station
% Johannes Boehm, Kamil Teke, 12.05.2009
% mod. 12.05.2009 Kamil Teke (describtion added)

 function [obs_per_stat,n_unk,T_] = stwisepar(per_stat,mjdstat,istat,mjd0,int)

mjd1 = min(mjdstat); % The time of the first scan in mjd [day]
mjd2 = max(mjdstat); % The time of the last scan in mjd [day]

t1 = (mjd1-mjd0)*24*60; % [minute] mjd0(midnight)--------mjd1(start)------(Session)------mjd2(end)-----
t2 = (mjd2-mjd0)*24*60; % [minute]

% determining the zwd estimation intervals and puting observation times in to them 
t10 = floor(t1/int.zwd)*int.zwd;
t20 = ceil(t2/int.zwd)*int.zwd;

T_.zwd = t10:int.zwd:t20; % The estimation intervals for zwd
n_unk.zwd = length(T_.zwd)-1; % The number of estimation intervals for zwd

% determining the clock estimation intervals and puting observation times in to them 
t11 = floor(t1/int.clk)*int.clk;
t22 = ceil(t2/int.clk)*int.clk;

T_.clk = t11:int.clk:t22; % The estimation intervals for CLK
n_unk.clk = length(T_.clk)-1; % The number of estimation intervals for zwd

% determining the east gradients estimation intervals and puting observation times in to them 
t12 = floor(t1/int.egr)*int.egr;
t23 = ceil(t2/int.egr)*int.egr;

T_.egr = t12:int.egr:t23; % The estimation intervals for east gradients
n_unk.egr = length(T_.egr)-1; % The number of estimation intervals for east gradients

% determining the north gradients estimation intervals and puting observation times in to them 
t13 = floor(t1/int.ngr)*int.ngr;
t24 = ceil(t2/int.ngr)*int.ngr;

T_.ngr = t13:int.ngr:t24; % The estimation intervals for north gradients
n_unk.ngr = length(T_.ngr)-1; % The number of estimation intervals for north gradients

% determining the coordinate estimation intervals and puting observation times in to them 
t14 = floor(t1/int.xyz)*int.xyz;
t25 = ceil(t2/int.xyz)*int.xyz;

T_.xyz = t14:int.xyz:t25; % The estimation intervals for coordinates
n_unk.xyz = length(T_.xyz)-1; % The number of estimation intervals for coordinates

obs_per_stat.minute = (per_stat(istat).mjd-mjd0)*24*60;

minute = obs_per_stat.minute;

obs_per_stat.no_int_zwd=zeros(1,length(minute));
for inter = 1:n_unk.zwd
    if inter == n_unk.zwd
    tst = (minute >= T_.zwd(inter)).*(minute <= T_.zwd(inter+1));
    else
    tst = (minute >= T_.zwd(inter)).*(minute < T_.zwd(inter+1));
    end
    obs_per_stat.no_int_zwd = obs_per_stat.no_int_zwd + tst*inter; % number of observation to the corresponding source in the 'inter' estimation interval
    n_unk.stm_zwd(inter) = sum(tst); % number of observation to the corresponding source in the 'inter' estimation interval
end

obs_per_stat.no_int_clk=zeros(1,length(minute));
for inter = 1:n_unk.clk
    if inter == n_unk.clk
    tst = (minute >= T_.clk(inter)).*(minute <= T_.clk(inter+1));
    else
    tst = (minute >= T_.clk(inter)).*(minute < T_.clk(inter+1));
    end
    obs_per_stat.no_int_clk = obs_per_stat.no_int_clk + tst*inter; % number of observation to the corresponding source in the 'inter' estimation interval
    n_unk.stm_clk(inter) = sum(tst); % number of observation to the corresponding source in the 'inter' estimation interval
end

obs_per_stat.no_int_ngr=zeros(1,length(minute));
for inter = 1:n_unk.ngr
    if inter == n_unk.ngr
    tst = (minute >= T_.ngr(inter)).*(minute <= T_.ngr(inter+1));
    else
    tst = (minute >= T_.ngr(inter)).*(minute < T_.ngr(inter+1));
    end
    obs_per_stat.no_int_ngr = obs_per_stat.no_int_ngr + tst*inter; % number of observation to the corresponding source in the 'inter' estimation interval
    n_unk.stm_ngr(inter) = sum(tst); % number of observation to the corresponding source in the 'inter' estimation interval
end

obs_per_stat.no_int_egr=zeros(1,length(minute));
for inter = 1:n_unk.egr
    if inter == n_unk.egr
    tst = (minute >= T_.egr(inter)).*(minute <= T_.egr(inter+1));
    else
    tst = (minute >= T_.egr(inter)).*(minute < T_.egr(inter+1));
    end
    obs_per_stat.no_int_egr = obs_per_stat.no_int_egr + tst*inter; % number of observation to the corresponding source in the 'inter' estimation interval
    n_unk.stm_egr(inter) = sum(tst); % number of observation to the corresponding source in the 'inter' estimation interval
end

obs_per_stat.no_int_xyz=zeros(1,length(minute));
for inter = 1:n_unk.xyz
    if inter == n_unk.xyz
    tst = (minute >= T_.xyz(inter)).*(minute <= T_.xyz(inter+1));
    else
    tst = (minute >= T_.xyz(inter)).*(minute < T_.xyz(inter+1));
    end
    obs_per_stat.no_int_xyz = obs_per_stat.no_int_xyz + tst*inter; % number of observation to the corresponding source in the 'inter' estimation interval
    n_unk.stm_xyz(inter) = sum(tst); % number of observation to the corresponding source in the 'inter' estimation interval
end

obs_per_stat.first=per_stat(istat).first;

%obs_per_stat.total = length(obs_per_stat.first);

obs_per_stat.mf = per_stat(istat).mf;
obs_per_stat.az = per_stat(istat).az;
obs_per_stat.zd = per_stat(istat).zd;

obs_per_stat.dx = per_stat(istat).dx;
obs_per_stat.dy = per_stat(istat).dy;
obs_per_stat.dz = per_stat(istat).dz;

obs_per_stat.dAO = per_stat(istat).dAO;

obs_per_stat.drg = per_stat(istat).drg;

if isfield(per_stat,'pAcr')
	obs_per_stat.pAcr = per_stat(istat).pAcr;
	obs_per_stat.pAce = per_stat(istat).pAce;
    obs_per_stat.pAcn = per_stat(istat).pAcn;
    obs_per_stat.pAsr = per_stat(istat).pAsr;
    obs_per_stat.pAse = per_stat(istat).pAse;
    obs_per_stat.pAsn = per_stat(istat).pAsn;
end
