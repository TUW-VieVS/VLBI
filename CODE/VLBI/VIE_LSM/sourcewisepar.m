% ************************************************************************
%   Description:
%   function to rearrange the data from scan based to source based
%
%   Reference:
%
%   Input:
%       'opt'               structure array         (for info. /DOC/opt.doc)
%       'scan'              structure array         (for info. /DOC/scan.doc)
%       'ntim'              (1,1)                   number of scans in a session
%       'mjdsource'  (num. of obs. to source,1)     epochs of the observations to the source in mjd
%       'isou'              (1,1)                   ith source
%       'mjd0'              (1,1)                   mjd of the 0:00 UTC of the first day of the session
%       'int_sou'           (1,optional)            estimation intervals of the source coordinates (RA-DE)
%
%   Output:
%       'per_source'       structure array     (for info. /DOC/per_stat.doc)
%       'n_unk'            structure array     number of estimation intervals for clocks, zwd, ngr, egr, xyz, sou
%       'T_'               structure array     estimation epochs for clocks, zwd, ngr, egr, xyz, sou
%
%   External calls:
%
%   Coded for VieVS:
%   25 Aug 2009 by Kamil Teke
%
%   Revision:
%   06 Dec 2009 by Kamil Teke: header added
%   03 May 2017 by A. Hellerschmied: Now compatible with satellite observations 
% ************************************************************************
function [per_source,n_unk,T_] = sourcewisepar(opt,scan,ntim,mjdsource,isou,mjd0,int_sou)

mjd1 = min(mjdsource); % The time of the first scan in mjd [day]
mjd2 = max(mjdsource); % The time of the last scan in mjd [day]

t1 = (mjd1-mjd0)*24*60; % [minute] mjd0(midnight)--------mjd1(start)------(Session)------mjd2(end)-----
t2 = (mjd2-mjd0)*24*60; % [minute]

% determining the source coor. estimation intervals and puting observation times in to them
t10 = floor(t1/int_sou)*int_sou;
t20 = ceil(t2/int_sou)*int_sou;

T_.source = t10:int_sou:t20; % The estimation intervals for source coor.
n_unk.sou = length(T_.source)-1; % The number of estimation intervals for source coor.

per_source.name = opt.source(isou).name; % The name of the source

i = 0; k = 0;
for itim = 1:ntim % number of scans per session
    for iobs = 1:scan(itim).nobs % number of observations per scan
        i = i + 1; % i : i. observation in the session
        iso = scan(itim).iso;
        if strcmp(scan(itim).obs_type, 'q')
            if iso == isou
                k = k + 1; % k : k. observation of the specific source
                per_source.mjd(k) = scan(itim).mjd;  % time of the observations per source [day]
                per_source.nob(k) = i; % row number of observation in o-c per source
                per_source.ra(k) = scan(itim).obs(iobs).psou(1);  % [cm/mas]
                per_source.de(k) = scan(itim).obs(iobs).psou(2);  % [cm/mas]
                per_source.minute(k) = (scan(itim).mjd - mjd0)*24*60; % The times of scans in minutes per source
                for h1 = 1 : n_unk.sou % n_unk.sou : number of source coor. estimation intervals in a session
                    if h1 == n_unk.sou
                        tst = (per_source.minute(k) >= T_.source(h1) && per_source.minute(k) <= T_.source(h1+1));
                    else
                        tst = (per_source.minute(k) >= T_.source(h1) && per_source.minute(k) < T_.source(h1+1));
                    end
                    if tst
                        per_source.no_int_rade(k) = h1; % h1. estimation interval that contains k. observation of the source
                    else
                    end
                end
            end
        end
    end
end

per_source.total = k;
