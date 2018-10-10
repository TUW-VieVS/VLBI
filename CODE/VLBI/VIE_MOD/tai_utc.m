% gives the difference TAI-UTC for time in MJD
% in accordance with IERS Bulletin C

function [tmu,varargout]=tai_utc(mjd,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: mjd >= 41317 (1. Jan. 1972)[days, UTC]
% Output: TAI-UTC [sec]

% last leap second: 1.Jan 2009 (Bul. C 37)

% last update: Feb. 2009, L. Plank 
%              May  2011, L. Plank: out of date actualized
%              Jul  2012, L. Plank: new leap second
%              Sep  2013, L. Plank: out of date actualized
%              Mar  2014, L. Plank: out of date actualized
%              Nov  2014, L. Plank: out of date actualized
%              Mar  2015, A. Hellerschmied: new leap second (2015-07-01), out of date actualized
%              2016-06-01, A. Hellerschmied: mjdmax changed to 2016-11-01
%			   2016-08-08, A. Girdiuk: varargin and varargout are added and should be used in couple only
%										varargin could be everything, for example, [tmu,zeit]=tai_utc(mjd,1) or [tmu,zeit]=tai_utc(mjd,'A')
%              2017-01-09, A. Hellerschmied: - mjdmax changed to 2017-06-01
%                                            - new leap second added (2017-01-01: + 37 sec)
%                                            - function revised, simplified and vectorized
%              2017-09-11, L. McCallum: mjdmax changed to 2017-12-31
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% no announce for leap second
mjdmax = modjuldat(2017,12,31,0);

% init.:
tmu = zeros(length(mjd), 1);

% display warning if time of interest is after Jan. 2016
if mjd > mjdmax
    disp('+++ please check for new leap second TAI-UTC tai_utc.m +++')
end

% for update add mjd of new leap second (1st value)
leap_sec_epochs = [57754, 57204,56109,54832,53736,51179,50630,50083,49534,49169,48804,48257,...
        47892,47161,46247,45516,45151,44786,44239,43874,43509,...
        43144,42778,42413,42048,41683,41499,41317];
    
% for update add new amount of tai-utc (1st value)    
sec = [37,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,...
       16,15,14,13,12,11,10];
   
for i = 1 : length(mjd)
    tmp     = sec(mjd(i) >= leap_sec_epochs);
    tmu(i)  = tmp(1);
end

% Optional output:
if ~isempty(varargin)
	varargout{1} = leap_sec_epochs;
end
