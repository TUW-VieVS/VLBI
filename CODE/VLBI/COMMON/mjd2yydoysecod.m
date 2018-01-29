% This function calculates the year (yr), the day of year (doy) and the 
% 5-digit seconds of day (secod) as required in SINEX output from modified 
% julian date.
%
% input
%       mjd     [n,1] or [1,n]       number
%       
% output
%       yyDoySecod  [n,3]       matrix
%          -yr      (:,1)       double      year (eg. 2008)     %%%% 2digits of year (eg. 2008=>08)
%          -doy     (:,2)       double      day of year [0, 365]
%          -secod   (:,3)       double      seconds of day    
%
%
% created
%   Matthias Madzak, 2 August 2010
%

function [yyDoySecod]=mjd2yydoysecod(mjd)

mjd=double(mjd);

% make size=[n,1] if size=[1,n]
if size(mjd, 2)>size(mjd, 1)
    mjd=mjd';
end

% preallocating current date matrix
curDate=zeros(length(mjd), 6);

% get date from mjd
[curDate(:,1), curDate(:,2), curDate(:,3), curDate(:,4), curDate(:,5), curDate(:,6)] = mjd2date(mjd);

% calculate day of year and seconds of day
doyWithYear=datenum(curDate);
%doyWithYear=datenum([curDate(:,1),zeros(length(mjd),2),curDate(:,4:6)]);
doyOfYear=datenum([curDate(:,1), repmat([0,0,0,0,0], size(curDate,1),1)]); % This is actually the serial date nr of start of year (month, day=0)
doy=floor(doyWithYear-doyOfYear);  % day of year of current time
% secod(:,1)=curDate(:,4)*60*60+curDate(:,5)*60+curDate(:,6);       % seconds of day of current time
secod=curDate(:,4)*60*60+curDate(:,5)*60+curDate(:,6);       % seconds of day of current time  (corr. by hana 18May11) 

yy=curDate(:,1);

% get 2 digits of year(s)
%yy=num2str(curDate(:,1));
%yy=yy(:, 3:end);

% make one output argument (matrix)
yyDoySecod=[yy, doy, secod];

