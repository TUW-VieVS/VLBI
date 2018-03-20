function y=date2mjd(ep);
%DATE2MJD Converts given epoch into modified Julian day.
% y=date2mjd(ep) computes the modified Julian date corresponding to the 
% epoch ep = [year month day hour min seconds]; if only [year month day]
% are passed the epoch is 00:00:00 on that day.
%
% If ep is an nx3 or nx6 array, each row of ep represents an epoch
% which is converted; the output argument is then an nx1 vector
% of modified Julian days.
%
% Note: if the year is lower than 100, it is assumed to be 19yy if yy>=80
%       and 20xx if yy<80.
%
% See also mjd2date, mjd2week, week2doy, week2date, week2mjd, doy2week.

% Created:  2006-12-07 AW
% Modified: 2006-12-14 AW       Handle multiple epoch input
%           2007-01-22 AW       Handling of year < 100

% (1) check input
if nargin < 1                    % epoch defined?
   error 'Too few arguments.'
end

h = size(ep);                    % is it nx3 or nx6?
if ~isnumeric(ep) || length(h) > 2 || ~ismember(h(2),[3 6])
   error 'ep must be an nx6 or nx3 array, n>=1.'
end

if h(2)==3                  % if needed, append 00:00:00
   ep(:,4:6)=0;
end

i = find(ep(:,1)<100);
if ~isempty(i)
   h = ep(i,1);
   j = find(h < 80);
   if ~isempty(j)
      h(j)=h(j)+100;
   end
   ep(i,1)=h+1900;
end

% (2) compute mjd
i = find(ep(:,2)<=2);       % modify year and month as required
ep(i,1)=ep(i,1)-1;
ep(i,2)=ep(i,2)+12;

                           % compute fraction of day   
frcday = ((ep(:,6)/60 + ep(:,5))/60 + ep(:,4))/24;   

                           % convert to Julian day
jd = fix(365.25*ep(:,1)) + fix(30.6001*(ep(:,2)+1)) + ep(:,3) ...
     + frcday + 1720981.5;
     
y  = jd - 2400000.5;       % and convert to modified Julian day
% --- Done.
