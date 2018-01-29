% Purpose  
%   Convert [yy,mm,dd,hh,min,sec] to modified julian time.
% History  
%   2010-04-06   Jing SUN   Created
%


function [mjd] = tmjd(yy, mm, dd, hh, min, sec)

if (mm <= 2)
   mm = mm + 12;
   yy = yy - 1;
end

if ((yy <= 1582) && (mm <= 10) && (dd <= 4))   % if date is before Oct 4, 1582
   b = 0;
else                                           % if date is after Oct 4, 1582   
   b = 2 - floor(yy / 100) + floor(yy / 400);
end

jd = floor(365.25d0 * (yy + 4716)) + floor(30.6001d0 * (mm + 1)) + dd + b - 1524.5d0;

mjd = jd - 2400000.5d0;

mjd = mjd + (hh * 3600.0d0 + min * 60.0d0 + sec) / 86400.0d0;


