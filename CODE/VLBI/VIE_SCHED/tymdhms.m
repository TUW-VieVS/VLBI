% Purpose  
%   Convert modified julian date to [yy,mm,dd,hh,min,secf].
% History  
%   2010-04-06   Jing SUN   Created
%


function [yy, mm, dd, hh, min, secf] = tymdhms(mjd)

% hour, minute, second
hh=floor((mjd-floor(mjd))*24);
min=floor((((mjd-floor(mjd))*24)-hh)*60);
secf=(((((mjd-floor(mjd))*24)-hh)*60)-min)*60; %float

% year, month, day
utcp1 = mjd + 2400000.5;
jd = floor(utcp1);
if (hh < 12)
    jd = jd + 1;
end
lx = jd + 68569;
nx = floor(4 * lx / 146097);
lx = lx - floor((146097 * nx + 3) / 4);                                                 
yy = floor(4000 * (lx + 1) / 1461001);                                              
lx = lx - floor(1461 * yy / 4) + 31;                                                  
mm = floor(80 * lx / 2447);                                                      
dd = lx - floor(2447 * mm / 80);                                                      
lx = floor(mm / 11);                                                              
mm = floor(mm + 2 - 12 * lx);                                                       
yy = floor(100 * (nx - 49) + yy + lx);        


