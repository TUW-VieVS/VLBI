% Purpose  
%   Check if the source is up at the station.        
% History  
%   2010-04-06   Jing SUN   Created
% 
% Changes:
%   2015-07-14, A. Hellerschmied: Bug-fix: sign of x-Angle (XYEW antenna mount) reversed
%   2016-05-11, M. Schartner: logical short-circuiting


function [lup] = zlup(mjd, az, el, ha, dc, staid, station, cutel)

% downtime
for idn = 1 : station(staid).downum  
    if ((mjd >= station(staid).downstart(idn)) && (mjd <= station(staid).downend(idn)))
        lup = false;
        return;
    end
end

% axis limit
if (strcmp(station(staid).axis, 'AZEL'))
    unaz = az;
    while (unaz < station(staid).lim11)
        unaz = unaz + 2 * pi;
    end
    while (unaz > station(staid).lim12)
        unaz = unaz - 2 * pi;
    end
    if ((unaz > station(staid).lim11) && (unaz < station(staid).lim12))
        lupaz = true;
    else
        lupaz = false;
    end
    if ((el > station(staid).lim21) && (el < station(staid).lim22))
        lupel = true;
    else 
        lupel = false;
    end
    lup = lupaz & lupel;
elseif (strcmp(station(staid).axis(1:4), 'HADC'))
    if ((ha > station(staid).lim11) && (ha < station(staid).lim12) && (dc > station(staid).lim21) && (dc < station(staid).lim22))
        lup = true;
    else 
        lup = false;
    end
elseif (strcmp(station(staid).axis(1:4), 'XYEW'))
    cel = cos(el);
    sel = sin(el);
    caz = cos(az);
    saz = sin(az);
    % x30 = -atan2(cel*caz,sel);
    x30 = atan2(cel*caz,sel);
    y30 = asin(cel*saz);
    if ((x30 > station(staid).lim11) && (x30 < station(staid).lim12) && (y30 > station(staid).lim21) && (y30 < station(staid).lim22))
        lup = true;
    else 
        lup = false;
    end
end      

% cutoff elevation
if (el > cutel) 
    lupc = true;
else
    lupc = false;
end
lup = lup & lupc;

% mask
hmasknum = station(staid).hmasknum;
hmask    = station(staid).hmask;
if ((lup == true) && (hmasknum > 0))
    if (mod(hmasknum,2) ~= 0)
        hmasktype = 1;   % step functions
    else
        hmasktype = 2;   % line segments
    end
    for i = 1 : floor(hmasknum/2)   
        ib = i * 2;
        if ((az >= hmask(ib-1)) && (az <= hmask(ib+1)))
            if (hmasktype == 1)
                elmask = hmask(ib);
            elseif (hmasktype == 2)
                elmask = ((hmask(ib+2) - hmask(ib)) / (hmask(ib+1) - hmask(ib-1))) * (az - hmask(ib-1)) + hmask(ib);
            end
            break;
        end
    end
    lupm = (el > elmask);
    lup = lup & lupm;
end


