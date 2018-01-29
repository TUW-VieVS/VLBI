% Purpose  
%   Calculate slewing time.        
% History  
%   2010-04-06   Jing SUN   Created
% 
% Changes:
%   2015-07-14, A. Hellerschmied: Bug-fix: sign of x-Angle (XYEW antenna mount) reversed
%   2016-07-11, M. Schartner: bugfix, now uses acceleration from station struct
%   2016-11-04, M. Schartner: new varargin parameter:
% 						      If varagin is 'directSlew' azel Antennas are forced to slew  
%							  to this exact azimuth not considering fastest way. 
function [st, unaz] = sslew(station, staobs, az, el, ha, dc, staid, PARA, varargin)

if (strcmp(station(staid).axis(1:4), 'AZEL'))
    % az slewing time
    if isequal(varargin,{'directSlew'}) 
        unaz = az;
        dfaz = abs(az-staobs(staid).az);
    else
        [unaz, dfaz] = sdfaz(station(staid).lim11, station(staid).lim12, staobs(staid).az, az);
    end
    
    if isempty(unaz)
        st = 999999;
        unaz = az;
    else
        
        [staz] = ssva2t(dfaz*180/pi, station(staid).rate1, station(staid).acc1);
        staz = staz + station(staid).c1;   % [sec]    
        % el slewing time
        dfel = abs(staobs(staid).el - el) * 180 / pi;
        [stel] = ssva2t(dfel, station(staid).rate2, station(staid).acc2);
        stel = stel + station(staid).c2;   % [sec]
        % slewing time
        st = ceil(max(staz, stel));
    end
elseif (strcmp(station(staid).axis(1:4), 'HADC'))
    % ha slewing time
    hanow  = staobs(staid).ha;
    hanew  = ha;
    dfha = abs(hanew - hanow);
    [stha] = ssva2t(dfha*180/pi, station(staid).rate1, station(staid).acc1);  % Calculate the slew time considering acceleration and deceleration [sec]. 
    stha = stha + station(staid).c1;   % [sec]
    % dc slewing time
    dcnow  = staobs(staid).dc;
    dcnew  = dc;
    dfdc = abs(dcnew - dcnow);
    [stdc] = ssva2t(dfdc*180/pi, station(staid).rate2, station(staid).acc2);
    stdc = stdc + station(staid).c2;   % [sec] 
    % slewing time
    st = ceil(max(stha, stdc));
    unaz = az; 
    
elseif (strcmp(station(staid).axis(1:4), 'XYEW'))
    % xnow, ynow
    cel = cos(staobs(staid).el);
    sel = sin(staobs(staid).el);
    caz = cos(staobs(staid).az);
    saz = sin(staobs(staid).az);
    % xnow = -atan2(cel*caz,sel);
    xnow = atan2(cel*caz,sel);
    ynow = asin(cel*saz);
    % xnew, ynew
    cel = cos(el);
    sel = sin(el);
    caz = cos(az);
    saz = sin(az);
    % xnew = -atan2(cel*caz,sel);
    xnew = atan2(cel*caz,sel);
    ynew = asin(cel*saz);
    % x mount slewing time
    dfx = abs(xnew - xnow);
    [stx] = ssva2t(dfx*180/pi, station(staid).rate1, station(staid).acc1);
    stx = stx + station(staid).c1;   % [sec] 
    % y mount slewing time
    dfy = abs(ynew - ynow);
    [sty] = ssva2t(dfy*180/pi, station(staid).rate2, station(staid).acc2);
    sty = sty + station(staid).c2;   % [sec] 
    % slewing time
    st = ceil(max(stx, sty));
    unaz = az; 
end