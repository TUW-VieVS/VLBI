% Purpose  
%   Check source position during the observation.        
%   return:
%       flag: 1 means subcon is ok
%       flag: 2 source not visible at start of the observation
%       flag: 3 source not visible at end of the observation
%       flag: 4 source not visible during the observation
%       flag: 5 problems with cable wrap consistency for AZEL antenna
%       
%       staid: station ID where the problem occures
%
% History  
%   2016-06-30 M. Schartner created
%   2016-11-03 M. Schartner: singleCheckSource now saves az el ha dc in subcon

function [flag,staid,subcon] = singleCheckSource(source, station, subcon, PARA, jpl)
    
stepsec = 5;   % [sec]

% load eph
mjdnum = 0;
for iscan = 1 : subcon.nscan
    for ista = 1 : subcon.scan(iscan).nsta
        mjdnum = mjdnum + 1;
        mjdsn1(mjdnum) = subcon.scan(iscan).startmjd;
        mjdnum = mjdnum + 1;
        mjdsn1(mjdnum) = subcon.scan(iscan).sta(ista).endmjd;
    end
end

[x, y] = sort(mjdsn1);
for i = 1 : mjdnum
    mjdsn2(i) = mjdsn1(y(i));
end 
leap = tai_utc(mjdsn2);   
for i = 1 : mjdnum
    tt(i,1) = mjdsn2(i) + (32.184 + leap(i,1)) / 86400;   
end
[ephem] = get_vEarth(tt, jpl);

% check source position
mjdn1 = 0;
for iscan = 1 : subcon.nscan
    srcid = subcon.scan(iscan).srcid;
    startmjd = subcon.scan(iscan).startmjd;
    ra = source(srcid).ra;
    de = source(srcid).de;
    for ista = 1 : subcon.scan(iscan).nsta 
        staid  = subcon.scan(iscan).sta(ista).staid;
        unazs  = subcon.scan(iscan).sta(ista).az;
        xtrs   = station(staid).xyz(1); 
        ytrs   = station(staid).xyz(2); 
        lon    = station(staid).llh(1); 
        lat    = station(staid).llh(2);
        endmjd = subcon.scan(iscan).sta(ista).endmjd; 
        % start of the observation
        mjdn1 = mjdn1 + 1;
        for i = 1 : mjdnum
            if (y(i) == mjdn1)
                mjdn2 = i;
                break;
            end
        end
        [azs, el, ha, dc] = zazel_r(startmjd, xtrs, ytrs, lon, lat, ra, de, mjdn2, ephem);
        [lups] = zlup(startmjd, azs, el, ha, dc, staid, station, PARA.MIN_CUTEL);
        if ~lups
            flag = 2;
            return
        end
        
        diffAz = subcon.scan(iscan).sta(ista).az-azs;
        diffAzRound = round(diffAz/(2*pi));
        unazs = azs+diffAzRound*(2*pi);
        
        subcon.scan(iscan).sta(ista).az = unazs;
        subcon.scan(iscan).sta(ista).el = el;
        subcon.scan(iscan).sta(ista).ha = ha;
        subcon.scan(iscan).sta(ista).dc = dc;

        
        % end of the observation
        mjdn1 = mjdn1 + 1;
        for i = 1 : mjdnum
            if (y(i) == mjdn1)
                mjdn2 = i;
                break;
            end
        end
        [aze, el, ha, dc] = zazel_r(endmjd, xtrs, ytrs, lon, lat, ra, de, mjdn2, ephem);
        [lupe] = zlup(endmjd, aze, el, ha, dc, staid, station, PARA.MIN_CUTEL); 
        if ~lupe
            flag = 3;
            return
        end
        
        % during the observation
        mjd = startmjd + stepsec/86400;  
        while (mjd <= endmjd)
            [az, el, ha, dc] = zazel_s(mjd, lon, lat, ra, de);
            [lup] = zlup(mjd, az, el, ha, dc, staid, station, PARA.MIN_CUTEL); 
            if ~lup
                flag = 4;
                return
            end
            mjd = mjd + stepsec/86400;
        end
        
        % check cable wrap consistency for AZEL antenna
        if strcmp(station(staid).axis(1:4), 'AZEL')
            [~, dfaz] = sdfaz(station(staid).lim11, station(staid).lim12, unazs, aze);
            if (dfaz > pi)
                flag = 5;
                return
            end
        end
        
    end
end
flag = 1;
staid = 0;