% Purpose  
%   Calculate scan with particular sources.
% History  
%   2010-04-06   Jing SUN   Created
%   2016-07-11   Matthias Schartner: checks subcon with rigorous model for az el ha dc
%   2016-09-07   Matthias Schartner: changes for singleCheckSource
%   2016-11-03   Matthias Schartner: singleCheckSource now saves az el ha dc in subcon


function [subcon] = spsource(srcid, station, twin, staobs, source, srcobs, obsmode, PARA, jpl)

% initialization
stanum = length(station);
srcnum = length(source);
subcontmp4.nscan=0;

% endmjdmax
endmjdmax = 0.0;
for ista = 1 : stanum
    if (staobs(ista).endmjd > endmjdmax)
        endmjdmax = staobs(ista).endmjd;
    end
end

% check the visibility at stations
downstanum = 0;
downstasn(1:stanum) = 0; 
ra = source(srcid).ra;
de = source(srcid).de;
for ista = 1 : stanum
    endmjd = staobs(ista).endmjd;
    [az, el, ha, dc] = zazel_s(endmjd, station(ista).llh(1), station(ista).llh(2), ra, de);
    [lup] = zlup(endmjd, az, el, ha, dc, ista, station, PARA.MIN_CUTEL);
    if (lup == 1)
        downstanum = downstanum + 1;
        downstasn(downstanum) = ista;
    end
end    

% loop
if (downstanum < 2)
    return;
end
subcontmp.nscan = 1;
subcontmp.scan(1).srcid    = srcid;
subcontmp.scan(1).startmjd = 0.0;
subcontmp.scan(1).nsta     = downstanum;
for ista = 1 : subcontmp.scan(1).nsta
    subcontmp.scan(1).sta(ista).staid = downstasn(ista);
    subcontmp.scan(1).sta(ista).slewtime = 0.0;
    subcontmp.scan(1).sta(ista).startmjd = 0.0;
    subcontmp.scan(1).sta(ista).az       = 0.0;
    subcontmp.scan(1).sta(ista).el       = 0.0;
    subcontmp.scan(1).sta(ista).ha       = 0.0;
    subcontmp.scan(1).sta(ista).dc       = 0.0;
    subcontmp.scan(1).sta(ista).duration = 0.0;
    subcontmp.scan(1).sta(ista).endmjd   = 0.0;
end
[subcontmp1] = sstartmjd(source, station, twin, staobs, subcontmp, PARA);
if (subcontmp1.nscan == 0)
    return;
end
[subcon] = sduration(source, station, twin, obsmode, subcontmp1, PARA);
if (subcon.nscan  == 0)
    return;
end

[flag, staid,subcon] = singleCheckSource(source, station, subcon, PARA, jpl);
while flag ~= 1
    subcon.scan.nsta = subcon.scan.nsta-1;
    boolStaId = find(staid==[subcon.scan.sta.staid]);
    subcon.scan.sta(boolStaId)=[];
    [subcon] = sduration(source, station, twin, obsmode, subcon, PARA);
    if (subcon.nscan  == 0)
        return;
    end
    [flag, staid,subcon] = singleCheckSource(source, station, subcon, PARA, jpl);
end


