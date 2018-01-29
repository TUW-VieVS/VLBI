% Purpose
%   Calculate scan with fill-in mode.
% History
%   2010-04-06   Jing SUN   Created
%   2016-06-17 M. Schartner improved speed
%   2016-06-27 M. Schartner bugfix
%   2016-07-04 M. Schartner: added singleCheckSource
%	2016-07-11 M. Schartner: improved fill-in mode
%   2016-09-06 M. Schartner: changes for sched_manually
%   2016-09-07 M. Schartner: changes for singleCheckSource
%   2016-10-19 M. Schartner: PARA.PARA.MIN_STANUM_FI
%   2016-11-03 M. Schartner: singleCheckSource now saves az el ha dc in subcon
%   2017-05-29: M. Schartner: Bugfix

function [subcon_f2] = sscanfi2(station, twin, staobs, source, srcobs, obsmode, endmjdmaxfi2, PARA, jpl)

% initialization
stanum = length(station);
srcnum = length(source);
subcon_f2.nscan = 0;

% endmjdmax
endmjdmax = 0.0;
for ista = 1 : stanum
    if (staobs(ista).endmjd > endmjdmax)
        endmjdmax = staobs(ista).endmjd;
    end
end

% Station lon lat
llh = [station.llh];
lon = llh(1:3:end);
lat = llh(2:3:end);

% which Station ends bevore observation ends?
waitsec = (endmjdmaxfi2-[staobs.endmjd])*86400;
bool = waitsec>PARA.MIN_SCAN + 10 + PARA.SOURCE + PARA.TAPETM + PARA.IDLE + PARA.CALIBRATION;

if sum(bool)<PARA.MIN_STANUM_FI
    return
end

for isrc = 1 : srcnum
    srcsta(isrc).downstanum = 0;
    srcsta(isrc).downstasn(1:stanum) = 0;
    % check the visibility at stations
    ra = source(isrc).ra;
    de = source(isrc).de;
    endmjd = [staobs.endmjd];
    [az, el, ha, dc] = zazel_s(endmjd, lon, lat, ra, de);
    for staid = 1 : stanum
        if bool(staid)
            [lup] = zlup(endmjd(staid), az(staid), el(staid), ha(staid), dc, staid, station, PARA.MIN_CUTEL);
            if lup
                srcsta(isrc).downstanum = srcsta(isrc).downstanum + 1;
                srcsta(isrc).downstasn(srcsta(isrc).downstanum) = staid;
            end
        end
    end
end

% loop
num = 0;
for isrc = 1 : srcnum

    % check if the source is repeated in a short interval
    if (endmjdmaxfi2-srcobs(isrc).obsmjd)*1440 < PARA.MIN_SRCRP
        continue;
    end

    if (srcsta(isrc).downstanum < 2)
        continue;
    end

    if source(isrc).star
        continue;
    end

    subcontmp.nscan = 1;
    subcontmp.scan(1).srcid    = isrc;
    subcontmp.scan(1).startmjd = 0.0;
    subcontmp.scan(1).nsta     = srcsta(isrc).downstanum;
    newnum = 0;
    for ista = 1 : subcontmp.scan(1).nsta
        subcontmp.scan(1).sta(ista).staid = srcsta(isrc).downstasn(ista);
        [az, el, ha, dc] = zazel_s(staobs(subcontmp.scan(1).sta(ista).staid).endmjd, station(subcontmp.scan(1).sta(ista).staid).llh(1), station(subcontmp.scan(1).sta(ista).staid).llh(2), source(subcontmp.scan(1).srcid).ra, source(subcontmp.scan(1).srcid).de);
        [cat] = zcoverage(az, el);
        ifnewcat = 1;
        for icat = 1 : staobs(subcontmp.scan(1).sta(ista).staid).ncat
            if (staobs(subcontmp.scan(1).sta(ista).staid).catsn(icat)==cat)
                ifnewcat = 0;
                break;
            end
        end
        newnum = newnum + ifnewcat;
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
        continue;
    end
    [subcontmp2] = sduration(source, station, twin, obsmode, subcontmp1, PARA);
    if (subcontmp2.nscan  == 0)
        continue;
    end
    subcontmp4 = subcontmp2;

    % check for source structure study
    iforsi = 1;
    if (PARA.FORSI == 1)
        for iscan = 1 : subcontmp4.nscan
            [obs_ok(iscan), grout(iscan)]= calc_soustruc_pscan(subcontmp4.scan(iscan),source,station);
            subcontmp4.scan(iscan).grout = grout;
        end
        iforsi = sum(obs_ok(1:subcontmp4.nscan));
    end
    if (iforsi == 0)
        continue;
    end
    % check tmpendmjdmax
    while true
        tmpend = 0.0;
        for ista = 1 : subcontmp4.scan(1).nsta
            if(subcontmp4.scan(1).sta(ista).endmjd > tmpend)
                tmpend = subcontmp4.scan(1).sta(ista).endmjd;
                badsta = ista;
            end
        end
        if (tmpend <= endmjdmaxfi2)  %%%

            if subcontmp4.scan.nsta>=PARA.MIN_STANUM_FI
                num = num + 1;
                subcon_all(num) = subcontmp4;
                new(num) = newnum;
                break;
            end

        end

        if subcontmp4.scan(1).nsta <=2
            break;
        end
        subcontmp6 = subcontmp4;
        subcontmp6.scan.startmjd = max([subcontmp6.scan.sta.startmjd]);
        subcontmp6.scan(1).nsta = subcontmp6.scan(1).nsta-1;
        subcontmp6.scan(1).sta(badsta)=[];

        [subcontmp7] = sduration(source, station, twin, obsmode, subcontmp6, PARA);
        if (subcontmp7.nscan == 0)
            break;
        end
        subcontmp4 = subcontmp7;
        % check for source structure study
        if (PARA.FORSI == 1)
            for iscan = 1 : subcontmp4.nscan
                [obs_ok(iscan), grout(iscan)]= calc_soustruc_pscan(subcontmp4.scan(iscan),source,station);
                subcontmp4.scan(iscan).grout = grout;
            end
        end
    end
end
if (num == 0)
    subcon_f2.nscan = 0;
    return;
end

[yy] = ssubsort(source, station, staobs, srcobs, subcon_all, PARA);

for i = 0:length(subcon_all)-1
    thisSubcon =  subcon_all(yy(end-i));
    [flag,~,thisSubcon] = singleCheckSource(source, station, thisSubcon, PARA, jpl);
    if flag == 1
        % ### Pick the best sub-con ###
        subcon_f2 = thisSubcon;
        break
    end
end
