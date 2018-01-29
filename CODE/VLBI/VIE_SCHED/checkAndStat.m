% Purpose  
%   checks the schedule station per station. 
%   also displays some statistics for easy comparison of schedules. 
% History  
%   2016.11.04 M. Schartner: created
%	2017.02.20 M. Schartner: now returns statistic structure

function [ stat ] = checkAndStat( sched,station,source,PARA,meanSky )

fprintf(PARA.fid_footer,'========================= start vie_check =========================\n');

allScans = [sched.scan];

allScanNr = [sched.nscan];

if sum(allScanNr) ~= length(allScans)
    error('Number of scans doesn''t match with scan number!')
end

%% get single observations
c = 1;
for i = 1:length(sched)
    for j = 1:sched(i).nscan
        allScans(c).schedNr = i;
        c = c+1;
    end
end

allStanum = [allScans.nsta];
allStart = [allScans.startmjd];
allScanNr = [allScans.schedNr];
allSource = [allScans.srcid];
obs = [allScans.sta];
c = 1;
for i = 1:length(allStanum)
    for j = 1:allStanum(i)
        obs(c).startscan=allStart(i);
        obs(c).srcid = allSource(i);
        obs(c).schedNr = allScanNr(i);
        c = c+1;
    end
end

nsta = max([obs.staid]);

%% check times station by station
errorNr = 0;
const = PARA.SOURCE + PARA.TAPETM + PARA.IDLE + PARA.CALIBRATION;
totalSlewTime = zeros(1,nsta);
totalScanTime = zeros(1,nsta);
totalConstTime = zeros(1,nsta);
totalNumberScans = zeros(1,nsta);

for ista = 1:nsta
    fprintf(PARA.fid_footer,'  Checking station: %s\n',station(ista).name);
    istaObs = obs([obs.staid] == ista);
    totalNumberScans(ista)=length(istaObs);
    staobs = struct();
    staobs(ista).endmjd   = PARA.startmjd;
    staobs(ista).az       = station(ista).lim11;
    staobs(ista).el       = station(ista).lim12;
    staobs(ista).ha       = station(ista).lim11;
    staobs(ista).dc       = station(ista).lim12;
    lon = station(ista).llh(1);
    lat = station(ista).llh(2);
    for iobs = 1:length(istaObs)
        thisObs = istaObs(iobs);
        
        mjd = thisObs.startscan;
        ra = source(thisObs.srcid).ra;
        de = source(thisObs.srcid).de;

%         [az, el, ha, dc] = zazel_s(mjd, lon, lat, ra, de);            
        az = thisObs.az;
        el = thisObs.el;
        ha = thisObs.ha;
        dc = thisObs.dc;
        
        %% check if source is up
        lup = zlup(mjd, mod(az,2*pi), el, ha, dc, ista, station, PARA.MIN_CUTEL);
        
        if ~lup
            errorNr = errorNr +1;
            fprintf(PARA.fid_footer,'ERROR NR %d:\n',errorNr);
            [year, month, day, hour, minute, second] = tymdhms(mjd);
            [datestr] = tdatestr(year, month, day, hour, minute, second);    
            fprintf(PARA.fid_footer,'Source %s (id: %d) is not visible for station %s at %s (sched nr: %d)\n',source(thisObs.srcid).name,thisObs.srcid,station(ista).name,datestr,thisObs.schedNr);
            fprintf(PARA.fid_footer,'----------------------------------------------------------------------------\n');
        end
        
        %% check scantime
        totalScanTime(ista) = totalScanTime(ista)+thisObs.duration;
        if thisObs.duration < (thisObs.endmjd-thisObs.startscan)
            errorNr = errorNr +1;
            fprintf(PARA.fid_footer,'ERROR NR %d:\n',errorNr);
            [year, month, day, hour, minute, second] = tymdhms(mjd);
            [datestr] = tdatestr(year, month, day, hour, minute, second);    
            fprintf(PARA.fid_footer,'Not enough scan time for %s to scan %s (id: %d) at %s (sched nr: %d)\n',station(ista).name,source(thisObs.srcid).name,thisObs.srcid,datestr,thisObs.schedNr);
            fprintf(PARA.fid_footer,'    Scan duration:  %.1f\n',thisObs.duration);
            fprintf(PARA.fid_footer,'    Available time: %.1f\n',(thisObs.endmjd-thisObs.startscan)*86400);
            fprintf(PARA.fid_footer,'----------------------------------------------------------------------------\n');
        end
        
        %% check scanstart and stationstart
        if thisObs.startscan<thisObs.startmjd
            errorNr = errorNr +1;
            fprintf(PARA.fid_footer,'ERROR NR %d:\n',errorNr);
            fprintf(PARA.fid_footer,'Station starttime is later than scan starttime (VieVS bug) (sched nr: %d)\n',thisObs.schedNr);
            [year, month, day, hour, minute, second] = tymdhms(thisObs.startmjd);
            [datestr] = tdatestr(year, month, day, hour, minute, second);    
            fprintf(PARA.fid_footer,'    Station start time:  %s \n',datestr);

            [year, month, day, hour, minute, second] = tymdhms(thisObs.startscan);
            [datestr] = tdatestr(year, month, day, hour, minute, second);    
            fprintf(PARA.fid_footer,'    Scan start time:     %s \n',datestr);
            fprintf(PARA.fid_footer,'----------------------------------------------------------------------------\n');
        end
        
        %% check slewtime + constants 
        [st,unaz] = sslew(station, staobs, az, el, ha, dc, ista, PARA, 'directSlew');

        if iobs>1
            totalSlewTime(ista) = totalSlewTime(ista)+st;
            totalConstTime(ista) = totalConstTime(ista)+const;            
            if thisObs.startscan < (staobs(ista).endmjd+(st+const-1)/86400)    % 1 second tolerance 
                errorNr = errorNr +1;
                fprintf(PARA.fid_footer,'ERROR NR %d:\n',errorNr);
                fprintf(PARA.fid_footer,'Not enough slew time for %s to scan %s (between sched nr: %d-%d)\n',station(ista).name,source(thisObs.srcid).name,istaObs(iobs-1).schedNr,thisObs.schedNr);

                [year, month, day, hour, minute, second] = tymdhms(staobs(ista).endmjd);
                [datestr] = tdatestr(year, month, day, hour, minute, second);    
                fprintf(PARA.fid_footer,'    End of last scan:  %s (source: %s id: %d)\n',datestr,source(istaObs(iobs-1).srcid).name,istaObs(iobs-1).srcid);

                [year, month, day, hour, minute, second] = tymdhms(thisObs.startscan);
                [datestr] = tdatestr(year, month, day, hour, minute, second);    
                fprintf(PARA.fid_footer,'    Start of new scan: %s (source: %s id: %d)\n',datestr,source(thisObs.srcid).name,thisObs.srcid);

                fprintf(PARA.fid_footer,'    Slewtime:       %6.1f [sec] \n',st);
                fprintf(PARA.fid_footer,'    Constants:      %6.1f [sec] \n',const);
                fprintf(PARA.fid_footer,'    Available time: %6.1f [sec] \n',(thisObs.startscan-staobs(ista).endmjd)*86400);
                fprintf(PARA.fid_footer,'    late on source: %6.1f [sec] \n',(st+const)-(thisObs.startscan-staobs(ista).endmjd)*86400);
                fprintf(PARA.fid_footer,'----------------------------------------------------------------------------\n');
            end
        end
        %% update for next scan
        [azEndmjd, el, ha, dc] = zazel_s(thisObs.endmjd, lon, lat, ra, de);  
        staobs(ista).endmjd = thisObs.startscan+thisObs.duration/86400;
        diffAz = az-azEndmjd;
        diffAzRound = round(diffAz/(2*pi));
        azEndmjd = azEndmjd+diffAzRound*(2*pi);
%         fprintf(PARA.fid_footer,'unaz: %.4f azEndmjd %.4f\n',unaz,azEndmjd)
        staobs(ista).az = azEndmjd;
        staobs(ista).el = el;
        staobs(ista).ha = ha;
        staobs(ista).dc = dc;
    end
end

if errorNr == 0
    fprintf(PARA.fid_footer,'    No problems found!\n');
end

%%
totalTime = (max([thisObs.endmjd])-PARA.startmjd)*86400;
stat = struct();
stat.pSlewTime = totalSlewTime/totalTime*100;
stat.pScanTime = totalScanTime/totalTime*100;
stat.pConstTime = totalConstTime/totalTime*100;
stat.pIdleTime = (totalTime-totalSlewTime-totalScanTime-totalConstTime)/totalTime*100;

fprintf(PARA.fid_footer,'.------------------.');
for ista = 1:nsta+2
    fprintf(PARA.fid_footer,'----------.');
end
fprintf(PARA.fid_footer,'\n|    statistics    |');
fprintf(PARA.fid_footer,'   mean   |');
fprintf(PARA.fid_footer,'    std   |');
for ista = 1:nsta
    fprintf(PARA.fid_footer,' %8s |',station(ista).name);
end
fprintf(PARA.fid_footer,'\n|------------------|');
for ista = 1:nsta+2
    fprintf(PARA.fid_footer,'----------|');
end

fprintf(PARA.fid_footer,'\n| Scan Time        |');
fprintf(PARA.fid_footer,' %6.2f %% |',mean(stat.pScanTime));
fprintf(PARA.fid_footer,' %6.2f %% |',std(stat.pScanTime));
fprintf(PARA.fid_footer,' %6.2f %% |',stat.pScanTime);

fprintf(PARA.fid_footer,'\n| Slew Time        |');
fprintf(PARA.fid_footer,' %6.2f %% |',mean(stat.pSlewTime));
fprintf(PARA.fid_footer,' %6.2f %% |',std(stat.pSlewTime));
fprintf(PARA.fid_footer,' %6.2f %% |',stat.pSlewTime);

fprintf(PARA.fid_footer,'\n| Idle Time        |');
fprintf(PARA.fid_footer,' %6.2f %% |',mean(stat.pIdleTime));
fprintf(PARA.fid_footer,' %6.2f %% |',std(stat.pIdleTime));
fprintf(PARA.fid_footer,' %6.2f %% |',stat.pIdleTime);

fprintf(PARA.fid_footer,'\n| Calibration Time |');
fprintf(PARA.fid_footer,' %6.2f %% |',mean(stat.pConstTime));
fprintf(PARA.fid_footer,' %6.2f %% |',std(stat.pConstTime));
fprintf(PARA.fid_footer,' %6.2f %% |',stat.pConstTime);

fprintf(PARA.fid_footer,'\n|------------------|');
for ista = 1:nsta+2
    fprintf(PARA.fid_footer,'----------|');
end
stat.totalNumberScans = totalNumberScans;
fprintf(PARA.fid_footer,'\n| Number of Scans  |');
fprintf(PARA.fid_footer,' %6.1f   |',mean(totalNumberScans));
fprintf(PARA.fid_footer,' %6.1f   |',std(totalNumberScans));
fprintf(PARA.fid_footer,' %4d     |',totalNumberScans);

stat.scansPerH = totalNumberScans/(totalTime/3600);
fprintf(PARA.fid_footer,'\n| Scans per hour   |');
fprintf(PARA.fid_footer,' %6.1f   |',mean(totalNumberScans)/(totalTime/3600));
fprintf(PARA.fid_footer,' %6.1f   |',std(totalNumberScans)/(totalTime/3600));
fprintf(PARA.fid_footer,' %6.1f   |',totalNumberScans/(totalTime/3600));
fprintf(PARA.fid_footer,'\n|------------------|');

for ista = 1:nsta+2
    fprintf(PARA.fid_footer,'----------|');
end

stat.meanSky = meanSky;
fprintf(PARA.fid_footer,'\n| Sky Coverage     |');
fprintf(PARA.fid_footer,' %6.2f   |',mean(meanSky));
fprintf(PARA.fid_footer,' %6.2f   |',std(meanSky));
fprintf(PARA.fid_footer,' %6.2f   |',meanSky);
fprintf(PARA.fid_footer,'\n''------------------''');
for ista = 1:nsta+2
    fprintf(PARA.fid_footer,'----------''');
end


stat.nobs = sum(([allScans.nsta].*([allScans.nsta]-1))/2);
stat.nscan = length(allScans);
fprintf(PARA.fid_footer,'\n\ntotal number of scans:        %7d\n',stat.nscan);
fprintf(PARA.fid_footer,'total number of observations: %7d\n\n',stat.nobs);

for i = 1:max([sched.nscan])
    fprintf(PARA.fid_footer,'Number of %d-scan subnets: %4d\n',i,sum([sched.nscan]==i));
end
fprintf(PARA.fid_footer,'\n');

for i = 2:max([allScans.nsta])
    fprintf(PARA.fid_footer,'Number of %d-station scans: %4d\n',i,sum([allScans.nsta]==i));
end
fprintf(PARA.fid_footer,'\n');

% fprintf(PARA.fid_footer,'\n\n=========================== statistics ============================\n');





