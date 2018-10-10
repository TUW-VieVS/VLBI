% Displays statistics for STAR mode
% 2016-07-11 M. Schartner: created

function [  ] = disp_stat_starmode( source,sched )

starSourceID = find([source.star]);
for i = 1:length(starSourceID)
    starObs(i) = struct('nscan',0,'sched',[],'duration',[],'time',[]);
end
allScan = [sched.scan];
scanSourceID = [allScan.srcid];


for isched = 1:length(sched)
    thisSched = sched(isched);
    for inscan = 1:thisSched.nscan
        thisScan = thisSched.scan(inscan);
        
        for iSTAR = 1:length(starSourceID)
            if starSourceID(iSTAR)==thisScan.srcid
                starObs(iSTAR).nscan = starObs(iSTAR).nscan+1;
                starObs(iSTAR).sched(starObs(iSTAR).nscan) = isched;
                starObs(iSTAR).duration(starObs(iSTAR).nscan) = max([thisScan.sta.duration]);
                starObs(iSTAR).time(starObs(iSTAR).nscan) =thisScan.startmjd;
                starObs(iSTAR).nsta(starObs(iSTAR).nscan) =thisScan.nsta;
            else
                continue
            end
        end
    end
end

fprintf('\n*** Observed STAR Sources: ***\n')
for i = 1:length(starObs)
    if starObs(i).nscan>0
        fprintf('Source %s: %2d scan(s):\n',source(starSourceID(i)).name,starObs(i).nscan);
        for j = 1:starObs(i).nscan
        	[year, month, day, hour, minute, second] = tymdhms(starObs(i).time(j));
            [datestr] = tdatestr(year, month, day, hour, minute, second);    
            fprintf('    sched: %3d | starttime: %s | duration:%4d [sec] | stations: %d\n',starObs(i).sched(j),datestr,starObs(i).duration(j),starObs(i).nsta(j))
        end
        
    end
end
fprintf('\n')

end

