% This function removes all scans that observe one of the badSources and
% tries to squeeze in another scan. 
%
% CREATED: 21.02.17 - Matthias Schartner
%
% CHANGELOG: 


function [ sched ] = optimizer_quick_fillin( sched, station, twin, source, obsmode, PARA, goodSourceID, badSourceID )

fprintf(PARA.fid_body,'\nStart to check forfilled optimization condition:\n');
fprintf(PARA.fid_body,'==============================================\n');
goodSources = source(goodSourceID);
jpl = load_jpl_Earth('jpl_421');

newSched = struct('nscan',[],'scan',[]);
isched = 1;
% for isched=1:length(sched)
while 1
    if isched>length(sched)
        break
    end
    flag_bad_scan = 0;
    
    for iscan = 1:length(sched(isched).nscan)
        if any(sched(isched).scan(iscan).srcid == badSourceID)
            flag_bad_scan = 1;
        end
    end
    
    if flag_bad_scan
        fprintf(PARA.fid_body,'---------- removed sched number %4d --------- \n',isched);
        %% get sourceID of all previouse scans inside PARA.MINSRCRP
        thisMJD = sched(isched).scan(iscan).startmjd;

        sourceIds_around = [];
        counter = 1;
        while 1
            if isched-counter<1
                break
            end
            sched_before = sched(isched-counter);
            for iscan_before = 1:length(sched_before.nscan)

                mjd_before = sched_before.scan(iscan_before).startmjd;
                if ((thisMJD-mjd_before)*24*60)<PARA.MIN_SRCRP
                    idx = find(sched_before.scan(iscan_before).srcid == goodSourceID);
                    if ~isempty(idx)
                        sourceIds_around(end+1) = idx;
                    end
                end
            end

            if (thisMJD-mjd_before)*24*60>PARA.MIN_SRCRP
                break
            end
            counter = counter+1;
        end


        counter = 1;
        while 1
            if isched+counter>length(sched)
                break
            end
            sched_after = sched(isched+counter);
            for iscan_after = 1:length(sched_after.nscan)

                mjd_after = sched_after.scan(iscan_after).startmjd;
                if ((mjd_after-thisMJD)*24*60)<PARA.MIN_SRCRP
                    idx = find(sched_after.scan(iscan_after).srcid == goodSourceID);
                    if ~isempty(idx)
                        sourceIds_around(end+1) = idx;
                    end
                end
            end

            if (mjd_after-thisMJD)*24*60>PARA.MIN_SRCRP
                break
            end
            counter = counter+1;
       end


        %% create dummy srcobs
        srcobs = struct('obsmjd',0,'nscan',0);
        for i = 1:length(goodSourceID)
            if any(i == sourceIds_around)
                srcobs(i).obsmjd = thisMJD-(PARA.MIN_SRCRP)/1440/2;
                srcobs(i).nscan = 1;
            else
                srcobs(i).obsmjd = 0;
                srcobs(i).nscan = 0;
            end
        end

        %% create dummy staobs
        staobs = struct('ifirst',[],'endmjd',[],'az',[],'el',[],'ha',[],'dc',[],'ncat',[],'cattn',[],'catsn',[]);

        counter = 1;

        found_sta = false(size(station));
        mjd_before = zeros(size(station))+Inf;
        while 1
            if isched-counter<1
                break
            end

            sched_before = sched(isched-counter);
            for iscan_before = 1:sched_before.nscan
                scan_before = sched_before.scan(iscan_before);
                for ista = 1:scan_before.nsta
                    this_staid = scan_before.sta(ista).staid;
                    if found_sta(this_staid)==0

                        found_sta(this_staid)=1;

                        staobs(this_staid).ifirst = 0;

                        mjd = scan_before.sta(ista).endmjd;
                        mjd_before(ista)=mjd;
                        
                        ra = source(scan_before.srcid).ra;
                        de = source(scan_before.srcid).de;
                        lon = station(this_staid).llh(1);
                        lat = station(this_staid).llh(2);

                        [azEndmjd, el, ha, dc] = zazel_s(mjd, lon, lat, ra, de);
                        azStartmjd = scan_before.sta(ista).az;
                        diffAz = azStartmjd-azEndmjd;
                        diffAzRound = round(diffAz/(2*pi));
                        az = azEndmjd+diffAzRound*(2*pi);

                        staobs(this_staid).endmjd = mjd;
                        staobs(this_staid).az = az;
                        staobs(this_staid).el = el;
                        staobs(this_staid).ha = ha;
                        staobs(this_staid).dc = dc;
                        staobs(ista).ncat = 0;
                        staobs(ista).cattn = [];
                        staobs(ista).catsn = [];
                    end
                end
            end
            if all(found_sta)
                break
            end
            counter = counter+1;
        end
        
        
        
        % check if all stations have observed before
        if ~all(found_sta)
            for i = 1:length(staobs)
                if(found_sta(i)==0)
                    staobs(ista).ifirst = 1;
                    staobs(ista).endmjd = 0;
                    staobs(ista).az = 0;
                    staobs(ista).el = 0;
                    staobs(ista).ha = 0;
                    staobs(ista).dc = 0;
                    staobs(ista).ncat = 0;
                    staobs(ista).cattn = [];
                    staobs(ista).catsn = [];
                end
            end
        end
        
        %% get the next required az el ha dz and mjd for each station 
        counter = 1;
        
        found_sta = false(size(station));
        az_end = zeros(size(station));
        el_end = zeros(size(station));
        ha_end = zeros(size(station));
        dc_end = zeros(size(station));
        endTime = zeros(size(station));
        waitsec = zeros(size(station));
        
        while 1
            if isched+counter>length(sched)
                break
            end

            sched_after = sched(isched+counter);
            for iscan_after = 1:sched_after.nscan
                scan_after = sched_after.scan(iscan_after);
                for ista = 1:scan_after.nsta
                    this_staid = scan_after.sta(ista).staid;
                    if found_sta(this_staid)==0
                        
                        found_sta(this_staid)=1;
                        
                        az_end(this_staid) = scan_after.sta(ista).az;
                        el_end(this_staid) = scan_after.sta(ista).el;
                        ha_end(this_staid) = scan_after.sta(ista).ha;
                        dc_end(this_staid) = scan_after.sta(ista).dc;
                        endTime(this_staid) = scan_after.startmjd;
                        waitsec(this_staid) = floor((endTime(this_staid)-staobs(this_staid).endmjd)*86400);
                    end
                end
            end
            if all(found_sta)
                break
            end
            counter = counter+1;
        end

        endOBS.az = az_end;
        endOBS.el = el_end;
        endOBS.ha = ha_end;
        endOBS.dc = dc_end;
        endOBS.endTime = endTime;
        endOBS.waitsec = waitsec;
        
        %% remove schedule in sched
        sched(isched)=[];
        isched = isched-1;
        %% squeeze in new scan
        counter = 0;
        while true
            [thisSubcon] = sscanfi1(station, twin, staobs, goodSources, srcobs, obsmode, [], PARA, jpl, endOBS);
            if (thisSubcon.nscan > 0)
                if counter == 0
                    fprintf(PARA.fid_body,'    new scan:\n');
                end
                [staobs, srcobs] = supdate(station, goodSources, thisSubcon, staobs, srcobs, [], PARA);
                thisSubcon.scan.srcid = goodSourceID(thisSubcon.scan(1).srcid);
                sched = [sched(1:isched) thisSubcon sched(isched+1:end)];
                isched = isched+1;
                newSched(end+1) = thisSubcon;
                counter = counter+1;
                
            elseif (thisSubcon.nscan == 0) && counter == 0
                fprintf(PARA.fid_body,'    no other possible scan found\n');
                break;
            else
                break;
            end
        end
         
        fprintf(PARA.fid_body,'==============================================\n');
    else
        newSched(end+1) = sched(isched);
    end
    isched = isched+1;
end
newSched = newSched(2:end);
end

