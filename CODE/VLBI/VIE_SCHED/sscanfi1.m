% Purpose  
%   Calculate scan with fill-in mode.
% History  
%   2016-07-01 M. Schartner: created
%   2016-07-04 M. Schartner: added singleCheckSource
%	2016-07-11 M. Schartner: improved fill-in mode
%   2016-09-06 M. Schartner: changes for sched_manually
%   2016-09-07 M. Schartner: changes for singleCheckSource
%   2016-10-19 M. Schartner: PARA.MIN_STANUM_FI
%   2016-10-27 M. Schartner: bugfix, removeing of addsec, rewrote big parts of this function 
%   2016-11-03 M. Schartner: singleCheckSource now saves az el ha dc in subcon
%   2017-02-24 M. Schartner: new varargin input parameter for optimizer_quick_fillin 

function [subcon_f1] = sscanfi1(station, twin, staobs, source, srcobs, obsmode, subcon, PARA, jpl, varargin)

% initialization
stanum = length(station);
srcnum = length(source);
subcon_f1.nscan = 0;
minScanTime = PARA.MIN_SCAN;
calib = PARA.SOURCE + PARA.TAPETM + PARA.IDLE + PARA.CALIBRATION;
waitsec = zeros(stanum,1);
endTime = zeros(stanum,1);
raEnd = zeros(stanum,1);
deEnd = zeros(stanum,1);
azEnd = zeros(stanum,1);
elEnd = zeros(stanum,1);
haEnd = zeros(stanum,1);
dcEnd = zeros(stanum,1);
% calc waittime and starttime
boolParticipate = false(stanum,1);


if isempty(varargin)
    minEndTime = zeros(subcon.nscan,1);
    for iscan=1 : subcon.nscan
        minEndTime(iscan) = subcon.scan(iscan).startmjd + min([subcon.scan(iscan).sta.duration])/86400;
        for ista = 1:subcon.scan(iscan).nsta
            idx = subcon.scan(iscan).sta(ista).staid;
            endTime(idx) = subcon.scan(iscan).startmjd;
            srcid = subcon.scan(iscan).srcid;
            raEnd(idx) = source(srcid).ra;
            deEnd(idx) = source(srcid).de;
            azEnd(idx) = subcon.scan(iscan).sta(ista).az;
            elEnd(idx) = subcon.scan(iscan).sta(ista).el;
            haEnd(idx) = subcon.scan(iscan).sta(ista).ha;
            dcEnd(idx) = subcon.scan(iscan).sta(ista).dc;
            waitsec(idx) = floor((subcon.scan(iscan).startmjd-staobs(ista).endmjd)*86400);
            boolParticipate(idx) = true;
        end
    end

    if ~all(boolParticipate)
        endTime(~boolParticipate)=min(minEndTime);
        waitsec(~boolParticipate)=floor((endTime(~boolParticipate)-[staobs(~boolParticipate).endmjd]')*86400);
    end
    
else
    azEnd = varargin{1}.az;
    elEnd = varargin{1}.el;
    haEnd = varargin{1}.ha;
    dcEnd = varargin{1}.dc;
    endTime = varargin{1}.endTime;
    waitsec = varargin{1}.waitsec;
end
boolSta = waitsec>minScanTime+2*calib+ 30; %assumption: min 30sec slewtime

nSubcons = 0;
% allSubcons = struct('nscan',[],'scan',[]);
if sum(boolSta)>2
    llh = [station.llh];
    lon = llh(1:3:end);
    lat = llh(2:3:end);
    for isrc = 1 : srcnum
        
        % check if the source is repeated in a short interval
        if (min(endTime)-srcobs(isrc).obsmjd)*1440 < PARA.MIN_SRCRP
            continue;
        end
        if ~isempty(subcon)
            if any(isrc == [subcon.scan.srcid] )
                continue;
            end
        end
        if source(isrc).star
            continue;
        end
        downstanum = 0;
        downstasn(1:stanum) = 0;  
        % check the visibility at stations
        ra = source(isrc).ra;
        de = source(isrc).de;
        endmjd = [staobs.endmjd];
        [az, el, ha, dc] = zazel_s(endmjd, lon, lat, ra, de);
        for staid = 1 : stanum
            if boolSta(staid)
                [lup] = zlup(endmjd(staid), az(staid), el(staid), ha(staid), dc, staid, station, PARA.MIN_CUTEL);
                if lup
                    downstanum = downstanum + 1;
                    downstasn(downstanum) = staid;
                end
            end
        end    
        
        % is source up for more than 2 stations 
        if (downstanum < PARA.MIN_STANUM_FI)
            continue;
        end
        
        % create Subconfiguration with 1 scan
        subcontmp = [];
        subcontmp.nscan = 1;
        subcontmp.scan(1).srcid    = isrc;
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
            continue;
        end     
        
        scanStartTime = zeros(stanum,1);
        timeLeftAfterScan = zeros(stanum,1);
        for ista = 1:subcontmp1.scan.nsta
            scanStartTime(subcontmp1.scan.sta(ista).staid) = subcontmp1.scan.startmjd;
            timeLeftAfterScan(subcontmp1.scan.sta(ista).staid) = (endTime(subcontmp1.scan.sta(ista).staid)-scanStartTime(subcontmp1.scan.sta(ista).staid))*86400;
        end
        

        % is enough time to scan and at least 10sec slew?
        bool = timeLeftAfterScan > (minScanTime + 10) ;
        if sum(bool) < PARA.MIN_STANUM_FI
            continue;
        end
        
        % if a scan got removed:
        if subcontmp1.scan.nsta ~= sum(bool)
            del = 0;
            for i=1:subcontmp1.scan.nsta
                if ~bool(subcontmp1.scan.sta(i-del).staid)
                    subcontmp1.scan.sta(i-del) = [];
                    del = del+1;
                end
            end
            subcontmp1.scan.nsta = sum(bool);
            [subcontmp1] = sstartmjd(source, station, twin, staobs, subcontmp1, PARA);
        end
        
        % calc scan duration
        [subcontmp2] = sduration(source, station, twin, obsmode, subcontmp1, PARA);
        if (subcontmp2.nscan  == 0)
            continue;
        end
        
        scanEndTime = zeros(stanum,1);
        timeLeftAfterScan = zeros(stanum,1);
        for ista = 1:subcontmp2.scan.nsta
            scanEndTime(subcontmp2.scan.sta(ista).staid) = subcontmp2.scan.startmjd+subcontmp2.scan.sta(ista).duration/86400;
            timeLeftAfterScan(subcontmp2.scan.sta(ista).staid) = ...
                (endTime(subcontmp2.scan.sta(ista).staid)-...
                scanEndTime(subcontmp2.scan.sta(ista).staid))*86400;
        end

        bool = timeLeftAfterScan > 10 ;
        if sum(bool) < PARA.MIN_STANUM_FI
            continue;
        end
        % if a scan got removed:
        if subcontmp2.scan.nsta ~= sum(bool)
            subcontmp2.scan.nsta = sum(bool);
            del = 0;
            for i=1:subcontmp2.scan.nsta
                if ~bool(subcontmp2.scan.sta(i-del).staid)
                    subcontmp2.scan.sta(i-del) = [];
                    del = del+1;
                end
            end
            subcontmp2 = sstartmjd(source, station, twin, staobs, subcontmp2, PARA);
            if (subcontmp2.nscan  == 0)
                continue;
            end
            subcontmp2 = sduration(source, station, twin, obsmode, subcontmp2, PARA);
            if (subcontmp2.nscan  == 0)
                continue;
            end
        end
                
        % + slew to new source
        staobs_temp = supdate(station, source, subcontmp2, staobs, srcobs, 0, PARA, 'noOutput');
        
        % calc slewtime to other source
%         staid = [subcontmp2.scan.sta.staid]; 
        slewtime2 = zeros(stanum,1);
        for ista = 1:subcontmp2.scan.nsta
            staid = subcontmp2.scan.sta(ista).staid;
            slewtime2(staid) = sslew(station, staobs_temp, azEnd(staid), elEnd(staid), haEnd(staid), dcEnd(staid), staid, PARA, 'directSlew');
        end
        
        bool = timeLeftAfterScan > slewtime2 + calib;
        for ista = 1:stanum
            if ~any(ista == [subcontmp2.scan.sta.staid]) && bool(ista) == 1
                bool(ista)=0;
            end
        end
        
        if sum(bool) < PARA.MIN_STANUM_FI
            continue;
        end
        
        if ~isequal(find(bool)',[subcontmp2.scan.sta.staid])
            del = 0;
            for i=1:subcontmp2.scan.nsta
                if ~bool(subcontmp2.scan.sta(i-del).staid)
                    subcontmp2.scan.sta(i-del) = [];
                    del = del+1;
                end
            end
            subcontmp2.scan.nsta = sum(bool);
            
            subcontmp2 = sstartmjd(source, station, twin, staobs, subcontmp2, PARA);
            
            if (subcontmp2.nscan  == 0)
                continue;
            end

            subcontmp2 = sduration(source, station, twin, obsmode, subcontmp2, PARA);
            if (subcontmp2.nscan  == 0)
                continue;
            end

        end
        
        % source structure test
        iforsi = 1;
        if (PARA.FORSI == 1)
            for iscan = 1 : subcontmp2.nscan
                [obs_ok(iscan), grout(iscan)]= calc_soustruc_pscan(subcontmp2.scan(iscan),source,station);
                subcontmp2.scan(iscan).grout = grout;
            end
            iforsi = sum(obs_ok(1:subcontmp2.nscan));
        end
        if (iforsi == 0)
            continue;
        end
                
        % final adding
        if subcontmp2.scan.nsta>=PARA.MIN_STANUM_FI
            nSubcons = nSubcons+1;
            subcon_all(nSubcons)=subcontmp2;
        end
    end
    
    % take the best scan
    if nSubcons>0
        [yy] = ssubsort(source, station, staobs, srcobs, subcon_all, PARA);
        
        for i = 0:length(subcon_all)-1
            thisSubcon =  subcon_all(yy(end-i));
            [flag,~,thisSubcon] = singleCheckSource(source, station, thisSubcon, PARA, jpl);
            if flag == 1
                % ### Pick the best sub-con ###
                subcon_f1 = thisSubcon;
                break
            end
        end        
    end
end
