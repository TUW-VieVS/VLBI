function [  ] = checkSkd( fname )
clc, clear, close all;
format compact
file = '/home/mschartn/shares/home/Axel_scheduling/R1912/VieSchedpp/r1912_v086.skd';

fid = fopen(file);

obsmode.samprate = 16;
obsmode.bits = 2;
obsmode.chanumband = [10 6];

while ~feof(fid)
    fline = fgetl(fid);

    if strcmp(strtrim(fline),'$PARAM')
        fprintf('read param ... ');
        PARA = readParam(fid);
        fprintf('done\n');
        break;
   end
end

frewind(fid);
while ~feof(fid)
    fline = fgetl(fid);

    if strcmp(strtrim(fline),'$SOURCES')
        fprintf('read sources ... ');
        source = readSource(fid);
        fprintf('done\n');
        break;
   end
end

frewind(fid);
while ~feof(fid)
    fline = fgetl(fid);

    if strcmp(strtrim(fline),'$FLUX')
        fprintf('read fluxInfo ... ');
        source = readFlux(fid,source);
        fprintf('done\n');
        break;
   end
end

frewind(fid);
while ~feof(fid)
    fline = fgetl(fid);
    if strcmp(strtrim(fline),'$STATIONS')
        fprintf('read stations ... ');
        station = readStation(fid);
        fprintf('done\n');
        break;
    end
end

frewind(fid);
while ~feof(fid)
    fline = fgetl(fid);
    if strcmp(strtrim(fline),'$SKED')
        fprintf('read scans, checking cable wrap and source visibility ...\n');
        sched = readSkd(fid,source,station,PARA);
        fprintf('done\n');
        break;
    end
end


fprintf('checking slew times and displaying statistics ...\n')
[meanSky,h] = checkSkyCoverage( sched,PARA,station );
stat = checkAndStat( sched,station,source,PARA,meanSky );
fprintf('done\n')

twin = struct();
twin.num = 0;
% ongs(source, station, sched, NGS_name, twin);

fprintf('checking scan durations ...\n')
twin.num = 0;
checkScanDuration(sched,station,twin,source,obsmode,PARA)
fprintf('done\n')

end

function[PARA] = readParam(fid)
    PARA.MIN_CUTEL = deg2rad(5);
    PARA.MULTISCHED = 0;
    PARA.disp_sky_coverage = 0;
    PARA.fid_footer = 1;
    PARA.MAX_BANDNUM = 2;
    PARA.MAX_FLUXPARA = 30;
    PARA.MAX_SEFDPARA = 4;
    PARA.STARMODE = 0;
    PARA.WAVEL = [0.034896 3.8];
    
    while 1
        fline = fgetl(fid);
        if strcmp(fline(1),'$')
            break
        end
        if strcmp(fline(1),'*')
            continue;
        end
       split = strsplit(fline);
        
        bool = strcmp(split,'START');
        if any(bool)
            val = find(bool)+1;
            
            date = split{val};
            date = date(1:end);
            yy = str2double(date(1:2));
            d1 = datenum(sprintf('01-Jan-20%02d 00:00:00',yy));

            ddd = str2double(date(3:5));
            hh = str2double(date(6:7));
            mm = str2double(date(8:9));
            ss = str2double(date(10:11));
            d2 = d1+ddd-1+hh/24+mm/1440+ss/86400;
            startTime = datetime(d2,'ConvertFrom','datenum');
            dVector = datevec(startTime);

            PARA.startmjd = mjuliandate(dVector);

        end
        
        bool = strcmp(split,'END');
        if any(bool)
            val = find(bool)+1;
            
            date = split{val};
            date = date(1:end);
            yy = str2double(date(1:2));
            d1 = datenum(sprintf('01-Jan-20%02d 00:00:00',yy));

            ddd = str2double(date(3:5));
            hh = str2double(date(6:7));
            mm = str2double(date(8:9));
            ss = str2double(date(10:11));
            d2 = d1+ddd-1+hh/24+mm/1440+ss/86400;
            startTime = datetime(d2,'ConvertFrom','datenum');
            dVector = datevec(startTime);

            PARA.endmjd = mjuliandate(dVector);
        end

        bool = strcmp(split,'SOURCE');
        if any(bool)
            val = find(bool)+1;
            PARA.SOURCE = str2double(split{val});
        end
        
        bool = strcmp(split,'TAPETM');
        if any(bool)
            val = find(bool)+1;
            PARA.TAPETM = str2double(split{val});
        end
        
        bool = strcmp(split,'IDLE');
        if any(bool)
            val = find(bool)+1;
            PARA.IDLE = str2double(split{val});
        end
        
        bool = strcmp(split,'CALIBRATION');
        if any(bool)
            val = find(bool)+1;
            PARA.CALIBRATION = str2double(split{val});
        end

        bool = strcmp(split,'CORSYNCH');
        if any(bool)
            val = find(bool)+1;
            PARA.CORSYNCH = str2double(split{val});
        end

        bool = strcmp(split,'MAXSCAN');
        if any(bool)
            val = find(bool)+1;
            PARA.MAX_SCAN = str2double(split{val});
        end
        
        bool = strcmp(split,'MINSCAN');
        if any(bool)
            val = find(bool)+1;
            PARA.MIN_SCAN = str2double(split{val});
        end
    end

end

function [station] = readStation(fid)

station = struct();
counter = 1;
while 1
    fline = fgetl(fid);
    if strcmp(fline(1),'$')
        break
    end
    if strcmp(fline(1),'*')
        continue;
    end
    if strcmp(fline(1),'A')
        C = textscan(fline,'%s %s %s %s %f %f %f %f %f %f %f %f %f %s %s %s %s');
        station(counter).name = C{3}{1};
        station(counter).id = C{2}{1};
        station(counter).po = upper(C{15});
        station(counter).eq = upper(C{16});
        station(counter).axis = C{4}{1};
        
        station(counter).rate1 = C{6}/60;
        station(counter).c1 = C{7};
        station(counter).acc1 = station(counter).rate1;
        station(counter).acc1 = station(counter).rate1;
        station(counter).rate2 = C{10}/60;
        station(counter).c2 = C{11};
        station(counter).acc2 = station(counter).rate2;
        station(counter).acc2 = station(counter).rate2;

        station(counter).lim11 = deg2rad(C{8});
        station(counter).lim12 = deg2rad(C{9});
        station(counter).lim21 = deg2rad(C{12});
        station(counter).lim22 = deg2rad(C{13});
        
        station(counter).hmasknum = 0;
        station(counter).hmask = 0;
        overlapping= (station(counter).lim12-station(counter).lim11)-2*pi;
        if(overlapping>0)
            station(counter).azc2 = station(counter).lim12;
            station(counter).azc1 = station(counter).lim12-overlapping;
            station(counter).azn2 = station(counter).lim12-overlapping;
            station(counter).azn1 = station(counter).lim11+overlapping;
            station(counter).azw2 = station(counter).lim11+overlapping;
            station(counter).azw1 = station(counter).lim11;
        else
            station(counter).azc2 = station(counter).lim12;
            station(counter).azc1 = station(counter).lim12;
            station(counter).azn2 = station(counter).lim12;
            station(counter).azn1 = station(counter).lim11;
            station(counter).azw2 = station(counter).lim11;
            station(counter).azw1 = station(counter).lim11;
        end
        station(counter).downum = 0;
        station(counter).hmasknum = 0;
        counter = counter+1;
    end


    if strcmp(fline(1),'P')
        C = textscan(fline,'%s %s %s %f %f %f %f %f %f %[^\n]');
        for i = 1:counter-1
            if strcmp(station(i).name,C{3})
                station(i).xyz = [C{4} C{5} C{6}];
                [station(i).llh(2),station(i).llh(1),station(i).llh(3)] = xyz2ell([C{4} C{5} C{6}]);
                break;
            end
        end
    end
    
    if strcmp(fline(1),'T')
        split = strsplit(fline);
        for i = 1:counter-1
            if strcmpi(station(i).eq,split{2})
                station(i).sefdpara = zeros(2,4);
                station(i).sefdpara(1,1) = str2double(split{7});
                station(i).sefdpara(2,1) = str2double(split{9});
                if(length(split)>=13)
                    if(strcmp(split{10},'X'))
                        station(i).sefdpara(1,2:4) = str2double(split(11:13));
                    end
                    if(strcmp(split{10},'S'))
                        station(i).sefdpara(2,2:4) = str2double(split(11:13));
                    end
                end
                if(length(split)>=17)
                    if(strcmp(split{14},'X'))
                        station(i).sefdpara(1,2:4) = str2double(split(15:17));
                    end
                    if(strcmp(split{14},'S'))
                        station(i).sefdpara(2,2:4) = str2double(split(15:17));
                    end
                end
                
                station(i).minsnr = [20    15];
                break;
            end
        end
    end
    
    if(strcmp(fline(1),'H'))
        data = strsplit(fline);
        id = upper(data(2));
        mask = str2double(data(3:end));
        mask(isnan(mask)) = [];
        mask = deg2rad(mask);
        station(strcmp(id,[station.po])).hmask = mask;
        station(strcmp(id,[station.po])).hmasknum = length(mask);
    end
end

if length({station.id})~=length(unique({station.id}))
    error('No unique id for stations')
end
end

function [source] = readSource(fid)
    source = struct();
    counter = 1;
    while 1
        fline = fgetl(fid);
        if strcmp(fline(1),'*')
            continue;
        end

        if strcmp(fline(1),'$')
            break
        end
        C = textscan(fline,'%s %s %f %f %f %f %f %f %[^\n]');
        source(counter).name = C{1}{1};
        if strcmp(C{2},'$')
            source(counter).commoname = '        ';
        else
            source(counter).commoname = C{2}{1};
            source(counter).name = C{2}{1};
        end
        source(counter).ra = (C{3} + C{4}/60 + C{5}/3600) * 15 * pi / 180;
        neg = C{6}<0;
        source(counter).de = deg2rad(abs(C{6}) + C{7}/60 + C{8}/3600);
        if(neg)
            source(counter).de = -1*source(counter).de;
        end
        source(counter).info = C{9}{1};
        source(counter).star = 0;
        source(counter).fluxpara = zeros(2,30);

        source(counter).ifp = 0; 
        source(counter).interested=0;
        source(counter).star = 0;
        
        counter = counter+1;
    end
end

function [source] = readFlux(fid,source)

    while 1
        fline = fgetl(fid);
        if strcmp(fline(1),'*')
            continue;
        end

        if all(fline == -1) || strcmp(fline(1),'$')
            break
        end
        split = strsplit(fline);
        fluxname = split{1};
        fluxband = split{2};
        fluxtype = split{3};
        fluxpara = str2double(split(4:end));
        fluxpara(isnan(fluxpara)) = [];
        
        idx = find(strcmp(fluxname,{source.name}));
        if(isempty(idx))
            idx = find(strcmp(fluxname,{source.name}));
        end
        
        
        if(fluxband == 'X')
            if(fluxtype == 'B')
                source(idx(1)).fluxpara(1,1:length(fluxpara)) = fluxpara;
                source(idx(1)).fluxpara(1,end) = 1;
            else
                if (source(idx(1)).fluxpara(1,1) < 1.0d-3)
                    source(idx(1)).fluxpara(1,1:6) = fluxpara;
                elseif (source(idx(1)).fluxpara(1,7) < 1.0d-3)
                    source(idx(1)).fluxpara(1,7:12) = fluxpara;
                elseif (source(idx(1)).fluxpara(1,13) < 1.0d-3)
                    source(idx(1)).fluxpara(1,13:18) = fluxpara;
                end 
                source(idx(1)).fluxpara(1,end) = 2;
            end
        else
            if(fluxtype == 'B')
                source(idx(1)).fluxpara(2,1:length(fluxpara)) = fluxpara;
                source(idx(1)).fluxpara(2,end) = 1;
            else
                if (source(idx(1)).fluxpara(2,1) < 1.0d-3)
                    source(idx(1)).fluxpara(2,1:6) = fluxpara;
                elseif (source(idx(1)).fluxpara(2,7) < 1.0d-3)
                    source(idx(1)).fluxpara(2,7:12) = fluxpara;
                elseif (source(idx(1)).fluxpara(2,13) < 1.0d-3)
                    source(idx(1)).fluxpara(2,13:18) = fluxpara;
                end 
                source(idx(1)).fluxpara(2,end) = 2;
            end
        end
        
        
        
    end

end

function [sched] = readSkd(fid,source,station,PARA)
    jpl = load_jpl_Earth('jpl_421');
    

    sched = struct();
    counter = 1;
    while 1
        fline = fgetl(fid);
        if strcmp(fline(1),'*')
            continue;
        end
        if strcmp(fline(1),'$')
            break
        end
        C = strsplit(strtrim(fline));
        srcname = C{1};
        date = C{5};
%                 maxdur = C{6};
        stat = C{10};
        ids = stat(1:2:end);
        cflags = stat(2:2:end);
        nsta = length(ids);
        duration = zeros(nsta,1);
        for i=1:nsta
            duration(nsta-i+1) = str2double(C{end-i+1});
        end

        sched(counter).nscan = 1;

        scan = struct();
        yy = str2double(date(1:2));
        d1 = datenum(sprintf('01-Jan-20%02d 00:00:00',yy));

        ddd = str2double(date(3:5));
        hh = str2double(date(6:7));
        mm = str2double(date(8:9));
        ss = str2double(date(10:11));
        d2 = d1+ddd-1+hh/24+mm/1440+ss/86400;
        startTime = datetime(d2,'ConvertFrom','datenum');
        dVector = datevec(startTime);

        scan.startmjd = mjuliandate(dVector);
        srcid = find(strcmp(srcname,{source.name}));
        if(isempty(srcid))
            srcid = find(strcmp(srcname,{source.commoname}));
        end
        scan.srcid = srcid;
        ra = source(scan.srcid).ra;

        de = source(scan.srcid).de;
        scan.nsta = nsta;

        sta = struct();
        for i = 1:scan.nsta
            id = ids(i);
            sta(i).staid = find(strcmp(id,{station.id}));
            sta(i).startmjd = scan.startmjd;
            thisStation = station(sta(i).staid);
            lat = thisStation.llh(2);
            lon = thisStation.llh(1);

            [az, el, ha, dc] = zazel_s(scan.startmjd, lon, lat, ra, de);
            sta(i).az = az;
            sta(i).el = el;
            sta(i).ha = ha;
            sta(i).dc = dc;
            sta(i).slewtime = 0;
            
            sta(i).endmjd = scan.startmjd + duration(i)/86400;
            sta(i).duration = ceil((sta(i).endmjd-sta(i).startmjd)*86400);
        end
        scan.sta = sta;

        subcon = struct();
        
        subcon.nscan =1;
        subcon.scan = scan;
        [flag,stid,subcon] = singleCheckSource(source, station, subcon, PARA, jpl);
        
        if~(flag == 1 || flag == 5)
            warning('scan %d source %s %s station %s error flag in singleCheckSource',counter,source(srcid).name,source(srcid).commoname, station(stid).name)
        end
        scan = subcon.scan;
        for i = 1:scan.nsta
            az = scan.sta(i).az;
            thisStation = station(scan.sta(i).staid);
            if ~strcmp(thisStation.axis,'AZEL')
                continue
            end
            if(cflags(i)=='C')
                while(az<thisStation.azc2)
                    az = az+2*pi;
                end
                while(az>thisStation.azc1)
                    az = az-2*pi;
                end
                az = az+2*pi;
                if(az<thisStation.azc1 || az>thisStation.azc2)
                    maz = mean([thisStation.azc1 thisStation.azc2]);
                    azu = az;
                    azl = az-2*pi;
                    if(abs(azu-maz) < abs(azl-maz))
                        az = azu;
                    else
                        az = azl;
                    end
                    azdiff = min([abs(az-thisStation.azn1) abs(az-thisStation.azn2)]);
                    warning('scan %d source %s %s station %s azimuth is in false cable wrap! off by: %f deg',counter,source(srcid).name,source(srcid).commoname, thisStation.name, rad2deg(azdiff))
                end
                
            elseif(cflags(i) == 'W')
                while(az<thisStation.azw2)
                    az = az+2*pi;
                end
                while(az>thisStation.azw1)
                    az = az-2*pi;
                end
                az = az+2*pi;
                if(az<thisStation.azw1 || az>thisStation.azw2)
                    maz = mean([thisStation.azw1 thisStation.azw2]);
                    azu = az;
                    azl = az-2*pi;
                    if(abs(azu-maz) < abs(azl-maz))
                        az = azu;
                    else
                        az = azl;
                    end
                    azdiff = min([abs(az-thisStation.azn1) abs(az-thisStation.azn2)]);
                    warning('scan %d source %s %s station %s azimuth is in false cable wrap! off by: %f deg',counter,source(srcid).name,source(srcid).commoname, thisStation.name, rad2deg(azdiff))
                end
            else
                while(az<thisStation.azn2)
                    az = az+2*pi;
                end
                while(az>thisStation.azn1)
                    az = az-2*pi;
                end
                az = az+2*pi;
                if(az<thisStation.azn1 || az>thisStation.azn2)
                    maz = mean([thisStation.azn1 thisStation.azn2]);
                    azu = az;
                    azl = az-2*pi;
                    if(abs(azu-maz) < abs(azl-maz))
                        az = azu;
                    else
                        az = azl;
                    end
                    azdiff = min([abs(az-thisStation.azn1) abs(az-thisStation.azn2)]);
                    warning('scan %d source %s %s station %s azimuth is in false cable wrap! off by: %f deg',counter,source(srcid).name,source(srcid).commoname, thisStation.name, rad2deg(azdiff))
                end
            end
            scan.sta(i).az = az;
        end
        
        sched(counter).scan = scan;
        counter = counter+1;
    end
end

function [] = checkScanDuration(sched,station,twin,source,obsmode,PARA)

for i=1:length(sched)
    
    subcon = struct();
    subcon.nscan = 1;
    subcon.scan = sched(i).scan;
    
    [subcon_s,~] = sduration(source, station, twin, obsmode, subcon, PARA);
    if(subcon_s.nscan == 0)
            warning('scan %d source %s %s no valid scan because of scan duration',...
            i,source(subcon.scan.srcid).name,source(subcon.scan.srcid).commoname);
            continue;
    end
    for j = 1:subcon.scan.nsta
        duration1 = subcon.scan.sta(j).duration;
        try
        duration2 = subcon_s.scan.sta(j).duration;
        catch
            fprintf('test');
        end
        if(subcon.scan.nsta ~= subcon_s.scan.nsta)
            warning('scan %d source %s %s at least one station has larger scan duration',...
            i,source(subcon.scan.srcid).name,source(subcon.scan.srcid).commoname);
            break;
        end
        
        if(duration1<duration2)
            warning('scan %d source %s %s station %s not enough scan duration!\nrequired: %d\nscheduled: %d',...
            i,source(subcon.scan.srcid).name,source(subcon.scan.srcid).commoname,station(subcon.scan.sta(j).staid).name,duration2, duration1);
            break
        end
        if(duration1>duration2+5 && duration1>21)
%             fprintf('scan %d source %s %s station %s too long!\nrequired: %d\nscheduled: %d\n',...
%             i,source(subcon.scan.srcid).name,source(subcon.scan.srcid).commoname,station(subcon.scan.sta(j).staid).name,duration2, duration1);
        end
    end
    
end

end