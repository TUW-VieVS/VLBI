% Purpose  
%   Write $SKED (scheduled observations).
% History  
%   2010-04-06   Jing SUN   Created

% CHANGES:
% - 2015-07-22, Lucia Plank: revision to use all required information from the GSFC catalog files 


function osked(source, station, obsmode, sched, fid_skd, PARA)

%
num = 0;
for isched = 1 : length(sched)
    for iscan = 1 : sched(isched).nscan
        num = num + 1;
        startmjd(num) = sched(isched).scan(iscan).startmjd;
        startmjdsn(num,1) = isched;
        startmjdsn(num,2) = iscan;
    end
end
[x, y] = sort(startmjd(1:num));
%
for i = 1 : num
    isched = startmjdsn(y(i),1);
    iscan = startmjdsn(y(i),2);
    srcid = sched(isched).scan(iscan).srcid;
    if (source(srcid).commoname(1) ~= ' ')
        srcname(1:8) = source(srcid).commoname(1:8);
    else
        srcname(1:8) = source(srcid).name(1:8);
    end
    startmjd = sched(isched).scan(iscan).startmjd;   %%%
    [year, month, day, hh, min, secf] = tymdhms(startmjd);
    [cdays] = tdays(year, month, day);
    dsec = round(secf);
    min(dsec==60)=min(dsec==60)+1;
    dsec(dsec==60)=0;
    hh(min==60)=hh(min==60)+1;
    min(min==60)=0;
    stadur = 0.0;
    for ista = 1 : sched(isched).scan(iscan).nsta
        if (sched(isched).scan(iscan).sta(ista).duration > stadur)
            stadur = sched(isched).scan(iscan).sta(ista).duration; 
        end   
    end  
    fprintf(fid_skd, '%s %3d %s PREOB  ', srcname,PARA.CALIBRATION,char(obsmode.sx));
    fprintf(fid_skd, '%s%s%02d%02d%02d', num2str(year-2000),cdays,hh,min,dsec);
    fprintf(fid_skd, '%8d MIDOB%8d POSTOB ', stadur, PARA.IDLE);
    for ista = 1 : sched(isched).scan(iscan).nsta
        fprintf(fid_skd, '%s', station(sched(isched).scan(iscan).sta(ista).staid).id);
        if (strcmp(station(sched(isched).scan(iscan).sta(ista).staid).axis, 'HADC') | strcmp(station(sched(isched).scan(iscan).sta(ista).staid).axis, 'XYEW'))
            fprintf(fid_skd, '-');
        elseif (strcmp(station(sched(isched).scan(iscan).sta(ista).staid).axis, 'AZEL'))
            az = sched(isched).scan(iscan).sta(ista).az;
            if (az>=station(sched(isched).scan(iscan).sta(ista).staid).azn1)&(az<=station(sched(isched).scan(iscan).sta(ista).staid).azn2)
                fprintf(fid_skd, '-');
            elseif (az>=station(sched(isched).scan(iscan).sta(ista).staid).azc1)&(az<=station(sched(isched).scan(iscan).sta(ista).staid).azc2)
                fprintf(fid_skd, 'C');
            elseif (az>=station(sched(isched).scan(iscan).sta(ista).staid).azw1)&(az<=station(sched(isched).scan(iscan).sta(ista).staid).azw2)
                fprintf(fid_skd, 'W');           
            end
        end
    end
    for ista = 1 : sched(isched).scan(iscan).nsta
       fprintf(fid_skd, ' 1F000000');
    end
    fprintf(fid_skd, ' YYNN');
    for ista = 1 : sched(isched).scan(iscan).nsta
       fprintf(fid_skd, '%6d', sched(isched).scan(iscan).sta(ista).duration);
    end
    fprintf(fid_skd, '\n');
end


