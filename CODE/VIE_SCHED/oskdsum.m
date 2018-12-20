% Purpose  
%   Write output file in SKD-sum format.
% History  
%   2010-04-06   Jing SUN   Created
%   2016-06-20 	 M. Schartner changes to improve speed


function oskdsum(source, station, obsmode, sched, fn_skdsum, PARA)

% open output file
fid_skdsum = fopen(fn_skdsum, 'w'); 
% 
fprintf(fid_skdsum, 'Session %s\n', PARA.EXPER);
%
[yy, mm, dd, hh, min, sec] = tymdhms(PARA.startmjd);
[cdays] = tdays(yy, mm, dd);
[datestr] = tdatestr(yy, mm, dd, hh, min, sec);
fprintf(fid_skdsum, 'Start %s (%4d-%s)\n', datestr, yy, cdays);
[yy, mm, dd, hh, min, sec] = tymdhms(PARA.endmjd);
[cdays] = tdays(yy, mm, dd);
[datestr] = tdatestr(yy, mm, dd, hh, min, sec);
fprintf(fid_skdsum, 'End   %s (%4d-%s)\n', datestr, yy, cdays);
%
stanum = length(station);
fprintf(fid_skdsum, 'Stations    %d\n', stanum);
for ista = 1 : stanum
    fprintf(fid_skdsum, '%s  %s  %s\n', station(ista).po, station(ista).name(1:8), station(ista).id);
end
fprintf(fid_skdsum, '\n');
%
fprintf(fid_skdsum, '   SKED Summary from file ./%s.skd for experiment %s\n', lower(PARA.EXPER), PARA.EXPER);  
fprintf(fid_skdsum, '\n');
fprintf(fid_skdsum, '      (all scans with at least one subnet station)\n');  
fprintf(fid_skdsum, '\n');
%
fprintf(fid_skdsum, '           4 chars/hour\n');
fprintf(fid_skdsum, '  SOURCE | 0           3           6           9          12          15          18          21          | #SCANS #OBS #Obs/bl\n');
srcnum = length(source);
srcsn_s(1:srcnum) = 0;
srcnum_s = 0;
for isched = 1 : length(sched)
    for iscan = 1 : sched(isched).nscan
        srcid = sched(isched).scan(iscan).srcid;
        if (size(find(srcsn_s(1:srcnum)==srcid),2) == 0)
            srcnum_s = srcnum_s + 1 ;
            srcsn_s(srcnum_s) = srcid;   
        end
    end
end
for isrc = 1 : srcnum_s
    rasort(isrc) = source(srcsn_s(isrc)).ra; 
end
[x, y] = sort(rasort);
for isrc = 1 : srcnum_s
    srcid = srcsn_s(y(isrc));
    fprintf(fid_skdsum, ' %s|', source(srcid).name);
    ctmp(1:96) = ' ';
    scanum1(isrc) = 0;
    obsnum1(isrc) = 0;
    for isched = 1 : length(sched)
        for iscan = 1 : sched(isched).nscan
            if (sched(isched).scan(iscan).srcid == srcid)
                [year, month, day, hour, minute, second] = tymdhms(sched(isched).scan(iscan).startmjd);
                cdex = hour * 4;
                if(minute < 15)
                    cdex = cdex + 1;
                elseif (minute>=15)&(minute<30)
                    cdex = cdex + 2;
                elseif (minute>=30)&(minute<45)
                    cdex = cdex + 3;
                elseif(minute>=45)&(minute<60)
                    cdex = cdex + 4;
                end
                ctmp(cdex) = 'x';
                scanum1(isrc) = scanum1(isrc) + 1;
                obsnum1(isrc) = obsnum1(isrc) + (sched(isched).scan(iscan).nsta*(sched(isched).scan(iscan).nsta-1))/2;
            end
        end
    end
    fprintf(fid_skdsum, '%s|', ctmp(1:96));
    fprintf(fid_skdsum, ' %5d', scanum1(isrc));
    fprintf(fid_skdsum, ' %5d', obsnum1(isrc));
    fprintf(fid_skdsum, ' %5.1f\n', 0.0);
end
fprintf(fid_skdsum, ' Total scans, obs:                                                                                          %5d %5d\n', sum(scanum1(1:srcnum_s)), sum(obsnum1(1:srcnum_s)));
fprintf(fid_skdsum, '\n');
%
fprintf(fid_skdsum, '           4 chars/hour\n');
fprintf(fid_skdsum, '  STATN  | 0           3           6           9          12          15          18          21          | #SCANS #OBS %%OBS\n');
for i = 1 : stanum
    fprintf(fid_skdsum, ' %s|', station(i).name(1:8));
    ctmp(1:96) = ' ';
    scanum2(i) = 0;
    obsnum2(i) = 0;
    for isched = 1 : length(sched)
        for iscan = 1 : sched(isched).nscan
            for ista = 1 : sched(isched).scan(iscan).nsta
                if (sched(isched).scan(iscan).sta(ista).staid == i)
                    [year, month, day, hour, minute, second] = tymdhms(sched(isched).scan(iscan).startmjd);
                    cdex = hour * 4;
                    if(minute < 15)
                        cdex = cdex + 1;
                    elseif (minute>=15)&(minute<30)
                        cdex = cdex + 2;
                    elseif (minute>=30)&(minute<45)
                        cdex = cdex + 3;
                    elseif(minute>=45)&(minute<60)
                        cdex = cdex + 4;
                    end
                    ctmp(cdex) = 'x';
                    scanum2(i) = scanum2(i) + 1;
                    obsnum2(i) = obsnum2(i) + (sched(isched).scan(iscan).nsta*(sched(isched).scan(iscan).nsta-1))/2;
                    break;
                end
            end
        end
    end
    fprintf(fid_skdsum, '%s|', ctmp(1:96));
    fprintf(fid_skdsum, ' %5d', scanum2(i));
    fprintf(fid_skdsum, ' %5d', obsnum2(i));
    fprintf(fid_skdsum, ' %5.1f\n', 100*obsnum2(i)/sum(obsnum1(1:srcnum_s)));
end
fprintf(fid_skdsum, '\n');

num = 0;
for isched = 1 : length(sched)
    for iscan = 1 : sched(isched).nscan
        thisScan = sched(isched).scan(iscan);
        idsource = sched(isched).scan(iscan).srcid;
        for sta1 = 1:thisScan.nsta-1;
            id1 = sched(isched).scan(iscan).sta(sta1).staid;
            for sta2 = sta1+1:thisScan.nsta
                id2 = sched(isched).scan(iscan).sta(sta2).staid;
                num = num+1;
                blsrc(num,1) = id1;
                blsrc(num,2) = id2;
                blsrc(num,3)  = idsource;
            end
        end
    end
end
tmp = blsrc(:,1:2);
blsrc(:,1:2) = sort(tmp,2);
[C,ia,ic] = unique(blsrc,'rows');

xtmp1 = length(ic);
xtmp2 = 0;
for i=1:max(ic)
    xtmp2 = xtmp2 + numel(ic(ic==i))^2;
end
num = length(ia);
M = mode(ic);
most = C(M,:);
sta1 = most(1);
sta2 = most(2);
src  = most(3);
maxva = numel(ic(ic==M));
minva = 10000;
for i=1:max(ic)
    tmp = numel(ic(ic==i));
    if tmp==1
        minva = 1;
        break
    end
    if tmp<minva
        minva=tmp;
    end
end


% for i1 = 1 : stanum-1
%     for i2 = i1+1 : stanum
%         for isched = 1 : length(sched)
%             for iscan = 1 : sched(isched).nscan
%                 ifuse1 = 0;
%                 ifuse2 = 0;
%                 for ista = 1 : sched(isched).scan(iscan).nsta
%                     if (sched(isched).scan(iscan).sta(ista).staid==i1)
%                         ifuse1 = 1;
%                     end
%                     if (sched(isched).scan(iscan).sta(ista).staid==i2)
%                         ifuse2 = 1;
%                     end
%                 end
%                 if ((ifuse1==0) | (ifuse2==0))
%                     continue;
%                 end
%                 ifuse0 = 0;
%                 for i = 1 : num
%                     if ((sched(isched).scan(iscan).srcid == blsrc(i).src) & (((i1==blsrc(i).sta1)&(i2==blsrc(i).sta2))|((i1==blsrc(i).sta2)&(i2==blsrc(i).sta1))) )
%                         ifuse0 = 1;
%                         blsrc(i).obsnum = blsrc(i).obsnum + 1;
%                         break;
%                     end
%                 end
%                 if (ifuse0 == 0)
%                     num = num + 1;
%                     blsrc(num).sta1 = i1;
%                     blsrc(num).sta2 = i2;
%                     blsrc(num).src  = sched(isched).scan(iscan).srcid;
%                     blsrc(num).obsnum = 1;
%                 end
%             end
%         end  
%     end
% end
% xtmp1 = 0;
% xtmp2 = 0;
% minva = 999999;
% maxva = 0;
% 
% for i = 1 : num
%     xtmp1 = xtmp1 + blsrc(i).obsnum;
%     xtmp2 = xtmp2 + blsrc(i).obsnum * blsrc(i).obsnum;
%     if (blsrc(i).obsnum < minva)
%         minva = blsrc(i).obsnum;
%     end
%     if (blsrc(i).obsnum > maxva)
%         maxva = blsrc(i).obsnum;
%         sta1 = blsrc(i).sta1;
%         sta2 = blsrc(i).sta2;
%         src  = blsrc(i).src;
%     end
% end


fprintf(fid_skdsum, ' Average number of obs. per baseline per source(normalized by up-time) =   %5.1f\n', (xtmp1/num));
fprintf(fid_skdsum, ' Min = %6d', minva);
fprintf(fid_skdsum, '  Max = %6d', maxva);
fprintf(fid_skdsum, ' (Baseline %s-%s on %s )', station(sta1).po, station(sta2).po, source(src).name(1:8));
fprintf(fid_skdsum, '  RMS = %5.1f\n', sqrt(xtmp2/num));
fprintf(fid_skdsum, '\n');
fprintf(fid_skdsum, ' Total time: %9.1f minutes (%6.1f hours)\n',  (PARA.endmjd-PARA.startmjd)*1440, (PARA.endmjd-PARA.startmjd)*24);
fprintf(fid_skdsum, '\n');
%
for i = 1 : ceil(stanum/6)
    if (i == 1)
        fprintf(fid_skdsum, ' Key:     ');
    else
        fprintf(fid_skdsum, '          ');
    end
    if (((i-1)*6+1)<=stanum)
        fprintf(fid_skdsum, '%s=%s   ', station((i-1)*6+1).po, station((i-1)*6+1).name(1:8));
    end
    if (((i-1)*6+2)<=stanum)
        fprintf(fid_skdsum, '%s=%s   ', station((i-1)*6+2).po, station((i-1)*6+2).name(1:8));
    end
    if (((i-1)*6+3)<=stanum)
        fprintf(fid_skdsum, '%s=%s   ', station((i-1)*6+3).po, station((i-1)*6+3).name(1:8));
    end
    if (((i-1)*6+4)<=stanum)
        fprintf(fid_skdsum, '%s=%s   ', station((i-1)*6+4).po, station((i-1)*6+4).name(1:8));
    end
    if (((i-1)*6+5)<=stanum)
        fprintf(fid_skdsum, '%s=%s   ', station((i-1)*6+5).po, station((i-1)*6+5).name(1:8));
    end
    if (((i-1)*6+6)<=stanum)
        fprintf(fid_skdsum, '%s=%s   ', station((i-1)*6+6).po, station((i-1)*6+6).name(1:8));
    end
    fprintf(fid_skdsum, '\n');
end
fprintf(fid_skdsum, '                   ');                   
for i = 1 : stanum
    fprintf(fid_skdsum, '%s   ',station(i).po);
end
fprintf(fid_skdsum, 'Avg\n');
for i = 1 : stanum
    obs(i)  = 0.0;
    cal(i)  = 0.0;
    slew(i) = 0.0;
    idle(i) = 0.0;
    for isched = 1 : length(sched)
        for iscan = 1 : sched(isched).nscan
            for ista = 1 : sched(isched).scan(iscan).nsta
                if (sched(isched).scan(iscan).sta(ista).staid == i)
                    obs(i)  = obs(i) + sched(isched).scan(iscan).sta(ista).duration;
                    cal(i)  = cal(i) + PARA.CALIBRATION + PARA.SOURCE + PARA.TAPETM + PARA.IDLE;
                    slew(i) = slew(i) + sched(isched).scan(iscan).sta(ista).slewtime;
                    idle(i) = idle(i) + (sched(isched).scan(iscan).startmjd - sched(isched).scan(iscan).sta(ista).startmjd)*86400;
                end
            end
        end
    end
end
fprintf(fid_skdsum, ' %% obs. time:      ');
avg = 0;
for i = 1 : stanum
    fprintf(fid_skdsum, '%2d   ', round(100*obs(i)/(PARA.duration*3600)));
    avg = avg + 100*obs(i)/(PARA.duration*3600);
end
fprintf(fid_skdsum, '%2d\n', round(avg/stanum));
fprintf(fid_skdsum, ' %% cal. time:      ');
avg = 0;
for i = 1 : stanum
    fprintf(fid_skdsum, '%2d   ', round(100*cal(i)/(PARA.duration*3600)));
    avg = avg + 100*cal(i)/(PARA.duration*3600);
end
fprintf(fid_skdsum, '%2d\n', round(avg/stanum));
fprintf(fid_skdsum, ' %% slew time:      ');
avg = 0;
for i = 1 : stanum
    fprintf(fid_skdsum, '%2d   ', round(100*slew(i)/(PARA.duration*3600)));
    avg = avg + 100*slew(i)/(PARA.duration*3600);
end
fprintf(fid_skdsum, '%2d\n', round(avg/stanum));
fprintf(fid_skdsum, ' %% idle time:      ');
avg = 0;
for i = 1 : stanum
    fprintf(fid_skdsum, '%2d   ', round(100*idle(i)/(PARA.duration*3600)));
    avg = avg + 100*idle(i)/(PARA.duration*3600);
end
fprintf(fid_skdsum, '%2d\n', round(avg/stanum));
fprintf(fid_skdsum, ' total # scans:  ');
avg = 0;
for i = 1 : stanum
    fprintf(fid_skdsum, '%4d ', scanum2(i));
    avg = avg + scanum2(i);
end
fprintf(fid_skdsum, '%4d\n', round(avg/stanum));
fprintf(fid_skdsum, ' # scans/hour :  ');
avg = 0;
for i = 1 : stanum
    fprintf(fid_skdsum, '%4d ', ceil(scanum2(i)/((PARA.endmjd-PARA.startmjd)*24)));
    avg = avg + scanum2(i)/((PARA.endmjd-PARA.startmjd)*24);
end
fprintf(fid_skdsum, '%4d\n', ceil(avg/stanum));
fprintf(fid_skdsum, ' Avg scan (sec): ');
avg = 0;
for i = 1 : stanum
    fprintf(fid_skdsum, '%4d ', round(obs(i)/scanum2(i)));
    avg = avg + obs(i)/scanum2(i);
end
fprintf(fid_skdsum, '%4d\n', round(avg/stanum));
fprintf(fid_skdsum, ' Total GBytes:   ');
avg = 0;
for i = 1 : stanum
    fprintf(fid_skdsum, '%4d ', ceil(obs(i)*(obsmode.samprate*obsmode.channelnum*obsmode.bits/1024)/8));
    avg = avg + ceil(obs(i)*(obsmode.samprate*obsmode.channelnum*obsmode.bits/1024)/8);
end
fprintf(fid_skdsum, '%4d\n', round(avg/stanum));
fprintf(fid_skdsum, '\n');
%
fprintf(fid_skdsum, '      # OF OBSERVATIONS BY BASELINE\n');
fprintf(fid_skdsum, '  |  ');
for i = 1 : stanum
    fprintf(fid_skdsum, '%s   ',  station(i).po);
end
fprintf(fid_skdsum, 'StnTotal\n');
fprintf(fid_skdsum, '-----');
for i = 1 : stanum
    fprintf(fid_skdsum, '-----');
end
fprintf(fid_skdsum, '--------\n');
for i = 1 : stanum-1
    fprintf(fid_skdsum, '%s|',  station(i).po);
    for j = 1 : i
        fprintf(fid_skdsum, '     ');
    end
    for j = i+1 : stanum
        num = 0;
        for isched = 1 : length(sched)
            for iscan = 1 : sched(isched).nscan
                ifuse1 = 0;
                ifuse2 = 0;
                for ista = 1 : sched(isched).scan(iscan).nsta
                    if (sched(isched).scan(iscan).sta(ista).staid == i)
                        ifuse1 = 1;
                    end
                    if (sched(isched).scan(iscan).sta(ista).staid == j)
                        ifuse2 = 1;
                    end
                end
                if ((ifuse1==0) | (ifuse2==0))
                    continue;
                end
                num = num + 1;
            end
        end
        fprintf(fid_skdsum, '%4d ', num);
    end
    fprintf(fid_skdsum, '%7d\n', obsnum2(i));
end
fprintf(fid_skdsum, '%s|',  station(stanum).po);
for j = 1 : stanum
    fprintf(fid_skdsum, '     ');
end
fprintf(fid_skdsum, '%7d\n', obsnum2(stanum));
fprintf(fid_skdsum, '\n');
num(1:stanum) = 0;
for isched = 1 : length(sched)
    for iscan = 1 : sched(isched).nscan
        num(sched(isched).scan(iscan).nsta) = num(sched(isched).scan(iscan).nsta) + 1;
    end
end
for i = 2 : stanum
    fprintf(fid_skdsum, ' Number of  %s-station scans: %d\n', num2str(i), num(i));
end
fprintf(fid_skdsum, '\n');
fprintf(fid_skdsum, ' Total # of scans, observations:   %d     %d\n', sum(scanum1(1:srcnum_s)), sum(obsnum1(1:srcnum_s)));

% number of scans in a subnet
scanum(1:10) = 0;
for isched = 1 : length(sched)
    scanum(sched(isched).nscan) = scanum(sched(isched).nscan) + 1;
end
for i = 1 : 10
    fprintf(fid_skdsum, ' Number of  %s-scan subnets: %d\n', num2str(i), scanum(i));
end

% close output file
fclose(fid_skdsum);


