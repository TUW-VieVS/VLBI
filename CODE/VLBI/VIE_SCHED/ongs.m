% Purpose  
%   Write output file in NGS format.
% History  
%   2010-04-06   Jing SUN   Created
%   2016-09-06 Matthias Schartner: uses old output format in tdatestr.m 
%   2017-02-17 Matthias Schartner: filepath in header is now in an extra
%                                  row, (bugfix for long pathnames)


function ongs(source, station, sched, fn_ngs, twin)

% open output file
fid_ngs = fopen(fn_ngs, 'w'); 

% Header Card

if length(fn_ngs)<70
     fprintf(fid_ngs, 'DATA IN NGS FORMAT FROM DATABASE: \n');   
     fprintf(fid_ngs, '    %s\n', fn_ngs);
end

fprintf(fid_ngs, 'Observed delays and rates in card #2, modified errors in card #9\n');

% Site Cards
stanum = length(station);
for ista = 1 : stanum
    iftwin2 = 0;
    for i = 1 : twin.num
        if (ista == twin.sn(i,2))
            iftwin2 = 1;
            break;
        end
    end
    if(iftwin2 == 1)
        continue;
    end
    fprintf(fid_ngs, '%s  %15.5f%15.5f%15.5f %s%10.5f\n', station(ista).name(1:8), station(ista).xyz(1:3), station(ista).axis, station(ista).offset);
end
fprintf(fid_ngs, '$END\n');

% Radio source position cards
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
    srcid = srcsn_s(isrc);
    fprintf(fid_ngs, '%s  %02d %02d %12.6f %s%02d %02d %12.6f\n', source(srcid).name(1:8), source(srcid).rah, source(srcid).ram, source(srcid).ras, source(srcid).sign, source(srcid).ded, source(srcid).dem, source(srcid).des);         
end
fprintf(fid_ngs, '$END\n');

% Auxiliary parameters
fprintf(fid_ngs, '  .8212990000000D+04           GR PH\n');
fprintf(fid_ngs, '$END\n');

% Data Cards
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

obsnum = 0;
for i = 1 : num
    isched = startmjdsn(y(i),1);
    iscan = startmjdsn(y(i),2);
    srcid = sched(isched).scan(iscan).srcid;
    startmjd = sched(isched).scan(iscan).startmjd;   %%%
    [year, month, day, hour, minute, second] = tymdhms(startmjd);
    [datestr] = tdatestr(year, month, day, hour, minute, second,1);
    for s1 = 1 : sched(isched).scan(iscan).nsta-1
        staid1 = sched(isched).scan(iscan).sta(s1).staid;
        for s2 = (s1+1) : sched(isched).scan(iscan).nsta
            staid2 = sched(isched).scan(iscan).sta(s2).staid;
            obsnum = obsnum + 1;
            sta1name = station(staid1).name;
            sta2name = station(staid2).name;
            if (twin.num > 0)    % twin telescopes use same station name 
                for itwin = 1 : twin.num
                    if (staid1 == twin.sn(itwin,2)) 
                        sta1name = station(twin.sn(itwin,1)).name;
                    end
                    if (staid2 == twin.sn(itwin,2)) 
                        sta2name = station(twin.sn(itwin,1)).name;
                    end
                end
            end
            fprintf(fid_ngs, '%s  %s  %s %s          %8d%s\n', sta1name, sta2name, source(srcid).name, datestr, obsnum, '01');
            fprintf(fid_ngs, '    0000000.00000000    .00000  -000000.0000000000    .00000 0      I '); fprintf(fid_ngs, '%8d%s\n', obsnum,'02');
            fprintf(fid_ngs, '    .00000    .00000    .00000    .00000   0.000000000000000        0.'); fprintf(fid_ngs, '%8d%s\n', obsnum,'03');
            fprintf(fid_ngs, '       .00   .0       .00   .0       .00   .0       .00   .0          '); fprintf(fid_ngs, '%8d%s\n', obsnum,'04');
            fprintf(fid_ngs, '   -.00000   -.00000    .00000    .00000    .00000    .00000          '); fprintf(fid_ngs, '%8d%s\n', obsnum,'05');
            fprintf(fid_ngs, '     0.000    00.000   000.000   000.000    00.000    00.000 0 0      '); fprintf(fid_ngs, '%8d%s\n', obsnum,'06');
            fprintf(fid_ngs, '        0.0000000000    .00000        -.0000000000    .00000  0       '); fprintf(fid_ngs, '%8d%s\n', obsnum,'08');
            fprintf(fid_ngs, '          0.00000000    .00000        0.0000000000    .00000 0      I '); fprintf(fid_ngs, '%8d%s\n', obsnum,'09');
        end
    end
end

% close output file
fclose(fid_ngs);


