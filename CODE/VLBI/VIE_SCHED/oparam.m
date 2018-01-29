% Purpose  
%   Write $PARAM (parameters used by sked and durdg).
% History  
%   2010-04-06   Jing SUN       Created
%   2015-03-19   David MAYER    Parameter SCHEDULER added/ SOFTWARE VERSION
%                               updated
%	2016-11-04 Matthias Schartner Bugfix


function oparam(fid_skd, station, srcnum_s, PARA)

%
fprintf(fid_skd, 'DESCRIPTION %s\n', PARA.DESCRIPTION);
fprintf(fid_skd, 'SCHEDULING_SOFTWARE VIE_SCHED\n');
p = path;
str_pointer = strfind(p, 'VIE_SCHED');
vie_sched_version_str = p(str_pointer(1) : (str_pointer(1) + length('VIE_SCHED_V') + 1));
fprintf(fid_skd, 'SOFTWARE_VERSION %s\n', vie_sched_version_str);
t = clock;
[tstr1] = sprintf('%4d/%02d/%02d', round(t(1)),round(t(2)),round(t(3)));
[tstr2] = sprintf('%02d:%02d:%02d', round(t(4)),round(t(5)),round(t(6)));
fprintf(fid_skd, 'SCHEDULE_CREATE_DATE %s %s\n', tstr1, tstr2);
[year, month, day, hour, minute, secf] = tymdhms(PARA.startmjd);
dsec = round(secf);
minute(dsec==60)=minute(dsec==60)+1;
dsec(dsec==60)=0;
hour(minute==60)=hour(minute==60)+1;
minute(minute==60)=0;
[cdays] = tdays(year, month, day);
[tsstr] = sprintf('%4d%s%02d%02d%02d', year,cdays,hour,minute,dsec);
[year, month, day, hour, minute, secf] = tymdhms(PARA.endmjd);
dsec = round(secf);
minute(dsec==60)=minute(dsec==60)+1;
dsec(dsec==60)=0;
hour(minute==60)=hour(minute==60)+1;
minute(minute==60)=0;
[cdays] = tdays(year, month, day);
[testr] = sprintf('%4d%s%02d%02d%02d', year,cdays,hour,minute,dsec);
fprintf(fid_skd, 'SCHEDULER  %s   CORRELATOR %s START %s END   %s\n',PARA.SCHEDULER, PARA.CORRELATOR, tsstr, testr);
%
fprintf(fid_skd, 'CALIBRATION%6d', PARA.CALIBRATION);
fprintf(fid_skd, ' CORSYNCH%9d',   PARA.CORSYNCH);
fprintf(fid_skd, ' DURATION%9d\n', PARA.SCANDURA);
%
fprintf(fid_skd, 'EARLY%12d', 0);
fprintf(fid_skd, ' IDLE%13d', PARA.IDLE);
fprintf(fid_skd, ' LOOKAHEAD%8d\n', 20);
%
fprintf(fid_skd, 'MAXSCAN%10d', PARA.MAX_SCAN);
fprintf(fid_skd, ' MINSCAN%10d', PARA.MIN_SCAN);
fprintf(fid_skd, ' MINIMUM%10d\n', 0);
%
fprintf(fid_skd, 'MIDTP%12d', 0);
fprintf(fid_skd, ' MODULAR%10d', 1);
fprintf(fid_skd, ' MODSCAN%10d', 10);
fprintf(fid_skd, ' PARITY%11d\n', 100);
%
fprintf(fid_skd, 'SETUP%12d', PARA.CALIBRATION);
fprintf(fid_skd, ' SOURCE%11d', PARA.SOURCE);
fprintf(fid_skd, ' TAPETM%11d', PARA.TAPETM);
fprintf(fid_skd, ' WIDTH%12d\n', 0);
%
fprintf(fid_skd, 'CONFIRM         %s', 'Y');
fprintf(fid_skd, ' VSCAN           %s\n', 'Y');
%
fprintf(fid_skd, 'DEBUG           %s', 'N');
fprintf(fid_skd, '   KEEP_LOG      %s', 'N');
fprintf(fid_skd, ' VERBOSE         %s\n', 'N');
%
fprintf(fid_skd, 'PRFLAG       %s', 'YYNN');
fprintf(fid_skd, ' SNR          %s\n', 'AUTO');
%
fprintf(fid_skd, 'FREQUENCY   SX');
fprintf(fid_skd, ' PREOB      PREOB  MIDOB     MIDOB  POSTOB     POSTOB\n');
%
fprintf(fid_skd, 'ELEVATION _  %4.1f\n', PARA.MIN_CUTEL*180/pi);

% TAPE
stanum = length(station);
fprintf(fid_skd, 'TAPE_MOTION _  START&STOP\n');
fprintf(fid_skd, 'TAPE_TYPE ');
for i = 1 : stanum
    fprintf(fid_skd, '%s %s ', station(i).po, deblank(station(i).equip(1:8)));
end
fprintf(fid_skd, '\n');
fprintf(fid_skd, 'TAPE_ALLOCATION _  SCHEDULED\n');
% SNR
num = 0;
blnum = stanum*(stanum-1)/2;
fprintf(fid_skd, 'SNR ');
for i = 1 : (stanum-1)
    for j = (i+1) : stanum
        fprintf(fid_skd, '%s-%s X %4d ', station(i).po, station(j).po, min(station(i).minsnr(1),station(j).minsnr(1)));
        fprintf(fid_skd, '%s-%s S %4d ', station(i).po, station(j).po, min(station(i).minsnr(2),station(j).minsnr(2)));
        num = num + 1;
        if ((rem(num,3) == 0) & (num ~= blnum))
            fprintf(fid_skd, '\n');
            fprintf(fid_skd, 'SNR ');
        end
    end
end
fprintf(fid_skd, '\n');
fprintf(fid_skd, 'SNR MARGIN X   0 MARGIN S   0\n');
fprintf(fid_skd, 'SNR AST_MARGIN X   0 AST_MARGIN S   0\n');
% SCAN
num = 0;
fprintf(fid_skd, 'SCAN   ');
for isrc = 1 : srcnum_s
    fprintf(fid_skd, '%3d %3d ', isrc, PARA.SCANDURA);
    num = num + 1;
    if ((rem(num,8)==0) & (num ~= srcnum_s))
        fprintf(fid_skd, '\n');
        fprintf(fid_skd, 'SCAN   ');
    end
end
fprintf(fid_skd, '\n');


