% Purpose  
%   Write output file in SKD format.
% History  
%   2010-04-06   Jing SUN   Created
%   2015-03-19   David MAYER   down time added
%   2015-07-22   Lucia PLANK ocodes revised
%   2016-05-10   Matthias SCHARTNER bugfix
%   2016-07-11   Matthias Schartner: oscource requires INFILE.sourcestar


function oskd(source, station, obsmode, sched, INFILE, fn_skd, PARA)

% number of observed sources
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

% open output file
fid_skd = fopen(fn_skd, 'w'); 

% write output file
%------$EXPER (experiment title)
fprintf(fid_skd, '$EXPER %s\n', PARA.EXPER);

%------$PARAM (parameters used by sked and durdg)
fprintf(fid_skd, '$PARAM\n');
oparam(fid_skd, station, srcnum_s, PARA);

%------$OP (automatic scheduling options)
fprintf(fid_skd, '$OP\n');
optimize(station, srcnum_s, fid_skd);

%------auxiliary
fprintf(fid_skd, '$MAJOR\n');
fprintf(fid_skd, 'Subnet ');
for ista = 1 : length(station)
    fprintf(fid_skd, '%s', station(ista).po);
end
fprintf(fid_skd, '\n');
fprintf(fid_skd, 'SkyCov             No\n');
fprintf(fid_skd, 'AllBlGood          No\n');
fprintf(fid_skd, 'MaxAngle       180.00\n');
fprintf(fid_skd, 'MinAngle        15.00\n');
fprintf(fid_skd, 'MinBetween         %d\n', PARA.MIN_SRCRP);
fprintf(fid_skd, 'MinSunDist      %5.2f\n', PARA.MIN_SUNDIST*180/pi);
fprintf(fid_skd, 'MaxSlewTime       300\n');
fprintf(fid_skd, 'TimeWindow          1.00\n');
fprintf(fid_skd, 'MinSubNetSize       2\n');
fprintf(fid_skd, 'NumSubNet           1\n');
fprintf(fid_skd, 'Best               60\n');
fprintf(fid_skd, 'FillIn            No\n');
fprintf(fid_skd, 'FillMinSub          3\n');
fprintf(fid_skd, 'FillMinTime       120\n');
fprintf(fid_skd, 'FillBest           80\n');
fprintf(fid_skd, 'Add_ps             30.0\n');
fprintf(fid_skd, 'SNRWts            Yes\n');
fprintf(fid_skd, '$MINOR\n');
fprintf(fid_skd, 'Astro           No   Abs          1.00\n');
fprintf(fid_skd, 'BegScan         No   Rel          1.00\n');
fprintf(fid_skd, 'EndScan         No   Rel          1.00\n');
fprintf(fid_skd, 'LowDec          No   Abs          1.00\n');
fprintf(fid_skd, 'NumLoEl         No   Abs          0.00        0.00\n');
fprintf(fid_skd, 'NumRiseSet      No   Abs          1.00\n');
fprintf(fid_skd, 'NumObs          No   Rel          1.00\n');
fprintf(fid_skd, 'SkyCov          No   Rel          1.00\n');
fprintf(fid_skd, 'SrcEvn          No   Abs          1.00 NONE\n');
fprintf(fid_skd, 'SrcWt           No   Abs          0.00\n');
fprintf(fid_skd, 'StatEvn         No   Abs          1.00 NONE\n');
fprintf(fid_skd, 'StatIdle        No   Abs          1.00\n');
fprintf(fid_skd, 'StatWt          No   Abs          1.00\n');
fprintf(fid_skd, 'TimeVar         No   Abs          1.00\n');
fprintf(fid_skd, '$ASTROMETRIC\n');
fprintf(fid_skd, '$SRCWT\n');
fprintf(fid_skd, '$STATWT\n');
fprintf(fid_skd, '$BROADBAND\n');
fprintf(fid_skd, '$DOWNTIME\n');
for i_down = 1:length([station.downum])
    for ii_down = 1:station(i_down).downum
        fprintf(fid_skd, '%s ',station(i_down).po); % station code
        %down time begin
        [down_start_DOY]=mjd2yydoysecod(station(i_down).downstart(ii_down));
        [down_start_year, down_start_month, down_start_day, down_start_hour, down_start_minu, down_start_sec] = mjd2date(station(i_down).downstart(ii_down));
        fprintf(fid_skd, '%04d-%03d-%02d:%02d:%02d ',down_start_year, down_start_DOY(2), down_start_hour, down_start_minu, round(down_start_sec));
        %down time end
        [down_start_DOY]=mjd2yydoysecod(station(i_down).downend(ii_down));
        [down_start_year, down_start_month, down_start_day, down_start_hour, down_start_minu, down_start_sec] = mjd2date(station(i_down).downend(ii_down));
        fprintf(fid_skd, '%04d-%03d-%02d:%02d:%02d\n',down_start_year, down_start_DOY(2), down_start_hour, down_start_minu, round(down_start_sec));
    end
end
fprintf(fid_skd, '$CATALOGS_USED\n');
fprintf(fid_skd, 'SOURCE    unknown     unknown\n');
fprintf(fid_skd, 'FLUX      unknown     unknown\n');
fprintf(fid_skd, 'ANTENNA   unknown     unknown\n');
fprintf(fid_skd, 'POSITION  unknown     unknown\n');
fprintf(fid_skd, 'EQUIP     unknown     unknown\n');
fprintf(fid_skd, 'MASK      unknown     unknown\n');
fprintf(fid_skd, 'MODES     unknown     unknown\n');
fprintf(fid_skd, 'FREQ      unknown     unknown\n');
fprintf(fid_skd, 'REC       unknown     unknown\n');
fprintf(fid_skd, 'RX        unknown     unknown\n');
fprintf(fid_skd, 'LOIF      unknown     unknown\n');
fprintf(fid_skd, 'TRACKS    unknown     unknown\n');
fprintf(fid_skd, 'HDPOS     unknown     unknown\n');

%------$SOURCES (list of sources for this experiment)
fprintf(fid_skd, '$SOURCES\n');
osource(source, sched, INFILE.source, INFILE.sourcestar, fid_skd);

%------$STATIONS (list of stations in this experiment)
fprintf(fid_skd, '$STATIONS\n');
stanum = length(station);
for ista = 1 : stanum
    staname(ista,1:8) = station(ista).name(1:8); 
end
% the A lines : antenna limits and rates
oantenna(staname, INFILE.antenna, fid_skd);
% the P lines : station position iformation
oposition(staname, INFILE.position, fid_skd);
% the T lines : station datea acquisition terminal information
staeqid = {station.eq};                         % added station id
oequip(staname, staeqid, INFILE.equip, fid_skd, PARA);
% the H lines : horizon mask
omask(staname, INFILE.mask, fid_skd);

station=orec(station,INFILE.rec,obsmode);
% get frequency information from freq.cat
obsmode=ofreq(obsmode,INFILE.freq); 

%------$CODES (frequency sequences and station LOs)
fprintf(fid_skd, '$CODES\n');
[obsmode] = ocodes(station, obsmode, INFILE, fid_skd, PARA);

%------$SKED (scheduled observations)
fprintf(fid_skd, '$SKED\n');
osked(source, station, obsmode, sched, fid_skd, PARA);

%------$HEAD (tape recorder head positions)
fprintf(fid_skd, '$HEAD\n');
ohead(station, obsmode, INFILE, fid_skd);

%------$FLUX (flux densities for each source)
fprintf(fid_skd, '$FLUX\n');
oflux(source, sched, INFILE.flux, fid_skd, PARA);

%------$PROC (station procedures)

% close output file
fclose(fid_skd);


