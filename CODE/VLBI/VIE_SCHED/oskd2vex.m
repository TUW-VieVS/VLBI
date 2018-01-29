% Purpose  
%   Write output file in VEX format from skd file.
% History  
%   2013-05-06   Jing SUN   Created
%   2015-07-23   Lucia PLANK major revision 
%   2015-07-27   Lucia PLANK tracks section
%   2015-10-02   David MAYER fixed bug in $MODE block



function oskd2vex(fn_sked, fn_vex)

CodesNum = 14;
PARA.MAX_BANDNUM = 2;
PARA.MAX_SEFDPARA = 4;
%
% open .skd file
fidskd = fopen(fn_sked, 'r'); 
if (fidskd < 0)
    error('    Error: open skd file %s !\n', fn_sked);
end
% open .vex file
fid = fopen(fn_vex, 'w'); 
if (fid < 0)
    error('    Error: open vex file %s !\n', fn_vex);
end
%
% reading $STATIONS information from skd file
[station] = iskdstation(fidskd);
stanum = length(station);
%
% 
fprintf(fid, 'VEX_rev = 1.5;\n');
frewind(fidskd);
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if ((linelength>16) & strcmp(line(1:16),'SOFTWARE_VERSION')) 
        [tmps, count, errmsg, nextindex] = sscanf(line(17:linelength), '%s', 1);
        fprintf(fid, '*  scheduling_software version %s\n', tmps);
        clear tmps;
        break;
    end  
end
fprintf(fid, '$GLOBAL;\n');
frewind(fidskd);
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if ((linelength>6) & strcmp(line(1:6),'$EXPER')) 
        [EXPER, count, errmsg, nextindex] = sscanf(line(7:linelength), '%s', 1);
        fprintf(fid, '    ref $EXPER = %s;\n', EXPER);
        fprintf(fid, '    ref $SCHEDULING_PARAMS = VIE_SCHED_PARAMS;\n');
        fprintf(fid, '$EXPER;\n');
        fprintf(fid, '  def %s;\n', EXPER);
        fprintf(fid, '    exper_name = %s;\n', EXPER);
        break;
    end
end
frewind(fidskd);
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if ((linelength>11) & strcmp(line(1:11),'DESCRIPTION')) 
        [tmps, count, errmsg, nextindex] = sscanf(line(12:linelength), '%s', 1);
        fprintf(fid, '    exper_description = %s;\n', tmps);
        clear tmps;
        break;
    end 
end
frewind(fidskd);
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if ((linelength>9) & strcmp(line(1:9),'SCHEDULER')) 
        [tmps, count, errmsg, nextindex] = sscanf(line(10:linelength), '%s', 1);
        fprintf(fid, '    PI_name = %s;\n', tmps);
        clear tmps;
        index = 10 + nextindex - 1;
        [tmps, count, errmsg, nextindex] = sscanf(line(index:linelength), '%s', 1);
        index = index + nextindex - 1;
        [tmps, count, errmsg, nextindex] = sscanf(line(index:linelength), '%s', 1);
        fprintf(fid, '    target_correlator = %s;\n', tmps);
        clear tmps;
        break;
    end
end
fprintf(fid, '  enddef;\n');
fprintf(fid, '*-----------------------   end $EXPER              ----------------------*\n');
%
% % $MODE
fprintf(fid, '*----------------------- begin $MODE               ----------------------*\n');
fprintf(fid, '$MODE;\n');
frewind(fidskd);
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if ((linelength<6) | ~strcmp(line(1:6),'$CODES')) 
        continue;
    end
    while ~feof(fidskd)
        line = fgetl(fidskd);
        linelength = length(line);
        if (strcmp(line(1:2),'F '))
            [tmps1, count, errmsg, nextindex] = sscanf(line(3:linelength), '%s', 1);
            index = 3 + nextindex - 1;
            [tmps2, count, errmsg, nextindex] = sscanf(line(index:linelength), '%s', 1);
            fprintf(fid, '  def %s.%s;\n', tmps1,tmps2);
            MODENAME = strcat(tmps1,'-',tmps2);
            clear tmps1;
            clear tmps2;
            break;
        end
    end
end
fprintf(fid, '    ref $FREQ = %s01', MODENAME);
for ista = 1 : stanum
    fprintf(fid, ':%s', station(ista).po(1:2));
end
fprintf(fid, ';\n');


% read CODES C section
nchan = 0;
frewind(fidskd);
ncodes=1;
stat_inblock = 0;
while  ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if (linelength>=6) && (strcmp(line(1:6),'$CODES'))
        line = fgetl(fidskd);
        skdcodes(ncodes).stations=line(15:end);
    end
    
    if (linelength>=5) && (strcmp(line(1:5),'C SX '))
        temp = '';
        temp = regexp(skdcodes(ncodes).stations,' ','split'); stat_inblock = stat_inblock + length(temp); %get number of stations in this block
        for ifreq=1:CodesNum
        nchan = nchan + 1;
        readlin=textscan(line, '%s %s %s %s %s %s %s %s %s');
        tmps7 = char(readlin{9});
        skdcodes(ncodes).lin(nchan,:)=readlin;
        skdcodes(ncodes).channel(nchan).ul='U';
        if  isempty(strfind(tmps7,',,'))
            nchan = nchan + 1;
            skdcodes(ncodes).channel(nchan).ul='L';
            skdcodes(ncodes).lin(nchan,:)=readlin;
        end
        clear tmps7;
        if ifreq==14
            nchan=0;
            ncodes=ncodes+1;
        end
        line=fgetl(fidskd);
        if strcmp(line(1:2),'F ')
           skdcodes(ncodes).stations=line(15:end);
           skdcodes(ncodes).lin = {};
        end
        end  
        if stat_inblock == stanum
            break
        end
    end
end



frewind(fidskd);
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if ((linelength<6) | ~strcmp(line(1:6),'$CODES')) 
        continue;
    end
    while ~feof(fidskd)
        line = fgetl(fidskd);
        linelength = length(line);
        if (line(1) == '$')
            break;
        end
        if ~strcmp(line(1:2),'L ')
            continue;
        end
        [tmps, count, errmsg, nextindex] = sscanf(line(3:linelength), '%s', 1);
        for ista = 1 : stanum
            if strcmp(line(3), station(ista).id)
                for ic=1:CodesNum
                    readlsec=textscan(line,'%s %s %s %s %s %s %s %s');
                    station(ista).bbc(ic) = readlsec{5};
                    station(ista).if(ic)  = readlsec{6};
                    line = fgetl(fidskd);
                end
            end
        end
    end
end
bbcnum = 1;
bbc(1).bn = station(1).bbc;
bbc(1).sn = strcat(':',station(1).po);
for ista = 2 : stanum
    ifnew = 1;
    for i = 1 : bbcnum
        if strcmp(station(ista).bbc, bbc(i).bn)
            bbc(i).sn = strcat(bbc(i).sn,':',station(ista).po);    
            ifnew = 0;
        end
    end
    if (ifnew == 1)
        bbcnum = bbcnum + 1;
        bbc(bbcnum).bn = station(ista).bbc;
        bbc(bbcnum).sn = strcat(':',station(ista).po);
    end
end
for ibbc = 1 : bbcnum
    fprintf(fid, '    ref $BBC = %s%02d%s;\n', MODENAME,ibbc,bbc(ibbc).sn); 
end
ifrenum = 1;
ifre(1).bn = station(1).if;
ifre(1).bbc = station(1).bbc;
ifre(1).sn = strcat(':',station(1).po);

for ista = 2 : stanum
    ifnew = 1;
    for i = 1 : ifrenum
        if strcmp(station(ista).if, ifre(i).bn)
            ifre(i).sn = strcat(ifre(i).sn,':',station(ista).po);    
            ifnew = 0;
        end
    end
    if (ifnew == 1)
        ifrenum = ifrenum  + 1;
        ifre(ifrenum).bbc = station(ista).bbc;
        ifre(ifrenum).bn = station(ista).if;
        ifre(ifrenum).sn = strcat(':',station(ista).po);
    end
end
for iif = 1 : ifrenum
    fprintf(fid, '    ref $IF = %s%02d%s;\n', MODENAME,iif,ifre(iif).sn); 
end

for it=1:length(skdcodes)
    defnam=char(skdcodes(it).lin{1,7});
    fprintf(fid, '    ref $TRACKS = %s_2f_1b-SX0%d',defnam(1:end-2),it);
    l=textscan(skdcodes(it).stations,'%s ');
    lst=l{1};
    for is=1:length(lst)
        ind=strmatch(lst{is},{station.name});
        if isempty(ind)
            fprintf(1,'WARNING!: Station %s not found in skd file!!!',lst{is});
        end
        fprintf(fid, ':%s', station(ind).po);
    end
    fprintf(fid, ';\n'); clear l lst
end

for ista = 1 : stanum
    station(ista).head = ' ';
end
frewind(fidskd);
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if ((linelength<5) | ~strcmp(line(1:5),'$HEAD')) 
        continue;
    end
    while ~feof(fidskd)
        line = fgetl(fidskd);
        linelength = length(line);
        if (line(1) == '$')
            break;
        end
        for ista = 1 : stanum
            if strcmp(line(1), station(ista).id)
                station(ista).head = strcat(station(ista).head, line(8:linelength));
                break;
            end
        end
    end
end
headnum = 1;
head(1).bn = station(1).head;
head(1).sn = strcat(':',station(1).po);
for ista = 2 : stanum
    ifnew = 1;
    for i = 1 : headnum
        if strcmp(station(ista).head, head(i).bn)
            head(i).sn = strcat(head(i).sn,':',station(ista).po);    
            ifnew = 0;
        end
    end
    if (ifnew == 1)
        headnum = headnum + 1;
        head(headnum).bn = station(ista).head;
        head(headnum).sn = strcat(':',station(ista).po);
    end
end
for ih = 1 : headnum
    fprintf(fid, '    ref $HEAD_POS = Mk341-SX%02dS%02d%s;\n', ih,(2*ih-1),head(ih).sn);
end
for ih = 1 : headnum
    fprintf(fid, '    ref $PASS_ORDER = Mk341-SX%02dS%02d%s;\n', ih,(2*ih-1),head(ih).sn);
end
fprintf(fid, '    ref $ROLL = NO_ROLL');
for ista = 1 : stanum
    fprintf(fid, ':%s', station(ista).po);
end
fprintf(fid, ';\n');
fprintf(fid, '    ref $PHASE_CAL_DETECT = Standard');
for ista = 1 : stanum
    fprintf(fid, ':%s', station(ista).po);
end
fprintf(fid, ';\n');
fprintf(fid, '    ref $TRACKS = Mark4_format');
for ista = 1 : stanum
    if ~strcmp(station(ista).equip,'MARK5B')
        fprintf(fid, ':%s', station(ista).po);
    end
end
fprintf(fid, ';\n');
fprintf(fid, '    ref $TRACKS = Mark5B_format');
for ista = 1 : stanum
    if strcmp(station(ista).equip,'MARK5B')
        fprintf(fid, ':%s', station(ista).po);
    end
end
fprintf(fid, ';\n');
fprintf(fid, '  enddef;\n');
fprintf(fid, '*-----------------------   end $MODE               ----------------------*\n');
%
% $SCHEDULING_PARAMS
fprintf(fid, '*----------------------- begin $SCHEDULING_PARAMS  ----------------------*\n');
fprintf(fid, '$SCHEDULING_PARAMS;\n');
fprintf(fid, '  def VIE_SCHED_PARAMS;\n');
fprintf(fid, 'start_literal( );\n');
% $PARAM
frewind(fidskd);
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if ((linelength>=6) & strcmp(line(1:6),'$PARAM')) 
        fprintf(fid, '%s\n', line(1:linelength));
        break;
    end
end
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if (line(1)=='$')
        break;
    end
    fprintf(fid, '%s\n', line(1:linelength));
end
% $OP
frewind(fidskd);
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if ((linelength>=3) & strcmp(line(1:3),'$OP')) 
        fprintf(fid, '%s\n', line(1:linelength));
        break;
    end
end
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if (line(1)=='$')
        break;
    end
    fprintf(fid, '%s\n', line(1:linelength));
end
% $MAJOR
frewind(fidskd);
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if ((linelength>=6) & strcmp(line(1:6),'$MAJOR')) 
        fprintf(fid, '%s\n', line(1:linelength));
        break;
    end
end
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if (line(1)=='$')
        break;
    end
    fprintf(fid, '%s\n', line(1:linelength));
end
% $MINOR
frewind(fidskd);
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if ((linelength>=6) & strcmp(line(1:6),'$MINOR')) 
        fprintf(fid, '%s\n', line(1:linelength));
        break;
    end
end
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if (line(1)=='$')
        break;
    end
    fprintf(fid, '%s\n', line(1:linelength));
end
% $ASTROMETRIC
frewind(fidskd);
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if ((linelength>=12) & strcmp(line(1:12),'$ASTROMETRIC')) 
        fprintf(fid, '%s\n', line(1:linelength));
        break;
    end
end
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if (line(1)=='$')
        break;
    end
    fprintf(fid, '%s\n', line(1:linelength));
end
% $SRCWT
frewind(fidskd);
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if ((linelength>=6) & strcmp(line(1:6),'$SRCWT')) 
        fprintf(fid, '%s\n', line(1:linelength));
        break;
    end
end
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if (line(1)=='$')
        break;
    end
    fprintf(fid, '%s\n', line(1:linelength));
end
% $STATWT
frewind(fidskd);
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if ((linelength>=7) & strcmp(line(1:7),'$STATWT')) 
        fprintf(fid, '%s\n', line(1:linelength));
        break;
    end
end
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if (line(1)=='$')
        break;
    end
    fprintf(fid, '%s\n', line(1:linelength));
end
%
fprintf(fid, '$BROADBAND\n');
% $DOWNTIME
frewind(fidskd);
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if ((linelength>=9) & strcmp(line(1:9),'$DOWNTIME')) 
        fprintf(fid, '%s\n', line(1:linelength));
        break;
    end
end
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if (line(1)=='$')
        break;
    end
    fprintf(fid, '%s\n', line(1:linelength));
end
%
fprintf(fid, '$STAT_SEFD\n');
for ista = 1 : stanum
    fprintf(fid, '%s', station(ista).name);
    fprintf(fid, ' X  ', station(ista).name);
    fprintf(fid, '%d  ', station(ista).sefdpara(1,1));
    for j = 2 : PARA.MAX_SEFDPARA
        if (station(ista).sefdpara(1,j) > 1d-3)
            fprintf(fid, '%f ', station(ista).sefdpara(1,j));
        end
    end
    fprintf(fid, '\n');
    fprintf(fid, '%s', station(ista).name);
    fprintf(fid, ' S  ', station(ista).name);
    fprintf(fid, '%d  ', station(ista).sefdpara(2,1));
    for j = 2 : PARA.MAX_SEFDPARA
        if (station(ista).sefdpara(2,j) > 1d-3)
            fprintf(fid, '%f ', station(ista).sefdpara(2,j));
        end
    end
    fprintf(fid, '\n');
end
fprintf(fid, '$CATALOGS_USED\n');
fprintf(fid, 'SOURCE    unknown     unknown\n');
fprintf(fid, 'FLUX      unknown     unknown\n');
fprintf(fid, 'ANTENNA   unknown     unknown\n');
fprintf(fid, 'EQUIP     unknown     unknown\n');
fprintf(fid, 'MASK      unknown     unknown\n');
fprintf(fid, 'MODES     unknown     unknown\n');
fprintf(fid, 'FREQ      unknown     unknown\n');
fprintf(fid, 'REC       unknown     unknown\n');
fprintf(fid, 'RX        unknown     unknown\n');
fprintf(fid, 'LOIF      unknown     unknown\n');
fprintf(fid, 'TRACKS    unknown     unknown\n');
fprintf(fid, 'HDPOS     unknown     unknown\n');
fprintf(fid, '$FLUX\n');
frewind(fidskd);
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if ((linelength>=5) & strcmp(line(1:5),'$FLUX')) 
        break;
    end
end
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if (strcmp(line(1),'$') | strcmp(line(1),' ')) 
        break;
    end
    fprintf(fid, '%s\n', line(1:linelength));
end
fprintf(fid, '$END\n');
fprintf(fid, 'end_literal( );\n');
fprintf(fid, '  enddef;\n');
fprintf(fid, '*-----------------------   end $SCHEDULING_PARAMS  ----------------------*\n');
%
% % $STATION
fprintf(fid, '*----------------------- begin $STATION            ----------------------*\n');
fprintf(fid, '$STATION;\n');
for ista = 1 : stanum
    fprintf(fid, '  def %s;\n', station(ista).po);
    fprintf(fid, '    ref $SITE = %s;\n', deblank(station(ista).name));
    fprintf(fid, '    ref $ANTENNA = %s;\n', deblank(station(ista).antname));
    fprintf(fid, '    ref $DAS = %s_rack;\n', station(ista).rack);
    fprintf(fid, '    ref $DAS = %s_%s;\n', station(ista).po,deblank(station(ista).eq));
    fprintf(fid, '    ref $DAS = %s_recorder;\n', station(ista).equip);
    fprintf(fid, '    ref $DAS = thick_tape;\n');
    fprintf(fid, '    ref $DAS = high_density;\n');
    fprintf(fid, '  enddef;\n');
end
fprintf(fid, '*-----------------------   end $STATION            ----------------------*\n');
%
% $ANTENNA
fprintf(fid, '*----------------------- begin $ANTENNA            ----------------------*\n');
fprintf(fid, '$ANTENNA;\n');
for ista = 1 : stanum
    fprintf(fid, '  def %s;\n', deblank(station(ista).antname));
    fprintf(fid, '    antenna_diam = %6.2f m;\n', station(ista).diam);
    if strcmp(station(ista).axis(1:4),'AZEL')
        axis1 = 'az';
        axis2 = 'el';
    elseif strcmp(station(ista).axis(1:4),'HADC')
        axis1 = 'ha';
        axis2 = 'dec';
    elseif strcmp(station(ista).axis(1:4),'XYEW')
        axis1 = 'x';
        axis2 = 'yew';
    end
    fprintf(fid, '    axis_type = %s : %s;\n', deblank(axis1),deblank(axis2));   
    fprintf(fid, '    axis_offset = %10.5f m;\n', station(ista).offset);
    fprintf(fid, '    antenna_motion = %s : %5.1f deg/min : %5d sec;\n', axis1, station(ista).rate1, station(ista).c1);
    fprintf(fid, '    antenna_motion = %s : %5.1f deg/min : %5d sec;\n', axis2, station(ista).rate2, station(ista).c2);
    fprintf(fid, '    pointing_sector = &n : %s : %6.1f deg : %6.1f deg : %s : %6.1f deg : %6.1f deg;\n', axis1,station(ista).lim11,station(ista).lim12,axis2,station(ista).lim21,station(ista).lim22);
%     fprintf(fid, '    pointing_sector = &ccw : %s : %6.1f deg : %6.1f deg : %s : %6.1f deg : %6.1f deg;\n', lower(station(ista).axis(1:2)),station(ista).azw1,station(ista).azw2,lower(station(ista).axis(3:4)),station(ista).lim21,station(ista).lim22);
%     fprintf(fid, '    pointing_sector = &n : %s : %6.1f deg : %6.1f deg : %s : %6.1f deg : %6.1f deg;\n', lower(station(ista).axis(1:2)),station(ista).azn1,station(ista).azn2,lower(station(ista).axis(3:4)),station(ista).lim21,station(ista).lim22);
%     fprintf(fid, '    pointing_sector = &cw : %s : %6.1f deg : %6.1f deg : %s : %6.1f deg : %6.1f deg;\n', lower(station(ista).axis(1:2)),station(ista).azc1,station(ista).azc2,lower(station(ista).axis(3:4)),station(ista).lim21,station(ista).lim22);
    clear axis1;
    clear axis2;
    fprintf(fid, '  enddef;\n');
end
fprintf(fid, '*-----------------------   end $ANTENNA            ----------------------*\n');
%
% $BBC
fprintf(fid, '*----------------------- begin $BBC                ----------------------*\n');
fprintf(fid, '$BBC;\n');
for ibbc = 1 : bbcnum
    fprintf(fid, '  def %s%02d;\n', MODENAME,ibbc); 
    nc = length(bbc(ibbc).bn)/CodesNum;
    for i = 1 : CodesNum
        if (nc == 2)
            fprintf(fid, '    BBC_assign = &BBC%02d : %02d : &IF_%s;\n', i,i,char(bbc(ibbc).bn{i}));
        elseif (nc == 1)
            fprintf(fid, '    BBC_assign = &BBC%02d : %02d : &IF_%s;\n', i,i,char(bbc(ibbc).bn{i}));
        end
    end
    fprintf(fid, '  enddef;\n');
end
fprintf(fid, '*-----------------------   end $BBC                ----------------------*\n');
%
% $DAS
fprintf(fid, '*----------------------- begin $DAS                ----------------------*\n');
fprintf(fid, '$DAS;\n');
fprintf(fid, '  def %s_rack;\n', station(1).rack);
fprintf(fid, '    electronics_rack_type = %s;\n', station(1).rack);
fprintf(fid, '  enddef;\n');
racknum = 1;
rackname(1).name = station(1).rack;
for ista = 2 : stanum
    ifnew = 1;
    for i = 1 : racknum
        if strcmp(station(ista).rack, rackname(i).name)
            ifnew = 0;
            break;
        end
    end
    if (ifnew == 1)
        fprintf(fid, '  def %s_rack;\n', station(ista).rack);
        fprintf(fid, '    electronics_rack_type = %s;\n', station(ista).rack);
        fprintf(fid, '  enddef;\n');
        racknum = racknum + 1;
        rackname(racknum).name = '';
        rackname(racknum).name = strcat(rackname(racknum).name,station(ista).rack);       
    end
end
fprintf(fid, '  def %s_recorder;\n', station(1).equip);
fprintf(fid, '    record_transport_type = %s;\n', station(1).equip);
fprintf(fid, '  enddef;\n');
recordernum = 1;
recordername(1).name = station(1).equip;
for ista = 2 : stanum
    ifnew = 1;
    for i = 1 : recordernum
        if strcmp(station(ista).equip, recordername(i).name)
            ifnew = 0;
        end
    end
    if (ifnew == 1)
        fprintf(fid, '  def %s_recorder;\n', station(ista).equip);
        fprintf(fid, '    record_transport_type = %s;\n', station(ista).equip);
        fprintf(fid, '  enddef;\n');
        recordernum = recordernum + 1;
        recordername(recordernum).name = '';
        recordername(recordernum).name = strcat(recordername(recordernum).name,station(ista).equip);
    end
end
for ista = 1 : stanum
    fprintf(fid, '  def %s_%s;\n', station(ista).po, deblank(station(ista).eq));
    fprintf(fid, '    recording_system_ID = %s;\n', deblank(station(ista).eq));
    fprintf(fid, '  enddef;\n');
end
fprintf(fid, '  def 1_recorder;\n');
fprintf(fid, '    number_drives = 1;\n');
fprintf(fid, '  enddef;\n');
fprintf(fid, '  def 2_recorder;\n');
fprintf(fid, '    number_drives = 2;\n');
fprintf(fid, '  enddef;\n');
fprintf(fid, '  def low_density;\n');
fprintf(fid, '    record_density = 33333 bpi;\n');
fprintf(fid, '  enddef;\n');
fprintf(fid, '  def high_density;\n');
fprintf(fid, '    record_density = 56250 bpi;\n');
fprintf(fid, '  enddef;\n');
fprintf(fid, '  def thick_tape;\n');
fprintf(fid, '    tape_length = 8800 ft;\n');
fprintf(fid, '  enddef;\n');
fprintf(fid, '  def thin_tape;\n');
fprintf(fid, '    tape_length = 17400 ft;\n');
fprintf(fid, '  enddef;\n');
fprintf(fid, '*-----------------------   end $DAS                ----------------------*\n');
%
% $FREQ
fprintf(fid, '*----------------------- begin $FREQ               ----------------------*\n');
fprintf(fid, '$FREQ;\n');
fprintf(fid, '  def %s01;\n', MODENAME);
 for ifreq=1:length(skdcodes(1).lin)
      readlin=skdcodes(1).lin(ifreq,:);
      tmps1 = char(readlin{3});
      tmps2 = char(readlin{4});
      tmps3 = char(readlin{8});
      tmps4 = char(readlin{6});
      ul= skdcodes(1).channel(ifreq).ul;
      fprintf(fid, '    chan_def = &%s : %s MHz : %s : %6.3f MHz : &CH%02d : &BBC%02d : &U_cal;\n', tmps1,tmps2,ul,str2num(tmps3),ifreq,str2num(tmps4));
      clear tmps1 tmps2 tmps4 tmps3;
end
frewind(fidskd);
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if (linelength>=5) & (strcmp(line(1:5),'R SX ')) 
        [tmpr] = sscanf(line(6:linelength),'%f',1);
        fprintf(fid, '    sample_rate = %4.1f Ms/sec;\n', tmpr);
        clear tmpr;
        break;
    end
end
fprintf(fid, '  enddef;\n');
fprintf(fid, '*-----------------------   end $FREQ               ----------------------*\n');
%
% $HEAD_POS
fprintf(fid, '*----------------------- begin $HEAD_POS           ----------------------*\n');
fprintf(fid, '$HEAD_POS;\n');
for ih = 1 : headnum
    fprintf(fid, '  def Mk341-SX%02dS%02d;\n', ih,(2*ih-1));
    [x1] = find(head(ih).bn=='(');
    [x2] = find(head(ih).bn==')');
    for j = 1 : length(x1)
        fprintf(fid, '    headstack_pos = %2d : %5s um;\n', j, head(ih).bn(x1(j)+1:x2(j)-1));
    end
    fprintf(fid, '  enddef;\n');
    clear x1;
    clear x2;
end
fprintf(fid, '*-----------------------   end $HEAD_POS           ----------------------*\n');
%
% $IF
fprintf(fid, '*----------------------- begin $IF                 ----------------------*\n');
fprintf(fid, '$IF;\n');
for iif = 1 : ifrenum
    fprintf(fid, '  def %s%02d;\n', MODENAME,iif);
     nnc = ' ';
    for j = 1 : CodesNum
        if ~strcmp(char(ifre(iif).bbc{j}), nnc)
            nnc = char(ifre(iif).bbc{j});
            fprintf(fid, '    if_def = &IF_%s : %s : R : %7.1f MHz : U : 1 MHz : 0 Hz;\n', nnc,nnc,str2num(char(ifre(iif).bn{j})));
        end
    end 
    fprintf(fid, '  enddef;\n');
end
fprintf(fid, '*-----------------------   end $IF                 ----------------------*\n');
%
% $PHASE_CAL_DETECT
fprintf(fid, '*----------------------- begin $PHASE_CAL_DETECT   ----------------------*\n');
fprintf(fid, '$PHASE_CAL_DETECT;\n');
fprintf(fid, '  def Standard;\n');
fprintf(fid, '    phase_cal_detect = &U_cal : 1;\n');
fprintf(fid, '  enddef;\n');
fprintf(fid, '*-----------------------   end $PHASE_CAL_DETECT   ----------------------*\n');
%
% $PASS_ORDER
fprintf(fid, '*----------------------- begin $PASS_ORDER         ----------------------*\n');
fprintf(fid, '$PASS_ORDER;\n');
for ih = 1 : headnum
    fprintf(fid, '  def Mk341-SX%02dS%02d;\n', ih,(2*ih-1));
    if strfind(head(ih).bn,'E1')>0
        fprintf(fid, '    pass_order =   1A :   2A :   3A :   4A :   5A :   6A :   7A :   8A :   9A :  10A :  11A :  12A :  13A :  14A;\n');
    else
        fprintf(fid, '    pass_order =   1A :   2A :   3A :   4A :   5A :   6A :   7A :   8A :   9A :  10A :  11A :  12A;\n');
    end
    fprintf(fid, '  enddef;\n');
end
fprintf(fid, '*-----------------------   end $PASS_ORDER         ----------------------*\n');
%
% $ROLL
fprintf(fid, '*----------------------- begin $ROLL               ----------------------*\n');
fprintf(fid, '$ROLL;\n');
fprintf(fid, '  def NO_ROLL;\n');
fprintf(fid, '    roll = off;\n');
fprintf(fid, '  enddef;\n');
fprintf(fid, '*-----------------------   end $ROLL               ----------------------*\n');
%
% $SCHED
fprintf(fid, '*----------------------- begin $SCHED              ----------------------*\n');
fprintf(fid, '$SCHED;\n');
%
frewind(fidskd);
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if (linelength>=5) & (strcmp(line(1:5),'$SKED')) 
        break;
    end
end
nsched = 0;
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if (line(1)=='$')
        break;
    end
    nsched = nsched + 1;
    scanstart(nsched,1:9) = line(24:32);  
    scanstart(nsched,10) = ' ';
end
for i = 1 : nsched-1
    if (scanstart(i,10)~=' ')
        continue;
    end
    jsc = nsched;
    for j = i+1 : nsched
        if strcmp(scanstart(i,1:9), scanstart(j,1:9))
            scanstart(i,10) = 'a';
            scanstart(j,10) = 'b';
            jsc = j;
            break;
        end
    end
    jsd = nsched;
    for j = jsc+1 : nsched
        if strcmp(scanstart(i,1:9), scanstart(j,1:9))
            scanstart(j,10) = 'c';
            jsd = j;
            break;
        end
    end
    for j = jsd+1 : nsched
        if strcmp(scanstart(i,1:9), scanstart(j,1:9))
            scanstart(j,10) = 'd';
            break;
        end
    end
end
%
frewind(fidskd);
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if (linelength>=5) & (strcmp(line(1:5),'$SKED')) 
        break;
    end
end
dandt0=' ';
nsched = 0;
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if (line(1)=='$')
        break;
    end
    nsched = nsched + 1;
    srcname = line(1:8);
    dandt = line(24:34);
    [stasn, count, errmsg, nextindex] = sscanf(line(65:linelength),'%s',1);
    nsta = length(stasn)/2;
    index = 65 + nextindex - 1;
    [tmps, count, errmsg, nextindex] = sscanf(line(index:linelength), '%s', (nsta+1));
    index = index + nextindex - 1;
    [tmpd, count, errmsg, nextindex] = sscanf(line(index:linelength), '%d', nsta);
    fprintf(fid, '  scan %s;\n', deblank(strcat(dandt(3:5),'-',dandt(6:7),dandt(8:9),scanstart(nsched,10))));
    fprintf(fid, '    start = %4dy%sd%sh%sm%ss;\n', (2000+str2num(dandt(1:2))),dandt(3:5),dandt(6:7),dandt(8:9),dandt(10:11));
    fprintf(fid, '    mode = GEOSX.SX;\n');
    fprintf(fid, '    source = %s;\n', deblank(srcname));
    for i = 1 : stanum
        for ista = 1 : nsta
            staid = stasn((ista-1)*2+1);
            if(station(i).id == staid)
                po = station(i).po;
                cable = stasn((ista-1)*2+2);
                if (cable == '-')
                    ca = 'n';
                elseif (cable == 'C')
                    ca = 'cw';
                elseif (cable == 'W')
                    ca = 'ccw';
                end
                fprintf(fid, '    station = %s :    0 sec : %5d sec :     0 ft : 1A : &%s : 1;\n', po,tmpd(ista),ca);
                break;
            end
        end  
    end
    fprintf(fid, '  endscan;\n');    
end
fprintf(fid, '*-----------------------   end $SCHED              ----------------------*\n');
%
% $SITES
fprintf(fid, '*----------------------- begin $SITES              ----------------------*\n');
fprintf(fid, '$SITE;\n');
for ista = 1 : stanum
    fprintf(fid, '  def %s;\n', deblank(station(ista).name));
    fprintf(fid, '    site_type = fixed;\n');
    fprintf(fid, '    site_name = %s;\n', deblank(station(ista).name));
    fprintf(fid, '    site_ID = %s;\n', station(ista).po);
    fprintf(fid, '    site_position = %12.3f m : %12.3f m : %12.3f m;\n', station(ista).xyz(1),station(ista).xyz(2),station(ista).xyz(3));
    if (station(ista).hmasknum > 0)
        fprintf(fid, '    horizon_map_az = ');
        fprintf(fid, '%5.1f deg', station(ista).hmask(1));
        for i = 3 : station(ista).hmasknum
            if (mod(i,2)==1)
                fprintf(fid, ' : %5.1f', station(ista).hmask(i));
            end
        end
        fprintf(fid, ';\n');
        fprintf(fid, '    horizon_map_el = ');
        fprintf(fid, '%5.1f deg', station(ista).hmask(2));
        for i = 3 : station(ista).hmasknum
            if (mod(i,2)==0)
                fprintf(fid, ' : %5.1f', station(ista).hmask(i));
            end
        end
        fprintf(fid, ';\n');
    end
    fprintf(fid, '    occupation_code = %08d;\n', round(station(ista).occ));
    fprintf(fid, '  enddef;\n');
end
fprintf(fid, '*-----------------------   end $SITES              ----------------------*\n');
%
% $SOURCE
fprintf(fid, '*----------------------- begin $SOURCE             ----------------------*\n');
fprintf(fid, '$SOURCE;\n');
frewind(fidskd);
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if (linelength>=8) & (strcmp(line(1:8),'$SOURCES')) 
        break;
    end
end
while ~feof(fidskd)
    line = fgetl(fidskd);
    linelength = length(line);
    if (line(1)=='$')
        break;
    end
    if (line(11)=='$')
        sourcename = line(2:9);
    else
        sourcename = line(11:18);
    end
    tmpr = sscanf(line(19:linelength),'%f',6);
    fprintf(fid, '  def %s;\n', deblank(sourcename));
    fprintf(fid, '    source_type = star;\n');
    fprintf(fid, '    source_name = %s;\n', deblank(sourcename));
    fprintf(fid, '    IAU_name = %s;\n', line(2:9));
    fprintf(fid, '    ra = %02dh%02dm%08.5fs;\n', tmpr(1),tmpr(2),roundn(tmpr(3),-5));
    if (length(strfind(line(19:linelength),'-'))>0)
        fprintf(fid, '    dec = -%02dd%02d''%07.4f";\n', abs(tmpr(4)),tmpr(5),roundn(tmpr(6),-4));
    else
        fprintf(fid, '    dec = %02dd%02d''%07.4f";\n', tmpr(4),tmpr(5),roundn(tmpr(6),-4));
    end
    fprintf(fid, '    ref_coord_frame = J2000;\n');
    fprintf(fid, '  enddef;\n');   
end
fprintf(fid, '*-----------------------   end $SOURCE             ----------------------*\n');
%
% $TRACKS
fprintf(fid, '*----------------------- begin $TRACKS             ----------------------*\n');
fprintf(fid, '$TRACKS;\n');
for it=1:length(skdcodes)
    defnam=char(skdcodes(it).lin{1,7});
    fprintf(fid, '  def %s_2f_1b-SX0%d;\n',defnam(1:end-2),it);
    for ich=1:length(skdcodes(it).lin)
        pass=char(skdcodes(it).lin{ich,9});
        if strcmp(pass(1:3),'101')
            npass=2;
        else
            npass=1;
        end 
         o=textscan(pass,'%d( %d, %d, %d, %d)');
         upper=[o{2},o{4}]; lower=[o{3},o{5}];
         if isempty(o{3})
             o=textscan(pass,'%d( %d,, %d)'); 
             upper=[o{2},o{3}]; 
         end
        if  strcmp(defnam(end),'1') % 1:1 fanout
            if strcmp(skdcodes(it).channel(ich).ul,'U')
                fprintf(fid, '    fanout_def =   : &CH%02d : sign : %d : %02d\n',ich,npass,upper(1)+3);
                fprintf(fid, '    fanout_def =   : &CH%02d : mag  : %d : %02d\n',ich,npass,upper(2)+3);
            else
                fprintf(fid, '    fanout_def =   : &CH%02d : sign : %d : %02d\n',ich,npass,lower(1)+3);
                fprintf(fid, '    fanout_def =   : &CH%02d : mag  : %d : %02d\n',ich,npass,lower(2)+3);
            end     
        else
            if strcmp(skdcodes(it).channel(ich).ul,'U')
                fprintf(fid, '    fanout_def =   : &CH%02d : sign : %d : %02d : %02d;\n',ich,npass,upper(1)+3,upper(1)+5);
                fprintf(fid, '    fanout_def =   : &CH%02d : mag  : %d : %02d : %02d;\n',ich,npass,upper(2)+3,upper(2)+5);
            else
                fprintf(fid, '    fanout_def =   : &CH%02d : sign : %d : %02d : %02d;\n',ich,npass,lower(1)+3,lower(1)+5);
                fprintf(fid, '    fanout_def =   : &CH%02d : mag  : %d : %02d : %02d;\n',ich,npass,lower(2)+3,lower(2)+5);
            end     
        end
    end
        fprintf(fid, '  enddef;\n');
end
        fprintf(fid, '  def Mark4_format;\n');
        fprintf(fid, '    track_frame_format = Mark4;\n');
        fprintf(fid, '  enddef;\n');
        fprintf(fid, '  def Mark5B_format;\n');
        fprintf(fid, '    track_frame_format = Mark5B;\n');
        fprintf(fid, '  enddef;\n');
 fprintf(fid, '*-----------------------   end $TRACKS             ----------------------*\n');
% close output file
fclose(fidskd);
fclose(fid);


