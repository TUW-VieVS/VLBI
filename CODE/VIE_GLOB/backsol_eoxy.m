% works only if the estimation intervals for all 5 eop are identical

function backsol_eoxy(pathGS,path_outglob,DIROUT,DIRIN)


pthL2 = pathGS.path_in; % ../DATA/LEVEL2/

curDate=clock;          % current date and time
fidOffic=fopen([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/eoxy_' DIRIN '.txt'],'wt');

fprintf(fidOffic,'# Created on %02.0f.%02.0f.%04.0f at %02.0f:%02.0f:%02.0f local time\n', curDate(3), curDate(2), curDate(1), curDate(4), curDate(5), curDate(6));
fprintf(fidOffic,'# Earth orientation parameters from the VLBI global solution (session-wise reduced parameters) \n');
fprintf(fidOffic,'# VLBI VieVS solution                                                                          \n');
fprintf(fidOffic,'# Analysis center: VIE (TU Wien, Austria)                                                      \n');
fprintf(fidOffic,'#                                                                                              \n');
fprintf(fidOffic,'# Time argument: UTC                                                                           \n');
fprintf(fidOffic,'# Nutation angles are wrt IAU 2006 nutation/precession                                         \n');
fprintf(fidOffic,'#                                                                                              \n');
fprintf(fidOffic,'# EOXY file contains the series of estimates of the Earth orientation                          \n');
fprintf(fidOffic,'# parameters obtained from processing VLBI experiments.                                        \n');
fprintf(fidOffic,'#                                                                                              \n');
fprintf(fidOffic,'# Format:                                                                                      \n');
fprintf(fidOffic,'#                                                                                              \n');
fprintf(fidOffic,'# The line which starts from # is considered as a comment                                      \n');
fprintf(fidOffic,'#                                                                                              \n');
fprintf(fidOffic,'# Field Columns Format Units   Meaning                                                         \n');
fprintf(fidOffic,'#                                                                                              \n');
fprintf(fidOffic,'#   1   2-13    F12.6  days    Modified Julian date of the UTC time tag                        \n');
fprintf(fidOffic,'#   2   15-23   F10.6  arcsec  The estimate of X pole coordinate                               \n');
fprintf(fidOffic,'#   3   25-33   F10.6  arcsec  The estimate of Y pole coordinate                               \n');
fprintf(fidOffic,'#   4   35-44   F10.7  sec     The UT1-UTC function                                            \n');
fprintf(fidOffic,'#   5   46-54   F10.6  mas     Celestial Pole Offset dX with                                   \n');
fprintf(fidOffic,'#                              respect to IAU 2006 nutation/precession              \n');
fprintf(fidOffic,'#   6   56-64   F10.6  mas     Celestial Pole Offset dY with                                   \n');
fprintf(fidOffic,'#                              respect to IAU 2006 nutation/precession              \n');
fprintf(fidOffic,'#   7   66-74   F9.6   arcsec  Formal uncertainty of X pole coordinate                         \n');
fprintf(fidOffic,'#   8   76-84   F9.6   arcsec  Formal uncertainty of Y pole coordinate                         \n');
fprintf(fidOffic,'#   9   86-95   F10.7  sec     Formal uncertainty of UT1-UTC function                          \n');
fprintf(fidOffic,'#  10   97-105  F9.6   mas     Formal uncertainty of nutation dX                               \n');
fprintf(fidOffic,'#  11   107-115 F9.6   mas     Formal uncertainty of nutation dY                               \n');
fprintf(fidOffic,'#  12   117-123 F7.2   psec    Weighted root mean square of postfit residual of                \n');
fprintf(fidOffic,'#                              the solution                                           \n');
fprintf(fidOffic,'#  13   125-130 F6.4   --      Correlation between the estimates of X-pole                     \n');
fprintf(fidOffic,'#                              positions and Y-pole position                       \n');
fprintf(fidOffic,'#  14   132-137 F6.4   --      Correlation between the estimates of X-pole                     \n');
fprintf(fidOffic,'#                              positions and UT1-TAI angle                         \n');
fprintf(fidOffic,'#  15   139-144 F6.4   --      Correlation between the estimates of Y-pole                     \n');
fprintf(fidOffic,'#                              positions and UT1-TAI angle                         \n');
fprintf(fidOffic,'#  16   146-151 F6.4   --      Correlation between the estimates of nutation                   \n');
fprintf(fidOffic,'#                              dX and nutation dY                                  \n');
fprintf(fidOffic,'#  17   153-158 I6     --      Number of used observations in the session                      \n');
fprintf(fidOffic,'#  18   160-165 A6     --      IVS session code                                                \n');
fprintf(fidOffic,'#  19   167-171 F5.2   hours   Session duration                                                \n');
fprintf(fidOffic,'#  20   173-181 F9.6   asc/day Rate of X pole coordinate \n');
fprintf(fidOffic,'#  21   183-191 F9.6   asc/day Rate of Y pole coordinate \n');
fprintf(fidOffic,'#  22   193-202 F10.7  sec     Length of day \n');
fprintf(fidOffic,'#  23   204-212 F9.6   mas/day Rate of dX \n');
fprintf(fidOffic,'#  24   214-222 F9.6   mas/day Rate of dY \n');
fprintf(fidOffic,'#  25   224-232 F9.6   asc/day Formal uncertainty of X pole coordinate rate \n');
fprintf(fidOffic,'#  26   234-242 F9.6   asc/day Formal uncertainty of Y pole coordinate rate \n');
fprintf(fidOffic,'#  27   244-253 F10.7  sec     Formal uncertainty of length of day \n');
fprintf(fidOffic,'#  28   255-256 F9.6   mas/day Formal uncertainty of dX rate \n');
fprintf(fidOffic,'#  29   258-259 F9.6   mas/day Formal uncertainty of dY rate \n');
fprintf(fidOffic,'#  30   261-324 A64    --      Network ID. The alphabetically ordered sequence \n');
fprintf(fidOffic,'#                              of two-letter IVS station identifiers. Only those \n');
fprintf(fidOffic,'#                              stations which provided observations used in the \n');
fprintf(fidOffic,'#                              solution are listed. The station names are defined \n');
fprintf(fidOffic,'#                              in the IVS document ivscontrol/ns-codes.txt \n');
fprintf(fidOffic,'# \n');
fprintf(fidOffic,'# If a given parameter was not estimated a filler, -0, is placed. \n');
fprintf(fidOffic,'# \n');



writeFormat=' %12.6f%10.6f%10.6f%11.7f%10.6f%10.6f%10.6f%10.6f%11.7f%10.6f%10.6f  %6.0f%7.0f%7.0f%7.0f%7.0f %6.0f %6s %5.2f    %6.0f%10.0f%11.0f%10.0f%10.0f%10.0f%10.0f%11.0f%3.0f%3.0f %s\n';



%%
EOPs = ['xpol'; 'ypol';'dut1'; 'dX  '; 'dY  '];
formatEOP = '%s %f %f %f %s %f %f';
[eop.param, eop.mjd, eop.val, eop.mx, eop.ses, eop.sesstart, eop.sesend]=textread([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/eop_' DIRIN '.txt'], formatEOP, 'commentstyle','matlab');
%val = mjd, val, std



eop.d1= eop.mjd- (eop.sesstart - 0.5); % within 12h from midnight
eop.d2= (eop.sesend + 0.5) - eop.mjd;

eop.use = logical (eop.d1 > 0 & eop.d2 >0);

eop.param(~eop.use)=[];
eop.mjd(~eop.use)=[];
eop.val(~eop.use)=[];
eop.mx(~eop.use)=[];
eop.ses(~eop.use)=[];
eop.sesstart(~eop.use)=[];
eop.sesend(~eop.use)=[];
eop.d1(~eop.use)=[];
eop.d2(~eop.use)=[];
eop.use(~eop.use)=[];

eop.dXdYapr = zeros(length(eop.mjd),1);


[sesuni, idsesuni] = unique(eop.ses);
mjdsesuni = eop.sesstart(idsesuni);


% get session code from master.txt

for i=1:size(sesuni,1)
    yrsall(i,:) = str2num(sesuni{i}(1:2));
end

yrs = unique(yrsall);


for i = 1:length(yrs)
    [sesname_y, sescode_y] = read_master(yrs(i));
    
    if i == 1
        sesname = sesname_y;
        sescode = sescode_y;
    else 
        sesname = [sesname; sesname_y];
        sescode = [sescode; sescode_y];
    end
    
end



idxpol = strcmp(eop.param,'xpol');
idypol = strcmp(eop.param,'ypol');
iddut1 = strcmp(eop.param,'dut1');
iddX = strcmp(eop.param,'dX');
iddY = strcmp(eop.param,'dY');


for i = 1:size(sesuni,1)
        
    load([pthL2 DIRIN '/' sesuni{i} '_par_glob.mat']);
   

       
    if i==1
        load(glob2.opt.trf{1})
    end
    
    
    
    
    idses = strcmp(eop.ses, sesuni{i});
    
    idxpolses = idxpol&idses;    
    mjdi = eop.mjd(idxpolses);
    idp=[];
    for j = 1:length(mjdi)
        idp = [idp; find(glob2.eopapr.mjd == mjdi(j))];
    end
    apr=[];
    apr = glob2.eopapr.xp(idp);   
    eop.apr(idxpolses,1) =apr; %as
    
    
    
    idypolses = idypol&idses;    
    mjdi = eop.mjd(idypolses);
    idp=[];
    for j = 1:length(mjdi)
        idp = [idp; find(glob2.eopapr.mjd == mjdi(j))];
    end
    apr=[];
    apr = glob2.eopapr.yp(idp);   
    eop.apr(idypolses,1) =apr; %as
    
   
    
    iddut1ses = iddut1&idses;    
    mjdi = eop.mjd(iddut1ses);
    idp=[];
    for j = 1:length(mjdi)
        idp = [idp; find(glob2.eopapr.mjd == mjdi(j))];
    end
    apr=[];
    apr = glob2.eopapr.ut1(idp);   
    eop.apr(iddut1ses,1) =apr; %s
   
    
    
    iddXses = iddX&idses;    
    mjdi = eop.mjd(iddXses);
    idp=[];
    for j = 1:length(mjdi)
        idp = [idp; find(glob2.eopapr.mjd == mjdi(j))];
    end
    apr=[];
    apr = glob2.eopapr.dX(idp);   
    eop.apr(iddXses,1) =apr; %as

    
    iddYses = iddY&idses;    
    mjdi = eop.mjd(iddYses);
    idp=[];
    for j = 1:length(mjdi)
        idp = [idp; find(glob2.eopapr.mjd == mjdi(j))];
    end
    apr=[];
    apr = glob2.eopapr.dY(idp);   
    eop.apr(iddYses,1) =apr; %as
    
    
    % stations in the session
    ant = {glob2.opt.stat.name};
    ancode=[''];
    for j=1:length(ant)
       idst=strcmp({superstations.name},ant(j)) ;
       ancode(j,:) = superstations(idst).code(2:3);
        
    end
     ancods =sort(ancode);
     eop.ancode(idxpolses,1) = {ancods}; % save only by x pole
     eop.numobs(idxpolses,1) = glob2.opt.total_obs;
     
     idsescode = strcmp(sesname,sesuni{i}(1:9));
     if sum(idsescode)==1
        eop.sescode(idxpolses,1) = sescode(idsescode);
     else
        eop.sescode(idxpolses,1) = {'------'}; 
     end
     

end




idxpol = strcmp(eop.param,'xpol');
idypol = strcmp(eop.param,'ypol');
iddut1 = strcmp(eop.param,'dut1');
iddX = strcmp(eop.param,'dX');
iddY = strcmp(eop.param,'dY');


[mjdsessort, idsessort]= sort(mjdsesuni);

for i = 1:length(idsessort) % session-wise
    sesa = sesuni(idsessort(i),:);
    idses = [];
    idses = strcmp(eop.ses, sesa);  
    
    idxpolses = idxpol&idses; % xp   
    mjdi = eop.mjd(idxpolses);
    
    outeop=[];
    outeop(:,1) = eop.val(idxpolses)./1000 + eop.apr(idxpolses); % as
    outeop(:,6) = eop.mx(idxpolses)./1000;
    
    
    idypolses = idypol&idses; % yp       
    outeop(:,2) = eop.val(idypolses)./1000 + eop.apr(idypolses); % as
    outeop(:,7) = eop.mx(idypolses)./1000;
   
    iddut1ses = iddut1&idses;       
    outeop(:,3) = eop.val(iddut1ses)./1000 + eop.apr(iddut1ses); % s
    outeop(:,8) = eop.mx(iddut1ses)./1000;
    
    iddXses = iddX&idses;       
    outeop(:,4) = eop.val(iddXses) + eop.apr(iddXses).*1000; % mas
    outeop(:,9) = eop.mx(iddXses);
    
    iddYses = iddY&idses;   
    outeop(:,5) = eop.val(iddYses) + eop.apr(iddYses).*1000; % mas
    outeop(:,10) = eop.mx(iddYses);
        
    
    sesdur = (eop.sesend(idxpolses) - eop.sesstart(idxpolses)).*24; %h
    ancode = eop.ancode(idxpolses);
    numobs = eop.numobs(idxpolses);
    
    sescode = eop.sescode{idxpolses};

    
    for j = 1:length(mjdi) % mjds in the session
        an = (join(string(ancode{j})));
        idsp = isspace(an);
        anout = char(an);
        anout(idsp) = [];
        
        
        fprintf(fidOffic,writeFormat, mjdi(j), outeop(j,:), -0, -0, -0, -0, -0, numobs(j), sescode, sesdur(j), -0, -0, -0, -0, -0,    -0, -0, -0, -0, -0,   anout );
        
    end


end
fclose(fidOffic);



