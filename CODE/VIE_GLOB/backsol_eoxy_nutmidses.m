
function backsol_eoxy_nutmidses(pathGS,path_outglob,DIROUT,DIRIN)

% clear all
% load('../DATA/GLOB/pathGS.mat')
% path_outglob = '../OUT/GLOB/';
% DIROUT = pathGS.out;
% DIRIN = pathGS.L2;

pthL2 = pathGS.path_in; % ../DATA/LEVEL2/

curDate=clock;          % current date and time
% writeFormat='%12.6f%10.7f%10.7f%11.8f%7.4f%7.4f%10.7f%10.7f%11.8f%9.4f%9.4f %8.2f %9.0f%10.0f%10.0f%10.0f %6.0f %8s %7.2f %10.0f%11.0f%12.0f%8.0f%8.0f%11.0f%11.0f%12.0f%8.0f%8.0f     %s\n';
writeFormat='%12.6f%11.7f%11.7f%12.8f%8.4f%8.4f%11.7f%11.7f%12.8f%10.4f%10.4f %8.2f %9.0f%10.0f%10.0f%10.0f %6.0f %16s %7.2f %10.0f%11.0f%12.0f%8.0f%8.0f%11.0f%11.0f%12.0f%8.0f%8.0f     %s\n';

formatEOP = '%s %f %f %f %s %f %f';


for ifile = 1:2

    if ifile ==1
        [eop.param, eop.mjd, eop.val, eop.mx, eop.ses, eop.sesstart, eop.sesend]=textread([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/eop_' DIRIN '.txt'], formatEOP, 'commentstyle','matlab');
        
        % write 2 midnights into the .eoxy
        d1= eop.mjd- (eop.sesstart - 0.5); % within 12h from midnight
        d2= (eop.sesend + 0.5) - eop.mjd;
        eop.use = logical (d1 > 0 & d2 >0);
        fidOffic=fopen([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/eoxy_' DIRIN '.txt'],'wt');
    else
        eop =[];
        [eop.param, eop.mjd, eop.val, eop.mx, eop.ses, eop.sesstart, eop.sesend]=textread([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/eop_' DIRIN '.txt'], formatEOP, 'commentstyle','matlab');
        
        % write only 1 midnight which is inside the session, or which is more close to the
        % start or end if the session did not run over 0 UT

        mjdu = round((eop.sesstart+eop.sesend)./2);
        eop.use = logical(eop.mjd==mjdu);
        fidOffic=fopen([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/eoxy_1offset_' DIRIN '.txt'],'wt');
    end

    
    
    eop.param(~eop.use)=[];
    eop.mjd(~eop.use)=[];
    eop.val(~eop.use)=[];
    eop.mx(~eop.use)=[];
    eop.ses(~eop.use)=[];
    eop.sesstart(~eop.use)=[];
    eop.sesend(~eop.use)=[];
    eop.use(~eop.use)=[];

    
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
            headxpolcon = glob2.opt.xpol.coef; %mas
            headypolcon = glob2.opt.ypol.coef; %mas
            headut1con = glob2.opt.dut1.coef; %mas
            headdXcon = glob2.opt.nutdx.coef *1000; %uas
            headdYcon = glob2.opt.nutdx.coef *1000; %uas
            headtrf = glob2.opt.trf{2};
            headcrf = glob2.opt.crf{2};
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
        iddYses = iddY&idses;        
        
        mjd1 = eop.sesstart(iddYses);
        mjd2 = eop.sesend(iddYses);
        MJDmid = mjd1 + (mjd2-mjd1)/2;
        
        
        if contains(glob2.eopapr.interp,'lagrange')
            eop.apr(iddXses,1) = lagint4v(glob2.eopapr.mjd, glob2.eopapr.dX, MJDmid(1));
            eop.apr(iddYses,1) = lagint4v(glob2.eopapr.mjd, glob2.eopapr.dY, MJDmid(1));
            eop.mjd(iddXses,1) = MJDmid;
            eop.mjd(iddYses,1) = MJDmid;        
        elseif contains(glob2.eopapr.interp,'linear')
            eop.apr(iddXses,1) = interp1(glob2.eopapr.mjd, glob2.eopapr.dX, MJDmid(1),'linear','extrap');
            eop.apr(iddYses,1) = interp1(glob2.eopapr.mjd, glob2.eopapr.dY, MJDmid(1),'linear','extrap');
            eop.mjd(iddXses,1) = MJDmid;
            eop.mjd(iddYses,1) = MJDmid;        
    
        end
    
        
        % stations in the session
        ant = {glob2.opt.stat.name};
        ancode=[''];
        ancods='';
        for j=1:length(ant)
           idst=strcmp({superstations.name},ant(j)) ;
           ancode(j,:) = [superstations(idst).code(2:3) '-'];
            
        end
         ancods =sortrows(ancode);
         eop.ancode(idxpolses,1) = {ancods}; % save only by x pole
         eop.numobs(idxpolses,1) = glob2.opt.total_obs;
         
         idsescode = strcmp(sesname,sesuni{i}(1:9));
         if sum(idsescode)==1
            eop.sescode(idxpolses,1) = sescode(idsescode);
         else
            eop.sescode(idxpolses,1) = {glob2.opt.session_name}; 
         end
         
         
         
    
    end
    
    
    [dastart(1), dastart(2), dastart(3)]=mjd2date(min(eop.sesstart));
    [daend(1), daend(2), daend(3)]=mjd2date(max(eop.sesend));
    
    %%
    
    fprintf(fidOffic,'%%=IVS-EOP 3.0 VIE %04.0f-%02.0f-%02.0fT%02.0f:%02.0f:%02.0f VIE %04.0f-%02.0f-%02.0fT%02.0f:%02.0f:%02.0f %04.0f-%02.0f-%02.0fT%02.0f:%02.0f:%02.0f  UTC R\n',curDate(1), curDate(2), curDate(3), curDate(4), curDate(5), curDate(6), dastart(1), dastart(2), dastart(3),0,0,0, daend(1), daend(2), daend(3),0,0,0);
    fprintf(fidOffic,'+HEADER\n');
    fprintf(fidOffic,'GENERATION_TIME		%04.0f-%02.0f-%02.0fT%02.0f:%02.0f:%02.0f\n', curDate(1), curDate(2), curDate(3), curDate(4), curDate(5), curDate(6));
    fprintf(fidOffic,'DATA_START		%04.0f-%02.0f-%02.0fT%02.0f:%02.0f:%02.0f\n',dastart(1), dastart(2), dastart(3),0,0,0);
    fprintf(fidOffic,'DATA_END		%04.0f-%02.0f-%02.0fT%02.0f:%02.0f:%02.0f\n',daend(1), daend(2), daend(3),0,0,0);
    fprintf(fidOffic,'DESCRIPTION		Earth orientation parameters from the VLBI global solution (session-wise reduced parameters) \n');
    fprintf(fidOffic,'ANALYSIS_CENTER		VIE (TU Wien, Austria) \n');
    fprintf(fidOffic,'CONTACT			VIE Analysis Center (vlbi_contact@tuwien.ac.at) \n');
    fprintf(fidOffic,'SOFTWARE		VieVS v3.2 \n');
    fprintf(fidOffic,'TECHNIQUE		VLBI \n');
    fprintf(fidOffic,'NUTATION_TYPE		CIO-BASED  \n');
    fprintf(fidOffic,'ROTATION_TYPE		UT1-UTC_LOD \n');
    fprintf(fidOffic,'TRF_APRIORI		%s \n', headtrf);
    fprintf(fidOffic,'CRF_APRIORI		%s \n', headcrf);
    fprintf(fidOffic,'EOP_SUBDAILY		DESAI-SIBOIS \n');
    %fprintf(fidOffic,'EOP_APRIORI		finals2000A.all \n');
    fprintf(fidOffic,['EOP_ESTIMATED		XPOL_BSP_1	%6.3f mas \n'], headxpolcon);
    fprintf(fidOffic,['EOP_ESTIMATED		YPOL_BSP_1	%6.3f mas \n'], headypolcon);
    fprintf(fidOffic,['EOP_ESTIMATED		DUT1_BSP_1	%6.3f mas \n'], headut1con);
    fprintf(fidOffic,['EOP_ESTIMATED		DX		%6.3f uas \n'], headdXcon);
    fprintf(fidOffic,['EOP_ESTIMATED		DY		%6.3f uas \n'], headdYcon);
    fprintf(fidOffic,'-HEADER\n');
    fprintf(fidOffic,'+DATA\n');
    fprintf(fidOffic,'# All fields are in free format separated by blanks.\n');
    fprintf(fidOffic,'# If a given parameter was not estimated a filler, NaN, is placed. \n');
    fprintf(fidOffic,'#     1           2        3          4        5      6       7         8         9          10       11       12       13        14        15        16       17     18       19       20         21          22       23      24        25        26          27        28      29          30                                                               31 \n'); 
    fprintf(fidOffic,'#   epoch       xPol      yPol      dUT1      dX     dY     sig_xP    sig_yP    sig_UT     sig_dX   sig_dY    wRMS   cor_xPyP  cor_xPUT  cor_yPUT  cor_dXdY   nObs   sessID   span      xPolR      yPolR       LOD      dXR     dYR    sig_xPR    sig_yPR     sig_LOD   sig_dXR sig_dYR     network                                                         comments \n'); 
    fprintf(fidOffic,'#   [MJD]       [as]      [as]       [s]     [mas]  [mas]    [as]      [as]      [s]        [mas]    [mas]    [ps]     [-]       [-]       [-]       [-]      [-]     [-]      [h]     [as/d]     [as/d]       [s]    [mas/d] [mas/d]   [as/d]     [as/d]       [s]     [mas/d] [mas/d]       [-]                                                              [-] \n'); 
    
    %%
    
    
    
    
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
        
        sesWRMS = NaN;
        MJDmid = eop.sesstart(idxpolses) + sesdur/2/24;
    
        outeop_nutNaN = outeop;
        outeop_nutNaN(:,4) = NaN;
        outeop_nutNaN(:,5) = NaN;
        outeop_nutNaN(:,9) = NaN;
        outeop_nutNaN(:,10) = NaN;
         
         
         
        for j = 1:length(mjdi) % mjds in the session
            an = (join(string(ancode{j})));
            idsp = isspace(an);
            anout = char(an);
            anout(idsp) = [];
            anout(end) = [];
            
            fprintf(fidOffic,writeFormat, mjdi(j), outeop_nutNaN(j,:), sesWRMS, NaN, NaN, NaN, NaN, numobs(j), sescode, sesdur(j),  NaN, NaN, NaN, NaN,  NaN, NaN, NaN, NaN, NaN, NaN,   anout );
            
        end   
        mjdi = eop.mjd(iddXses);
        %dXdY - middle of the session
        fprintf(fidOffic,writeFormat, mjdi(1), NaN, NaN, NaN, outeop(1,4:5), NaN, NaN, NaN, outeop(1,9:10), sesWRMS, NaN, NaN, NaN, NaN, numobs(1), sescode, sesdur(1),  NaN, NaN, NaN, NaN,  NaN, NaN, NaN, NaN, NaN, NaN,   anout );
            
    
    end
    
    
    fprintf(fidOffic,'-DATA\n');
    fprintf(fidOffic,'%%IVS-EOP 3.0 END\n');
    fclose(fidOffic);

end



