% Function for the back solution of the sessionwise reduced parameters,
% can be run after vie_glob

% The resulting .txt files will be saved in /BACKWARD_SOLUTION/DIROUT/



% Coded for VieVS by Hana Spicakova
% 1 Aug 2011

% Revised on 
% 11 Feb 2013 by Hana Krásná: compability with version '2.1'
% 24 Jul 2013 by Hana Krásná: changed output for sources
% 16 Jan 2018 by Hana Krásná: relative paths changed to starting directory WORK, script changed to a function; backward_solution(DIRIN, DIROUT), DIRIN and DIROUT according vie_glob
%%
 
function backward_solution(DIRIN, DIROUT)

% Which parameters do you want to estimate?
solbck.eop=1;
solbck.zwd=1;
solbck.tgr=1;
solbck.ant=1;
solbck.sou=1;



%%
path_outglob = '../OUT/GLOB/';
path_Nb = [path_outglob 'BACKWARD_DATA/' DIROUT '/'];

load([path_outglob '_ESTIMATES/' DIROUT '/globsol_' DIRIN '.mat'])
load([path_Nb 'parsplit_' DIRIN])  %parsplit
load([path_Nb 'paths_' DIRIN])  %paths
   

lse = size(globsol.sessions,2);
ncond=globsol.nr_of_conditions;
Q11 = globsol.Q(1:end-ncond,1:end-ncond);
x1 = globsol.x;


for ise= lse:-1: 1
    % counter
    if mod(ise,100)==0
        fprintf(' processing session %4.0f \n',ise)
    end 
    load([path_Nb 'b2_' DIRIN '_' num2str(ise)])  %b2
    load([path_Nb 'N22_' DIRIN '_' num2str(ise)]) %N22
    load([path_Nb 'N21_' DIRIN '_' num2str(ise)]) %N21
    
    redpos = parsplit(ise).redpos;
    redcoor = parsplit(ise).redcoor;
    redsou = parsplit(ise).redsou;
    
    %% constraints   
        
    if ~isempty(redcoor) 
        
        [xx idco]=intersect(redpos,redcoor);
        clear xx

        nco=length(redcoor);

        H=zeros(nco,length(redpos));
        for i=1:nco
           H(i,idco(i))=1;
        end

        cons(1:nco) = 1/4; % weight     
        Ph=diag(cons);    
        HTPH =  H'*Ph*H;    
        N22 = N22 + HTPH;

        w(1:nco)=0; %0 mas - ideal walue       
        b2(idco) = b2(idco) + Ph*w';
        
        clear Ph w H HTPH cons
    end
        
        %%
    if ~isempty(redsou) 
        
        [xx idsr]=intersect(redpos,redsou);
        clear xx

        nsr=length(redsou);

        H=zeros(nsr,length(redpos));
        for i=1:nsr
           H(i,idsr(i))=1;
        end

        cons(1:nsr) = 1/4; % weight     
        Ph=diag(cons);    
        HTPH =  H'*Ph*H;    
        N22 = N22 + HTPH;

        w(1:nsr)=0; %0 mas - ideal walue
        b2(idsr) = b2(idsr) + Ph*w';

        clear Ph w H HTPH cons

    end
    
    %%
    
    invN22 = inv(N22);

    x2 = invN22*b2 - invN22*N21*x1;
    % Covariance matrix
    Q22 = invN22 + invN22*N21*Q11*N21'*invN22;
    varpar=full(sqrt(diag(Q22)));

    
    bckdsol(ise).x(:,1)=x2; %values
    bckdsol(ise).x(:,2)=varpar; %standard deviations
    bckdsol(ise).x(:,3)=parsplit(ise).redpos'; %columns in the old orig. N matrix
    bckdsol(ise).Q = Q22;
    
    clear x2 varpar Q22 b2 N22 N21 invN22
end


if exist([path_outglob 'BACKWARD_SOLUTION/' DIROUT])~=7
   mkdir([path_outglob 'BACKWARD_SOLUTION/' DIROUT]);
end


%  save([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/bckdsol'],'bckdsol');
%%
% load([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/bckdsol'])

if solbck.eop==1
    fidEOP = fopen([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/eop_' DIRIN '.txt'],'wt');
    fprintf(fidEOP,'\n %% parameter   mjd    value    standard deviation\n');
    fprintf(fidEOP,' %% Units are [mas] for xpol, ypol, dX, dY;  [ms] for dut1\n\n');
end
    
if solbck.zwd==1
    fidZWD = fopen([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/zwd_' DIRIN '.txt'],'wt');
    fprintf(fidZWD,'\n %% parameter   station    mjd    value    standard deviation\n');
    fprintf(fidZWD,' %% Units are [cm] \n\n');
    formatZWD = '%c%c%c   %c%c%c%c%c%c%c%c      %10.5f    %8.3f    %8.4f  %s \n';
end    

if solbck.tgr==1
    fidTGR = fopen([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/tgr_' DIRIN '.txt'],'wt');
    fprintf(fidTGR,'\n %% parameter   station    mjd    value    standard deviation\n');
    fprintf(fidTGR,' %% Units are [cm] \n\n');
    formatTGR = '%c%c%c   %c%c%c%c%c%c%c%c      %10.5f    %8.3f    %8.4f  %s \n';
end

if solbck.ant==1
    fidANT = fopen([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/ant_' DIRIN '.txt'],'wt');
    fprintf(fidANT,'\n %% parameter   station    mjd    value    standard deviation\n');
    fprintf(fidANT,' %% Units are [cm] \n\n');
    formatANT = '%c%c%c%c%c   %c%c%c%c%c%c%c%c      %10.5f    %8.3f    %8.4f  %s \n';
end

if solbck.sou==1
    fidSOU = fopen([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/sou_' DIRIN '.txt'],'wt');
    fprintf(fidSOU,'\n %% parameter   source    mjd    value    standard deviation\n');
    fprintf(fidSOU,' %% Units are [mas] \n\n');
    formatSOU = '%c%c%c   %c%c%c%c%c%c%c%c      %10.5f    %8.3f    %8.4f  \n';
    
    
    
    curDate=clock;          % current date and time
    fidSOUCAT = fopen([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/sou_catalogue_' DIRIN '.txt'],'wt');
    fprintf(fidSOUCAT, 'Created on %02.0f.%02.0f.%04.0f at %02.0f:%02.0f:%02.0f local time\n', curDate(3), curDate(2), curDate(1), curDate(4), curDate(5), curDate(6));
    fprintf(fidSOUCAT,'------------------------------------------------------------------------------------------------------------------------------------------------\n');
    fprintf(fidSOUCAT,'IVS name  IERS name  Right Ascension   Declination         Uncertainty       Corr.   Mean    First   Last    Nb    Nb     Nb    Database name \n');
    fprintf(fidSOUCAT,'                        J2000.0         J2000.0           R.A.      Dec.     RA-Dc   MJD     MJD     MJD     sess. obs.   obs.     \n');
    fprintf(fidSOUCAT,'                       h  m  s           o  m  as           s         as              of observation span          delay  rate                 \n');      
    fprintf(fidSOUCAT,'------------------------------------------------------------------------------------------------------------------------------------------------\n');
    formatSOUCAT = ' %c%c%c%c%c%c%c%c %c%c%c%c%c%c%c%c %02s %02s %011.8f %s%02s %02s %010.7f %10.8f %9.7f %7.4f %7.1f %7.1f %7.1f  %4s %6.0f %6s %s \n';
end


for ise= lse:-1:  1
    if mod(ise,100)==0
        fprintf(' processing session %4.0f \n',ise) 
    end 

    load([path_Nb 'parGS_' DIRIN '_' num2str(ise)]) %parGS
    oldcol={parGS.oldcol};
    

    load([pathGS.path_in pathGS.L2 '/' globsol.sessions{ise} '_par_glob.mat']); % compatible with version 2.1
        
    for i=1:size(bckdsol(ise).x,1)
        for j=1:length(oldcol)
            if iscell(oldcol{j})==0
                k=(find((oldcol{j})==bckdsol(ise).x(i,3)));
            end
            if ~isempty(k)
                pname=parGS(j).name;
                
                % Output for EOP
                if solbck.eop==1
                    if strcmp(pname,'xpol') || strcmp(pname,'ypol') || strcmp(pname,'dut1') || strcmp(pname,'dX') || strcmp(pname,'dY')
                        formatEOP = '%c%c%c%c    %10.5f    %8.5f   %8.6f  %s \n';
                        if strcmp(pname,'dX') || strcmp(pname,'dY'); formatEOP = '%c%c      %10.5f    %8.5f   %8.6f  %s \n';end
                        mjd=parGS(j).mjd(k);
                        value=bckdsol(ise).x(i,1);
                        stdev=bckdsol(ise).x(i,2);
                        if strcmp(pname,'dut1'); value = value/15; stdev = stdev/15;  end % mas --> ms for UT
                        fprintf(fidEOP, formatEOP, pname, mjd, value , stdev , globsol.sessions{ise});
                        clear mjd value stdev
                    end
                end
                
                % Output for zwd
                if solbck.zwd==1
                    if strcmp(pname,'zwd')
                        % find name of the antenna
                        oldcolZWD={glob2.x.zwd.col};
                        for jj=1:length(oldcolZWD)
                            kk=(find((oldcolZWD{jj})==bckdsol(ise).x(i,3)));
                            if ~isempty(kk)
                                antenna = glob2.x.antenna(jj).name;
                            end
                        end
                        mjd=parGS(j).mjd(k);
                        value=bckdsol(ise).x(i,1);
                        stdev=bckdsol(ise).x(i,2);
                        fprintf(fidZWD, formatZWD, pname, antenna, mjd, value, stdev, globsol.sessions{ise});
                        clear mjd value stdev antenna oldcolZWD
                    end
                end
                
                % Output for troposheric gradients
                if solbck.tgr==1
                    if strcmp(pname,'ngr') || strcmp(pname,'egr')
                        if strcmp(pname,'ngr')
                            % find name of the antenna
                            oldcolTGR={glob2.x.ngr.col};
                        elseif strcmp(pname,'egr')
                            oldcolTGR={glob2.x.egr.col};
                        end   
                        for jj=1:length(oldcolTGR)
                            kk=(find((oldcolTGR{jj})==bckdsol(ise).x(i,3)));
                            if ~isempty(kk)
                                antenna = glob2.x.antenna(jj).name;
                            end
                        end
                        mjd=parGS(j).mjd(k);
                        value=bckdsol(ise).x(i,1);
                        stdev=bckdsol(ise).x(i,2);
                        fprintf(fidTGR, formatTGR, pname, antenna, mjd, value, stdev, globsol.sessions{ise});
                        clear mjd value stdev antenna oldcolTGR
                    end
                end     
                    
                    
                % Output for antennas
                if solbck.ant==1
                    if strcmp(pname,'ant_x') || strcmp(pname,'ant_y') || strcmp(pname,'ant_z')
                        % find name of the antenna
                        if strcmp(pname,'ant_x')
                            oldcolANT={glob2.x.coorx.col};
                        elseif strcmp(pname,'ant_y')
                            oldcolANT={glob2.x.coory.col};
                        elseif strcmp(pname,'ant_z')
                            oldcolANT={glob2.x.coorz.col};
                        end    
                        for jj=1:length(oldcolANT)
                            kk=(find((oldcolANT{jj})==bckdsol(ise).x(i,3)));
                            if ~isempty(kk)
                                antenna = glob2.x.antenna(jj).name;
                            end
                        end
                        
                        mjd=glob2.opt.midnight; % same for all stations, one estimate per session
                        value=bckdsol(ise).x(i,1);
                        stdev=bckdsol(ise).x(i,2);
                        fprintf(fidANT, formatANT, pname, antenna, mjd, value, stdev, globsol.sessions{ise});
                        clear mjd value stdev antenna oldcolANT
                    end
                end
                
                % Output for sources
                if solbck.sou==1
                   
                    
                    if strcmp(pname,'sra') || strcmp(pname,'sde') 
                        % find name of the source
                        if strcmp(pname, 'sra')
                            oldcolSOU=glob2.x.col_soura;
                        elseif strcmp(pname, 'sde')
                            oldcolSOU=glob2.x.col_soude;
                        end
                        for jj=1:length(oldcolSOU)
                            kk=(find((oldcolSOU(jj))==bckdsol(ise).x(i,3)));
                            if ~isempty(kk)
                                source = glob2.x.source(jj).name;
                            end
                        end

                        mjd=glob2.opt.midnight;
                        value=bckdsol(ise).x(i,1);
                        stdev=bckdsol(ise).x(i,2);
                        fprintf(fidSOU, formatSOU, pname, source, mjd, value, stdev);
                        clear mjd value stdev source oldcolSOU
                    end
                    
                    
                    if strcmp(pname,'sra') || strcmp(pname,'sde') 
                        % find name of the source
                        if strcmp(pname, 'sra')
                            oldcolSOURA=glob2.x.col_soura; % all sources in the session
                            oldcolSOUDe=glob2.x.col_soude;
                            for jj=1:length(oldcolSOURA)
                                kk=(find((oldcolSOURA(jj))==bckdsol(ise).x(i,3)));
                                if ~isempty(kk)
                                    kkRA = oldcolSOURA(jj);
                                    kkDe = oldcolSOUDe(jj);
                                    source_ivs = glob2.x.source(jj).name;
                                    source_iers = glob2.x.source(jj).IERSname;
                                    numobs =  glob2.opt.source(jj).total_obs;
                                    
                                   %the a priori coordinates
                                    RA_apr = glob2.opt.source(jj).ra; %rad
                                    De_apr = glob2.opt.source(jj).de; %rad
                                end
                            end

                            mjd=(glob2.opt.first_scan+glob2.opt.last_scan)/2;

                            ilineRA = find(kkRA == bckdsol(ise).x(:,3)); %(= i)
                            ilineDe = find(kkDe == bckdsol(ise).x(:,3));

                            dRA=bckdsol(ise).x(ilineRA,1); %mas
                            stdevRA=bckdsol(ise).x(ilineRA,2); %mas
                            dDe=bckdsol(ise).x(ilineDe,1); %mas
                            stdevDe=bckdsol(ise).x(ilineDe,2); %mas
                            
                            
                            RA_est=RA_apr + dRA/1000/3600/180*pi; %[rad]
                            RA_est_h = RA_est/pi*12; %[h]
                            RA_h = floor(RA_est_h); %[h]
                            RA_min=floor((RA_est_h-RA_h)*60); %[min]
                            RA_sec=(((RA_est_h-RA_h)*60)-RA_min)*60; %[sec]
                            
                            De_est=De_apr + dDe/1000/3600/180*pi; %[rad]
                            De_est_deg=De_est/pi*180; %[deg]
                            De_deg=fix(De_est_deg); %[deg]
                            De_min1=(De_est_deg-De_deg)*60; %[min]
                            id=find(De_min1<0);
                            De_min1(id)=De_min1(id)*(-1);
                            De_min=floor(De_min1);
                            De_sec=(De_min1-De_min)*60; %[sec]

                            res= [RA_h RA_min RA_sec De_deg De_min De_sec De_est_deg]; % results

                            if res(7)<0
                                sign = '-';
                            else
                                sign = ' ';
                            end
                            
                            errRA = stdevRA/15/1000; %[s]
                            errDe = stdevDe/1000; %[as]
                            
                            
                            % correlations
                            Q = bckdsol(ise).Q;
                            corrRADe= Q(ilineRA,ilineDe) / sqrt(Q(ilineRA,ilineRA) * Q(ilineDe,ilineDe) );
                            
                            fprintf(fidSOUCAT, formatSOUCAT, source_ivs, source_iers,  num2str(res(1)), num2str(res(2)), res(3), sign,   num2str(abs(res(4))), num2str(res(5)), res(6) , errRA, errDe, corrRADe, mjd, mjd, mjd, '1', numobs, '0', globsol.sessions{ise});
                            clear mjd dRA dDe stdevRA stdevDe source_ivs source_iers numobs oldcolSOU kkRA kkDe RAapr Deapr sign res corrRADe ilineRA ilineDe
                        end
                    end
                end
            end
        end
    end
    

end

if solbck.eop==1; fclose(fidEOP); end
if solbck.zwd==1; fclose(fidZWD); end
if solbck.tgr==1; fclose(fidTGR); end
if solbck.ant==1; fclose(fidANT); end
if solbck.sou==1; fclose(fidSOU); end
if solbck.sou==1; fclose(fidSOUCAT); end

%  save -v7.3 ([path_outglob 'BACKWARD_SOLUTION/' DIROUT '/bckdsol'],'bckdsol');

    


