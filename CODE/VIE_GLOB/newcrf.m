% ************************************************************************
%   Description:
%   This function creates a .txt file with new sources coordinates. 
%
%   Input:										
%       globsol             structure with estimates
%       paths               paths to directories
%                           
%   Output:                
%      'crf_*.txt'
%
%   External calls: 	
%      -               					    											
%       
%   Coded for VieVS: 
%   04 Aug 2010 by Hana Spicakova
%
%   Revision: 
%   18 Oct 2012 by Hana Krásná: Output coded in a format of an offical CRF
%               catalogue
%   24 Oct 2012 by Hana Krásná: changes according to change of the units of
%               RA in globsol
%   24 Mar 2013 by Hana Krasna: changed according to supersource file
%   23 Jul 2013 by Hana Krasna: slight change of the output format
%   03 Nov 2016 by Hana Krasna: slight change of the output format
%   10 Apr 2017 by David Mayer: added IVS format
%%

function newcrf(globsol,paths,crf)
%paths=pathGS % separate function


if exist([paths.path_out 'CRF/' paths.out])~=7
	mkdir([paths.path_out 'CRF/' paths.out]);
end

fid=fopen([paths.path_out 'CRF/' paths.out '/crf_' paths.L2 '.txt'],'wt');

fprintf(fid,'%% A priori catalogue of source positions used for the analysis: %s \n', globsol.source.apriori_cat{1});
fprintf(fid,'%% Estimates dRA, dDe from LEVEL2 data: %s \n\n\n',paths.L2);

fprintf(fid,'%% source         RA [h min sec]              De [° min arcsec ]    \n\n');

% RA
RA_apr=globsol.source.apriori_rade(:,1);  %[rad]
dRA = globsol.source.drade(:,1); %[ms]

RA_est=RA_apr + dRA/1000/3600/12*pi; %[rad]
RA_est_h = RA_est/pi*12; %[h]

RA_h = floor(RA_est_h); %[h]
RA_min=floor((RA_est_h-RA_h)*60); %[min]
RA_sec=(((RA_est_h-RA_h)*60)-RA_min)*60; %[sec]


% De
De_apr=globsol.source.apriori_rade(:,2); %[rad]
dDe = globsol.source.drade(:,2); %[mas]

De_est=De_apr + dDe/1000/3600/180*pi; %[rad]

De_est_deg=De_est/pi*180; %[deg]


De_deg=fix(De_est_deg); %[deg]
De_min1=(De_est_deg-De_deg)*60; %[min]
id=find(De_min1<0);
De_min1(id)=De_min1(id)*(-1);
De_min=floor(De_min1);

De_sec=(De_min1-De_min)*60; %[sec]

[refname_s,ind]=sortrows(globsol.source.refname.IERS);
a1 = char(refname_s);

res= [RA_h RA_min RA_sec De_deg De_min De_sec De_est_deg]; % results
sres=res(ind,:); % sorted results

numSou = size(globsol.source.refname.IERS,1);


for i=1:numSou  % acoording to the columns in N-matrix
    fprintf(fid,'%c%c%c%c%c%c%c%c        %2.0f  %2.0f  %11.8f         %3.0f  %2.0f  %10.7f     \n', a1(i,:), sres(i,1), sres(i,2), sres(i,3), sres(i,4), sres(i,5), sres(i,6));
end

fclose(fid);


%% Create an ASCII CRF Catalogue in an "official" format

for i=1:length(crf)
    crfnames(i,:)=crf(i).IERSname;
end

souactiv = globsol.source.souactiv;

Q = globsol.Q;
dradecol = globsol.source.dradecol;

for k=1:numSou

    sou(k).refnameIERS(1:8)=globsol.source.refname.IERS(k,:);
    sou(k).refnameIVS(1:8)=globsol.source.refname.IVS(k,:);
    if isempty(sou(k).refnameIVS) || strcmp(sou(k).refnameIVS, ' ')
        sou(k).refnameIVS = '        ';
    end
	
    
    idcrf=(find(strcmp(cellstr(crfnames),cellstr(globsol.source.refname.IERS(k,:)))==1));
    if ~isempty(idcrf)
        sou(k).designation=crf(idcrf).designation;
        if isempty(sou(k).designation)
            sou(k).designation='ICRF J               ';
        end
    else
        sou(k).designation='ICRF J               ';  
    end

    idDatum = (find(strcmp(cellstr(globsol.source.datum),cellstr(globsol.source.refname.IERS(k,:)))==1));
    if ~isempty(idDatum)
        sou(k).dat='D';
    else
        sou(k).dat=' ';  
    end
  
    
    idObsSou=find(souactiv(k,:)>0);
    soumjdall = souactiv(end,idObsSou);
    sou(k).first_last_mean_O = [min(soumjdall) max(soumjdall) (min(soumjdall)+max(soumjdall))/2] ;
    sou(k).Nexp = length(idObsSou);
    sou(k).Nobs = sum(souactiv(k,idObsSou));
   
    
    
    % correlation
    sou(k).correl = Q(dradecol(k,1),dradecol(k,2)) / sqrt(Q(dradecol(k,1),dradecol(k,1)) * Q(dradecol(k,2),dradecol(k,2)) );
    
    
end

mRADe = globsol.source.sigma_rade./1000; %[s, as]


% sort
RAdeg = res(:,1) + res(:,2)./60 + res(:,3)./3600;
[a,idsort]=sort(RAdeg);



% get current date and time
curDate=clock;          % current date and time

fidOffic=fopen([paths.path_out 'CRF/' paths.out '/crf_catalogue_' paths.L2 '.txt'],'wt');
fidOffic_IVSformat = fopen([paths.path_out 'CRF/' paths.out '/crf_catalogue_IVSformat_' paths.L2 '.txt'],'wt');

fprintf(fidOffic, 'Created on %02.0f.%02.0f.%04.0f at %02.0f:%02.0f:%02.0f local time\n', curDate(3), curDate(2), curDate(1), curDate(4), curDate(5), curDate(6));
fprintf(fidOffic_IVSformat, 'Created on %02.0f.%02.0f.%04.0f at %02.0f:%02.0f:%02.0f local time\n\n', curDate(3), curDate(2), curDate(1), curDate(4), curDate(5), curDate(6));

fprintf(fidOffic,'\n\n\n\n');
fprintf(fidOffic,'---------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf(fidOffic,'ICRF  Designation     IERS Des. Inf. Right Ascension   Declination         Uncertainty            Corr.   Mean     First    Last      Nb     Nb     IVS name\n');
fprintf(fidOffic,'                                         J2000.0         J2000.0           R.A.      Dec.         RA-Dc   MJD      MJD      MJD       sess.  obs.       \n');
fprintf(fidOffic,'                                      h  m  s           o  m  as           s         as                   of observation span               \n');
fprintf(fidOffic,'---------------------------------------------------------------------------------------------------------------------------------------------------\n');

fprintf(fidOffic_IVSformat,'Columns  Units   Meaning                                                            \n');
fprintf(fidOffic_IVSformat,'   1     --      IVS designation                                                    \n');
fprintf(fidOffic_IVSformat,'   2     --      IERS designation                                                   \n');
fprintf(fidOffic_IVSformat,' 3-5     hms     Right ascension in hour, minutes and seconds of time               \n');
fprintf(fidOffic_IVSformat,' 6-8     dms     Declination in degree, minutes and seconds of degree               \n');
fprintf(fidOffic_IVSformat,'   9     sec     Formal uncertainty of the right ascension in second of time        \n');
fprintf(fidOffic_IVSformat,'  10     arcsec  Formal uncertainty of the declination                              \n');
fprintf(fidOffic_IVSformat,'  11     --      Correlation between right ascension and declination                \n');
fprintf(fidOffic_IVSformat,'  12     mjd     Mean epoch of observation                                          \n');
fprintf(fidOffic_IVSformat,'  13     mjd     First epoch of observation                                         \n');
fprintf(fidOffic_IVSformat,'  14     mjd     Last epoch of observation                                          \n');
fprintf(fidOffic_IVSformat,'  15     --      Number of sessions                                                 \n');
fprintf(fidOffic_IVSformat,'  16     --      Number of delays                                                   \n');
fprintf(fidOffic_IVSformat,'  17     --      Number of delay rates                                              \n');
fprintf(fidOffic_IVSformat,'  18     --      Estimation flag: GLO for global or ARC for session-wise            \n\n');

writeFormat='%21s  %8s  %1s  %02s %02s %011.8f  %s%02s %02s %010.7f  %10.8f %9.7f  %6.3f  %7.1f %7.1f %7.1f  %5.0f %6.0f      %8s\n';
writeFormat_IVS='%8s %8s    %02s %02s %011.8f    %s%02s %02s %010.7f%15.8f%14.7f%8.4f%8.1f%8.1f%8.1f%6.0f%7.0f%7.0f %3s\n';

for i=1:numSou
    
    if res(idsort(i),7)<0
        sign = '-';
    else
        sign = ' ';
    end
    
    if sou(idsort(i)).refnameIERS(1,1:3) == 'VIE'
        sou(idsort(i)).refnameIERSnoVIE = '        ';
    else
        sou(idsort(i)).refnameIERSnoVIE = sou(idsort(i)).refnameIERS;
    end

    
    fprintf(fidOffic,writeFormat, sou(idsort(i)).designation, sou(idsort(i)).refnameIERSnoVIE, sou(idsort(i)).dat,...
       num2str(res(idsort(i),1)), num2str(res(idsort(i),2)), res(idsort(i),3),...
       sign, num2str(abs(res(idsort(i),4))), ...
       num2str(res(idsort(i),5)), res(idsort(i),6), mRADe(idsort(i),:),...
       sou(idsort(i)).correl, sou(idsort(i)).first_last_mean_O(:,3),...
       sou(idsort(i)).first_last_mean_O(:,1:2), sou(idsort(i)).Nexp , sou(idsort(i)).Nobs , sou(idsort(i)).refnameIVS );

    fprintf(fidOffic_IVSformat,writeFormat_IVS, sou(idsort(i)).refnameIVS, sou(idsort(i)).refnameIERS, ...
       num2str(res(idsort(i),1)), num2str(res(idsort(i),2)), res(idsort(i),3),...
       sign, num2str(abs(res(idsort(i),4))), ...
       num2str(res(idsort(i),5)), res(idsort(i),6), mRADe(idsort(i),:),...
       sou(idsort(i)).correl, sou(idsort(i)).first_last_mean_O(:,3),...
       sou(idsort(i)).first_last_mean_O(:,1:2), sou(idsort(i)).Nexp , sou(idsort(i)).Nobs, 0, 'GLO');
	   
end

fclose(fidOffic);
fclose(fidOffic_IVSformat);
