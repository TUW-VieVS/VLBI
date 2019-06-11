% This script reads in Goddard's source.cat as well as the most complete 
% ICRF2 list, vievsCrf, blokq.dat and makes one vievsCrf out of it.

% created by Matthias Madzak

% 20 March 2013 changed by Hana Krasna

clear all
close all
clc

% run from /VieVS/VLBI/CODE/MISC
pathF = '../../CRF/';

% filenames
icrf2nonVCSFilename=[pathF 'data/icrf2-non-vcs.dat'];
icrf2VCSFilename=[pathF 'data/icrf2-vcs-only.dat'];
vieCrf10aFilename=[pathF 'data/VieCRF10a.txt'];
sourceCatFilename=[pathF 'data/source.cat'];
blokqFilename=[pathF '../TRF/data/blokq.dat'];

translationFile=[pathF 'data/source_translation_IERS_IVS.dat'];

vievsCrfFilename=[pathF 'data/vievsCrf.txt'];


% load translation table
fid=fopen(translationFile);
transl=textscan(fid, ' %8c  %8c\n', 'commentstyle', '%'); % with %8s it crashed when a blank is within a source name, e.g. "3C 279"
fclose(fid);


% read source.cat ----------------
fid=fopen(sourceCatFilename);
sourceCatData=textscan(fid, '%8s %8s %f %f %f %f %f %f %f %f %11s', ...
    'commentstyle', '*', 'delimiter', '|');
fclose(fid);
% make cell string out of char array (this had problems - empty lines!)
%sourceCatData{11}=cellstr(sourceCatData{11});


% read icrf2 non VCS ---------------------
fid=fopen(icrf2nonVCSFilename);
icrf2nonVCSdata=textscan(fid, '%21c  %11c  %f %f %f %f %f %f %f %f %f %f %f %f %f %f', ...
    'headerlines', 23);
fclose(fid);


% read icrf2 VCS ---------------------
fid=fopen(icrf2VCSFilename);
icrf2VCSdata=textscan(fid, '%21c  %8s  %f %f %f  %f %f %f  %f %f  %f  %f %f %f %f %f', ...
    'headerlines', 20);
fclose(fid);


% read VieCRF10a -----------------
fid=fopen(vieCrf10aFilename);
vieCrfData=textscan(fid, '%8s %f %f %f %f %f %f', 'CommentStyle', '%',...
    'delimiter', '|');
fclose(fid);
%vieCrfData{1}=cellstr(vieCrfData{1});

% Separate the IERS and IVS with the internal translation table
% source_translation_IERS_IVS.dat
% if in VieCRF10a.dat was the IVS name, it is saved in column 8 and 
% in column 1 the name is overwritten with the IERS name 

for iSource=1:size(vieCrfData{1},1)
    curSouInVieCRF = strcmp(cellstr(transl{2}),deblank(vieCrfData{1}(iSource,1)));

    % if found -> write correct IERS name and IVS shift into column 8
    if sum(curSouInVieCRF)==1
        vieCrfData{8}(iSource,1) = vieCrfData{1}(iSource,:);
        vieCrfData{1}(iSource,1) = {transl{1}(find(curSouInVieCRF),1:8)};
    else
        vieCrfData{8}(iSource,1) ={'$'};
    end

end


% read blokq -----------------
fid=fopen(blokqFilename);

% preallocating
blokqData=cell(10000,8);
weAreAtSource=0;
blokqSourceInd=1;

while ~feof(fid)
    line=[fgetl(fid), '                                                        '];
        
    % if we are at source catalog
    if strcmp(line(1:31), '$$    SOURCE CATALOG goes here:')
        weAreAtSource=1;
    end
       
    % if '    ' are now at beginning of line - we have a source
    if weAreAtSource==1 && strcmp(line(1:4), '    ')
        % ivsname
        blokqData{blokqSourceInd,1}=line(5:12);   
        blokqData{blokqSourceInd,2}='';
        blokqData{blokqSourceInd,3}=str2double(line(15:16));
        blokqData{blokqSourceInd,4}=str2double(line(18:19));
        blokqData{blokqSourceInd,5}=str2double(line(21:29));
        blokqData{blokqSourceInd,6}=str2double(line(35:37));
        blokqData{blokqSourceInd,7}=str2double(line(39:40));
        blokqData{blokqSourceInd,8}=str2double(line(42:49));
        blokqData{blokqSourceInd,9}=line(58:end);
        blokqSourceInd=blokqSourceInd+1;
    end
end
        
fclose(fid);
blokqData(blokqSourceInd:end,:)=[];

% In blokq.dat the IERS and IVS (common) names are mixed together.
% Separate the IERS and IVS with the internal translation table
% source_translation_IERS_IVS.dat

% if in blokq.dat was the IVS name, it is saved in column 2 and 
% in column 1 the name is overwritten with the IERS name 

for iSource=1:size(blokqData,1)
    curSouInBlocq = strcmp(cellstr(transl{2}), deblank(blokqData{iSource,1}));

    % if found -> write correct IERS name and IVS shift into column 2
    if sum(curSouInBlocq)==1
        blokqData{iSource,2} = blokqData{iSource,1};
        blokqData{iSource,1} = transl{1}(find(curSouInBlocq),:);
    else
        blokqData{iSource,2} ={'$'};
    end

end

% in blockq.dat some sources are given several times. Check and take only
% those coordinates which appears as first
i=1;
for iSource=1:size(blokqData,1)
    curSouInBlocq = strcmp(blokqData(:,1),blokqData{iSource,1});
        if sum(curSouInBlocq)>1
            multip(i)={num2str(find(curSouInBlocq == 1))};
            i=i+1;
        end
end
DelMulti=[];
for iMulti=1:size(multip,2)
    a=[];
    a=find(strcmp(multip(iMulti),multip)==1);
    DelMulti = [DelMulti a(2:end)];
end

DelMulti=unique(DelMulti);
multip(DelMulti)=[];

Delblockq=[];
for i=1:length(multip)
   Delblockq = [Delblockq; str2num(multip{i})]; 
end

blokqData(Delblockq,:)=[];
        


% ====================
% write vievsCrf
% ====================
fid=fopen(vievsCrfFilename, 'w');

fprintf(fid, '* This is vievsCrf, the backup celestial reference frame for the\n* Vienna VLBI Software (VieVS).\n');
fprintf(fid, '* This file was originally computed from ICRF2 (icrf2-non-vcs.dat + icrf2-vcs-only.dat), SKED''s\n* source.cat, blokq.dat.\n*\n');
fprintf(fid, '*IERS-Name Common    hh mm ss.ssss    sdd mm as.sssss epoch  0.0     source    datum\n*\n');

% first: loop over icrf2, later add missing sources from source.cat file
% therefore: get logicals whether source.cat-source was already taken
VCSAlreadyTaken=zeros(size(icrf2VCSdata{1},1),1);
sourceCatAlreadyTaken=zeros(size(sourceCatData{1},1),1);
vieCrftxtAlreadyTaken=zeros(size(vieCrfData{1},1),1);
blokqDatAlreadyTaken=zeros(size(blokqData,1),1);

myFormat=' %-8s  %-8s  %02.0f %02.0f %09.6f  %+03.0f %02.0f %08.5f %6.1f %03.1f  %-15s %1.0f \n';

for iSource=1:length(icrf2nonVCSdata{1})
    
    % find current source (from icrf2) in source.cat (search for equal
    % IAU-Name == IERS Name)
    
    curSouInSourceCat=strcmp(sourceCatData{1},icrf2nonVCSdata{2}(iSource,1:8));% | strcmp(sourceCatData{2}, icrf2nonVCSdata{2}(iSource));
    curSouInVCS=      strcmp(icrf2VCSdata{2}, icrf2nonVCSdata{2}(iSource,1:8));   
    curSouInVieCrf=   strcmp(vieCrfData{1},   icrf2nonVCSdata{2}(iSource,1:8));
    curSouInBlokq=    strcmp(blokqData(:,1),  icrf2nonVCSdata{2}(iSource,1:8));
    
    % find IVS (current) name in the translation table 
    curSouInTransl = strcmp(cellstr(transl{1}), icrf2nonVCSdata{2}(iSource,1:8));
    if sum(curSouInTransl)==1
        curCommonName = char(transl{2}(find(curSouInTransl),:));
    else
        curCommonName='$';
    end
    
    
    % write data to vievsCrf
    fprintf(fid, myFormat, ...
        icrf2nonVCSdata{2}(iSource,1:8), curCommonName, ...
        icrf2nonVCSdata{3}(iSource), icrf2nonVCSdata{4}(iSource), icrf2nonVCSdata{5}(iSource), ...
        icrf2nonVCSdata{6}(iSource), icrf2nonVCSdata{7}(iSource), icrf2nonVCSdata{8}(iSource), ...
        2000, 0, 'ICRF2_non_VCS', 1);


    % make logical to one (for "already added to vievsCrf")
    % if found in just 
    sourceCatAlreadyTaken(curSouInSourceCat)=1;
    VCSAlreadyTaken(curSouInVCS)=1;
    vieCrftxtAlreadyTaken(curSouInVieCrf)   =1;
    blokqDatAlreadyTaken(curSouInBlokq)     =1;

end


% check if sources in only VCS have not been printed yet
% most probably all sources have to be added
if sum(VCSAlreadyTaken)<size(VCSAlreadyTaken,1)
    % get sources to still be added
    sourcesStillToPrint=find(VCSAlreadyTaken==0);
    
    for iSource=1:length(sourcesStillToPrint)
        
        % find IVS (current) name in the translation table 
        curSouInTransl = strcmp(cellstr(transl{1}), char(icrf2VCSdata{2}(iSource)));
        if sum(curSouInTransl)==1
            curCommonName = char(transl{2}(find(curSouInTransl),:));
        else
            curCommonName='$';
        end
    
        
        fprintf(fid, myFormat, ...
            char(icrf2VCSdata{2}(iSource)), curCommonName, ...
            icrf2VCSdata{3}(iSource), icrf2VCSdata{4}(iSource), icrf2VCSdata{5}(iSource), ...
            icrf2VCSdata{6}(iSource), icrf2VCSdata{7}(iSource), icrf2VCSdata{8}(iSource), ...
            2000, 0, 'ICRF2_VCS', 1);
 

        % also tell source.cat, that cur source has already been taken
        curSouInSourceCat=strcmp(sourceCatData{1}, icrf2VCSdata{2}{sourcesStillToPrint(iSource)});
        sourceCatAlreadyTaken(curSouInSourceCat)=1;
        
        
        % also tell "VieCRF10a", that cur source has already been taken
        curSouInVieCrf=strcmp(vieCrfData{1}, icrf2VCSdata{2}{sourcesStillToPrint(iSource)});
        vieCrftxtAlreadyTaken(curSouInVieCrf)=1;
        
        % "blocq"
        curSouInBlokq=strcmp(blokqData(:,1), icrf2VCSdata{2}{sourcesStillToPrint(iSource)});
        blokqDatAlreadyTaken(curSouInBlokq)=1;
        
    end
end


% check if sources in source.cat have not been printed yet
if sum(sourceCatAlreadyTaken)<size(sourceCatAlreadyTaken,1)
    % get sources to still be plotted
    sourcesStillToPrint=find(sourceCatAlreadyTaken==0);
    
    for iSource=1:length(sourcesStillToPrint)
        fprintf(fid, myFormat,...
            sourceCatData{1}{sourcesStillToPrint(iSource)},...
            sourceCatData{2}{sourcesStillToPrint(iSource)},...
            sourceCatData{3}(sourcesStillToPrint(iSource)),...
            sourceCatData{4}(sourcesStillToPrint(iSource)),...
            sourceCatData{5}(sourcesStillToPrint(iSource)),...
            sourceCatData{6}(sourcesStillToPrint(iSource)),...
            sourceCatData{7}(sourcesStillToPrint(iSource)),...
            sourceCatData{8}(sourcesStillToPrint(iSource)),...
            2000,0,...
            sourceCatData{11}{sourcesStillToPrint(iSource)},...
            1);
               
        % also tell "VieCRF10a", that cur source has already been taken
        curSouInVieCrf=strcmp(vieCrfData{1}, sourceCatData{1}{sourcesStillToPrint(iSource)}) | ...
                       strcmp(vieCrfData{1}, sourceCatData{2}{sourcesStillToPrint(iSource)});
        vieCrftxtAlreadyTaken(curSouInVieCrf)=1;
        
        % "blocq"
        curSouInBlokq=strcmp(blokqData(:,1), sourceCatData{1}{sourcesStillToPrint(iSource)}) | ...
                      strcmp(blokqData(:,1), sourceCatData{2}{sourcesStillToPrint(iSource)});
        blokqDatAlreadyTaken(curSouInBlokq)=1;
        
    end
end

% check if sources in vieCrf10a have not been added to vievsCrf
if sum(vieCrftxtAlreadyTaken)<size(vieCrftxtAlreadyTaken,1)
    sourcesStillToPrint=find(vieCrftxtAlreadyTaken==0);
    
    for iSource=1:length(sourcesStillToPrint)
        fprintf(fid, myFormat, ...
            char(vieCrfData{1}{sourcesStillToPrint(iSource)}),...
            char(vieCrfData{8}(sourcesStillToPrint(iSource))),...
            vieCrfData{2}(sourcesStillToPrint(iSource)),...
            vieCrfData{3}(sourcesStillToPrint(iSource)),...
            vieCrfData{4}(sourcesStillToPrint(iSource)),...
            vieCrfData{5}(sourcesStillToPrint(iSource)),...
            vieCrfData{6}(sourcesStillToPrint(iSource)),...
            vieCrfData{7}(sourcesStillToPrint(iSource)),...
            2000,0,...
            'VieCRF10a',...
            1);
        
        % also tell "blokq" that this source was already added to list
        curSouInBlokq=strcmp(blokqData(:,1), vieCrfData{1}{sourcesStillToPrint(iSource)});
        blokqDatAlreadyTaken(curSouInBlokq)=1;
            
    end
    
    
end

% check if sources in vieCrf10a have not been added to vievsCrf
if sum(blokqDatAlreadyTaken)<size(blokqDatAlreadyTaken,1)
    sourcesStillToPrint=find(blokqDatAlreadyTaken==0);
    
    for iSource=1:length(sourcesStillToPrint)
        fprintf(fid, myFormat,...
            char(blokqData{sourcesStillToPrint(iSource),1}),...
            char(blokqData{sourcesStillToPrint(iSource),2}),...
            blokqData{sourcesStillToPrint(iSource),3},...
            blokqData{sourcesStillToPrint(iSource),4},...
            blokqData{sourcesStillToPrint(iSource),5},...
            blokqData{sourcesStillToPrint(iSource),6},...
            blokqData{sourcesStillToPrint(iSource),7},...
            blokqData{sourcesStillToPrint(iSource),8},...
            2000,0,...
            'blokq.dat',...
            1);
    end
end
    
fclose(fid);


