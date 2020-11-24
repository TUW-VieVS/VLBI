% This function creates a supersource (.mat) file where  all source
% dependent information is stored
% 18.09.2012, Matthias Madzak

% 27 March 2013 changed by Hana Kr�sn�
% 29 April 2013 by Hana Kr�sn� bug corrected for negativ deglinations with zero degree
% 27 Jan 2014 by Hana Krasna: icrf2nonVCS added into supersource struct
% 14 Jul 2016 by David Mayer: added IVSnames and GSFC2015b catalogue to supersource struct
% 24 Sep 2016 by Hana Krasna: IVS_SrcNamesTable.txt instead of
% source_translation_IERS_IVS.dat, related changes to bookkeeping of IERS and
% IVS names, VieCRF13 instead of VieCRF10a, GUI updated, UserCRF can contain IERS
% or IVS names
% 09 Aug 2017 by Hana Krasna: related changes to update of vievsCrf.txt
% (only IVS names are stored there now). Automatic creation of dummy IERS
% names (VIE_xxxx) in case that there is no IERS name given in
% IVS_SrcNamesTable.txt
% 23 Aug 2017 by Hana Krasna: ismissing function removed as it is not
% compatible with older Matlab versions
%
%%
function [varargout]=mk_supersourceFile(inFiles, outFile, UserCrfNames)

% ==============
% 1. get inFiles
% ==============

ICRF2File=inFiles(1).name;
ICRF2FileVCS=inFiles(2).name;
ICRF3sxFile=inFiles(3).name;
VieCRF13File=inFiles(4).name;
vievsCrfFile=inFiles(5).name;
GSFC2015bFile =inFiles(6).name;
userCrfFile=inFiles(7).name;
IVSnamesFile =inFiles(8).name;

% =============
% 2. read files
% =============

% ----------------------------
% 2.1 vievsCrf.txt (main loop)
% ----------------------------
fprintf('Loading VieVSCrf.txt (BACKUP)\n\n')

fid=fopen(vievsCrfFile);
% get data from file
data=textscan(fid, '%8s %f %f %f %03s %f %f %f %f %19c %f', 'commentstyle', '*');
fclose(fid);


% get number of sources
nSources=size(data{1},1);
source(nSources)=struct;


for iSource=1:nSources
     curIVSname=[data{1}{iSource}, '          '];
    
    source(iSource).IVSname=curIVSname(1:8);
    source(iSource).IERSname='        ';
    
    % Right Ascension
    source(iSource).vievsCrf.ra= (data{2}(iSource) + data{3}(iSource)/60 + data{4}(iSource)/3600) * pi /12;
    
    % Declination 
    declchar = data{5}{iSource}; %(hana 29Apr13 - corrected for negative declinations -00deg)
    if strcmp(declchar(1), '-')
        signum=-1;
    else
        signum=1;
        if data{6}(iSource)<0 % consider format variant where sdd=0 and minus sign is written to mm
           signum=-1;
           data{6}(iSource)=data{6}(iSource)*signum;
        end
    end

    source(iSource).vievsCrf.de=signum*...
        (abs(str2double(declchar))+data{6}(iSource)/60+data{7}(iSource)/3600)*pi/180; % [rad]
    
    % epoch
    source(iSource).vievsCrf.epoch=data{8}(iSource);
    
    % motion
    source(iSource).vievsCrf.velocity=data{9}(iSource);
    
    % source (where from)
    source(iSource).vievsCrf.whereFrom=data{10}(iSource,:);
    
    % in_datum?
    if isnan(data{11}(iSource))
        source(iSource).vievsCrf.in_crf=1;
    else
        source(iSource).vievsCrf.in_crf=data{11}(iSource);
    end
    
end


% load IVS translation table from Goddard #############################
fid=fopen(IVSnamesFile);
IVSnames_raw=textscan(fid, '%8s %*2s %16s %*2s %10s %*2s %8s %16s %2f %*1s %2f %*1s %9.7f %*2s %1s %2f %*1s %2f %*1s %9.7f %*[^\n]', 'whitespace','', 'commentstyle', '#');
IVSnames.IVSname = IVSnames_raw{1,1};
IVSnames.designation = IVSnames_raw{1,2};
IVSnames.IERSname = IVSnames_raw{1,4};
for i=1:length(IVSnames.IERSname)
    if strcmp(IVSnames.IERSname(i),'-       ')
        IVSnames.IERSname(i)=IVSnames.IVSname(i);
    end
end
IVSnames.ra = (IVSnames_raw{1,6} + IVSnames_raw{1,7}./60 + IVSnames_raw{1,8}./3600).*pi/12;
IVSnames.de = (IVSnames_raw{1,10} + IVSnames_raw{1,11}./60 + IVSnames_raw{1,12}./3600).*pi/180;
IVSnames.de(strcmp(IVSnames_raw{1,9},'-')) = IVSnames.de(strcmp(IVSnames_raw{1,9},'-')).*(-1);
clear IVSnames_raw IVSnamesFile
fclose(fid);


% write info about translation table
fprintf('Translation table applied:\n');
iIERS=0;
for iSource=1:nSources
    curIVSname = ['        '];
    curIVSname = source(iSource).IVSname;

    % find source in translation table 
    curIVSnameInIERStable=~cellfun(@isempty, strfind(cellstr(IVSnames.IVSname), (curIVSname)));
    indCurSoInTransl=find(curIVSnameInIERStable);
    if ~isempty(indCurSoInTransl)
        source(iSource).IERSname=char(IVSnames.IERSname(indCurSoInTransl));  
        source(iSource).designation=['ICRF ' char(IVSnames.designation(indCurSoInTransl))];
    end
    
    % if IVS name differs from IERS name or is missing in the Goddard
    % translation table
    if strcmp(curIVSname,source(iSource).IERSname)==0
        if source(iSource).IERSname == '        '
            iIERS = iIERS+1;
            cIERS = num2str(iIERS,'%05.0f');
            
%             fprintf('No IERS name found in IVS_SrcNamesTable.txt for %s. id(%1.0f) \nAt the moment the IERS name VIE_%04.0f will be written to supersourcefile.\n\n',...
%                 curIVSname,iSource,iIERS);
            source(iSource).IERSname = ['VIE' cIERS];
        else
%             fprintf('source(%1.0f).IVSname: %s \nsource(%1.0f).IERSname: %c%c%c%c%c%c%c%c \n\n',...
%                 iSource, curIVSname, iSource, source(iSource).IERSname);
        end
    end
    
    
%     %% DELETE when VieVS automatic processing on server will use the VieVS Version 3.0 !!!
%     source(iSource).commonName = {source(iSource).IVSname};

end






fprintf('Loading gsf2015b_astro.sou.txt \n\n')

%load GSFC2015b catalog #############################
fid=fopen(GSFC2015bFile);

% GSFC2015b_raw=textscan(fid, '%*9s %8s %*1s %2f %*1s %2f %*1s %11.9f %*1s %1s %2f %*1s %2f %*1s %10.8f %*1s %10.9f %*1s %9.8f %*1s %1s %5f %*1s %7.1f %*1s %7.1f %*1s %7.1f %*1s %5f %*1s %6f %*1s %6f %*1s %3s%*[^\n]', 'whitespace','', 'commentstyle', '#', 'HeaderLines', 8);
% 
% GSFC2015b.IVSname = GSFC2015b_raw{1,1};
% GSFC2015b.ra = (GSFC2015b_raw{1,2} + GSFC2015b_raw{1,3}./60 + GSFC2015b_raw{1,4}./3600).*pi/12;
% GSFC2015b.de = (GSFC2015b_raw{1,6} + GSFC2015b_raw{1,7}./60 + GSFC2015b_raw{1,8}./3600).*pi/180;
% GSFC2015b.de(strcmp(GSFC2015b_raw{1,5},'-')) = GSFC2015b.de(strcmp(GSFC2015b_raw{1,5},'-')).*(-1);
% GSFC2015b.ra_sigma = GSFC2015b_raw{1,9}./3600*pi/12;
% GSFC2015b.de_sigma = GSFC2015b_raw{1,10}./3600*pi/180;
% GSFC2015b.corr = GSFC2015b_raw{1,12};
% GSFC2015b.corr(strcmp(GSFC2015b_raw{1,11},'-')) = GSFC2015b.corr(strcmp(GSFC2015b_raw{1,11},'-')).*(-1);
% GSFC2015b.meanMjd = GSFC2015b_raw{1,13};
% GSFC2015b.firstMjd = GSFC2015b_raw{1,14};
% GSFC2015b.lastMjd = GSFC2015b_raw{1,15};
% GSFC2015b.numberSess = GSFC2015b_raw{1,16};
% GSFC2015b.numberObs = GSFC2015b_raw{1,17};
% GSFC2015b.estType = GSFC2015b_raw{1,19};

GSFC2015b_raw=textscan(fid, '%*8c %8c  %f %f %11c  %3c %f %10c    %f %f %f    %f %f %f    %5c %6c %3c  %*[^\n]',  'commentstyle', '#','delimiter', ' ', 'HeaderLines', 8);

GSFC2015b.IVSname = GSFC2015b_raw{1,1};
GSFC2015b.ra = (GSFC2015b_raw{1,2} + GSFC2015b_raw{1,3}./60 + str2num(GSFC2015b_raw{1,4})./3600).*pi/12;
GSFC2015b.de = (abs(str2num(GSFC2015b_raw{1,5})) + GSFC2015b_raw{1,6}./60 + str2num(GSFC2015b_raw{1,7})./3600).*pi/180;
GSFC2015b.de(strcmp(cellstr(GSFC2015b_raw{5}(:,1)),'-')) = GSFC2015b.de(strcmp(cellstr(GSFC2015b_raw{5}(:,1)),'-')).*(-1);
GSFC2015b.ra_sigma = GSFC2015b_raw{1,8}./3600*pi/12;
GSFC2015b.de_sigma = GSFC2015b_raw{1,9}./3600*pi/180;
GSFC2015b.corr = GSFC2015b_raw{1,10};
GSFC2015b.meanMjd = GSFC2015b_raw{1,11};
GSFC2015b.firstMjd = GSFC2015b_raw{1,12};
GSFC2015b.lastMjd = GSFC2015b_raw{1,13};
GSFC2015b.numberSess = str2num(GSFC2015b_raw{1,14});
GSFC2015b.numberObs = str2num(GSFC2015b_raw{1,15});
GSFC2015b.estType = GSFC2015b_raw{1,16};

clear GSFC2015b_raw GSFC2015bFile
fclose(fid);

% write info about translation table
fprintf('Compare GSFC2015b with the VieVSCRF backup.\n');

    %##############################
for iSource = 1:size(GSFC2015b.IVSname,1)

    index_IVS = strcmp(GSFC2015b.IVSname(iSource,:), cellstr({source.IVSname}));
    if sum(index_IVS) == 1
        if (sqrt((GSFC2015b.ra(iSource) - source(index_IVS).vievsCrf.ra)^2+(GSFC2015b.de(iSource) - source(index_IVS).vievsCrf.de)^2)*180/pi)>10/3600
            fprintf('Check if source %s was matched correctly. Difference in RA is %5.1f [as] and in DEC %5.1f [as]. source(%1.0f) \n', GSFC2015b.IVSname(iSource,:),(GSFC2015b.ra(iSource) - source(index_IVS).vievsCrf.ra)*180/pi*3600,(GSFC2015b.de(iSource) - source(index_IVS).vievsCrf.de)*180/pi*3600, find(index_IVS));
        end
        source(index_IVS).GSFC2015b.ra = GSFC2015b.ra(iSource);
        source(index_IVS).GSFC2015b.de = GSFC2015b.de(iSource);
        source(index_IVS).GSFC2015b.ra_sigma = GSFC2015b.ra_sigma(iSource);
        source(index_IVS).GSFC2015b.de_sigma = GSFC2015b.de_sigma(iSource);
        source(index_IVS).GSFC2015b.corr = GSFC2015b.corr(iSource);
        source(index_IVS).GSFC2015b.meanMjd = GSFC2015b.meanMjd(iSource);
        source(index_IVS).GSFC2015b.firstMjd = GSFC2015b.firstMjd(iSource);
        source(index_IVS).GSFC2015b.lastMjd = GSFC2015b.lastMjd(iSource);
        source(index_IVS).GSFC2015b.numberSess = GSFC2015b.numberSess(iSource);
        source(index_IVS).GSFC2015b.numberObs = GSFC2015b.numberObs(iSource);
        source(index_IVS).GSFC2015b.estType = GSFC2015b.estType(iSource);
    else
        fprintf('\nSource %s from GSFC2015b not found in VieVSCRF or multiple entry. Solve the problem! Check vievsCrf.txt\n\n', GSFC2015b.IVSname(iSource,:))
    end
end

    %##############################

    
fprintf('\nLoading icrf2-non-vcs.dat\n\n')
fid=fopen(ICRF2File);
data=textscan(fid, '%21c  %11c  %f %f %f %03s %f %f %f %f %f %f %f %f %f %f', ...
    'headerlines', 23, 'delimiter', '\n');
fclose(fid);



% match coordinates from this catalogue
for k=1:size(data{1},1)
 
    % in icrf2-non-vcs only iers names are always given
    antennaFound=strcmp(data{2}(k,1:8), {source.IERSname});
    
    % find index where to put the new data
    if sum(antennaFound)>0
        if sum(antennaFound)>1 % should not happen
            fprintf('More than one source with name %-8s found in struct.\n', data{2}(k,1:8));
        end
        indexInSourceStruct=find(antennaFound);
    else
        fprintf('Source %-8s in ICRF2 was not found in vievsCrf - append to struct.\n', data{2}(k,1:8)); % should not happen
        indexInSourceStruct=size(source,2)+1;
        source(indexInSourceStruct).IERSname=data{2}(k,1:8); % this name exists when the source is found (so is only written here)
    end
    
    % since there might be more sources found with the same IAU_name - go
    % through all of them and add data
    % --> THIS SHOULD NOT HAPPEN!!! (hana)
    for iFoundIndices=1:length(indexInSourceStruct)
        % one of the found sources
        curIndInSourceStruct=indexInSourceStruct(iFoundIndices);
        % names
%         source(curIndInSourceStruct).designation=data{1}(k,:);
        % defining source or not
        if strcmp(data{2}(k,11), 'D')
            source(curIndInSourceStruct).icrf2.defining=1;
            source(curIndInSourceStruct).icrf2nonVCS.defining=1;
        else
            source(curIndInSourceStruct).icrf2.defining=0;
            source(curIndInSourceStruct).icrf2nonVCS.defining=0;
        end

        
        % Save the coordinates in icrf2 and to icrf2nonVCS
        
        % Right Ascension
        source(curIndInSourceStruct).icrf2.ra=(data{3}(k)+data{4}(k)/60+...
            data{5}(k)/3600)*pi/12;                                 % [rad]
        source(curIndInSourceStruct).icrf2nonVCS.ra=(data{3}(k)+data{4}(k)/60+...
            data{5}(k)/3600)*pi/12;                                 % [rad]


        % Declination 
        declchar = data{6}{k}; %(hana 29Apr13 - corrected for negative declinations -00deg)
        if strcmp(declchar(1), '-')
            signum=-1;
        else
            signum=1;
        end
        source(curIndInSourceStruct).icrf2.de=signum*...
            (abs(str2double(declchar))+data{7}(k)/60+data{8}(k)/3600)*pi/180; % [rad]
        source(curIndInSourceStruct).icrf2nonVCS.de=signum*...
            (abs(str2double(declchar))+data{7}(k)/60+data{8}(k)/3600)*pi/180; % [rad]

        
        % sigmas
        source(curIndInSourceStruct).icrf2.ra_sigma=data{9}(k)/3600*pi/12;  % [rad]
        source(curIndInSourceStruct).icrf2.de_sigma=data{10}(k)/3600*pi/180; % [rad]
        source(curIndInSourceStruct).icrf2nonVCS.ra_sigma=data{9}(k)/3600*pi/12;  % [rad]
        source(curIndInSourceStruct).icrf2nonVCS.de_sigma=data{10}(k)/3600*pi/180; % [rad]

        % correlation
        source(curIndInSourceStruct).icrf2.corr=data{11}(k);
        source(curIndInSourceStruct).icrf2nonVCS.corr=data{11}(k);

        % observation span
        source(curIndInSourceStruct).icrf2.meanMjd=data{12}(k);
        source(curIndInSourceStruct).icrf2.firstMjd=data{13}(k);
        source(curIndInSourceStruct).icrf2.lastMjd=data{14}(k);
        source(curIndInSourceStruct).icrf2nonVCS.meanMjd=data{12}(k);
        source(curIndInSourceStruct).icrf2nonVCS.firstMjd=data{13}(k);
        source(curIndInSourceStruct).icrf2nonVCS.lastMjd=data{14}(k);

        % number of sources/observations
        source(curIndInSourceStruct).icrf2.numberSess=data{15}(k);
        source(curIndInSourceStruct).icrf2.numberObs=data{16}(k);
        source(curIndInSourceStruct).icrf2nonVCS.numberSess=data{15}(k);
        source(curIndInSourceStruct).icrf2nonVCS.numberObs=data{16}(k);

    end
end


fprintf('Loading icrf2-vcs-only.dat\n\n')
fid=fopen(ICRF2FileVCS);
    data=textscan(fid, '%21c  %8s  %f %f %f  %03s %f %f  %f %f  %f  %f %f %f %f %f', ...
                  'headerlines', 20);
fclose(fid);

% match the coordinates from this catalogue
for k=1:size(data{1},1)
    
    % try to find source name
    % in icrf2-vcs-only only iers names are always given
    antennaFound=strcmp(data{2}(k), {source.IERSname});
    
    % find index where to put the new data
    if sum(antennaFound)>0
        if sum(antennaFound)>1
            fprintf('More than one source with name %-8s found in struct.\n', char((data{2}(k)))); % should not happen
        end
        indexInSourceStruct=find(antennaFound);
    else
        fprintf('Source %-8s in ICRF2 was not found in vievsCrf - append to struct.\n', char((data{2}(k)))); % should not happen
        indexInSourceStruct=size(source,2)+1;
        source(indexInSourceStruct).IERSname=data{2}(k); % this name exists when the source is found (so is only written here)
    end
    
    % since there might be more sources found with the same IAU_name - go
    % through all of them and add data
    % --> THIS SHOULD NOT HAPPEN!!! (hana)
    for iFoundIndices=1:length(indexInSourceStruct)
        % one of the found sources
        curIndInSourceStruct=indexInSourceStruct(iFoundIndices);
        
        % names
%         source(curIndInSourceStruct).designation=data{1}(k,:);

        % VCS source is never defining
        source(curIndInSourceStruct).icrf2.defining=0;

        % Right Ascension
        source(curIndInSourceStruct).icrf2.ra=(data{3}(k)+data{4}(k)/60+...
            data{5}(k)/3600)*pi/12;                                 % [rad]


        % Declination 
        declchar = data{6}{k}; %(hana 29Apr13 - corrected for negative declinations -00deg)
        if strcmp(declchar(1), '-')
            signum=-1;
        else
            signum=1;
        end
        source(curIndInSourceStruct).icrf2.de=signum*...
            (abs(str2double(declchar))+data{7}(k)/60+data{8}(k)/3600)*pi/180; % [rad]

        % sigmas
        source(curIndInSourceStruct).icrf2.ra_sigma=data{9}(k)/3600*pi/12;  % [rad]
        source(curIndInSourceStruct).icrf2.de_sigma=data{10}(k)/3600*pi/180; % [rad]

        % correlation
        source(curIndInSourceStruct).icrf2.corr=data{11}(k);

        % observation span
        source(curIndInSourceStruct).icrf2.meanMjd=data{12}(k);
        source(curIndInSourceStruct).icrf2.firstMjd=data{13}(k);
        source(curIndInSourceStruct).icrf2.lastMjd=data{14}(k);

        % number of sources/observations
        source(curIndInSourceStruct).icrf2.numberSess=data{15}(k);
        source(curIndInSourceStruct).icrf2.numberObs=data{16}(k);
    end
end



fprintf('Loading ICRF3sx\n\n')
fid=fopen(ICRF3sxFile);
data=textscan(fid, '%21c  %11c  %f %f %f %03s %f %f %f %f %f %f %f %f %f %f %f', ...
    'headerlines', 22, 'delimiter', '\n');

fclose(fid);

for i = 1: length(data{1})
    index_bool = strcmp(data{2}(i,1:8), {source.IERSname}');
    if ~sum(index_bool)
        index_bool = strcmp(data{2}(i,1:8), {source.IVSname}');
        if ~sum(index_bool) 
            fprintf(1,'Source %s is in ICRF3 but was not found in supersource file. It is not written into the supersource file.\n', data{2}(i,1:8));
            continue
        end
    end
    
    RA = deg2rad((data{3}(i) + data{4}(i)/60 + data{5}(i)/3600)*15); %[rad]
    if strncmp('-', data{6}(i),1)
        DE = deg2rad((str2double(data{6}(i)) - data{7}(i)/60 - data{8}(i)/3600)); %[rad]
    else
        DE = deg2rad((str2double(data{6}(i)) + data{7}(i)/60 + data{8}(i)/3600)); %[rad]
    end
    
    ang_sep = sqrt(((RA-source(index_bool).vievsCrf.ra)*cos(DE))^2 + (DE-source(index_bool).vievsCrf.de)^2);
    if ang_sep > deg2rad(10/3600) % check if sources is less than 10 as from vievsCrf --> if yes than it is a wrong match
        fprintf(1,'Possible wrong match (angular separation of %4.0f as) for source %s it will still be written into the supersource file. Please check by hand\n', rad2deg(ang_sep)*3600, data{2}(i,1:8));
    end
    
    if strcmp('D', data{2}(i,end))
        source(index_bool).icrf3sx.defining = 1;
    else
        source(index_bool).icrf3sx.defining = 0;
    end
    source(index_bool).icrf3sx.ra = RA; %[rad]
    source(index_bool).icrf3sx.de = DE; %[rad]
    source(index_bool).icrf3sx.ra_sigma = deg2rad(data{9}(i)*15/3600); %[rad]
    source(index_bool).icrf3sx.de_sigma = deg2rad(data{10}(i)/3600); %[rad]
    source(index_bool).icrf3sx.corr = data{11}(i); 
    source(index_bool).icrf3sx.meanMjd = data{12}(i); 
    source(index_bool).icrf3sx.firstMjd = data{13}(i); 
    source(index_bool).icrf3sx.lastMjd = data{14}(i); 
    source(index_bool).icrf3sx.numberSess = data{15}(i); 
    source(index_bool).icrf3sx.numberObs = data{16}(i); 
    
end



fprintf('Loading VieCRF13.txt\n\n')
fid=fopen(VieCRF13File);
data=textscan(fid, '%21s  %11c  %f %f %f %03s %f %f %f %f %f %f %f %f %f %f', 'headerlines', 10, 'delimiter', '\n');
fclose(fid);



% match coordinates from this catalogue
for k=1:size(data{1},1)
 
    % in VieCRF only iers names are always given
    antennaFound=strcmp(data{2}(k,1:8), {source.IERSname});
    indexInSourceStruct=find(antennaFound);
    
    % since there might be more sources found with the same IAU_name - go
    % through all of them and add data
    % --> THIS SHOULD NOT HAPPEN!!! (hana)
    for iFoundIndices=1:length(indexInSourceStruct)
        % one of the found sources
        curIndInSourceStruct=indexInSourceStruct(iFoundIndices);
        % defining source or not
        if strcmp(data{2}(k,11), 'D')
            source(curIndInSourceStruct).VieCRF13.defining=1;
        else
            source(curIndInSourceStruct).VieCRF13.defining=0;
        end

        
        % Right Ascension
        source(curIndInSourceStruct).VieCRF13.ra=(data{3}(k)+data{4}(k)/60+...
            data{5}(k)/3600)*pi/12;                                 % [rad]
        % Declination 
        declchar = data{6}{k}; 
        if strcmp(declchar(1), '-')
            signum=-1;
        else
            signum=1;
        end
        source(curIndInSourceStruct).VieCRF13.de=signum*...
            (abs(str2double(declchar))+data{7}(k)/60+data{8}(k)/3600)*pi/180; % [rad]

        
        % sigmas
        source(curIndInSourceStruct).VieCRF13.ra_sigma=data{9}(k)/3600*pi/12;  % [rad]
        source(curIndInSourceStruct).VieCRF13.de_sigma=data{10}(k)/3600*pi/180; % [rad]

        % correlation
        source(curIndInSourceStruct).VieCRF13.corr=data{11}(k);

        % observation span
        source(curIndInSourceStruct).VieCRF13.meanMjd=data{12}(k);
        source(curIndInSourceStruct).VieCRF13.firstMjd=data{13}(k);
        source(curIndInSourceStruct).VieCRF13.lastMjd=data{14}(k);

        % number of sources/observations
        source(curIndInSourceStruct).VieCRF13.numberSess=data{15}(k);
        source(curIndInSourceStruct).VieCRF13.numberObs=data{16}(k);

    end
end



% ------------
%  User CRF
% ------------
%
% current format of the User CRF:
%0002+200         0   4  35.75830819          20  19  42.3178689     
%0002-478         0   4  35.65550076         -47  36  19.6037040     
%0003+380         0   5  57.17539180          38  20  15.1490143 
%2356+385        23  59  33.18080106          38  50  42.3183135     
%2357-318        23  59  35.49154622         -31  33  43.8245431   

if ~isempty(userCrfFile)
    fprintf('Loading %s (user crf file)\n\n', userCrfFile)

    fid=fopen(userCrfFile);
        data=textscan(fid, '%8c %f %f %f %03s %f %f', 'commentstyle', '%'); %if using 8s instead of 8c problems with blank spaces in the source name occure
    fclose(fid);
    
    if UserCrfNames.IERS
        for i=1: size(data{1},1)
            sourceFound=strcmp(data{1}(i,1:8), {source.IERSname});
            indexInSourceStructs=find(sourceFound);
            if sum(sourceFound)>0
                indexInSourceStruct=indexInSourceStructs(1);
            else
                indexInSourceStruct=length(source)+1;
                source(indexInSourceStruct).IERSname=data{1}(i,1:8);
                fprintf('Source %s not found in supersource file (IERS names). Appended at the end.\n', data{1}(i,1:8))
                
                % Try to find the new source in the Goddard translation
                % table, if not found write the same IERS and IVS name.
                curIERSname = source(indexInSourceStruct).IERSname;
                % find source in translation table 
                curIERSnameInIERStable=~cellfun(@isempty, strfind(cellstr(IVSnames.IERSname), strtrim(curIERSname)));
                indCurSoInTransl=find(curIERSnameInIERStable);
                if ~isempty(indCurSoInTransl)
                    source(indexInSourceStruct).IVSname=char(IVSnames.IVSname(indCurSoInTransl));  
                    fprintf('IVS name %s matched from the translation table.\n',source(indexInSourceStruct).IVSname)
                else
                    source(indexInSourceStruct).IVSname=data{1}(i,1:8);
                    fprintf('Source %s is stored under the same IERS and IVS name.\n',data{1}(i,1:8))
                end
            end
            source(indexInSourceStruct).userCrf.ra=(data{2}(i)+data{3}(i)/60+data{4}(i)/3600)*pi/12; % [rad]

            % Declination 
            declchar = data{5}{i};
            if strcmp(declchar(1), '-')
                signum=-1;
            else
                signum=1;
            end
            source(indexInSourceStruct).userCrf.de=signum*(abs(str2double(declchar))+data{6}(i)/60+data{7}(i)/3600)*pi/180; % [rad]
        end
        
    elseif UserCrfNames.IVS  
        for i=1:size(data{1},1)
            sourceFound=strcmp(data{1}(i,1:8), {source.IVSname});
            indexInSourceStructs=find(sourceFound);
            if sum(sourceFound)>0
                indexInSourceStruct=indexInSourceStructs(1);
            else
                indexInSourceStruct=length(source)+1;
                source(indexInSourceStruct).IVSname=data{1}(i,1:8);
                fprintf('Source %s not found in supersource file (IVS names). Appended at the end.\n', data{1}(i,1:8))
                
                % Try to find the new source in the Goddard translation
                % table, if not found write the same IERS and IVS name.
                curIVSname = source(indexInSourceStruct).IVSname;
                % find source in translation table 
                curIVSnameInIVStable=~cellfun(@isempty, strfind(cellstr(IVSnames.IVSname), strtrim(curIVSname)));
                indCurSoInTransl=find(curIVSnameInIVStable);
                if ~isempty(indCurSoInTransl)
                    source(indexInSourceStruct).IERSname=char(IVSnames.IERSname(indCurSoInTransl));  
                    fprintf('IERS name %s matched from the translation table.\n',source(indexInSourceStruct).IERSname)
                else
                    source(indexInSourceStruct).IERSname=data{1}(i,1:8);
                    fprintf('Source %s is stored under the same IERS and IVS name.\n',data{1}(i,1:8))
                end
            end
            source(indexInSourceStruct).userCrf.ra=(data{2}(i)+data{3}(i)/60+data{4}(i)/3600)*pi/12; % [rad]

            % Declination 
            declchar = data{5}{i};
            if strcmp(declchar(1), '-')
                signum=-1;
            else
                signum=1;
            end
            source(indexInSourceStruct).userCrf.de=signum*(abs(str2double(declchar))+data{6}(i)/60+data{7}(i)/3600)*pi/180; % [rad]
        end
    end
else
    fprintf('\nNo user-CRF selected.\n')
end


% % last loop over all sources which have no vievsCrf (backup) coords -> add
% % them to vievsCrf
% sourceNoVievsCrfCoords=cellfun(@isempty, {source.vievsCrf});
% if sum(sourceNoVievsCrfCoords)>0
%     % --> we have sources with no vievsCrf coords
%     % get indices of those
%     sourceInd=find(sourceNoVievsCrfCoords);
%     
%     for iSource=1:length(sourceInd)
%         % if we have icrf2 coords -> use them
%         if ~isempty(source(sourceInd(iSource)).icrf2)
%             curRa=source(sourceInd(iSource)).icrf2.ra;
%             curDe=source(sourceInd(iSource)).icrf2.de;
%             
%         elseif ~isempty(source(sourceInd(iSource)).icrf1Ext2)
%             curRa=source(sourceInd(iSource)).icrf1Ext2.ra;
%             curDe=source(sourceInd(iSource)).icrf1Ext2.de;
%         elseif ~isempty(source(sourceInd(iSource)).VieCRF13)
%             curRa=source(sourceInd(iSource)).VieCRF13.ra;
%             curDe=source(sourceInd(iSource)).VieCRF13.de;
%         else
%             fprintf('Error: No vievsCrf, but also not other CRF?!\n')
%         end
%         source(sourceInd(iSource)).vievsCrf.ra=curRa;
%         source(sourceInd(iSource)).vievsCrf.de=curDe;
%     end
% end


% Check if the IERS names and  IVS names in the supersource.mat are unique
fprintf('\nFinal check - check if the IERS and IVS names in the supersource file are unique\n\n');

nSource=length({source.IERSname});
uniIERS=length(unique(cellstr(char(source.IERSname))));
uniIVS=length(unique(cellstr(char(source.IVSname))));

if nSource~=uniIERS || nSource~=uniIVS
    fprintf('!!!Multiple entries found:\n');

    b=[];
    [~,b]=unique(cellstr(char(source.IERSname)));
    IndMultiIERS=setdiff(1:nSource,b);
    for i=1:length(IndMultiIERS)
        if ~isempty(source(IndMultiIERS(i)).IERSname)
            fprintf('More than one source with the IERS name %-8s found in the supersource.mat! Please correct it!\n',source(IndMultiIERS(i)).IERSname)
        end
    end
    
    b=[];
    [~,b]=unique(cellstr(char(source.IVSname)));
    IndMultiIVS=setdiff(1:nSource,b);
    
    for i=1:length(IndMultiIVS)
        if ~isempty(source(IndMultiIVS(i)).IVSname)
            fprintf('More than one source with the IVS name %-8s found in the supersource.mat! Please correct it!\n',source(IndMultiIVS(i)).IVSname)
        end
    end
else
    fprintf('OK, the IERS and IVS names in the supersource file are unique\n\n');
end



% =============================================
% Save to disk and define output of function
% =============================================

fprintf('\nSaving supersource struct\n\n');

% save to output mat file
save(outFile, 'source');

% make temporary message
message=sprintf('Supersource file was successfully written to %s', outFile);

% if output is wanted = error/successful message
if nargout>0
    varargout{1}=message;
else 
    % printf message to screen
    fprintf(message);
end


