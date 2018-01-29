% check which sources in NGS file are missing in the supersource file
% output is written in a txt file and can be copied in vievsCrf.txt
% Hana Krasna, Jan 03, 2018


clear all
% load supersource file
load('../supersource.mat')

% load NGS file
ngsfile = '17APR23KV_N006';
ngsfil = ['../../DATA/NGS/2017/' ngsfile];

%%

fid_ngs = fopen(ngsfil,'r');

    if fid_ngs ~= -1
        wholeFile = textscan(fid_ngs,'%s','delimiter', '\n', 'whitespace', '');
        wholeFile = wholeFile{1};
        idx_line = 1;
        nlines = length(wholeFile);
    else
        error(' couln''t find NGS-file %s',ngsfil)
    end

    aa=strcmp(wholeFile,'$END');
    idEND = find(aa==1);
    first_line = idEND(1);
    last_line = idEND(2);
fclose(fid_ngs);

numsou = last_line-first_line-1;


for i=1:numsou
    sname(i,:) = wholeFile{i+first_line}(1:8);
        
    RA(i)=(str2double(wholeFile{i+first_line}(11:12)) + ...
          str2double(wholeFile{i+first_line}(14:15))/60 + ...
          str2double(wholeFile{i+first_line}(20:28))/3600) /12*pi;
      
    desi=wholeFile{i+first_line}(30);
    si=1;
    if strcmp(desi,'-')
        si = -1;
    end
    
    De(i)=si*(str2double(wholeFile{i+first_line}(31:32)) + ...
        str2double(wholeFile{i+first_line}(34:35))/60 + ...
        str2double(wholeFile{i+first_line}(40:48))/3600) /180*pi;
    
    
    
    
    RAh(i,:) = wholeFile{i+first_line}(11:12);
    RAm(i,:) = wholeFile{i+first_line}(14:15);
    RAs(i,:) = wholeFile{i+first_line}(20:28);

    Desig(i) = wholeFile{i+first_line}(30);
    
    Deg(i,:) = wholeFile{i+first_line}(31:32);
    Dem(i,:) = wholeFile{i+first_line}(34:35);
    Des(i,:) = wholeFile{i+first_line}(40:48);
end

k=0;
for i=1:numsou
    iso=sum(strcmp(sname(i,:), cellstr({source.IVSname})));
%      ii=find(strcmp(sname(i,:), cellstr({source.IVSname})))

    if iso==0
        k=k+1;
        newso(k)=i;
    end

end




fid=fopen(['missing_sources_' ngsfile '.txt'],'wt');
for i=1:length(newso)
    k=1;
    if strcmp(Desig(newso(i)),'-')
        k=-1;
    end
    
    fprintf(fid,'%c%c%c%c%c%c%c%c      %c%c %c%c %c%c%c%c%c%c%c%c%c  %3.0f %c%c  %c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c %c%c%c%c%c%c%c%c%c%c%c%c%c%c    %c\n', ...
        sname(newso(i),:), RAh(newso(i),:), RAm(newso(i),:), RAs(newso(i),:), k*str2double(Deg(newso(i),:)), Dem(newso(i),:), Des(newso(i),:),' 2000.0 0.0', ngsfile, '0' );
end
fclose(fid);


