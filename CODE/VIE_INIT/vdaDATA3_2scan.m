function anttab=vdaDATA3_2scan(DATA3)


% STATION information - DATA.3
ind=contains(DATA3,'DATA.3 AIR_TEMP');
ourblock=DATA3(ind);
anttab=zeros(size(ourblock,1),5); % change nst1obs+nst2obs
for i=1:size(ourblock,1)
    input_str = ourblock{i};
    seg = split(input_str);
    
    anttab(i,1) = str2double(seg(3));
    anttab(i,2) = str2double(seg(4));
    
    D=find(char(seg(7))=='D');
    seg{7}(1,D)='e';
    
    anttab(i,3) = str2double(seg(7));
end
id999=find(anttab(i,3)>-900);
anttab(id999,3) = anttab(id999,3)-273.15; % K -> Celsius


ind=contains(DATA3,'DATA.3 ATM_PRES');
ourblock=DATA3(ind);
for i=1:size(ourblock,1)
    input_str = ourblock{i};
    seg = split(input_str);
    
%     anttab(i,1) = str2double(seg(3));
%     anttab(i,2) = str2double(seg(4));
    
    D=find(char(seg(7))=='D');
    seg{7}(1,D)='e';
    
    anttab(i,4) = str2double(seg(7));
end
anttab(:,4) = anttab(:,4)./100; % Pa -> hPa

ind=contains(DATA3,'DATA.3 REL_HUMD');
ourblock=DATA3(ind);
for i=1:size(ourblock,1)
    input_str = ourblock{i};
    seg = split(input_str);
    
%     anttab(i,1) = str2double(seg(3));
%     anttab(i,2) = str2double(seg(4));
    
    D=find(char(seg(7))=='D');
    seg{7}(1,D)='e';
    
    anttab(i,5) = str2double(seg(7));
end


% %  Check data value:
tdry=anttab(:,3);
relHum=anttab(:,5);
ee = 6.1078 .* exp((17.1 .* tdry) ./ (235 + tdry)) .* relHum; % formula by Magnus * relative humidity

id999=find(ee<-100);
ee(id999)=-999;
anttab(:,6)=ee;



ind=contains(DATA3,'DATA.3 CABL_DEL');
ourblock=DATA3(ind);
for i=1:size(ourblock,1)
    input_str = ourblock{i};
    seg = split(input_str);
    
%     anttab(i,1) = str2double(seg(3));
%     anttab(i,2) = str2double(seg(4));
    
    D=find(char(seg(7))=='D');
    seg{7}(1,D)='e';
    
    anttab(i,7) = str2double(seg(7))*10^9; %ns
end

ind=contains(DATA3,'DATA.3 MEANCABL');
ourblock=DATA3(ind);
for i=1:size(ourblock,1)
    input_str = ourblock{i};
    seg = split(input_str);
    id = ismember(anttab(:,2),str2double(seg(5)));    
    D=find(char(seg(7))=='D');
    seg{7}(1,D)='e';
    anttab(id,8) = str2double(seg(7))*10^9; %ns
end



ind=contains(DATA3,'DATA.3 CABL_SGN');
ourblock=DATA3(ind);
for i=1:size(ourblock,1)
    input_str = ourblock{i};
    seg = split(input_str);
    id = ismember(anttab(:,2),str2double(seg(5)));
    
    anttab(id,9) = str2double(seg(7));
end

anttab(:,10)=(anttab(:,7)+anttab(:,8)) .* anttab(:,9);
