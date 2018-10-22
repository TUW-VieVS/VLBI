% input: fstation -> SINEX file VLBI-DISCONT.txt from
% http://vlbi.geod.uni-bonn.de/IVS-AC/

% output: structure antbr
% antbr.name - 8 char antenna name
% antbr.name - id code of antenna
% antbr.breaks - vector with breaks in mjd


% Coded by Hana Spicakova
% 31 March 2010
%%

fstation = 'VLBI-DISCONT_tsuk.txt';  %VLBI antenna/breaks list
fid = fopen(fstation);
i=1;j=1;
while ~feof(fid)
    tline = fgetl(fid);
    if strcmp('+SITE/ID',cellstr(tline));
        tline=fgetl(fid);
        while strcmp('-SITE/ID',cellstr(tline))==0
            A=tline';
            if (A(1)~='*')
                id_ref(i,:)=A(2:5);
                name(i,:)=A(22:30);
                i=i+1;
            end
            tline=fgetl(fid);
        end
    end
    
    if strcmp('+SOLUTION/DISCONTINUITY',cellstr(tline));
        tline=fgetl(fid);
        while strcmp('-SOLUTION/DISCONTINUITY',cellstr(tline))==0
            A=tline';
            if (A(1)~='*')
                id_br(j,:)=A(2:5);
                br(j,:)=A(17:28);
                j=j+1;
            end
            tline=fgetl(fid);
        end
    end

end
fclose(fid);

for i=1:size(br,1)
    if br(i,1)=='7' || br(i,1)=='8' || br(i,1)=='9'
        jj=1900;
    elseif br(i,1)=='0' || br(i,1)=='1' || br(i,1)=='2'
        jj=2000;
    end
    yr=str2num(br(i,1:2)) + jj;
    doy=str2num(br(i,4:6));
    jd=2400000.5;
    if doy~=0
        jd=doy2jd(yr,doy);
    end
    br_mjd(i)=jd-2400000.5;
end



for k=1:size(name,1)
    antbr(k).name=name(k,:);
    antbr(k).id=id_ref(k,:);
    antbr(k).break=[];
end


for j = 1: size(id_br,1)
    for k = 1:length(antbr)
        if strcmp(id_br(j,:),antbr(k).id)
            lbr=length(antbr(k).break);
            antbr(k).break(lbr+1) = br_mjd(j);
        end
    end
end


save VLBI_discont antbr
