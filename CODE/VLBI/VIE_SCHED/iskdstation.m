% Purpose  
%   Read $STATIONS information from skd file.
% History  
%   2013-05-06   Jing SUN   Created
%   2015-07-24   Lucia PLANK strtrim added for station name


function [station] = iskdstation(fid)

stanum = 0;
% A 
frewind(fid);
while ~feof(fid)
    line = fgetl(fid);
    linelength = length(line);
    if ((linelength>=9) & strcmp(line(1:9),'$STATIONS')) 
        break;
    end
end
while ~feof(fid)
    line = fgetl(fid);
    linelength = length(line);
    if (line(1) == '$')
        break;
    end
    if strcmp(line(1:2), 'A ')
        stanum = stanum + 1;
        station(stanum).id   = line(4);
        station(stanum).antname = line(6:13);
        station(stanum).axis = line(15:18);
        [tmprd, count, errmsg, nextindex] = sscanf(line(20:linelength), '%f', 10);
        station(stanum).offset = tmprd(1);
        station(stanum).rate1  = tmprd(2);        
        station(stanum).c1     = tmprd(3);      
        station(stanum).lim11  = tmprd(4); 
        station(stanum).lim12  = tmprd(5); 
        station(stanum).rate2  = tmprd(6);    
        station(stanum).c2     = tmprd(7);     
        station(stanum).lim21  = tmprd(8);  
        station(stanum).lim22  = tmprd(9);  
        station(stanum).diam   = tmprd(10);    
        clear tmprd;
        station(stanum).po(1:2) = ' ';
        station(stanum).eq(1:3) = ' ';
        station(stanum).ms(1:2) = ' ';
        index = 20 + nextindex - 1;
        [tmppo, count, errmsg, nextindex] = sscanf(line(index:linelength), '%s', 1);
        station(stanum).po(1:2) = tmppo(1:2);
        clear tmppo;
        index = index + nextindex - 1;
        [tmpeq, count, errmsg, nextindex] = sscanf(line(index:linelength), '%s', 1);
        station(stanum).eq(1:length(tmpeq)) = tmpeq(1:length(tmpeq));
        clear tmpeq;
        index = index + nextindex - 1;
        [tmpms, count, errmsg, nextindex] = sscanf(line(index:linelength), '%s', 1);
        if (count == 1)
            station(stanum).ms(1:2) = tmpms(1:2);
        end
        clear tmpms;
    end
end
% cable wrap
for ista = 1 : stanum
    az1 = station(ista).lim11;  
    az2 = station(ista).lim12;  
    station(ista).aznp = (az2-az1)/2 + az1;
    % neutral range
    az = station(ista).aznp;
    while ((az>=az1)&(az<=az2))
        az = az - 0.001;
        if ((az+360)<=az2) 
            break;
        end  
    end
    station(ista).azn1 = az + 0.001;
    az = station(ista).aznp;
    while ((az>=az1)&(az<=az2))
        az = az + 0.001;
        if ((az-360)>=az1)
            break;
        end  
    end
    station(ista).azn2 = az - 0.001;
    if (station(ista).azn1 < station(ista).lim11)
        station(ista).azn1 = station(ista).lim11;
    end
    if (station(ista).azn2 > station(ista).lim12)
        station(ista).azn2 = station(ista).lim12;
    end
    % clockwise range
    station(ista).azc1 = station(ista).azn2;
    station(ista).azc2 = az2;
    % count-clockwise range
    station(ista).azw1 = az1;
    station(ista).azw2 = station(ista).azn1;
end
% P 
frewind(fid);
while ~feof(fid)
    line = fgetl(fid);
    linelength = length(line);
    if ((linelength>=9) & strcmp(line(1:9),'$STATIONS')) 
        break;
    end
end
while ~feof(fid)
    line = fgetl(fid);
    linelength = length(line);
    if (line(1) == '$')
        break;
    end
    if strcmp(line(1:2), 'P ')
        [tmprd, count, errmsg, nextindex] = sscanf(line(16:linelength), '%f', 4);
        for ista = 1 : stanum
            if strcmp(line(4:5),station(ista).po)
                station(ista).name = strtrim(line(7:14));
                station(ista).xyz(1) = tmprd(1);
                station(ista).xyz(2) = tmprd(2);
                station(ista).xyz(3) = tmprd(3);
                station(ista).occ    = tmprd(4);
                break;
            end
        end
        clear tmprd;
    end
end
% T 
for ista = 1 : stanum
    station(ista).sefdpara(1,1:4) = 0.0;
    station(ista).sefdpara(2,1:4) = 0.0;
end
frewind(fid);
while ~feof(fid)
    line = fgetl(fid);
    linelength = length(line);
    if ((linelength>=9) & strcmp(line(1:9),'$STATIONS')) 
        break;
    end
end
while ~feof(fid)
    line = fgetl(fid);
    linelength = length(line);
    if (line(1) == '$')
        break;
    end
    if strcmp(line(1:2), 'T ')
        po(1:3) = ' ';
        if (line(5) == ' ')
            po(1:2) = line(6:7);
        else
            po(1:3) = line(5:7);
        end
        for ista = 1 : stanum
            if strcmp(po,station(ista).eq)
                break;
            end
        end
        maxindex = 0;
        [il] = findstr(line(34:linelength), ' X ');
        if (length(il) > 0)
            [tmprd, count, errmsg, nextindex] = sscanf(line(34+il(1)-1+2:linelength),'%f', 1);
            station(ista).sefdpara(1,1) = tmprd;
            maxindex = 34+il(1)-1+2 + nextindex - 1 ;
        end
        if (length(il) == 2)
            [tmprds, count, errmsg, nextindex] = sscanf(line(34+il(2)-1+2:linelength), '%f', 3);
            station(ista).sefdpara(1,2:4) = tmprds(1:3);
            maxindex = 34+il(2)-1+2 + nextindex - 1 ;
        end
        [il] = findstr(line(34:linelength), ' S ');
        if (length(il) > 0)
            [tmprd, count, errmsg, nextindex] = sscanf(line(34+il(1)-1+2:linelength),'%f', 1);
            station(ista).sefdpara(2,1) = tmprd;
            if ((34+il(1)-1+2 + nextindex - 1) > maxindex)
                maxindex = 34+il(1)-1+2 + nextindex - 1 ;
            end
        end
        if (length(il) == 2)
            [tmprds, count, errmsg, nextindex] = sscanf(line(34+il(2)-1+2:linelength), '%f', 3);
            station(ista).sefdpara(2,2:4) = tmprds(1:3);
            if (34+il(2)-1+2 + nextindex - 1) > maxindex
                maxindex = 34+il(2)-1+2 + nextindex - 1 ;
            end
        end
        [tmps, count, errmsg, nextindex] = sscanf(line(maxindex:linelength), '%s', 1);
        station(ista).rack  = tmps;
        index = maxindex + nextindex - 1;
        [tmps] = sscanf(line(index:linelength), '%s', 1);
        station(ista).equip = tmps;
        clear tmprd;
        clear tmprds;
        clear tmps;
    end % if strcmp(line(1:2), 'T ')
end
% H 
frewind(fid);
while ~feof(fid)
    line = fgetl(fid);
    linelength = length(line);
    if ((linelength>=9) & strcmp(line(1:9),'$STATIONS')) 
        break;
    end
end
while ~feof(fid)
    line = fgetl(fid);
    linelength = length(line);
    if (line(1) == '$')
        break;
    end
    if strcmp(line(1:2), 'H ')
        [tmprds, count, errmsg, nextindex] = sscanf(line(5:linelength), '%f');
        for ista = 1 : stanum
            if strcmp(line(3:4),station(ista).ms)
                break;
            end
        end
        station(ista).hmasknum              = length(tmprds);
        station(ista).hmask(1:station(ista).hmasknum) = tmprds(1:length(tmprds)); 
        clear tmprds;
    end
end


