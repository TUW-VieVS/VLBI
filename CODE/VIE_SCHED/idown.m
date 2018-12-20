% Purpose  
%   Read down.txt file.
% History  
%   2010-04-06   Jing SUN   Created
%   


function [station] = idown(filename, station, PARA)

% open down.txt file
fid = fopen(filename, 'r');
if (fid < 0)
    fprintf(PARA.fid_header,'    no down.txt file !\n');
    return;
end

% read the downtime information
stanum = length(station);
for ista = 1 : stanum
    frewind(fid);
    while ~feof(fid)
        line = fgetl(fid);
        linelength = length(line);
        if ((linelength < 9) | (~strcmp(line(1), ' ')))   %%%%%
            continue;
        end
        if ~strcmp(line(2:9), station(ista).name(1:8))
            continue;
        end
        [tmpd] = sscanf(line(13:linelength), '%d', 12);
        station(ista).downum = station(ista).downum + 1;
        station(ista).downstart(station(ista).downum) = tmjd(tmpd(1),tmpd(2),tmpd(3),tmpd(4),tmpd(5),tmpd(6));
        station(ista).downend(station(ista).downum)   = tmjd(tmpd(7),tmpd(8),tmpd(9),tmpd(10),tmpd(11),tmpd(12));
        fprintf(PARA.fid_header,'    There is downtime information for station %s : %s\n', station(ista).name, line(13:linelength));
    end
end

% close file
fclose(fid);


