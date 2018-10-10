% Purpose  
%   Read psource.txt file.
% History  
%   2010-04-06   Jing SUN   Created
%   


function [source] = ipsource(filename, source, PARA)

% open psource.txt file
fid = fopen(filename, 'r');
if (fid < 0)
    fprintf(PARA.fid_header,'    no psource.txt file !\n');
    return;
end

% read the list of particular sources
fprintf(PARA.fid_header,'    list of particular sources : \n');
srcnum = length(source);
while ~feof(fid)
    line = fgetl(fid);
    linelength = length(line);
    if ((linelength < 9) | (~strcmp(line(1), ' ')))   %%%%%
        continue;
    end
    srcname = line(2:9);
    fprintf(PARA.fid_header,'    %s ', srcname(1:8));
    ifuse = 0;
    for isrc = 1 : srcnum
        if strcmp(srcname, source(isrc).name(1:8))
            source(isrc).ifp = 1;
            yy  = str2num(line(11:14));
            mm  = str2num(line(15:16));
            dd  = str2num(line(17:18));
            hh  = str2num(line(19:20));
            min = str2num(line(21:22));
            sec = str2num(line(23:24));
            source(isrc).pst = tmjd(yy,mm,dd,hh,min,sec);
            source(isrc).pdt = str2num(line(25:28)); %[min]
            ifuse = 1;
            break;
        end
    end
    if (ifuse == 1)
        fprintf(PARA.fid_header,'with complete information\n');
    elseif (ifuse == 0)
        fprintf(PARA.fid_header,'without complete information\n');
    end
end

% close file
fclose(fid);


