% Purpose  
%   Read station mask information from the mask catalog file.
% History  
%   2010-04-06   Jing SUN   Created
%   


function [station] = imask(filename, station, PARA)

% open mask catalog file
fid = fopen(filename, 'r');
if (fid < 0)
    fprintf(PARA.fid_header,'    no mask catalog file %s !\n', filename);
    return;
end

% read the mask information
stanum = length(station);
for ista = 1 : stanum
    hmasknum = 0; 
    hmask(1:PARA.MAX_HOR) = 0.0;
    frewind(fid);
    while ~feof(fid)
        line = fgetl(fid);
        linelength = length(line);
        if ((linelength < 2) | (~strcmp(line(1:2), ' H')))   %%%%%
            continue;
        end
        [hc, count, errmsg, nextindex] = sscanf(line(1:linelength), '%s', 1);
        index = nextindex;
        [staname, count, errmsg, nextindex] = sscanf(line(index:linelength), '%s', 1);
        index = index + nextindex - 1;
        [mstr, count, errmsg, nextindex] = sscanf(line(index:linelength), '%s', 1);
        index = index + nextindex - 1;
        [vector, count] = sscanf(line(index:linelength), '%f');
        if  (strcmp(staname, deblank(station(ista).name(1:8)))) %%%strcmp(mstr, station(ista).ms)
            hmask(1:count) = vector(1:count);
            hmasknum = count;  
            while ~feof(fid)
                line = fgetl(fid);
                linelength = length(line);
                if ((linelength < 2) | (~strcmp(line(1:2), ' -')))
                    break;
                end
                [vector, count] = sscanf(line(4:linelength), '%f');
                hmask(hmasknum+1 : hmasknum+count) = vector(1:count);
                hmasknum = hmasknum + count;  
            end
            station(ista).hmasknum = hmasknum;
            station(ista).hmask(1:hmasknum) = hmask(1:hmasknum) * pi / 180.0;   % deg-->rad  
            fprintf(PARA.fid_header,'    There is mask information for station %s : ', station(ista).name(1:8));
            for i = 1 : station(ista).hmasknum
                fprintf(PARA.fid_header,'%5.1f ', station(ista).hmask(i)*180/pi);
            end
            fprintf(PARA.fid_header,'\n');
        end
    end
end

% close file
fclose(fid);


