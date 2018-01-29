% Purpose  
%   Read station position from the position catalog file.
% History  
%   2010-04-06   Jing SUN   Created
%   


function [station] = iposition(filename, station, PARA)

% open position catalog file
fid = fopen(filename, 'r');
if (fid < 0)
    error('    no position catalog file %s !\n', filename);
end

% read the position information
stanum = length(station);
for ista = 1 : stanum
    frewind(fid);
    ifposition = 0;
    while ~feof(fid)
        line = fgetl(fid);
        linelength = length(line);
        if ((linelength < 12) | (~strcmp(line(1), ' ')))   %%%%%
            continue;
        end
        [tmps1, count, errmsg, nextindex] = sscanf(line(1:linelength), '%s', 1);
        [tmps2, count, errmsg, nextindex] = sscanf(line(nextindex:linelength), '%s', 1);
        if (strcmp(tmps2, deblank(station(ista).name(1:8)))) %%%strcmp(line(2:3), station(ista).po(1:2))
            ifposition = 1;
            station(ista).xyz = sscanf(line(14:linelength), '%f', 3);
            pos(1,1) = station(ista).xyz(1);
            pos(1,2) = station(ista).xyz(2);
            pos(1,3) = station(ista).xyz(3);
            [station(ista).llh(2),station(ista).llh(1),station(ista).llh(3)] = xyz2ell(pos);
            fprintf(PARA.fid_header,'    %s\n', line(1:linelength));
            break;
        end
    end
    if (ifposition == 0)
        error('    There is no position information for station %s \n', station(ista).name);
    end
end

% close file
fclose(fid);


