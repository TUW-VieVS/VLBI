% Purpose  
%   Read twin.txt file to get the twin/multiple antennas at a site.
% History  
%   2014-03-20   Jing SUN   Created
%   


function [twin] = itwin(filename, station)

% initialize
twin.num = 0;
twin.sn(1,1) = 0;
twin.sn(1,2) = 0;
    
% open twin.txt file
fid = fopen(filename, 'r');
if (fid < 0)
    fprintf('No twin.txt file !\n');
    return;
end

% read twin.txt file
while ~feof(fid)
    line = fgetl(fid);
    linelength = length(line);
    if ~strcmp(line(1), ' ') 
        continue;
    end
    strcmp1 = sscanf(line(2:9), '%s', 1);
    strcmp2 = sscanf(line(11:18), '%s', 1);
    staid1 = 0;
    staid2 = 0;
    for ista = 1 : length(station)
        if strcmp(deblank(station(ista).name),deblank(strcmp1))
            staid1 = ista;
        end
        if strcmp(deblank(station(ista).name),deblank(strcmp2))
            staid2 = ista;
        end
    end
    if (staid1 > 0) & (staid2 > 0)
        twin.num = twin.num + 1;
        twin.sn(twin.num,1) = staid1;
        twin.sn(twin.num,2) = staid2;
    end
end

% close file
fclose(fid);


