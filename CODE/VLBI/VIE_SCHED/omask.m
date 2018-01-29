% Purpose  
%   Write $STATIONS (list of stations in this experiment).
%          the H lines : horizon mask
% History  
%   2010-04-06   Jing SUN   Created
%   


function omask(staname, filename, fid_skd)

% open mask catalog file
fid = fopen(filename, 'r');

% check and write mask information
stanum = size(staname,1);
for ista = 1 : stanum
    frewind(fid);
    while ~feof(fid)
        line = fgetl(fid);
        linelength = length(line);
        if ((linelength < 2) | (~strcmp(line(1:2), ' H')))   %%%%%
            continue;
        end
        [hc, count, errmsg, nextindex] = sscanf(line(1:linelength), '%s', 1);
        index = nextindex;
        [sta, count, errmsg, nextindex] = sscanf(line(index+1:linelength), '%s', 1);
        index = index + nextindex;
        if ~strcmp(sta, deblank(staname(ista,1:8)))
            continue;
        end
        fprintf(fid_skd, 'H %s', line(index+1:linelength));  
        while ~feof(fid)
            line = fgetl(fid);
            linelength = length(line);
            if ((linelength < 2) | (~strcmp(line(1:2), ' -')))
                break;
            end
            fprintf(fid_skd, '%s', line(3:linelength)); 
        end
        fprintf(fid_skd, '\n');
    end
end

% close file
fclose(fid);


