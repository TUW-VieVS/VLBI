% Purpose  
%   Write $STATIONS (list of stations in this experiment).
%          the P lines : station position iformation
% History  
%   2010-04-06   Jing SUN   Created
%   


function oposition(staname, filename, fid_skd)

% open position catalog file
fid = fopen(filename, 'r');

% check and write position information
stanum = size(staname,1);
for ista = 1 : stanum
    frewind(fid);
    while ~feof(fid)
        line = fgetl(fid);
        linelength = length(line);
        [tmps1, count, errmsg, nextindex] = sscanf(line(1:linelength), '%s', 1);
        [tmps2, count, errmsg, nextindex] = sscanf(line(nextindex:linelength), '%s', 1);
        if ((linelength > 12) & (line(1) == ' ') & strcmp(tmps2,deblank(staname(ista,1:8))))
            fprintf(fid_skd, 'P %s\n', line(1:linelength));
        end
    end
end

% close file
fclose(fid);


