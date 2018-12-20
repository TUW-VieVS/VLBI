% Purpose  
%   Write $STATIONS (list of stations in this experiment).
%          the A lines : antenna limits and rates
% History  
%   2010-04-06   Jing SUN   Created
%   


function oantenna(staname, filename, fid_skd)

% open antenna catalog file
fid = fopen(filename, 'r');

% check and write antenna information
stanum = size(staname,1);
for ista = 1 : stanum
    frewind(fid);
    while ~feof(fid)
        line = fgetl(fid);
        linelength = length(line);
        if ((linelength > 76) & (line(1) == ' ') & strcmp(line(4:11),staname(ista,1:8)))
            fprintf(fid_skd, 'A %s\n', line(1:linelength));
        end
    end
end

% close file
fclose(fid);


