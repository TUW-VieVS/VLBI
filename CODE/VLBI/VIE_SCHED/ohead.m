% Purpose  
%   Write $HEAD (tape recorder head positions).
% History  
%   2010-04-06   Jing SUN   Created
%   2015-07-22   Lucia PLANK small revision


function ohead(station, obsmode, INFILE, fid_skd)

stanum = length(station);
for ista = 1 : stanum
    stahdpos(ista,1:8) = ' ';
end 
% open rec catalog file
fid = fopen(INFILE.rec, 'r');
while ~feof(fid)
    line = fgetl(fid);
    linelength = length(line);
    if ((linelength > 1) & (line(1) == ' '))
        [recname] = sscanf(line, '%s', 1);  
        if strcmp(recname, obsmode.recname)
            while true
                line = fgetl(fid);
                linelength = length(line);
                if ((linelength < 2) | (line(2) ~= '-'))
                    break;
                end
                rname(1:8) = line(10:17);
                hdposname = sscanf(line(21:linelength), '%s', 1);
                for ista = 1 : stanum
                    if strcmp(rname(1:8),station(ista).name(1:8))
                        stahdpos(ista,1:length(hdposname)) = hdposname;
                        break;
                    end
                end              
            end
        end
    end
end
fclose(fid);

% open hdpos catalog file
fid = fopen(INFILE.hdpos, 'r');
for ista = 1 : stanum
    frewind(fid);
    while ~feof(fid)
        line = fgetl(fid);
        linelength = length(line);
        if ((linelength > 2) & (line(1) == ' ') & (line(2) ~= '-'))
            [hname] = sscanf(line, '%s', 1);  
            if strcmp(hname, deblank(stahdpos(ista,:)))
                while true
                    line = fgetl(fid);
                    linelength = length(line);
                    if ((linelength < 2) | (line(2) ~= '-'))
                        break;
                    end
                    fprintf(fid_skd, '%s %s   %s\n', station(ista).id, char(obsmode.sx), line(6:linelength));
                end              
            end
        end
    end
end
fclose(fid);


