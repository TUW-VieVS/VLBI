% Purpose  
%   Read station equip information from the equip catalog file.
% History  
%   2010-04-06   Jing SUN   Created
%
% Changes
%   2013-11-14   Hellerschmied Andreas: Exclude SEFD check for satellite scheduling, read new parameter: station.equip_rack
%   

function [station] = iequip(filename, station, PARA)

% open equip catalog file
fid = fopen(filename, 'r');
if (fid < 0)
    error('    no equip catalog file %s !\n', filename);
end

% read the equip information
stanum = length(station);
for ista = 1 : stanum
    frewind(fid);
    while ~feof(fid)
        line = fgetl(fid);
        linelength = length(line);
        if ((linelength < 9) | (~strcmp(line(1), ' ')))   %%%%%
            continue;
        end
        eqstr = sscanf(line(12:14), '%s', 1);
        if (strcmp(line(2:9), station(ista).name(1:8))) %strcmp(eqstr, deblank(station(ista).eq))
            % read SEFD parameters 
            for ib = 1 : PARA.MAX_BANDNUM
                bandfind(1) = ' ';
                bandfind(2) = PARA.BAND(ib);
                bandfind(3) = ' ';
                [il] = findstr(line(10:linelength), bandfind);
                if (length(il) > 0)
                    station(ista).sefdpara(ib,1) = sscanf(line(10+il(1)-1+2:linelength),'%f', 1);
                    fprintf(PARA.fid_header,'    %s%s%d', station(ista).name(1:8), bandfind(1:3), station(ista).sefdpara(ib,1));
                end
                if (length(il) == 2)
                    station(ista).sefdpara(ib,2:4) = sscanf(line(10+il(2)-1+2:linelength), '%f', 3);
                    fprintf(PARA.fid_header,'     %10.5f%10.5f%10.5f', station(ista).sefdpara(ib,2:4));
                end
                fprintf(PARA.fid_header,'\n');
            end 
            % read equip name (recorder)
            ln = 0; 
            sn = '';
            for ill = linelength : -1 : 1
                if ((ln>0) && (line(ill) == ' '))
                    break;
                end
                if(line(ill) == ' ')
                    continue;
                else
                    ln = ln + 1;
                    sn(ln) = line(ill);
                end   
            end
            for ill = 1 : ln
                station(ista).equip(ill) = sn(ln-(ill-1));
            end
            
            % read equip name (rack)
            sn = '';
            ln1 = 0; 
            fprintf(PARA.fid_header,'%s',line(linelength - ln-2));
            for ill = (linelength - ln) : -1 : 1
                if ((ln1>0) && (line(ill) == ' '))
                    break;
                end
                if(line(ill) == ' ')
                    continue;
                else
                    ln1 = ln1 + 1;
                    sn(ln1) = line(ill);
                end   
            end
            for ill = 1 : ln1
                station(ista).equip_rack(ill) = sn(ln1-(ill-1));
            end
        end
    end
    % check SEFD parameters
    ifsefd = 1;
    for ib = 1 : PARA.MAX_BANDNUM
        if (station(ista).sefdpara(ib,1) < 1.0d-3)
            ifsefd = 0;
        end
    end
    if ((ifsefd == 0) && (PARA.OPTIMIZATION ~= 3) ) % Hellerschmied: ..&& (PARA.OPTIMIZATION ~= 3) added!
        error('    There is no enough equip information for station %s \n', station(ista).name);
    end
end

% close file
fclose(fid);


