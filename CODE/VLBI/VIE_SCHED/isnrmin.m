% Purpose  
%   Read snrmin.txt file.
% History  
%   2010-04-06   Jing SUN   Created
%   


function [station] = isnrmin(filename, station, PARA)

% open snrmin.txt file
fid = fopen(filename, 'r');
if (fid < 0)
    fprintf(PARA.fid_header,'    no snrmin.txt file !\n');
    return;
end

% read the minimum snr information
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
        [station(ista).minsnr(1:PARA.MAX_BANDNUM)] = sscanf(line(10:linelength), '%d', PARA.MAX_BANDNUM);
        fprintf(PARA.fid_header,'    There is new SNR_min information for station %s :  ', station(ista).name);
        for iband = 1 : PARA.MAX_BANDNUM
            fprintf(PARA.fid_header,'%d ', station(ista).minsnr(iband));
        end
        fprintf(PARA.fid_header,'\n');
    end
end

% close file
fclose(fid);


