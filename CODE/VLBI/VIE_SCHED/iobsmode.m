% Purpose  
%   Read observing mode information from the modes catalog file.
% History  
%   2010-04-06   Jing SUN   Created
%   2016-02-24   Lucia PLANK enabling modes with names <16
%   


function [obsmode] = iobsmode(obsmodename, INFILE, PARA)

% initialize obsmode information struct
obsmode.name = obsmodename(1:end); 
while length(obsmode.name)<16
        obsmode.name=[obsmode.name,' '];
end
obsmode.freqname   = ' ';
obsmode.bw         = 0.0;
obsmode.samprate   = 0.0;
obsmode.recname    = ' ';
obsmode.channelnum = 0.0;
obsmode.fanout     = 0.0;
obsmode.bits       = 0.0;
obsmode.chanumband(1:PARA.MAX_BANDNUM) = 0;

% open modes catalog file
fprintf(PARA.fid_header,'3.1 read modes file: %s \n', INFILE.modes);
fid = fopen(INFILE.modes, 'r');
if (fid < 0)
    error('    no modes catalog file %s !\n', INFILE.modes);
end

% read the modes information
ifmodes = 0;
while ~feof(fid)
    line = fgetl(fid);
    linelength = length(line);
    if ((linelength > 16) & strcmp(line(1:16),obsmode.name(1:16)))
        ifmodes = 1;       
        [obsmode.freqname, count, errmsg, nextindex] = sscanf(line(18:linelength), '%s', 1);
        index = 18 + nextindex - 1;      
        [obsmode.bw, count, errmsg, nextindex] = sscanf(line(index:linelength), '%f', 1);
        index = index + nextindex - 1;     
        [obsmode.samprate, count, errmsg, nextindex] = sscanf(line(index:linelength), '%f', 1);
        index = index + nextindex - 1;
        [obsmode.recname, count, errmsg, nextindex] = sscanf(line(index:linelength), '%s', 1);
        [il] = findstr(obsmode.recname, '-');    %%trk-chan-fan-bit%%
        obsmode.channelnum = sscanf(obsmode.recname((il(1)+1):(il(2)-1)), '%f', 1);
        obsmode.fanout     = sscanf(obsmode.recname((il(2)+1):(il(3)-1)), '%f', 1);
        obsmode.bits       = sscanf(obsmode.recname((il(3)+1):length(obsmode.recname)), '%f', 1);
        fprintf(PARA.fid_header,'    %s\n', line(1:linelength));
        break;
    end
end
if (ifmodes == 0)
    error('    There is no modes information for obsmode %s \n', obsmode.name(1:16));
end

% close file
fclose(fid);

% number of channel at each band
for iband = 1 : PARA.MAX_BANDNUM
    obsmode.chanumband(iband) = PARA.CHANUM(iband);
end


