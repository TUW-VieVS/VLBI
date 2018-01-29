% Purpose  
%   (1) Read input information (param, station, source, observing mode).
%   (2) Do initial calculations.
% History  
%   2010-04-06   Jing SUN   Created
%   ....
%   ..no information available..
%   ....
%   2015-02-06 A. Hellerschmied: Added option to delete source.mat vis GUI
%       option.
%   2015-03-04 A. Hellerschmied: Observation mode name now defined in
%       param.txt (not hard-coded any more).
%   2016-01-01 M. Schartner, L. Plank: STAR mode added


function [station, twin, tagalongname, source5, obsmode, srcat, catpair, PARA] = vie_sched_i(staname, INFILE, PARA)

% station number
stanum = size(staname, 1);

% read param.txt file
fprintf(PARA.fid_header,'0. read param.txt file\n');
[PARA] = iparam(INFILE.param, stanum, PARA);

% read station information
fprintf(PARA.fid_header,'1. read station information \n');
[station] = istation(staname, INFILE, PARA);

fprintf(PARA.fid_header,'Number of stations to be scheduled is %d\n', length(station));
% check the twin/multiple antennas at a site
[twin] = itwin(INFILE.twin, station);
if (twin.num ~= 0)
    for i = 1 : twin.num
        fprintf(PARA.fid_header,'Twin telescopes :  %s --- %s\n', station(twin.sn(i,1)).name, station(twin.sn(i,2)).name);
    end
end
% check the station for tag-along mode
[tagalongname] = itagalong(INFILE.tagalong,PARA);


% Delete source.mat (if it exists) if this option is selected in the GUI 
if PARA.USE_NEW_SOURCE_File
    if (exist([PARA.pfolder 'source.mat']) > 0)
        delete([PARA.pfolder 'source.mat']);
        fprintf(PARA.fid_header,' => /DATE/LEVEL5/source.mat deleted!\n');
    end
end


% read source information
if (exist([PARA.pfolder 'source.mat']) > 0)      % load source mat file
    fprintf(PARA.fid_header,'2. load source information (from source.mat) \n');
    load ([PARA.pfolder 'source.mat'], 'source');
else                                             % creat source mat file
    fprintf(PARA.fid_header,'2. read source information (from source.cat) and create source.mat \n');
    [source] = isource(INFILE, PARA);
end

% *** STAR mode ***
if PARA.STARMODE == 1
    source = isource_star(INFILE, PARA, source);
end


[source] = isweight(INFILE.sweight, source);     % David Mayer, 2014 Aug 13
source2 = source;
srcnum2 = length(source2);
fprintf(PARA.fid_header,'  number of sources with complete flux parameters is %d\n', srcnum2);
% read file psource.txt and check particular sources
fprintf(PARA.fid_header,'  read file: %s \n', INFILE.psource);
[source2] = ipsource(INFILE.psource, source2, PARA);
% check maximum flux
srcnum3 = 0;

for isrc = 1 : srcnum2
    % *** STAR mode *** skip for star source
    if source(isrc).star==1 
        srcnum3 = srcnum3 + 1;
        source3(srcnum3) = source2(isrc);
        continue;
    end
    
    for iband = 1 : PARA.MAX_BANDNUM
        if (round(source2(isrc).fluxpara(iband,PARA.MAX_FLUXPARA)) == 1)
            fluxmax(iband) = source2(isrc).fluxpara(iband,2);
        elseif (round(source2(isrc).fluxpara(iband,PARA.MAX_FLUXPARA)) == 2)
            fluxmax(iband) = source2(isrc).fluxpara(iband,1);
        end
    end
    if ((min(fluxmax(1:PARA.MAX_BANDNUM)) >= PARA.MIN_FLUX) | (source2(isrc).ifp==1))
        srcnum3 = srcnum3 + 1;
        source3(srcnum3) = source2(isrc);
    else
        fprintf(PARA.fid_header,'  weak source : %s\n', source2(isrc).name);
    end
end
fprintf(PARA.fid_header,'  number of weak sources is %d\n', (srcnum2-srcnum3));
% check if the source is too close to the sun at the start and end of session
srcnum4 = 0;
[sunras, sundes] = ssunpo(PARA.startmjd);
[sunrae, sundee] = ssunpo(PARA.endmjd);
for isrc = 1 : srcnum3
    srcra = source3(isrc).ra;
    srcde = source3(isrc).de;
    cd2 = cos(srcde);
    sd2 = sin(srcde);
    % compute distance of the source from the sun
    cra = cos(sunras - srcra);
    cd1 = cos(sundes);
    sd1 = sin(sundes);
    arg = cd1 * cd2 * cra + sd1 * sd2;
    arcrs = atan2(sqrt(1-arg*arg), arg);
    cra = cos(sunrae - srcra);
    cd1 = cos(sundee);
    sd1 = sin(sundee);
    arg = cd1 * cd2 * cra + sd1 * sd2;
    arcre = atan2(sqrt(1-arg*arg), arg);
    if (((arcrs>=PARA.MIN_SUNDIST) & (arcre>=PARA.MIN_SUNDIST)) | (source3(isrc).ifp==1))
        srcnum4 = srcnum4 + 1;
        source4(srcnum4) = source3(isrc); 
    else
        fprintf(PARA.fid_header,'  source close to the Sun : %s %4.1f %4.1f\n', source3(isrc).name, arcrs*180/pi, arcre*180/pi);
    end
end
fprintf(PARA.fid_header,'  number of sources close to the Sun is %d\n', (srcnum3-srcnum4));
% sort source by ra
for isrc = 1 : srcnum4
    ra(isrc) = source4(isrc).ra;
end
[x, y] = sort(ra);
for isrc = 1 : srcnum4
    source5(isrc) = source4(y(isrc));
end
srcnum = length(source5);
if (srcnum > 0)
    fprintf(PARA.fid_header,'Number of sources to be scheduled is %d\n', srcnum);
else
    error(PARA.fid_header,'No source to be scheduled !\n');
end

% read observing mode information
fprintf(PARA.fid_header,'3. read observing mode information \n');
[obsmode] = iobsmode(PARA.OBSMODE_NAME, INFILE, PARA);
fprintf(PARA.fid_header,'  number of channels is %d (%d', obsmode.channelnum, obsmode.chanumband(1));
for iband = 2 : PARA.MAX_BANDNUM
    fprintf(PARA.fid_header,' + %d', obsmode.chanumband(iband));
end
fprintf(PARA.fid_header,')\n');
fprintf(PARA.fid_header,'  channel bandwidth is %d [MHz]\n', obsmode.bw);
fprintf(PARA.fid_header,'  sample rate is %d [MHz]\n',       obsmode.samprate);
fprintf(PARA.fid_header,'  quantification is %d bits\n',     obsmode.bits);

% get the distribution of sources for source-based strategy
if (PARA.OPTIMIZATION == 1)
    fprintf(PARA.fid_header,'4. calculate the distribution of sources\n');
    [cat, srcat, catpair] = isrcseq(source5, PARA);
else
    srcat   = 0.0;
    catpair = 0.0;
end


