% Purpose
%   Write $CODES (frequency sequences and station LOs).
% History
%   2010-04-06   Jing SUN   Created
%   2014-10-29   Jamie McCallum new tracks mode for Mk4
%   2015-05-01   David Mayer --> bug fixed: the comparison of record types
%   is no longer case sensitive
%   2015-07-22   Lucia Plank --> major revision for a more proper usage of
%   the information provided in the catalog files (different modes allowed)
%   2015-09-11   David Mayer and Andreas Hellerschmied --> Bug fix station indices
%   2016-05-19   Matthias Schartner --> Bugfix for comments with whitespace
%   2016-09-05   Lucia Plank --> fanout from tracksfile


function [obsmode] = ocodes(station, obsmode, INFILE, fid_skd, PARA)

% number of stations
% stanum = length(station);
% Work out how many CODES need to be written out.

trackmodes=unique([station.trackmode]);
numcsec=size(trackmodes,2);


for i=1:numcsec
    % the F line
    cstnum=strcmp(trackmodes(i),[station.trackmode]);
    stnum=find(cstnum);
    
    fprintf(fid_skd, 'F');
    fprintf(fid_skd, ' %s', obsmode.freqname);
    fprintf(fid_skd, ' %s   ', [obsmode.sx{1}]);
    for ista = 1 : length(stnum)
        sta_id = stnum(ista);
        fprintf(fid_skd, ' %s', deblank(station(sta_id).name(1:8)));
    end
    fprintf(fid_skd, '\n');
    
    % get tracks information
    [tracks,fanout]=otrack(INFILE.tracks,trackmodes(i));
    
    for ic=1:length(obsmode.freq)
        sx     = obsmode.sx{1};
        subgp  = char(obsmode.freq{ic,1});
        freq   = char(obsmode.freq{ic,3});
        pcal   = char(obsmode.freq{ic,7});
        ch     = char(tracks{ic,1});
        if length(ch)==1
            ch=[ch,' '];
        end
        tmode  = [PARA.TRACKSMODE,'1:'];
        dum    = char(trackmodes(i));
        %fanout = dum(7);
        bw     = obsmode.bw;
        tracki = char(tracks{ic,2});
        fprintf(fid_skd, 'C %s %s  %s  %s   %s %s%s %7.2f %s\n', sx, subgp,freq,pcal,ch,tmode,fanout, bw, tracki);
    end
end
% the R line  : the sample rate
fprintf(fid_skd, 'R %s  %7.3f\n', char(obsmode.sx), obsmode.samprate);

% the L lines : the LO setups
% get the station loif mode from the rx.cat file
fid = fopen(INFILE.rx, 'r');
while ~feof(fid)
    line = fgetl(fid);
    linelength = length(line);
    if ((linelength > 1) & (line(1) == ' '))
        [rxname] = sscanf(line, '%s', 1);
        if strcmp(rxname, obsmode.rxname)
            while true
                line = fgetl(fid);
                linelength = length(line);
                if ((linelength < 2) || isempty(strfind(line(2:3),'-')) )
                    break;
                end
                rname(1:8) = line(10:17);
                loifname = sscanf(line(22:linelength), '%s', 1);
                for ista = 1 : length(station)
                    if strcmp(rname(1:8),station(ista).name(1:8))
                        station(ista).staloif = loifname;
                        break;
                    end
                end
            end
        end
    end
end
fclose(fid);

% get loif information per station and write to skd file
fid = fopen(INFILE.loif, 'r');
for ista = 1 : length(station)
    frewind(fid);
    while ~feof(fid)
        line = fgetl(fid);
        linelength = length(line);
        [loifname] = sscanf(line, '%s', 1);
        if strcmp(loifname, station(ista).staloif)
            while true
                line = fgetl(fid);
                linelength = length(line);
                if ((linelength < 2) | (line(2) ~= '-'))
                    break;
                end
                remain = line;
                loif=textscan(remain,'%s %s %s %s %s %s');
                
                %                 [str1, remain] = strtok(remain,' ');
                %                 [str2, remain] = strtok(remain,' ');
                %                 [str3, remain] = strtok(remain,' ');
                %                 [str4, remain] = strtok(remain,' ');
                %                 [str5, remain] = strtok(remain,' ');
                %                 [str6, remain] = strtok(remain,' ');
                %                 str55(1:10) = ' ';
                %                 str55(1:length(str5)) = str5;
                %                 str22(1:2) = ' ';
                %                 str22(1:length(str2)) = str2;
                %                 str33(1:2) = ' ';
                %                 str33(1:length(str3)) = str3;
                fprintf(fid_skd, 'L %s %s', station(ista).id, char(obsmode.sx));
                fprintf(fid_skd, ' %s  %s %s %s %s\n', char(loif{4}), char(loif{3}), char(loif{5}), char(loif{2}), char(loif{6}));
            end
        end
    end
end
fclose(fid);

% the B line  : barrel roll
fprintf(fid_skd, 'B SX\n');


