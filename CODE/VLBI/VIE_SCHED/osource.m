% Purpose  
%   Write $SOURCES (list of sources for this experiment).
% History  
%   2010-04-06   Jing SUN   Created
%   2016-06-11	 Matthia Schartner: now also searches source coordinates in source_star.cat if it is not in source.cat


function osource(source, sched, filename, starfilename, fid_skd)

% open source catalog file 
ifid = 1;
fid = fopen(filename, 'r');
fid2 = fopen(starfilename, 'r');
if (fid < 0)
    ifid = 0;
end

% write $SOURCES
srcnum = length(source);
srcsn_s(1:srcnum) = 0;
srcnum_s = 0;
for isched = 1 : length(sched)
    for iscan = 1 : sched(isched).nscan
        srcid = sched(isched).scan(iscan).srcid;
        if (size(find(srcsn_s(1:srcnum)==srcid),2) == 0)
            srcnum_s = srcnum_s + 1 ;
            srcsn_s(srcnum_s) = srcid;
            
            if (ifid == 0)
                fprintf(fid_skd, ' %s $         %02d %02d %9.6f     %s%02d %02d %8.5f %s\n', source(srcid).name(1:8), source(srcid).rah, source(srcid).ram, source(srcid).ras, source(srcid).sign, source(srcid).ded, source(srcid).dem, source(srcid).des, 'external file');
            elseif (ifid == 1)
                ifile = 0;
                frewind(fid);
                while ~feof(fid)
                    line = fgetl(fid);
                    linelength = length(line);
                    if ((linelength > 9) & strcmp(line(2:9), source(srcid).name(1:8)))
                        fprintf(fid_skd, '%s\n', line(1:linelength));
                        ifile = 1;
                        break;
                    end           
                end
                if ifile == 0
                    frewind(fid2);
                    while ~feof(fid2)
                        line = fgetl(fid2);
                        linelength = length(line);
                        if ((linelength > 9) & strcmp(line(2:9), source(srcid).name(1:8)))
                            fprintf(fid_skd, '%s\n', line(1:linelength));
                            ifile = 1;
                            break;
                        end           
                    end
                end
                if (ifile == 0)
                    fprintf(fid_skd, ' %s $         %02d %02d %9.6f     %s%02d %02d %8.5f %s\n', source(srcid).name(1:8), source(srcid).rah, source(srcid).ram, source(srcid).ras, source(srcid).sign, source(srcid).ded, source(srcid).dem, source(srcid).des, 'external file');
                end  
            end   %ifid  
            
        end
    end
end

% close file
if (fid >= 0)
    fclose(fid);
end


