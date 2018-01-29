% Purpose  
%   Write /VieVS/CRF/crf_sched.txt file.
% History  
%   2010-04-06   Jing SUN   Created
%    


function ocrf(source, sched, INFILE, PARA)

% open output file
fn_crf  = [PARA.pvievs 'CRF/' 'crf_sched.txt'];
fid_crf = fopen(fn_crf, 'w'); 

% open source.cat file
fid = fopen(INFILE.source, 'r');
if (fid < 0)
    fprintf('    no source catalog file %s !\n', INFILE.source);
    return;
end

% write crf
srcnum = length(source);
srcsn_s(1:srcnum) = 0;
srcnum_s = 0;
for isched = 1 : length(sched)
    for iscan = 1 : sched(isched).nscan
        srcid = sched(isched).scan(iscan).srcid;
        if (size(find(srcsn_s(1:srcnum)==srcid),2) == 0)
            srcnum_s = srcnum_s + 1 ;
            srcsn_s(srcnum_s) = srcid;
            frewind(fid);
            ifind = 0;
            while ~feof(fid)
                line = fgetl(fid);
                linelength = length(line);
                if ((linelength > 9) & strcmp(line(2:9), source(srcid).name(1:8)))
                    fprintf(fid_crf, '%s\n', line(1:linelength));
                    ifind = 1;
                    break;
                end           
            end  
            if(ifind == 0)
                fprintf('    no position for source %s in source.cat file !\n', source(srcid).name(1:8));
            end
        end
    end
end

% close file
fclose(fid_crf);
fclose(fid);


