% Purpose  
%   Write $FLUX (flux densities for each source).
% History  
%   2010-04-06   Jing SUN   Created
%   2016-05-20   Schartner M. different file reading and indexing instead of loops


function oflux(source, sched, filename, fid_skd, PARA)

% open flux catalog file 
ifid = 1;
fid = fopen(filename, 'r');
if (fid < 0)
    ifid = 0;
end

% change textreading to improve speed 
C = textscan(fid,'%s %[^\n]','Delimiter',' ','Commentstyle','*'); 
fluxSrcName = C{1};
i = strcmp(fluxSrcName,'');
fluxSrcName(i)=[];
fluxOthers = C{2};
fluxOthers(i)=[];
% close file
if (fid >= 0)
    fclose(fid);
end

% write $FLUX
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
                for ib = 1 : PARA.MAX_BANDNUM
                    fprintf(fid_skd, '%s %s ',source(srcid).name(1:8),PARA.BAND(ib));
                    if (round(source(srcid).fluxpara(ib,PARA.MAX_FLUXPARA)) == 1)
                         fprintf(fid_skd, 'B ');
                    elseif (round(source(srcid).fluxpara(ib,PARA.MAX_FLUXPARA)) == 2)
                        fprintf(fid_skd, 'M ');
                    end
                    fprintf(fid_skd,'%3.1f ',  source(srcid).fluxpara(ib,1));
                    fprintf(fid_skd,'%5.2f ',  source(srcid).fluxpara(ib,2));
                    fprintf(fid_skd,'%7.1f\n', source(srcid).fluxpara(ib,3));
                end
            elseif (ifid == 1)
                
                % get Source Name and Common Source Name
                srcName = source(srcid).name(1:8);
                srcCommName = source(srcid).commoname(1:8);
                % Compare Names with Flux Catalog
                idx1 = strcmp(fluxSrcName,strtrim(srcName));
                idx2 = strcmp(fluxSrcName,strtrim(srcCommName));
                idx = idx1 | idx2;
                % check if there is a match
                if any(idx)
                    wholeLine = [fluxSrcName,fluxOthers]';
                    fprintf(fid_skd,'%s %s\n',wholeLine{:,idx});
                else 
                    for ib = 1 : PARA.MAX_BANDNUM
                        fprintf(fid_skd, '%s %s ',source(srcid).name(1:8),PARA.BAND(ib));
                        if (round(source(srcid).fluxpara(ib,PARA.MAX_FLUXPARA)) == 1)
                            fprintf(fid_skd, 'B ');
                        elseif (round(source(srcid).fluxpara(ib,PARA.MAX_FLUXPARA)) == 2)
                            fprintf(fid_skd, 'M ');
                        end
                        fprintf(fid_skd, '%3.1f ',  source(srcid).fluxpara(ib,1));
                        fprintf(fid_skd, '%5.2f ',  source(srcid).fluxpara(ib,2));
                        fprintf(fid_skd, '%7.1f\n', source(srcid).fluxpara(ib,3));
                    end
                end  
            end   %ifid
            
        end
    end
end


end


