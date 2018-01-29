% Purpose  
%   Calculate the schedule.
%   original functions vie_sched_s1 and vie_sched_s2 by Jing SUN
% History  
%   2016-07-11   Matthias Schartner:  Created
%   2016-09-06   Matthias Schartner:  uses bestSubcon.m  
%   2016-09-07   Matthias Schartner:  changes for singleCheckSource
%   2016-09-08   Matthias Schartner:  bugfix for source based scheduling
%   2016-09-08   Matthias Schartner:  changes for new fillin mode 1

function [sched, obs] = vie_sched_src(station, twin, source, obsmode, srcat, catpair, PARA, varargin)
%% #### initialize ####
if isempty(varargin)                                    % normal
    % initialize 'staobs'
    stanum = length(station);
    for ista = 1 : stanum
        staobs(ista).ifirst   = 1;
        staobs(ista).endmjd   = PARA.startmjd;
        staobs(ista).az       = station(ista).lim11;
        staobs(ista).el       = station(ista).lim12;
        staobs(ista).ha       = station(ista).lim11;
        staobs(ista).dc       = station(ista).lim12;
        staobs(ista).ncat     = 0;
        staobs(ista).cattn(1) = PARA.startmjd;
        staobs(ista).catsn(1) = 0;
    end
    % initialize 'srcobs'
    srcnum = length(source);
    for isrc = 1 : srcnum
        srcobs(isrc).obsmjd = 0.0;
        srcobs(isrc).nscan = 0.0;
    end
    % grid information
    srcn(1:length(srcat)) = 0;
    icatpair = 0;
    % load jpl earth
    jpl = load_jpl_Earth('jpl_421');
else                                                    % manually
    obs = varargin{1};
    staobs = obs.staobs;
    stanum = length(station);
    srcobs = obs.srcobs;
    srcnum = length(source);
    srcn = obs.srcn;
    icatpair = obs.icatpair;
    jpl = varargin{2};
end

% twin
twin.num = 0;
% fringe check
fringesrcid = 0;
ifringcheck = 0;
for isrc = 1 : srcnum
    if strcmp(deblank(source(isrc).name), deblank(PARA.SRCFRINGE))
        fringesrcid = isrc;
        break;
    end
end

%% #### schedule the session ####
nsched = 0;
scanuminfo = 0;
nmainsched = 0;
while true   %%% nsched

    % for STARMODE
    if PARA.STARMODE==1 && mod(nmainsched,PARA.CADENCE)==0 && nmainsched~=0
        starsource = 1;
        fprintf(PARA.fid_body,'     *** Observing STAR source ***\n');
    else
        starsource = 0;
    end
    % standard scheduling mode
    
    %% ##### Calc all possible sub-cons for the end of the last scan #####
    switch PARA.OPTIMIZATION
        case 1 
            %% #### source-based strategy ####
            [subcon,icatpair, srcn] = bestSubcon( 'source',station,twin,source,srcat,catpair,staobs,srcobs,icatpair,srcn,obsmode,jpl,PARA );
            
        case 2  
            %% #### station-based strategy ####
            subconall = bestSubcon( 'station',station,twin,source,staobs,srcobs,obsmode,PARA,starsource );
            
            % check if no scan is scheduled, e.g. all stations are down
            if subconall(1).nscan == 0 && length(subconall)==1
                subcon.nscan = 0;
            else
                % #### 2nd sort = final sort! ####
                [yy] = ssubsort(source, station, staobs, srcobs, subconall, PARA);

                % ### check the best sub-con ###
                for i = 0:length(subconall)-1
                    thisSubcon =  subconall(yy(end-i));
                    [flag,~,thisSubcon] = singleCheckSource(source, station, thisSubcon, PARA, jpl);
                    if flag == 1
                        % ### Pick the best sub-con ###
                        subcon = thisSubcon;
                        break
                    end
                end
                if i == length(subconall)-1 && flag ~= 1
                    subcon = subconall(yy(end));
                end
                
            end
            clear subconall;
    end

    %% #### fill-in 1 scheduling mode ####   
    if ((PARA.FILLINMODE == 1) || (PARA.FILLINMODE == 12)) && subcon.nscan>0
        while true
            [subcon_f1] = sscanfi1(station, twin, staobs, source, srcobs, obsmode, subcon, PARA, jpl);
            if (subcon_f1.nscan > 0)
                fprintf(PARA.fid_body,'------------ fill in 1 scan -----------  (sched: %3d)\n',nsched+1);
                [staobs, srcobs, endmjdmax, scanuminfo] = supdate(station, source, subcon_f1, staobs, srcobs, scanuminfo, PARA);
                nsched = nsched + 1;
                sched(nsched) = subcon_f1;
            elseif (subcon_f1.nscan == 0)
                break;
            end
        end
    end

    %% #### update structures ####
    fprintf(PARA.fid_body,'------------ %d-scan subnet ------------  (sched: %3d)\n',subcon.nscan,nsched+1);
    
    [staobs, srcobs, endmjdmax, scanuminfo] = supdate(station, source, subcon, staobs, srcobs, scanuminfo, PARA);
    if (subcon.nscan > 0)
        nsched = nsched + 1;
        nmainsched = nmainsched+1;
        sched(nsched) = subcon;
    end
    if(endmjdmax > PARA.endmjd)
        obs = struct();
        obs.staobs = staobs;
        obs.srcobs = srcobs;
        obs.srcn = srcn;
        obs.icatpair = icatpair;
        obs.srcat = srcat;
        break;
    end
    
    %% #### fill-in 2 scheduling mode ####   
    endmjdmaxfi2 = endmjdmax;
    if ((PARA.FILLINMODE == 2) || (PARA.FILLINMODE == 12)) && subcon.nscan>0
        while true
            [subcon_f2] = sscanfi2(station, twin, staobs, source, srcobs, obsmode, endmjdmaxfi2, PARA, jpl);
            if (subcon_f2.nscan > 0)
                fprintf(PARA.fid_body,'------------ fill in 2 scan -----------  (sched: %3d)\n',nsched+1);
                [staobs, srcobs, endmjdmax, scanuminfo] = supdate(station, source, subcon_f2, staobs, srcobs, scanuminfo, PARA);
                nsched = nsched + 1;
                sched(nsched) = subcon_f2;
            elseif (subcon_f2.nscan == 0)
                break;
            end
        end
    end
    if(endmjdmax > PARA.endmjd)
        obs = struct();
        obs.staobs = staobs;
        obs.srcobs = srcobs;
        obs.srcn = srcn;
        obs.icatpair = icatpair;
        obs.srcat = srcat;
        break;
    end
    
    %% #### particular sources ####
    for isrc = 1 : srcnum
        if (source(isrc).ifp==1)&&(endmjdmax>source(isrc).pst)&&((endmjdmax-srcobs(isrc).obsmjd)*1440>source(isrc).pdt)
            [subcon] = spsource(isrc, station, twin, staobs, source, srcobs, obsmode, PARA, jpl);
            if (subcon.nscan > 0)
                fprintf(PARA.fid_body,'=======================================\n');
                if ((PARA.FILLINMODE == 1) || (PARA.FILLINMODE == 12))
                    addSec = zeros(stanum,1);
                    while true
                        [subcon_f1,addSec] = sscanfi1(station, twin, staobs, source, srcobs, obsmode, subcon, PARA, addSec, jpl);
                        if (subcon_f1.nscan > 0)
                            fprintf(PARA.fid_body,'------------ fill in 1 scan -----------  (sched: %3d)\n',nsched+1);
                            [staobs, srcobs, endmjdmax, scanuminfo] = supdate(station, source, subcon_f1, staobs, srcobs, scanuminfo, PARA);
                            nsched = nsched + 1;
                            sched(nsched) = subcon_f1;
                        elseif (subcon_f1.nscan == 0)
                            break;
                        end
                    end
                end   
                fprintf(PARA.fid_body,'       ### PARTICULAR SOURCE ###\n');
                fprintf(PARA.fid_body,'------------ 1-scan subnet ------------  (sched: %3d)\n',nsched+1);
                [staobs, srcobs, endmjdmax, scanuminfo] = supdate(station, source, subcon, staobs, srcobs, scanuminfo, PARA);
                nsched = nsched + 1;
                sched(nsched) = subcon;
                % fill-in mode
                endmjdmaxfi2 = endmjdmax;
                if ((PARA.FILLINMODE == 2) || (PARA.FILLINMODE == 12))
                    while true
                        [subcon_f2] = sscanfi2(station, twin, staobs, source, srcobs, obsmode, endmjdmaxfi2, PARA, jpl);
                        if (subcon_f2.nscan > 0)
                            fprintf(PARA.fid_body,'------------ fill in 2 scan -----------  (sched: %3d)\n',nsched+1);
                            [staobs, srcobs, endmjdmax, scanuminfo] = supdate(station, source, subcon_f2, staobs, srcobs, scanuminfo, PARA);
                            nsched = nsched + 1;
                            sched(nsched) = subcon_f2;
                        elseif (subcon_f2.nscan == 0)
                            break;
                        end
                    end
                end
            end   %if (subcon.nscan > 0) 
        end   %
    end   %for isrc = 1 : srcnum
    if(endmjdmax > PARA.endmjd)
        obs = struct();
        obs.staobs = staobs;
        obs.srcobs = srcobs;
        obs.srcn = srcn;
        obs.icatpair = icatpair;
        obs.srcat = srcat;
        break;
    end
    
    %% #### sources for fringe check ####
    if (fringesrcid ~= 0) && (ifringcheck == 0)
        downstanum = 0;
        downstasn(1:stanum) = 0; 
        for ista = 1 : stanum
            endmjd = staobs(ista).endmjd;
            [az, el, ha, dc] = zazel_s(endmjd, station(ista).llh(1), station(ista).llh(2), ra, de);
            [lup] = zlup(endmjd, az, el, ha, dc, ista, station, PARA.MIN_CUTEL);
            if (lup == 1)
                downstanum = downstanum + 1;
                downstasn(downstanum) = ista;
            end
        end   
        if (downstanum < stanum)
            continue;
        end
        subcontmp.nscan = 1;
        subcontmp.scan(1).srcid    = fringesrcid;
        subcontmp.scan(1).startmjd = 0.0;
        subcontmp.scan(1).nsta     = downstanum;
        for ista = 1 : subcontmp.scan(1).nsta
            subcontmp.scan(1).sta(ista).staid = downstasn(ista);
            subcontmp.scan(1).sta(ista).slewtime = 0.0;
            subcontmp.scan(1).sta(ista).startmjd = 0.0;
            subcontmp.scan(1).sta(ista).az       = 0.0;
            subcontmp.scan(1).sta(ista).el       = 0.0;
            subcontmp.scan(1).sta(ista).ha       = 0.0;
            subcontmp.scan(1).sta(ista).dc       = 0.0;
            subcontmp.scan(1).sta(ista).duration = 0.0;
            subcontmp.scan(1).sta(ista).endmjd   = 0.0;
        end
        [subcontmp1] = sstartmjd(source, station, twin, staobs, subcontmp, PARA);
        if (subcontmp1.scan(1).nsta < stanum)
            continue;
        end
        for ista = 1 : subcontmp1.scan(1).nsta
            subcontmp1.scan(1).sta(ista).duration = PARA.MIN_SCAN;
            subcontmp1.scan(1).sta(ista).endmjd   = subcontmp1.scan(1).startmjd + subcontmp1.scan(1).sta(ista).duration / 86400.0;
        end
        [staobs, srcobs, endmjdmax, scanuminfo] = supdate(station, source, subcontmp1, staobs, srcobs, scanuminfo, PARA);
        nsched = nsched + 1;
        sched(nsched) = subcontmp1;
        ifringcheck = 1;
    end
    %% #### update even if no scan is generated #### 
    if(subcon.nscan == 0)
        for ista = 1 : stanum
            staobs(ista).endmjd = staobs(ista).endmjd+1.0/1440;
        end
    end 
    fprintf(PARA.fid_body,'=======================================\n');
    clear subcon
end

