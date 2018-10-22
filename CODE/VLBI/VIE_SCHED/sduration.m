% Purpose  
%   Calculate duration (scan length).
% History  
%   2010-04-06   Jing SUN   Created
%   2016-05-10	Matthias Schartner: improve speed slightly.
%   2016-05-25  Matthias Schartner: STAR mode added
%   2016-06-14  Matthias Schartner: changed the way deleted stations and
%                                   scans are handled
%   2016-07-11  Matthias Schartner: bugfix
%   2016-09-06 M. Schartner: bugfix, actbl as 2nd return value

function [subcon_s,baselineDuration] = sduration(source, station, twin, obsmode, subcon, PARA)

subcon_s.nscan = 0;
baselineDuration = struct();
% calculate duration
for iscan = 1 : subcon.nscan
    % scan
    scantmp = subcon.scan(iscan);
    % source flux
    srcid = scantmp.srcid;
    ra = source(srcid).ra;
    de = source(srcid).de;
    fluxpara(1:PARA.MAX_BANDNUM,1:PARA.MAX_FLUXPARA) = source(srcid).fluxpara(1:PARA.MAX_BANDNUM,1:PARA.MAX_FLUXPARA);
    %
    badstanum   = 0;
    badstasn(1) = 0;
    while true
        ifsta(1:scantmp.nsta) = 0;
        % baseline
        blnum = 0;
        for ibl = 1 : (scantmp.nsta - 1)
            for jbl = (ibl + 1) : scantmp.nsta
                %
                if ((size(find(badstasn(1:badstanum)==ibl),2) > 0) || (size(find(badstasn(1:badstanum)==jbl),2) > 0))
                    continue;
                end
                % station sefdpara
                staid1 = scantmp.sta(ibl).staid;
                staid2 = scantmp.sta(jbl).staid;
                el1 = scantmp.sta(ibl).el;
                el2 = scantmp.sta(jbl).el;
                sefdpara1(1:PARA.MAX_BANDNUM,1:PARA.MAX_SEFDPARA) = station(staid1).sefdpara(1:PARA.MAX_BANDNUM,1:PARA.MAX_SEFDPARA); 
                sefdpara2(1:PARA.MAX_BANDNUM,1:PARA.MAX_SEFDPARA) = station(staid2).sefdpara(1:PARA.MAX_BANDNUM,1:PARA.MAX_SEFDPARA); 
                % check the twin/multiple antennas at a site (1 - same source observations)
                if ((twin.num > 0) && (PARA.TWINMODE == 1))
                    for itwin = 1 : twin.num
                         if (staid1 == twin.sn(itwin,1)) 
                            for iband = 1 : PARA.MAX_BANDNUM
                                sefdpara1(iband,1) = sefdpara1(iband,1)/2;
                            end
                         end
                        if (staid2 == twin.sn(itwin,1))
                            for iband = 1 : PARA.MAX_BANDNUM
                                sefdpara2(iband,1) = sefdpara2(iband,1)/2;
                            end
                        end
                    end
                end
                % baseline
                blnum = blnum + 1;
                actbl(blnum).staid1 = staid1;
                actbl(blnum).staid2 = staid2;
                % calculate obsflux
                blx = station(staid1).xyz(1)-station(staid2).xyz(1);
                bly = station(staid1).xyz(2)-station(staid2).xyz(2);
                blz = station(staid1).xyz(3)-station(staid2).xyz(3);
                [actbl(blnum).obsflux] = sobsflux(scantmp.startmjd, blx, bly, blz, ra, de, fluxpara, PARA);
                % calculate duration
                for iband = 1 : PARA.MAX_BANDNUM
                    % sefd
                    [sefdel1] = ssefdel(el1, sefdpara1(iband,1), sefdpara1(iband,2), sefdpara1(iband,3), sefdpara1(iband,4));
                    [sefdel2] = ssefdel(el2, sefdpara2(iband,1), sefdpara2(iband,2), sefdpara2(iband,3), sefdpara2(iband,4));
                    % duration 
                    minsnr = min(station(staid1).minsnr(iband),station(staid2).minsnr(iband));
                    anum = (1.75 * minsnr / actbl(blnum).obsflux(iband)) ^ 2;
                    anu1 = sefdel1 * sefdel2;
                    anu2 = obsmode.samprate * 1.0d6 * obsmode.chanumband(iband) * obsmode.bits;
                    actbl(blnum).duration(iband) = ceil(anum * (anu1 / anu2) + PARA.CORSYNCH);  %%%CORSYNCH
                    if (actbl(blnum).duration(iband) > PARA.MAX_SCAN && source(srcid).ifp==0 && source(srcid).star==0)
                        ifsta(ibl) = ifsta(ibl) + actbl(blnum).duration(iband);
                        ifsta(jbl) = ifsta(jbl) + actbl(blnum).duration(iband);                 
                    elseif (actbl(blnum).duration(iband) > 2000) && (source(srcid).ifp==1) %%%particular sources
                        ifsta(ibl) = ifsta(ibl) + 1;
                        ifsta(jbl) = ifsta(jbl) + 1;
                    elseif (actbl(blnum).duration(iband) > 10000) && (source(srcid).star==1) %%% STAR mode
                        ifsta(ibl) = ifsta(ibl) + 1;
                        ifsta(jbl) = ifsta(jbl) + 1;
                    end
                end
            end 
        end
        
        % check flux and duration at each station 
        [bads, badsta] = max(ifsta(1:scantmp.nsta));
        if (bads > 0)
            badstanum = badstanum + 1;
            badstasn(badstanum) = badsta;
            clear actbl
        else 
            break;
        end
    end
    
    if PARA.STARMODE==1
        strongantidx = find(strcmp(PARA.STRONGANT,{station.name}));  %%%%%%%%%%%%%%%% STARMODE
        strongant = find([scantmp.sta.staid]==strongantidx,1);
    end
        
    % calculate duration for each station
    if source(srcid).star == 0;
        for ista = 1 : scantmp.nsta
            if (size(find(badstasn(1:badstanum)==ista),2) > 0)
                continue;
            end
            staid = scantmp.sta(ista).staid;
            stadur = 0;
            for ibl = 1 : blnum
                staid1 = actbl(ibl).staid1;
                staid2 = actbl(ibl).staid2;
                if ((staid1 == staid) || (staid2 == staid))
                    for iband = 1 : PARA.MAX_BANDNUM
                        if (actbl(ibl).duration(iband) > stadur)
                            stadur = actbl(ibl).duration(iband);
                        end
                    end
                end
            end
            if (stadur < PARA.MIN_SCAN)
                scantmp.sta(ista).duration = PARA.MIN_SCAN;             
            else
                scantmp.sta(ista).duration = stadur;   
            end      
            scantmp.sta(ista).endmjd = scantmp.startmjd + scantmp.sta(ista).duration / 86400.0; %%%
        end
    else         
        for ista = 1 : scantmp.nsta
            if (size(find(badstasn(1:badstanum)==ista),2) > 0)
                continue;
            end
            staid = scantmp.sta(ista).staid;
            stadur = 0;
            for ibl = 1 : blnum
                staid1 = actbl(ibl).staid1;
                staid2 = actbl(ibl).staid2;
                if ~isempty(strongantidx)
                    if (staid1 == strongantidx || staid2 == strongantidx)
                        if ((staid1 == staid) || (staid2 == staid))
                            for iband = 1 : PARA.MAX_BANDNUM
                                if (actbl(ibl).duration(iband) > stadur)
                                    stadur = actbl(ibl).duration(iband);
                                end
                            end
                        end
                    end
                else
                    if ((staid1 == staid) || (staid2 == staid))
                        for iband = 1 : PARA.MAX_BANDNUM
                            if (actbl(ibl).duration(iband) > stadur)
                                stadur = actbl(ibl).duration(iband);
                            end
                        end
                    end
                end
            end
            if (stadur < PARA.MIN_SCAN)
                scantmp.sta(ista).duration = PARA.MIN_SCAN;             
            else
                scantmp.sta(ista).duration = stadur;   
            end      
            scantmp.sta(ista).endmjd = scantmp.startmjd + scantmp.sta(ista).duration / 86400.0; %%%
        end
    end
    
    % exclude the strong antenna for normal scans in STAR mode
    if PARA.STARMODE==1
        if source(srcid).star==1
            if isempty(strongant)
                lst = length([scantmp.sta.staid]);
                for ir = 1:lst
                    badstanum = badstanum+1;
                    badstasn(badstanum)=scantmp.sta(ir).staid;
                end
                continue
            end
        else
            if ~isempty(strongant)
                badstanum = badstanum+1;
                badstasn(badstanum) = strongant;
            end
        end
    end
    
    % write output variable
    subcon_s.nscan = subcon_s.nscan + 1;
    subcon_s.scan(subcon_s.nscan).srcid = scantmp.srcid;
    subcon_s.scan(subcon_s.nscan).startmjd = scantmp.startmjd;
    subcon_s.scan(subcon_s.nscan).nsta = 0;
    for ista = 1 : scantmp.nsta
        if (size(find(badstasn(1:badstanum)==ista),2) > 0)
            continue;
        end
        subcon_s.scan(subcon_s.nscan).nsta = subcon_s.scan(subcon_s.nscan).nsta + 1;
        subcon_s.scan(subcon_s.nscan).sta(subcon_s.scan(subcon_s.nscan).nsta).staid    = scantmp.sta(ista).staid;
        subcon_s.scan(subcon_s.nscan).sta(subcon_s.scan(subcon_s.nscan).nsta).slewtime = scantmp.sta(ista).slewtime;
        subcon_s.scan(subcon_s.nscan).sta(subcon_s.scan(subcon_s.nscan).nsta).startmjd = scantmp.sta(ista).startmjd; 
        subcon_s.scan(subcon_s.nscan).sta(subcon_s.scan(subcon_s.nscan).nsta).az       = scantmp.sta(ista).az;
        subcon_s.scan(subcon_s.nscan).sta(subcon_s.scan(subcon_s.nscan).nsta).el       = scantmp.sta(ista).el;
        subcon_s.scan(subcon_s.nscan).sta(subcon_s.scan(subcon_s.nscan).nsta).ha       = scantmp.sta(ista).ha;
        subcon_s.scan(subcon_s.nscan).sta(subcon_s.scan(subcon_s.nscan).nsta).dc       = scantmp.sta(ista).dc;
        subcon_s.scan(subcon_s.nscan).sta(subcon_s.scan(subcon_s.nscan).nsta).duration = scantmp.sta(ista).duration;
        subcon_s.scan(subcon_s.nscan).sta(subcon_s.scan(subcon_s.nscan).nsta).endmjd   = scantmp.sta(ista).endmjd;
    end
    if (subcon_s.scan(subcon_s.nscan).nsta < 2)
        subcon_s.scan(subcon_s.nscan) =[];
        subcon_s.nscan = subcon_s.nscan - 1;
    end
    if subcon_s.nscan >0 && exist('actbl','var')
        baselineDuration(iscan).actbl=actbl;
    else
        baselineDuration(iscan).actbl = [];
    end
    clear actbl
end


