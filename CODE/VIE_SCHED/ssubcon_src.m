% Purpose  
%   Calculate subconfiguration (a group of scans).
% History  
%   2010-04-06   Jing SUN   Created
%   


function [subcon2, icatpair, srcn] = ssubcon_src(station, twin, source, srcat, catpair, staobs, srcobs, icatpair, srcn, PARA)

% initialization
stanum = length(station);
srcnum = length(source);
gridnum = length(srcat);
pairnum = length(catpair);
% endmjdmax
endmjdmax = 0.0;
for ista = 1 : stanum
    if (staobs(ista).endmjd > endmjdmax)
        endmjdmax = staobs(ista).endmjd;
    end
end

% select the source  
srcseq(1:PARA.SRCNUM) = 0;
loopnum = 0;
while true  % catpair
    loopnum = loopnum + 1;
    if(loopnum > pairnum)
        break;
    end
    icatpair = icatpair + 1;   %%%Jing
    if (icatpair > pairnum)
        icatpair = icatpair - pairnum;
    end
    for i = 1 : PARA.SRCNUM  
        pair = catpair(icatpair).pair(i) - 1;
        for i1 = 1 : gridnum % grid
            % select grid
            pair = pair + 1;
            if (pair > gridnum)
                pair = pair - gridnum;
            end
            if (srcat(pair).num == 0)
                continue;
            end
            ifar = 1;
            ifsrcrep = 1;
            nsta = 0;
            iflux = 1;
            % select source in this grid
            for i2 = 1 : srcat(pair).num  % source in the grid
                srcn(pair) = srcn(pair) + 1;
                if (srcn(pair) > srcat(pair).num)
                    srcn(pair) = srcn(pair) - srcat(pair).num;
                end
                srcid = srcat(pair).src(srcn(pair));
                ifar = 1;
                ifsrcrep = 1;
                nsta = 0;
                iflux = 1;
                % check the angular distance
                if ((PARA.SRCNUM==2) & (i == 2) & srcseq(1)>0)
                    srcid1 = srcseq(1);
                    srcid2 = srcid;
                    ra1 = source(srcid1).ra;
                    de1 = source(srcid1).de;
                    ra2 = source(srcid2).ra;
                    de2 = source(srcid2).de;
                    arg = cos(de1) * cos(de2) * cos(ra1 - ra2) + sin(de1) * sin(de2);
                    srcd = acos(arg)*180/pi;
                    if (srcd < 120.0)
                        ifar = 0;
                    end
                end
                if ((PARA.SRCNUM==4) & (i >= 2))
                    srcid1 = srcseq(1);
                    ra1 = source(srcid1).ra;
                    de1 = source(srcid1).de;
                    if(i==2)
                        srcid2 = srcid;
                        ra2 = source(srcid2).ra;
                        de2 = source(srcid2).de;
                        arg = cos(de1) * cos(de2) * cos(ra1 - ra2) + sin(de1) * sin(de2);
                        srcd = acos(arg)*180/pi;
                        if (srcd < 60.0)
                            ifar = 0;
                        end
                    elseif (i==3)
                        srcid3 = srcid;
                        ra3 = source(srcid3).ra;
                        de3 = source(srcid3).de;
                        arg = cos(de1) * cos(de3) * cos(ra1 - ra3) + sin(de1) * sin(de3);
                        srcd = acos(arg)*180/pi;
                        if (srcd < 60.0)
                            ifar = 0;
                        end
                        arg = cos(de2) * cos(de3) * cos(ra2 - ra3) + sin(de2) * sin(de3);
                        srcd = acos(arg)*180/pi;
                        if (srcd < 60.0)
                            ifar = 0;
                        end
                    elseif (i==4)
                        srcid4 = srcid;
                        ra4 = source(srcid4).ra;
                        de4 = source(srcid4).de;
                        arg = cos(de1) * cos(de4) * cos(ra1 - ra4) + sin(de1) * sin(de4);
                        srcd = acos(arg)*180/pi;
                        if (srcd < 60.0)
                            ifar = 0;
                        end
                        arg = cos(de2) * cos(de4) * cos(ra2 - ra4) + sin(de2) * sin(de4);
                        srcd = acos(arg)*180/pi;
                        if (srcd < 60.0)
                            ifar = 0;
                        end
                        arg = cos(de3) * cos(de4) * cos(ra3 - ra4) + sin(de3) * sin(de4);
                        srcd = acos(arg)*180/pi;
                        if (srcd < 60.0)
                            ifar = 0;
                        end
                    end
                end
                if (ifar == 0)
                    continue;
                end
                % check if the source is repeated in a short interval
                if ((endmjdmax-srcobs(srcid).obsmjd) < PARA.MIN_SRCRP/1440)
                    ifsrcrep = 0;
                    continue;
                end
                % check the visibility at stations
                for ista = 1 : stanum
                    [az, el, ha, dc] = zazel_s(endmjdmax, station(ista).llh(1), station(ista).llh(2), source(srcid).ra, source(srcid).de);
                    [lup] = zlup(endmjdmax, az, el, ha, dc, ista, station, PARA.MIN_CUTEL);
                    if (lup == 1)
                        nsta = nsta + 1;
                        stasn(nsta) = ista;
                    end
                end
                if (nsta < 2)
                    continue;
                end
                % source projected flux
                if (source(srcid).ifp == 0) 
                    ra = source(srcid).ra;
                    de = source(srcid).de;
                    fluxpara(1:PARA.MAX_BANDNUM,1:PARA.MAX_FLUXPARA) = source(srcid).fluxpara(1:PARA.MAX_BANDNUM,1:PARA.MAX_FLUXPARA);
                    for ibl = 1 : (nsta - 1)
                        for jbl = (ibl + 1) : nsta   
                            staid1 = stasn(ibl);
                            staid2 = stasn(jbl);  
                            blx = station(staid1).xyz(1)-station(staid2).xyz(1);
                            bly = station(staid1).xyz(2)-station(staid2).xyz(2);
                            blz = station(staid1).xyz(3)-station(staid2).xyz(3);
                            [obsflux] = sobsflux(endmjdmax, blx, bly, blz, ra, de, fluxpara, PARA);
                            if (min(obsflux(1:PARA.MAX_BANDNUM)) < PARA.MIN_FLUX)  
                                iflux = 0;
                            end
                        end
                    end
                end
                if (iflux == 1)
                    break;
                end 
            end %for i2 = 1 : srcat(pair1).num
            if ((nsta >= 2) & (iflux == 1) & (ifsrcrep == 1) & (ifar == 1)) 
                break;
            end
        end %for i1 = 1 : gridnum
        if ((nsta >= 2) & (iflux == 1) & (ifsrcrep == 1) & (ifar == 1)) 
            srcseq(i) = srcid;
            if (i == 1)
                icatpair = pair;
            end
        end
    end %for i = 1 : PARA.SRCNUM     
    if(prod(srcseq(1:PARA.SRCNUM)) == 0)
        subcon2.nscan = 0;
        continue;
    end

    % calculate the subconfiguration  
    for ista = 1 : stanum
        for i = 1 : PARA.SRCNUM
            srcid = srcseq(i);
            [azs(ista,i), els(ista,i), ha, dc] = zazel_s(endmjdmax, station(ista).llh(1), station(ista).llh(2), source(srcid).ra, source(srcid).de);     
            [lup] = zlup(endmjdmax, azs(ista,i), els(ista,i), ha, dc, ista, station, PARA.MIN_CUTEL);
            [cat] = zcoverage(azs(ista,i), els(ista,i));
            ifnewcat(i) = 1;
            for icat = 1 : staobs(ista).ncat
                if (staobs(ista).catsn(icat)==cat)
                    ifnewcat(i) = 0;
                    break;
                end
            end
            if (lup == 0)
                els(ista,i) = 0.0;
                ifnewcat(i) = 0;
            end
        end
        %
        [x, y] = sort(els(ista,1:PARA.SRCNUM));
        stavis(ista) = y(PARA.SRCNUM);
        for i = PARA.SRCNUM : -1 : 1
            if (i==PARA.SRCNUM) & (x(i)<PARA.MIN_CUTEL)
                stavis(ista) = 0;
                break;
            end
            if (x(i)<PARA.MIN_CUTEL)
                break;
            end
            if (x(i)>PARA.MIN_CUTEL)&(ifnewcat(y(i))==1)
                stavis(ista) = y(i);   
                break;
            end
        end
    end
    
    % check the twin/multiple antennas at a site (1 - same source observations)
    if ((twin.num > 0) & (PARA.TWINMODE == 1))
        for itwin = 1 : twin.num
            staid1 = twin.sn(itwin,1);
            staid2 = twin.sn(itwin,2);
            stavis(staid2) = 0;
        end
    end
    % check the twin/multiple antennas at a site (2 - continuous observations)
    if ((twin.num > 0) & (PARA.TWINMODE == 2))
        for itwin = 1 : twin.num
            staid1 = twin.sn(itwin,1);
            staid2 = twin.sn(itwin,2);
            if (staobs(staid1).endmjd < staobs(staid2).endmjd)
                stavis(staid1) = stavis(staid1);
                stavis(staid2) = 0;
            else
                stavis(staid2) = stavis(staid1);
                stavis(staid1) = 0;  
            end
        end   
    end
    % check the twin/multiple antennas at a site (3 - multidirectional observations)
    if ((twin.num > 0) & (PARA.TWINMODE == 3))
        for itwin = 1 : twin.num
            staid1 = twin.sn(itwin,1);
            staid2 = twin.sn(itwin,2);
            sortel1=els(staid1,1:PARA.SRCNUM);
            [x1, y1] = sort(sortel1(1:PARA.SRCNUM));
            el1 = x1(PARA.SRCNUM);
            az1 = azs(staid1, y1(PARA.SRCNUM));
            el2 = x1(PARA.SRCNUM-1);
            az2 = azs(staid1, y1(PARA.SRCNUM-1));    
            ha=0; dc=0;
            if (el1 >= PARA.MIN_CUTEL) & (el2 >= PARA.MIN_CUTEL)
                % p1
                [st11, unaz] = sslew(station, staobs, az1, el1, ha, dc, staid1, PARA); 
                [st22, unaz] = sslew(station, staobs, az2, el2, ha, dc, staid2, PARA);
                st1 = max(st11, st22);
                % p2
                [st21, unaz] = sslew(station, staobs, az1, el1, ha, dc, staid2, PARA); 
                [st12, unaz] = sslew(station, staobs, az2, el2, ha, dc, staid1, PARA); 
                st2 = max(st21, st12);      
                %
                if (st1 < st2)
                    stavis(staid1) = y1(PARA.SRCNUM);
                    stavis(staid2) = y1(PARA.SRCNUM-1);
                else
                   stavis(staid1) = y1(PARA.SRCNUM-1);
                   stavis(staid2) = y1(PARA.SRCNUM); 
                end
            elseif (el1 >= PARA.MIN_CUTEL) & (el2 < PARA.MIN_CUTEL)
                if (staobs(staid1).endmjd < staobs(staid2).endmjd)
                    stavis(staid1) = y1(PARA.SRCNUM);
                    stavis(staid2) = 0;
                else
                    stavis(staid1) = 0;
                    stavis(staid2) = y1(PARA.SRCNUM);
                end
            end
        end
    end
        
    % subcon1
    subcon1.nscan = 0;
    for i = 1 : PARA.SRCNUM
        subcon1.nscan = subcon1.nscan + 1;
        subcon1.scan(subcon1.nscan).srcid = srcseq(i);
        subcon1.scan(subcon1.nscan).startmjd = 0.0;
        subcon1.scan(subcon1.nscan).nsta = 0;
        for ista = 1 : stanum
            if (stavis(ista) == i)
                subcon1.scan(subcon1.nscan).nsta = subcon1.scan(subcon1.nscan).nsta + 1;
                subcon1.scan(subcon1.nscan).sta(subcon1.scan(subcon1.nscan).nsta).staid    = ista;
                subcon1.scan(subcon1.nscan).sta(subcon1.scan(subcon1.nscan).nsta).slewtime = 0.0;
                subcon1.scan(subcon1.nscan).sta(subcon1.scan(subcon1.nscan).nsta).startmjd = 0.0;
                subcon1.scan(subcon1.nscan).sta(subcon1.scan(subcon1.nscan).nsta).az       = 0.0;
                subcon1.scan(subcon1.nscan).sta(subcon1.scan(subcon1.nscan).nsta).el       = 0.0;
                subcon1.scan(subcon1.nscan).sta(subcon1.scan(subcon1.nscan).nsta).ha       = 0.0;
                subcon1.scan(subcon1.nscan).sta(subcon1.scan(subcon1.nscan).nsta).dc       = 0.0;
                subcon1.scan(subcon1.nscan).sta(subcon1.scan(subcon1.nscan).nsta).duration = 0.0;
                subcon1.scan(subcon1.nscan).sta(subcon1.scan(subcon1.nscan).nsta).endmjd   = 0.0;
            end
        end
    end
    % check the number of station in a scan (>=2)
    subcon2.nscan = 0;
    for iscan = 1 : subcon1.nscan
        if (subcon1.scan(iscan).nsta >= 2)
            subcon2.nscan = subcon2.nscan + 1;
            subcon2.scan(subcon2.nscan).srcid = subcon1.scan(iscan).srcid;
            subcon2.scan(subcon2.nscan).startmjd = 0.0;
            subcon2.scan(subcon2.nscan).nsta = subcon1.scan(iscan).nsta;
            for ista = 1 : subcon1.scan(iscan).nsta
                subcon2.scan(subcon2.nscan).sta(ista).staid    = subcon1.scan(iscan).sta(ista).staid; 
                subcon2.scan(subcon2.nscan).sta(ista).slewtime = 0.0; 
                subcon2.scan(subcon2.nscan).sta(ista).startmjd = 0.0; 
                subcon2.scan(subcon2.nscan).sta(ista).az       = 0.0;  
                subcon2.scan(subcon2.nscan).sta(ista).el       = 0.0; 
                subcon2.scan(subcon2.nscan).sta(ista).ha       = 0.0; 
                subcon2.scan(subcon2.nscan).sta(ista).dc       = 0.0; 
                subcon2.scan(subcon2.nscan).sta(ista).duration = 0.0; 
                subcon2.scan(subcon2.nscan).sta(ista).endmjd   = 0.0;
            end
        end
    end
    
    % get subconfiguration and break
    if(subcon2.nscan > 0)
        break;
    end
end


