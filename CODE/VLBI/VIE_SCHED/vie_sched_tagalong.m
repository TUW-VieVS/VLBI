% Purpose  
%   Calculate the schedule (tag-along mode for special station).
% History  
%   2013-12-30   Jing SUN   Created
%   


function [sched] = vie_sched_tagalong(station1, station2, source, obsmode, sched, PARA)

% initialize 'staobs'
stanum = length(station2);
for itag = 1 : stanum
    staobs(itag).ifirst   = 1;
    staobs(itag).endmjd   = PARA.startmjd;
    staobs(itag).az       = station2(itag).lim11;
    staobs(itag).el       = station2(itag).lim12;
    staobs(itag).ha       = station2(itag).lim11;
    staobs(itag).dc       = station2(itag).lim12;
end

% loop sched
for isched = 1 : length(sched)
    for iscan = 1 : sched(isched).nscan
        ra = source(sched(isched).scan(iscan).srcid).ra;
        de = source(sched(isched).scan(iscan).srcid).de;
        fluxpara(1:PARA.MAX_BANDNUM,1:PARA.MAX_FLUXPARA) = source(sched(isched).scan(iscan).srcid).fluxpara(1:PARA.MAX_BANDNUM,1:PARA.MAX_FLUXPARA);
        for itag = 1 : stanum
            % check the visibility at the start time of scan
            [az, el, ha, dc] = zazel_s(sched(isched).scan(iscan).startmjd, station2(itag).llh(1), station2(itag).llh(2), ra, de);
            [lup] = zlup(sched(isched).scan(iscan).startmjd, az, el, ha, dc, itag, station2, PARA.MIN_CUTEL);
            if (lup == 0)
                continue;
            end
            % check the slewing time
            if (staobs(itag).ifirst==1) 
                slewtime = 0.0;
                unaz = az;
                startmjd = PARA.startmjd;
            else
                [slewtime, unaz] = sslew(station2, staobs, az, el, ha, dc, itag, PARA);
                startmjd = staobs(itag).endmjd + slewtime/86400.0 + PARA.SOURCE/86400 + PARA.TAPETM/86400 + PARA.IDLE/86400 + PARA.CALIBRATION/86400;  
            end
            if (startmjd > sched(isched).scan(iscan).startmjd)
                continue;
            end
%             % check the visibility at the arrival time of slewing
%             [az, el, ha, dc] = zazel_s(startmjd, station2(itag).llh(1), station2(itag).llh(2), ra, de);
%             [lup] = zlup(startmjd, az, el, ha, dc, itag, station2, PARA.MIN_CUTEL);
%             if (lup == 0)
%                 continue;
%             end
            % check the scan length
            maxscanlength = 0.0;
            for ista = 1 : sched(isched).scan(iscan).nsta
                if (sched(isched).scan(iscan).sta(ista).duration >= maxscanlength)
                    maxscanlength = sched(isched).scan(iscan).sta(ista).duration;
                end
            end
            for ista = 1 : sched(isched).scan(iscan).nsta
                sefdpara1(1:PARA.MAX_BANDNUM,1:PARA.MAX_SEFDPARA) = station1(sched(isched).scan(iscan).sta(ista).staid).sefdpara(1:PARA.MAX_BANDNUM,1:PARA.MAX_SEFDPARA);
                sefdpara2(1:PARA.MAX_BANDNUM,1:PARA.MAX_SEFDPARA) = station2(itag).sefdpara(1:PARA.MAX_BANDNUM,1:PARA.MAX_SEFDPARA); 
                % calculate obsflux
                blx = station2(itag).xyz(1)-station1(sched(isched).scan(iscan).sta(ista).staid).xyz(1);
                bly = station2(itag).xyz(2)-station1(sched(isched).scan(iscan).sta(ista).staid).xyz(2);
                blz = station2(itag).xyz(3)-station1(sched(isched).scan(iscan).sta(ista).staid).xyz(3);
                [actbl(ista).obsflux] = sobsflux(sched(isched).scan(iscan).startmjd, blx, bly, blz, ra, de, fluxpara, PARA);
                % calculate duration
                for iband = 1 : PARA.MAX_BANDNUM
                    % sefd
                    [sefdel1] = ssefdel(sched(isched).scan(iscan).sta(ista).el, sefdpara1(iband,1), sefdpara1(iband,2), sefdpara1(iband,3), sefdpara1(iband,4));
                    [sefdel2] = ssefdel(el, sefdpara2(iband,1), sefdpara2(iband,2), sefdpara2(iband,3), sefdpara2(iband,4));
                    % duration 
                    minsnr = min(station1(sched(isched).scan(iscan).sta(ista).staid).minsnr(iband),station2(itag).minsnr(iband));
                    anum = (1.75 * minsnr / actbl(ista).obsflux(iband)) ^ 2;
                    anu1 = sefdel1 * sefdel2;
                    anu2 = obsmode.samprate * 1.0d6 * obsmode.chanumband(iband) * obsmode.bits;
                    actbl(ista).duration(iband) = ceil(anum * (anu1 / anu2) + PARA.CORSYNCH);  %%%CORSYNCH
                end
            end
            duration = 0;
            for ista = 1 : sched(isched).scan(iscan).nsta
                for iband = 1 : PARA.MAX_BANDNUM
                    if (actbl(ista).duration(iband) >= duration)
                        duration = actbl(ista).duration(iband);
                    end
                end
            end
            if (duration > maxscanlength)
                continue;
            end
            % add station to sched structure
            sched(isched).scan(iscan).nsta = sched(isched).scan(iscan).nsta + 1;
            sched(isched).scan(iscan).sta(sched(isched).scan(iscan).nsta).staid = itag + length(station1); %%%
            sched(isched).scan(iscan).sta(sched(isched).scan(iscan).nsta).slewtime = slewtime;
            sched(isched).scan(iscan).sta(sched(isched).scan(iscan).nsta).startmjd = startmjd;
            sched(isched).scan(iscan).sta(sched(isched).scan(iscan).nsta).az = unaz;
            sched(isched).scan(iscan).sta(sched(isched).scan(iscan).nsta).el = el;
            sched(isched).scan(iscan).sta(sched(isched).scan(iscan).nsta).ha = ha;
            sched(isched).scan(iscan).sta(sched(isched).scan(iscan).nsta).dc = dc;
            if (duration < PARA.MIN_SCAN)
                sched(isched).scan(iscan).sta(sched(isched).scan(iscan).nsta).duration = PARA.MIN_SCAN;             
            else
                sched(isched).scan(iscan).sta(sched(isched).scan(iscan).nsta).duration = duration;   
            end 
            sched(isched).scan(iscan).sta(sched(isched).scan(iscan).nsta).endmjd = sched(isched).scan(iscan).startmjd + duration / 86400.0;
           % update staobs variable
           staobs(itag).ifirst = 0;
           staobs(itag).endmjd = sched(isched).scan(iscan).sta(sched(isched).scan(iscan).nsta).endmjd;
           staobs(itag).az     = sched(isched).scan(iscan).sta(sched(isched).scan(iscan).nsta).az;
           staobs(itag).el     = sched(isched).scan(iscan).sta(sched(isched).scan(iscan).nsta).el;
           staobs(itag).ha     = sched(isched).scan(iscan).sta(sched(isched).scan(iscan).nsta).ha;
           staobs(itag).dc     = sched(isched).scan(iscan).sta(sched(isched).scan(iscan).nsta).dc;            
        end    
    end
end


