% Purpose  
%   Calculate start time of the scan using iteration.        
% History  
%   2010-04-06   Jing SUN   Created
%   2016-05-10	Matthias Schartner: improve speed useing:
%                                   - logical short-circuiting
%                                   - whole structure copying (line 113)
%   2016-06-14  Matthias Schartner: - changed the way deleted stations and
%                                     scans are handled
%                                   - vectorised zazel_s.m
%	2016-09-13 Matthias Scharner:   - bugfix
%	2016-11-05 Matthias Scharner:   - bugfix


function [subcon_s] = sstartmjd(source, station, twin, staobs, subconin, PARA)

% initialization
for iscan = 1 : subconin.nscan
    subconin.scan(iscan).startmjd = 0.0;
    for ista = 1 : subconin.scan(iscan).nsta 
        subconin.scan(iscan).sta(ista).startmjd = 0.0;
    end
end

% iteration
subcon = subconin;
for i = 1 : 2
    % calculate slewing time of each station
    for iscan = 1 : subcon.nscan
        nsta = subcon.scan(iscan).nsta;
        srcid = subcon.scan(iscan).srcid;
        ra = source(srcid).ra;
        de = source(srcid).de;
        
        staid = [subcon.scan(iscan).sta.staid];
        llh = [station(staid).llh];
        lon = llh(1:3:end);
        lat = llh(2:3:end);
        
        bool = [subcon.scan(iscan).sta.startmjd]>1;
        mjd = zeros(size(bool));
        if any(bool)
            mjd(bool) = [subcon.scan(iscan).startmjd];
        end
        mjd(~bool) = [staobs(~bool).endmjd];
        [az, el, ha, dc] = zazel_s(mjd, lon, lat, ra, de);    
        for ista = 1:nsta
            lup = zlup(mjd(ista), az(ista), el(ista), ha(ista), dc, staid(ista), station, PARA.MIN_CUTEL);
            if lup
                [st, unaz] = sslew(station, staobs, az(ista), el(ista), ha(ista), dc, staid(ista), PARA);
            else
                st = 999999;
                unaz = az(ista);
            end
            if (staobs(staid(ista)).ifirst == 1)
                subcon.scan(iscan).sta(ista).slewtime = 0;      % [sec]
                subcon.scan(iscan).sta(ista).startmjd = PARA.startmjd;
            else
                subcon.scan(iscan).sta(ista).slewtime = st;     % [sec]
                subcon.scan(iscan).sta(ista).startmjd = staobs(staid(ista)).endmjd + st/86400.0 + PARA.SOURCE/86400 + PARA.TAPETM/86400 + PARA.IDLE/86400 + PARA.CALIBRATION/86400;  
            end
            subcon.scan(iscan).sta(ista).az = unaz;
            subcon.scan(iscan).sta(ista).el = el(ista);
            subcon.scan(iscan).sta(ista).ha = ha(ista);
            subcon.scan(iscan).sta(ista).dc = dc;  
        end
    end
    
    %David Mayer - check if down time was set at the beginning of the
    %session; If so set beginning of scan for particular station to
    %earliest start of other stations
    for iscan = 1 : subcon.nscan
        for ista = 1 : subcon.scan(iscan).nsta
            staid = subcon.scan(iscan).sta(ista).staid;
            if (staobs(staid).ifirst == 1) && (min(station(staid).downstart) == PARA.startmjd)
                subcon.scan(iscan).sta(ista).startmjd = 9999999; % ????????????????????????????????????????
                subcon.scan(iscan).sta(ista).startmjd = min([subcon.scan(iscan).sta.startmjd]);
                %disp('temp');
            end
        end
    end
    % check the station arrived late
    iflatenum   = 0;
    iflatesn(1) = 0;
    for iscan = 1 : subcon.nscan
        if (source(subcon.scan(iscan).srcid).ifp==1)
            continue;
        end
        for ista = 1 : subcon.scan(iscan).nsta
            smjd(ista) = subcon.scan(iscan).sta(ista).startmjd;
        end
        [x] = sort(smjd);
        % tref = median(x(1:subcon.scan(iscan).nsta));
        if (subcon.scan(iscan).nsta >= 2)
           tref = x(2);   
        else
           tref = x(1);
        end
        clear smjd;
        for ista = 1 : subcon.scan(iscan).nsta
            if((subcon.scan(iscan).sta(ista).startmjd-tref)*86400 > PARA.MAX_WAIT)
                iflatenum = iflatenum + 1;
                iflatesn(iflatenum) = subcon.scan(iscan).sta(ista).staid;
            end
        end
    end
    
    % check maximum slewing time and calculate start time of the scan
    subcon_s = [];
    subcon_s.nscan = 0;
    for iscan = 1 : subcon.nscan
        this_nscan = subcon_s.nscan+1;
        
        subcon_s.nscan = this_nscan;
        subcon_s.scan(this_nscan).srcid = subcon.scan(iscan).srcid;
        subcon_s.scan(this_nscan).startmjd = 0.0;
        subcon_s.scan(this_nscan).nsta = 0.0;
        for ista = 1 : subcon.scan(iscan).nsta
            staid = subcon.scan(iscan).sta(ista).staid;
            if ((subcon.scan(iscan).sta(ista).slewtime<=PARA.MAXSLEWTIME) && (size(find(iflatesn(1:iflatenum)==staid),2)== 0))
                this_station = subcon_s.scan(this_nscan).nsta + 1;
                
                subcon_s.scan(this_nscan).nsta = this_station;
                subcon_s.scan(this_nscan).sta(this_station)=subcon.scan(iscan).sta(ista);
                subcon_s.scan(this_nscan).sta(this_station).duration = 0.0;
                subcon_s.scan(this_nscan).sta(this_station).endmjd   = 0.0;
                
                if (subcon_s.scan(this_nscan).sta(this_station).startmjd > subcon_s.scan(this_nscan).startmjd)
                    subcon_s.scan(this_nscan).startmjd = subcon_s.scan(this_nscan).sta(this_station).startmjd;
                end
            end
        end
        if (subcon_s.scan(this_nscan).nsta < 2)
            subcon_s.scan(this_nscan) =[];
            subcon_s.nscan = subcon_s.nscan - 1;
        end
    end

%     iteration
%     if (i==3 && ~isequal(subcon, subcon_s))
%         fprintf('not equal')
%     end
    

    clear subcon;
    subcon = subcon_s;
        %%
%     for ii = 1: subcon_s.nscan
%         if subcon_s.scan(ii).nsta ~= size(subcon_s.scan(ii).sta,2)
%             fprintf('Problem!\n');
%         end
%     end

end


