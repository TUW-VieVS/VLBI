% Purpose  
%   Get rise and set of sources at stations.
% History  
%   2010-04-06   Jing SUN   Created
%


function [riseset] = iriseset(source, station, PARA)

% initialize the riseset
srcnum = length(source);
stanum = length(station);
for ista = 1 : stanum
    for isrc = 1 : srcnum
        riseset(ista,isrc,1:PARA.MAX_RISESET) = 0.0;   
    end
end
    
% start and end of the session
startmjd = PARA.startmjd;
endmjd   = PARA.endmjd;

% calculate rise/set
for ista = 1 : stanum
    lon = station(ista).llh(1); 
    lat = station(ista).llh(2);
    for isrc = 1 : srcnum
        ra = source(isrc).ra;
        de = source(isrc).de;
        % visibility at the start
        mjd = startmjd;
        [az, el, ha, dc] = zazel_s(mjd, lon, lat, ra, de);
        [upold] = zlup(mjd, az, el, ha, dc, ista, station, PARA.MIN_CUTEL);  
        % upold  
        if (upold)
            rsnum = 1;
            rstime(rsnum,1) = mjd;
            rstime(rsnum,2) = 0.0;
        else
            rsnum = 0;
        end
        % loop from startmjd to endmjd
        while (mjd <= endmjd)
            mjd = mjd + 1 / 24;
            [az, el, ha, dc] = zazel_s(mjd, lon, lat, ra, de);
            [upnow] = zlup(mjd, az, el, ha, dc, ista, station, PARA.MIN_CUTEL); 
            % upnow
            if (upnow == upold)
                continue;
            end
            mjd = mjd - 1 / 24;
            for im = 1 : 59
                mjd = mjd + 1 / 1440;     
                [az, el, ha, dc] = zazel_s(mjd, lon, lat, ra, de);
                [upnow] = zlup(mjd, az, el, ha, dc, ista, station, PARA.MIN_CUTEL);  
                % upnow
                if (upnow == upold)
                    continue;
                end
                if((upnow == true) & (upold  == false))
                    rsnum = rsnum + 1;
                    rstime(rsnum,1) = mjd;
                    rstime(rsnum,2) = 0.0;
                elseif((upnow == false) & (upold  == true))
                    rstime(rsnum,2) = mjd - 1 / 1440;
                end
                upold = upnow;
                break;                        
            end % for im = 1 : 59            
        end % while (mjd <= endmjd)          
        % never set at endmjd
        if ((rsnum > 0) & (rstime(rsnum,2) <= 1.0d-3))
            rstime(rsnum,2) = mjd;
        end     
        % save riseset variable
        for i = 1 : rsnum
            riseset(ista,isrc,(i-1)*2+1) = rstime(i,1);   % [minute]
            riseset(ista,isrc,(i-1)*2+2) = rstime(i,2);   % [minute]
        end
    end
end


