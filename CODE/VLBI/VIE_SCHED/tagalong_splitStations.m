% excludes the special station for tag-along mode
%
% CREATED: 21.02.17 - Matthias Schartner

function [ station_main, station_tagalong ] = tagalong_splitStations( station, tagalongname )

stanum = length(station);
stanum1 = 0;
stanum2 = 0;
station_main = [];
station_tagalong = [];
for ista = 1 : stanum
    iftag = 0;
    % Check, if tag-along station is in the station list:
    for iss = 1 : size(tagalongname,1)
        if strcmp(deblank(station(ista).name),deblank(tagalongname(iss,:)))
            iftag = 1;
            break;
        end
    end
    if (iftag == 0) % normal station
        stanum1 = stanum1 + 1;
        if stanum1 == 1
            station_main = station(ista);
        else
            station_main(stanum1) = station(ista);
        end
    else % tag-along station
        stanum2 = stanum2 + 1;
        if stanum2 == 1
            station_tagalong = station(ista);
        else
            station_tagalong(stanum2) = station(ista);
        end
    end
end


end

