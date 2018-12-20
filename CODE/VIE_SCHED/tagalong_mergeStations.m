% merges the special station to station structure for tag-along mode
%
% CREATED: 21.02.17 - Matthias Schartner

function [ station ] = tagalong_mergeStations( station_main, station_tagalong )

station = struct();
stanum = 0;
for ista = 1 : length(station_main)
    stanum = stanum + 1;
    if stanum == 1
        station = station_main(ista);
    else
        station(stanum) = station_main(ista);
    end
end
for ista = 1 : length(station_tagalong)
    stanum = stanum + 1;
    station(stanum) = station_tagalong(ista);
end

end

