% This function trims the schedule to the selected timespan.
%
% CREATED: 21.02.17 - Matthias Schartner
%
% CHANGELOG: 2017-03-10 Matthias Schartner bugfix
%            2017-05-29: M. Schartner: Bugfix 

function [ sched ] = sched_trim( sched, PARA )

while 1
    removedSomething = false;
    iscanNew = 1;

    for iscan = 1:sched(end).nscan
        if sched(end).scan(iscanNew).startmjd + max([sched(end).scan(iscanNew).sta.duration])/86400 > PARA.endmjd
            sched(end).scan(iscanNew) = [];
            sched(end).nscan = sched(end).nscan-1;
            removedSomething = true;
        else
            iscanNew = iscanNew+1;
        end
    end

    if sched(end).nscan <1
        sched(end) = [];
    end

    if ~removedSomething
        break
    end
end

end
