% Purpose  
%   Calculate sky coverage number at the station.        
% History  
%   2010-04-06   Jing SUN   created
%


function [cat] = zcoverage(az, el)

% az(0--2*pi)
while (az < 0)
    az = az + 2 * pi;
end
while (az > 2*pi)
    az = az - 2 * pi;
end

% sky coverage(1--13)
if (el >= (60*pi/180))
    cat = 1; 
elseif ((el >= (30*pi/180)) & (el < (60*pi/180))) 
    if (az < (90*pi/180))
        cat = 2;
    elseif ((az >= (90*pi/180)) & (az < (180*pi/180)))
        cat = 3;
    elseif ((az >= (180*pi/180)) & (az < (270*pi/180)))
        cat = 4;
    elseif ((az >= (270*pi/180)) & (az < (360*pi/180)))
        cat = 5;
    end
elseif (el < (30*pi/180))
    if (az < (45*pi/180))
        cat = 6;
    elseif ((az >= (45*pi/180)) & (az < (90*pi/180)))
        cat = 7;
    elseif ((az >= (90*pi/180)) & (az < (135*pi/180)))
        cat = 8;
    elseif ((az >= (135*pi/180)) & (az < (180*pi/180)))
        cat = 9;
    elseif ((az >= (180*pi/180)) & (az < (225*pi/180)))
        cat = 10;
    elseif ((az >= (225*pi/180)) & (az < (270*pi/180)))
        cat = 11;
    elseif ((az >= (270*pi/180)) & (az < (315*pi/180)))
        cat = 12;
    elseif ((az >= (315*pi/180)) & (az < (360*pi/180)))
        cat = 13;
    end
end


