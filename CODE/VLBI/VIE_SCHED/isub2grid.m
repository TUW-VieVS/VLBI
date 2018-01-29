% Purpose  
%   Analyze the 2-source subnet. 
% History  
%   2010-04-06   Jing SUN   Created
%   


function [catpair] = isub2grid(cat)

gridnum = length(cat);
for igrid1 = 1 : gridnum
    % rangeid1
    rangeid1 = igrid1;
    ra1 = cat(rangeid1).ramean;
    de1 = cat(rangeid1).demean;
    % rangeid2
    ra2 = ra1 + pi; 
    if (ra2 >= 2*pi)
        ra2 = ra2 - 2*pi;
    end
    de2 = -1*de1;
    for i = 1 : gridnum
        if ((ra2>=cat(i).ramin)&(ra2<cat(i).ramax)&(de2>=cat(i).demin)&(de2<cat(i).demax))
            rangeid2 = i;
            break;
        end
    end
    % catpair
    catpair(igrid1).pair(1) = rangeid1;
    catpair(igrid1).pair(2) = rangeid2;
end


