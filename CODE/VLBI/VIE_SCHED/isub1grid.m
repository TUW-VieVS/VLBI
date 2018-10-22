% Purpose  
%   Analyze the 1-source subnet. 
% History  
%   2010-04-06   Jing SUN   Created
%   


function [catpair] = isub1grid(cat)

gridnum = length(cat);
for igrid1 = 1 : gridnum
    catpair(igrid1).pair(1) = igrid1;
end


