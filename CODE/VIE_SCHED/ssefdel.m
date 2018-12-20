% Purpose  
%   Adjust the zenith SEFD by elevation.
% Input   
%   el        : elevation [rad]
%   sefd      : zenith SEFD
%   y, c0, c1 : SEFD parameters
% Output   
%   sefdel : adjusted SEFD
% History  
%   2010-04-06   Jing SUN   created
%


function [sefdel] = ssefdel(el, sefd, y, c0, c1)

if (y > 1.0d-3)   % adjust the zenith SEFD
    xel = 1.0 / ((sin(el)) ^ y);
    fac = c0 * (xel ^ 0) + c1 * (xel ^ 1);
    sefdel = sefd * fac;
    if (sefdel < sefd)   %%%
        sefdel = sefd;
    end
else              % use zenith value
    sefdel = sefd;
end


