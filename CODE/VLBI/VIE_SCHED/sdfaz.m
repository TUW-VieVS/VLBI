% Purpose  
%   Calculate the azimuth difference between the NOW and the NEW source 
%   positions, taking into account cable wrap.
% Input    
%   lim11, lim12                        : antenna lower and upper limits for axis az [rad]
%   unaznow                             : the unambiguous NOW source azimuth [rad]
%   aznew                               : the ambiguous NEW source azimuth [rad]
%   flag_take_alternative_cable_wrap    : if true => Select alternative cable wrap, instead of the nearest one
% Output   
%   dfaz    : azimuth difference [rad]
%   unaznew : the unambiguous NEW source azimuth [rad]
% History  
%   2010-04-06   Jing SUN   created
%   2016-05-11   Matthais SCHARTNER: different computation to improve speed
%   2016-11-02, A. Hellerschmied: Added possibity to select an alternative cable wrap section ("flag_take_alternative_cable_wrap")


function [unaznew, dfaz] = sdfaz(lim11, lim12, unaznow, aznew, varargin)

if nargin == 5
    flag_take_alternative_cable_wrap = varargin{1};
else
    flag_take_alternative_cable_wrap = false;
end

% Get the first valid az value withint the limits:
while (aznew > lim11)
    aznew = aznew - 2 * pi;
end
while (aznew < lim11)
    aznew = aznew + 2 * pi;
end

n = floor((lim12-aznew)/(2*pi))+1;
if n == 1
    unaznew = aznew;
    dfaz = abs(aznew - unaznow);
else % calculate the possibile ambiguity of aznew
    unaznew = ones(n,1)*aznew+(0:n-1)'*2*pi;
    dfaz = abs(unaznew-unaznow);    
    
    % Select cable wrap section with shortest slew path:
    if ~flag_take_alternative_cable_wrap
        [dfaz, id] = min(dfaz);
        unaznew = unaznew(id);
    else
        % Take alternative cable wrap
        [dfaz,b] = sort(dfaz);
        dfaz = dfaz(2);
        unaznew = unaznew(b(2));
    end
end
 
end
