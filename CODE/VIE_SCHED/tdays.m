% Purpose  
%   Convert [yy,mm,dd] to days of year in string format.
% History  
%   2010-04-06   Jing SUN   Created
%


function [cdays] = tdays(yy, mm, dd)

x1 = datenum(yy, 1, 1);
x2 = datenum(yy, mm, dd);
days = (x2 - x1) + 1;

if (days < 10)
    cdays = strcat('00', num2str(days));
elseif (days >= 10) & (days < 100)
    cdays = strcat('0', num2str(days));
elseif (days >= 100) 
    cdays = num2str(days);
end


