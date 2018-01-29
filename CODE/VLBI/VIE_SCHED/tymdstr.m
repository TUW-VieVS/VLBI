% Purpose  
%   Convert [yy,mm,dd] to character string.
% History  
%   2010-04-06   Jing SUN   Created
%   


function [cymd] = tymdstr(yy, mm, dd)

% year
cyy = num2str(yy);

% month
switch mm
    case 1
        cmm = 'JAN';
    case 2
        cmm = 'FEB';
    case 3
        cmm = 'MAR';
    case 4
        cmm = 'APR';
    case 5
        cmm = 'MAY';
    case 6 
        cmm = 'JUN';
    case 7
        cmm = 'JUL';
    case 8
        cmm = 'AUG';
    case 9
        cmm = 'SEP';
    case 10
        cmm = 'OCT';
    case 11
        cmm = 'NOV';
    case 12
        cmm = 'DEC';
    otherwise
        cmm = '***';
end

% day
if (dd < 10)
    cdd = strcat('0', num2str(dd));
else
    cdd = num2str(dd);
end

% cymd
cymd = strcat(cyy(3:4), cmm, cdd);


