% Purpose  
%   Calculate ra and dec of the sun at julian epoch mjd.  
% History  
%   2010-04-06   Jing SUN   created
%   2011-09-14   From DBS, based on AENA formulas.
%                GOOD TO <0.01 DEGREES FROM 1950 TO 2050 (1984 AENA, p. C24)
%


function [sunra, sunde] = ssunpo(mjd)

% NUMBER OF DAYS SINCE J2000.0 
days = mjd - 51544.5;  % [days]
% MEAN SOLAR LONGITUDE
slon = 280.460 + 0.9856474 * days; 
slon = mod(slon, 360);
if (slon < 1.0d-3)
    slon = slon + 360.0;
end
% MEAN ANOMALY OF THE SUN
sanom = 357.528 + 0.9856003 * days; 
sanom = sanom * pi / 180.0;
sanom = mod(sanom, 2*pi);
if (sanom < 1.0d-3)
    sanom = sanom + 2 * pi;
end
% ECLIPTIC LONGITUDE AND OBLIQUITY OF THE ECLIPTIC 
ecllon = slon + 1.915 * sin(sanom) + 0.020 * sin(2*sanom);
ecllon = ecllon * pi / 180.0; 
ecllon = mod(ecllon, 2*pi); 
quad = ecllon / (0.5*pi); 
iquad = floor(1 + quad); 
obliq = 23.439 - 0.0000004 * days; 
obliq = obliq * pi / 180.0; 
% RIGHT ASCENSION AND DECLINATION
% (RA IS IN SAME QUADRANT AS ECLIPTIC LONGITUDE)
sunra = atan(cos(obliq) * tan(ecllon)); 
if (iquad == 2) 
    sunra = sunra + pi; 
elseif (iquad == 3) 
    sunra = sunra + pi; 
elseif (iquad == 4) 
    sunra = sunra + 2*pi; 
end
sunde = asin(sin(obliq) * sin(ecllon)); 


