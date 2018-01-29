% checks if the used STAR Parameters in para.txt are valid. 
% 20.05.16 M. Schartner

function [ PARA ] = checkStarSetup( station, PARA )

cadence = PARA.CADENCE;
if ~isnumeric(cadence) || cadence<0
    error('Parameter CADENCE in param.txt is not valid. Pleas use a numerical input >1')
end

strongant = PARA.STRONGANT;
sname = {station.name};
bool = strcmp(strongant,sname);
if ~any(bool)
    error('Use a STRONGANT in para.txt that is part of your selected station network')
end

if PARA.FILLINMODE ~= 0
    warning('FILLIN-mode is not supportet for STAR-mode jet!\n')
    PARA.FILLINMODE = 0;
end

end

