% Purpose  
%   Read acceleartion information from the acceleration catalog file.
% History  
%   2016-06-30   Matthias Schartner   Created
%
% CHANGES
%   


function [station] = iacceleration(filename, station, PARA);
fid = fopen(filename);
C = textscan(fid,'%s %f %f','Commentstyle','*','MultipleDelimsAsOne',1);
fclose(fid);
name = C{1};
accAz = C{2};
accEl = C{3};

for ista = 1:length(station)
    bool = strcmpi(name,station(ista).name);
    if sum(bool)==1
        station(ista).acc1 = accAz(bool);
        station(ista).acc2 = accEl(bool);
    else
        station(ista).acc1 = PARA.RATE1A;
        station(ista).acc2 = PARA.RATE2A;
    end
    
end

end

