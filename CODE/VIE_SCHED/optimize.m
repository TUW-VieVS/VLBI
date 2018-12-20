% Purpose  
%   Write $OP (automatic scheduling options).
% History  
%   2010-04-06   Jing SUN   Created
%   


function optimize(station, srcnum_s, fid_skd)

fprintf(fid_skd, 'XP F YP F DUT F PSI F EPS F\n');
for ista = 1 : length(station)
    fprintf(fid_skd, '%s AOFF F ARAT F COFF F CRT1 F CRT2 F X F Y F Z F\n', station(ista).po);
end
num = 0;
for isrc = 1 : srcnum_s
    num = num + 1;
    fprintf(fid_skd, '%3d F  ', isrc);
    if (rem(num,10) == 0)
        fprintf(fid_skd, '\n');
    end
end
if (rem(num,10) ~= 0)
    fprintf(fid_skd, '\n');
end

fprintf(fid_skd, 'XP F YP F DUT F PSI F EPS F\n');
for ista = 1 : length(station)
    fprintf(fid_skd, '%s AOFF F ARAT F COFF F CRT1 F CRT2 F X F Y F Z F\n', station(ista).po);
end
num = 0;
for isrc = 1 : srcnum_s
    num = num + 1;
    fprintf(fid_skd, '%3d F  ', isrc);
    if (rem(num,10) == 0)
        fprintf(fid_skd, '\n');
    end
end
if (rem(num,10) ~= 0)
    fprintf(fid_skd, '\n');
end


