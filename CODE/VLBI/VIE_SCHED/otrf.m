% Purpose  
%   Write /VieVS/TRF/trf_sched.txt file.
% History  
%   2010-04-06   Jing SUN   Created
%    


function otrf(station, PARA)

% open output file
fn_trf  = [PARA.pvievs 'TRF/' 'trf_sched.txt'];
fid_trf = fopen(fn_trf, 'w'); 

% write trf
stanum = length(station);
for ista = 1 : stanum
    fprintf(fid_trf, '%s     %15.5f %15.5f %15.5f%s\n', station(ista).name(1:8), station(ista).xyz(1:3), '   0.0  0.0 0.0   51544      0  99999');
end

% close output file
fclose(fid_trf);


