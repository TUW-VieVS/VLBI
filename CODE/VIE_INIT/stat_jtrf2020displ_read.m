
function data = stat_jtrf2020displ_read(jtrf_file)

%jtrf2020path = '../TRF/data/JTRF/jtrf2020_defining_station_position_xyz_archive_vlbi';
jtrf2020path = '../TRF/data/JTRF/jtrf2020_u2022_station_position_xyz';
fil = [jtrf2020path '/' jtrf_file];
if exist(fil,'file')
    fid = fopen(fil);
    data = textscan(fid,'%f %f%f%f %f%f%f %f%f%f %f%f%f%f%f%f%f %s','CommentStyle','#');
    fclose(fid);
else
    data =cell(1,18);
end

% # COLUMN 1  ...................................... Time Tag [decimal year] 
% # COLUMN 2  ...................................... X Displ [m]     
% # COLUMN 3  ...................................... Y Displ [m]     
% # COLUMN 4  ...................................... Z Displ [m]     
% # COLUMN 5  ...................................... X Stdev [m]     
% # COLUMN 6  ...................................... Y Stdev [m]     
% # COLUMN 7  ...................................... Z Stdev [m]     
% # COLUMN 8  ...................................... XY Correlation  
% # COLUMN 9  ...................................... XZ Correlation  
% # COLUMN 10 ...................................... YZ Correlation  
% # COLUMN 11 ...................................... Time Tag [seconds past J2000] (Leap Second Correction Not Applied) 
% # COLUMN 12 ...................................... Time Tag [year]    
% # COLUMN 13 ...................................... Time Tag [month]   
% # COLUMN 14 ...................................... Time Tag [day]     
% # COLUMN 15 ...................................... Time Tag [hour]    
% # COLUMN 16 ...................................... Time Tag [minutes] 
% # COLUMN 17 ...................................... Time Tag [seconds] 
% # COLUMN 18 ...................................... Data Type Flag 
