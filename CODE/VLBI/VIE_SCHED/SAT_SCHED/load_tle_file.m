% #########################################################################
% #     load_tle_file
% #########################################################################
%
% DESCRIPTION
%   Loads TLE files.
%
% CREATED  
%   2015-07-20     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
%
% INPUT
%   - PARA                : cheduling parameter structure
%
%
% OUTPUT
%   - PARA                : cheduling parameter structure (updated with TLE data)
%   - error_code          : Error Code (0 = no erros occured)
%   - error_msg           : Error Message (empty, if no errors occured)
%
% CHANGES:
%

function [PARA, error_code, error_msg] = load_tle_file(PARA)



% ##### Init #####
error_code = 0;
error_msg = '';
i_sat = 0;


% % ##### Preallocation #####
% tle_data        = struct('sat', []);
% tle_data.sat    = struct('header_str', [], 'line_1_str', [], 'line_2_str', []);


% ##### open TLE file: ##### 
fid_tle = fopen([PARA.TLE_FILEPATH, PARA.TLE_FILENAME], 'r'); 
if (fid_tle == -1)
    error_msg = ['Can not open TLE file: ', PARA.TLE_FILEPATH, PARA.TLE_FILENAME];
    error_code = 1;
    return;
end



% ##### Read TLE file #####
while (~feof(fid_tle))
    
    % TLE dataset count:
    i_sat = i_sat + 1;
    
    % Get data (3 lines = one complete TLE dataset)
    PARA.tle_data.sat(i_sat).header_str  = fgets(fid_tle, 30); % Header line
    PARA.tle_data.sat(i_sat).line_1_str  = fgets(fid_tle, 71); % Line 1
    PARA.tle_data.sat(i_sat).line_2_str  = fgets(fid_tle, 71); % Line 2

end


% ##### close TLE file: #####
if ( (exist('fid_tle', 'var') ) && (fid_tle ~= -1) )
    fclose(fid_tle);
end


return;