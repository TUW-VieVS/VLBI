% #########################################################################
% #     read_tle
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
% - 2022-04-27, H. Wolf: changed it in order to use it in order to load tle
% data as a-priori orbit information
% - 2024-12-16, H. Wolf: read TLE lines without blanks

function [TLE, error_code, error_msg] = read_tle(satOrbitFilePath, satOrbitFileName)
    % definition of path of the file containing the TLE data of the satellite
    TLE.TLE_FILENAME = satOrbitFileName;
    TLE.TLE_FILEPATH = satOrbitFilePath;                                
    
    % SGP4 prediction parameter setup:
    TLE.SGP4_GRAV_CONST            = 72;  % WGS72
    TLE.SGP4_VERIFICATION_MODE     = 0;
    
    error_code = 0;
    error_msg = '';
    iSat = 0;
           
    % open TLE file
    fid_tle = fopen([satOrbitFilePath, satOrbitFileName], 'r'); 
    if (fid_tle == -1)
        error_msg = ['Can not open TLE file: ', TLE.TLE_FILEPATH, TLE.TLE_FILENAME];
        error_code = 1;
        return;
    end
    
    % read TLE file
    while (~feof(fid_tle))
        iSat = iSat + 1;
        TLE.tle_data.sat(iSat).header_str  = deblank(fgets(fid_tle)); % Header line
        TLE.tle_data.sat(iSat).line_1_str  = deblank(fgets(fid_tle)); % Line 1
        TLE.tle_data.sat(iSat).line_2_str  = deblank(fgets(fid_tle)); % Line 2
    end
    
    if ( (exist('fid_tle', 'var') ) && (fid_tle ~= -1) )
        fclose(fid_tle);
    end
end