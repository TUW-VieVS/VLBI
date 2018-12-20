% -------------------------------------------------------------------------
%
%                              read_satellite_frequ_file
%
%   This function reads a satelite frequency file (text-format).
%
%   Author: 
%       2013-11-12 : Andreas Hellerschmied (heller182@gmx.at)
%   
%   changes       :
%   - 2015-06-16, A. Hellerschmied: 2nd frequency added to "sat_frequ" structure
%           
%
%   inputs        :
%   - filepath_sat_requ : Filepath of satellite frequency file
%   - filename_sat_frequ: Filename of satellite frequency file
%     
%
%   outputs       :
%   - error_code        : Error Code (0 = no erros occured)
%   - error_msg         : Error Message (empty, if no errors occured)
%   - sat_frequ         : Satellite Frequency structure
%    
%
%   locals        :
% 
%
%   coupling      :
%   
%
%   references    :
%
%-------------------------------------------------------------------------

function [sat_frequ, error_code, error_msg] = read_satellite_frequ_file(filepath_sat_frequ, filename_sat_frequ)

    % init:
    error_code = 0;
    error_msg = '';
    
    % prealocation of "sat_frequ" struct
    sat_frequ = struct('sat', []);
    sat_frequ.sat = struct('name', [], 'frequ', []);
    
    % open satellite frequency file
    fid_sat_frequ = fopen([filepath_sat_frequ, filename_sat_frequ], 'r'); 
    if (fid_sat_frequ == -1)
        error_msg = 'Can not open satellite frequency setup File at the given filepath.';
        error_code = 1;
        return;
    end
    
    
   % Read in 3 columns: <sat name> <frequ 1> <frequ 2>
   [temp_data] = textscan(fid_sat_frequ, '%s %s %s', 'CommentStyle', '*', 'Delimiter', {'=', ';'});
   
   % Write data to struct:
   for i_sat = 1 : length(temp_data{1})
       sat_frequ.sat(i_sat).name = temp_data{1}{i_sat};
       sat_frequ.sat(i_sat).frequ_1 = str2num(temp_data{2}{i_sat});
       sat_frequ.sat(i_sat).frequ_2 = str2num(temp_data{3}{i_sat});
   end;
   

    % close satellite frequency file
    if ( (exist('fid_sat_frequ', 'var') ) && (fid_sat_frequ ~= -1) )
        fclose(fid_sat_frequ);
    end
    

return;

