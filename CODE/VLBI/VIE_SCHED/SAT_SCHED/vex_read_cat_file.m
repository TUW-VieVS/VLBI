% -------------------------------------------------------------------------
%
%                              vex_read_cat_file
%
%   This function reads a VEX catalog file (text-format).
%
%   Author: 
%       2013-11-17 : Andreas Hellerschmied (andreas.hellerschmied@geo.tuwien.ac.at)
%   
%   changes       :
%           
%
%   inputs        :
%   - filepath_vex_cat : Filepath of VEX catalog file
%   - filename_vex_cat : Filename of VEX catalog file
%   - vex_para         : VEX parameter structure
%   - station_name     : Name of associated station
%   - block_ref        : Desigantion of the associated VEX $BLOCK
%     
%
%   outputs       :
%   - error_code        : Error Code (0 = no erros occured)
%   - error_msg         : Error Message (empty, if no errors occured)
%   - vex_para          : VEX parameter structure
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

function [vex_para, error_code, error_msg] = vex_read_cat_file(filepath_vex_cat, filename_vex_cat, vex_para, station_name, block_ref)

    % init:
    error_code = 0;
    error_msg = '';

    flag_found_entry = 0;
    
    parameter_label = eval(['vex_para.', station_name, '.', block_ref, '.label']);
    
    % open .cat file
    fid_vex_cat = fopen([filepath_vex_cat, filename_vex_cat], 'r'); 
    if (fid_vex_cat == -1)
        error_msg = 'Can not open VEX Catalog File at the given filepath.';
        error_code = 1;
        return;
    end

    % Read Data from Catalog File:
    [temp_data] = textscan(fid_vex_cat, '%s %s', 'CommentStyle', '*', 'Delimiter', '=');
    
    for i_para = 1 : length(temp_data{1})
       
        temp_str1 = temp_data{1}{i_para};
        temp_str2 = temp_data{2}{i_para};
        
        length_temp_str1 = length(temp_str1);
        
        % Find seeked entry with the right label in cat-data:
        if (    (temp_str1(1) == '>')                                       && ...
                (strcmp(temp_str1(2 : length_temp_str1), parameter_label))          )
            
            flag_found_entry = 1;
            
        elseif (flag_found_entry == 1)
            
            % Next entry in cat-data => break => finished!
            if (    (temp_str1(1) == '>')                                       && ...
                    (~strcmp(temp_str1(2 : length_temp_str1), parameter_label))          )
                break;
            end
            
            % Write parameter data to struct:
            try
            eval(['vex_para.',station_name, '.', block_ref, '.', temp_str1, '.str= temp_str2;']);
            catch error
                disp('error');
            end
            
        end
       
    end % i_para = 1 : length(temp_data{1})

    % close vex Catalog file
    if ( (exist('fid_vex_cat', 'var') ) && (fid_vex_cat ~= -1) )
        fclose(fid_vex_cat);
    end
    

return;

