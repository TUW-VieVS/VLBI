% -------------------------------------------------------------------------
%
%                              read_parameter_and_cat_files
%
%   This function reads a VEX parameter file (text-format) and ralated VEX 
%   catalog files (*.cat).
%
%   Author: 
%       2013-11-07 : Andreas Hellerschmied (heller182@gmx.at)
%   
%   changes       :
%           
%
%   inputs        :
%   - filepath_vex_para : Filepath of VEX parameter file
%   - filename_vex_para : Filename of VEX parameter file
%   - filepath_vex_cat  : Filepath of VEX catalog files
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
%   - vex_read_cat_file : Reads VEX catalog file.
%   
%
%   references    :
%
%-------------------------------------------------------------------------

function [vex_para, error_code, error_msg] = read_parameter_and_cat_files(filepath_vex_para, filename_vex_para, filepath_vex_cat)

    % init:
    error_code = 0;
    error_msg = '';
    
    % prealocation of "vex_para" struct
    vex_para = struct();
    
    % ##### open vex-parameter file #####
    fid_vex_para = fopen([filepath_vex_para, filename_vex_para], 'r'); 
    if (fid_vex_para == -1)
        error_msg = 'Can not open VEX Parameter File at the given filepath.';
        error_code = 1;
        return;
    end
    
    % ##### Read Data from Parameter File: #####
   [temp_data] = textscan(fid_vex_para, '%s %s', 'CommentStyle', '*', 'Delimiter', '=');
   
   for i_para = 1 : length(temp_data{1})
       
       temp_str1 = temp_data{1}{i_para};
       temp_str2 = temp_data{2}{i_para};
       
       station_name = '';
       
       % ##### Read data from catalog file, if a link is found: #####
       if temp_str1(1) == '>'
           
           % ### Get name of VEX parameter block: ###
            % Init:
            block_ref = '';
            block_name = '';
            
            % Remove ">" from the string:
            temp_str1 = temp_str1(2 : length(temp_str1));
           
            % Write parameter data to struct:
            eval(['vex_para.', temp_str1, '.label = temp_str2;']);
            
            for i_char = 1 : length(temp_str1) 
                if (temp_str1(i_char) == '.')
                    i_temp = i_char;
                    break; 
                end
                station_name(i_char) = temp_str1(i_char);                
            end
           
            for i_char = (i_temp + 1) : length(temp_str1) 
                if (    (temp_str1(i_char) == '(')  || ...
                        (temp_str1(i_char) == '=')          )
                    break; 
                end
                block_name(i_char - i_temp) = temp_str1(i_char);                
            end
            
            for i_char = (i_temp + 1) : length(temp_str1) 
                if (temp_str1(i_char) == '=')
                    break; 
                end
                block_ref(i_char - i_temp) = temp_str1(i_char);                
            end
            
            % Prepare filename of vex cat. file:
            filename_vex_cat = ['vex_', block_name, '.cat'];
            
            % ### Read vex cat file and load data to "vex_para" structure ###
            [vex_para, error_code, error_msg] = vex_read_cat_file(filepath_vex_cat, filename_vex_cat, vex_para, station_name, block_ref);
            if error_code
                error_msg = ['vex_read_cat_file: ', error_msg];
               return; 
            end
            
       else
                    
           % Write parameter data to struct:
           eval(['vex_para.', temp_str1, ' = temp_str2;']);
        end
        
   end % i_para = 1 : length(temp_data{1})
   

    % ##### close vex parameter file #####
    if ( (exist('fid_vex_para', 'var') ) && (fid_vex_para ~= -1) )
        fclose(fid_vex_para);
    end
    

return;

