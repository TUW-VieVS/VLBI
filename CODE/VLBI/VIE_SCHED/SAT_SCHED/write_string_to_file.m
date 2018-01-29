% #########################################################################
% #     write_string_to_file
% #########################################################################
%
% DESCRIPTION
%   This function writes an input string to a defined file.
%
%
% CREATED  
%   2015-08-28     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
%
%
% INPUT
% - out_string              - schedule data structure
% - filepath_out            - Filepath for output file
% - filename_out            - Filename for output file
%
%
% OUTPUT
% - error_code              - Error Code (0 = no erros occured)
% - error_msg               - Error Message (empty, if no errors occured)
%
% CHANGES:
%
function [error_code, error_msg] = write_string_to_file(out_string, filepath_out, filename_out)
 
    % ##### Init #####
    error_code = 0;
    error_msg = '';
    
    % ##### Initial checks #####
    % Is output string empty?
    if isempty(out_string)
        error_msg = 'String is empty!';
        error_code = 2;
        return;
    end
    
    
    % ##### Open file #####
    fid = fopen([filepath_out, filename_out], 'w');
    if (fid == -1)
        error_msg = ['Can not create/open file at the defined filepath: ', filepath_out, filename_out];
        error_code = 2;
        return;
    end
    
    
    % ##### Write string to file #####
    fprintf(fid, '%s', out_string);
    
    
    % ##### Close file #####
    if ( (exist('fid', 'var') ) && (fid ~= -1) )
        fclose(fid);
    end

return;
