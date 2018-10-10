% -------------------------------------------------------------------------
%
%                              vex_header
%
%   Writes Header Lines to VEX File.
%   The VEX Header contains:
%       - Software Version Descriptions
%       - VEX Format definition
%
%   Author: 
%       2013-11-07 : Andreas Hellerschmied (heller182@gmx.at)
%   
%   changes       :
%   - 2015-06-16, A. Hellerschmied: Updated for VieVS 2.3
%           
%
%   inputs        :
%   - fid_vex           : File ID of VEX File
%   - vex_para          : Data from VEX paramer file, structure
%     
%
%   outputs       :
%   - error_code        : Error Code (0 = no erros occured)
%   - error_msg         : Error Message (empty, if no errors occured)
%    
%
%   locals        :
% 
%
%   coupling      :
%   
%
%   references    :
%   - VEX File Definition/Example, Rev 1.5b1, 30. Jan. 2002, available
%     online at: http://www.vlbi.org/vex/
%   - VEX Parameter Tables, Rev 1.5b1, 29. Jan. 2002, available
%     online at: http://www.vlbi.org/vex/
%
%-------------------------------------------------------------------------

function [error_code, error_msg] = vex_header(fid_vex, vex_para)

    % Init
    error_code = 0;
    error_msg = '';
    
    % VEX Format definition:
    try
        fprintf(fid_vex, 'VEX_rev = %s;\n', vex_para.VEX_rev);
    catch error
        error_code = 1;
        error_msg = 'VEX_rev not available';  
        return;
    end
    
    % VieVS Release name/number:
    try
        fprintf(fid_vex, '*    VieVS version: %s\n', vex_para.vievs_version);
    catch error
        error_code = 1;
        error_msg = 'VieVS version not available';  
        return;
    end
    
    % Additional Header notes, if available:
    try
        fprintf(fid_vex, '*    %s\n', vex_para.header_note);
    catch error

    end

    fprintf(fid_vex, '* ---------------------------------------------------\n');
    fprintf(fid_vex, '*\n');
    
return;

