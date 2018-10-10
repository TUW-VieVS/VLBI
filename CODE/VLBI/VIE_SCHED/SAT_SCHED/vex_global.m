% -------------------------------------------------------------------------
%
%                              vex_global
%
%   Writes the $GLOBAL Block to a VEX File.
%   The VEX $GLOBAL Block contains:
%       - 
%
%   Author: 
%       2013-11-12 : Andreas Hellerschmied (heller182@gmx.at)
%   
%   changes       :
%   - 2015-06-16, A. Hellerschmied: Updated for VieVS 2.3
%           
%
%   inputs        :
%   - fid_vex           : File ID of VEX File
%   - vex_para          : Data from VEX paramer file, structure
%   - sched_data        : scheduling data structure
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

function [error_code, error_msg] = vex_global(fid_vex, vex_para, sched_data)

    % Init
    error_code = 0;
    error_msg = '';
    
    fprintf(fid_vex, '$GLOBAL;\n');
    fprintf(fid_vex, '*\n');
    
    % Reference to $EXPER Block:
    try
        fprintf(fid_vex, 'ref $EXPER = %s;\n', sched_data.experiment_label);
    catch error
        error_code = 1;
        error_msg = 'Experiment Label not available';  
        return;
    end
    
    
    % Additional notes, if available:
    try
        fprintf(fid_vex, '*    %s;\n', vex_para.global_note);
    catch error

    end

    fprintf(fid_vex, '* ---------------------------------------------------\n');
    fprintf(fid_vex, '*\n');
    
return;

