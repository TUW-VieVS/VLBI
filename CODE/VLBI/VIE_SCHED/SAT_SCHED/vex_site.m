% -------------------------------------------------------------------------
%
%                              vex_site
%
%   Writes the $SITE Block to a VEX File.
%       - 
%
%   Author: 
%       2013-11-14 : Andreas Hellerschmied
%   
%   changes       :
%   - 2015-06-18, A. Hellerschmied: Updated for VieVS 2.3
%   - 2015-10-20, A. Hellerschmied: Possibility added to write combined VEX
%               files: Site blocks are written for all stations defined by the input
%               argument "stat_id_list".
%           
%
%   inputs        :
%   - fid_vex           : File ID of VEX File
%   - vex_para          : Data from VEX paramer file, structure
%   - sched_data        : scheduling data structure
%   - stat_id_list      : List of IDs of stations which should be considered
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

function [error_code, error_msg] = vex_site(fid_vex, vex_para, sched_data, stat_id_list)

    % ##### Init #####
    error_code = 0;
    error_msg = '';
    all_labels_str_array = {};
    
    
    % ##### Wite block #####
    if isempty(stat_id_list) % Error
        error_code = 2;
        error_msg = 'Station ID list is empty';  
        return;
    else
        
        fprintf(fid_vex, '$SITE;\n');
        fprintf(fid_vex, '*\n');
        
        % ##### loop over all stations in ID list ##### 
        for i_stat = 1 : length(stat_id_list)
            
            stat_id      = stat_id_list(i_stat);
            % station_name = sched_data.stat(stat_id).name; 
            
            % Define Site:
            try
                fprintf(fid_vex, 'def %s;\n', sched_data.stat(stat_id).site_label);
            catch error
                error_code = 1;
                error_msg = 'Station name not available';  
                return;
            end


            % Site_type:
            try
                fprintf(fid_vex, '    site_type = %s;\n', vex_para.site_type);
            catch error
                error_code = 1;
                error_msg = 'Site Type not available';  
                return;
            end


            % Site_name:
            try
                fprintf(fid_vex, '    site_name = %s;\n', sched_data.stat(stat_id).name);
            catch error
                error_code = 1;
                error_msg = 'Site Name not available';  
                return;
            end


            % site_ID:
            try
                fprintf(fid_vex, '    site_ID = %s;\n', sched_data.stat(stat_id).label);
            catch error
                error_code = 1;
                error_msg = 'Site ID not available';  
                return;
            end


            % Site_Position: 
            try
                fprintf(fid_vex, '    site_position = %6.4f m : %6.4f m : %6.4f m ;\n', sched_data.stat(stat_id).trf_x, sched_data.stat(stat_id).trf_y, sched_data.stat(stat_id).trf_z);
            catch error
                error_code = 1;
                error_msg = 'Site Position not available';  
                return;
            end


            % Site_Velocity:
            try
                fprintf(fid_vex, '    site_velocity = %6.4f m/yr : %6.4f m/yr : %6.4f m/yr;\n', sched_data.stat(stat_id).trf_vx, sched_data.stat(stat_id).trf_vy, sched_data.stat(stat_id).trf_vz);
            catch error
                error_code = 1;
                error_msg = 'Site Velocity not available';  
                return;
            end


            % Site_Position_Epoch:
            try
                fprintf(fid_vex, '    site_position_epoch = %s;\n', num2str(sched_data.stat(stat_id).trf_epoch));
            catch error
                error_code = 1;
                error_msg = 'Site Position Epoch not available';  
                return;
            end


            % Site_Position_Ref:
            try
                fprintf(fid_vex, '    site_position_ref = %s;\n', sched_data.stat(stat_id).site_position_ref);
            catch error
                fprintf(fid_vex, '*   site_position_ref = <= Not available!\n');
            end


            % Additional notes, if available:
            try
                fprintf(fid_vex, '*   %s\n', vex_para.site_note);
            catch error

            end

            fprintf(fid_vex, 'enddef;\n');
            fprintf(fid_vex, '*\n');
            
        end % for i_stat = 1 : length(stat_id_list)

        fprintf(fid_vex, '* ---------------------------------------------------\n');
        fprintf(fid_vex, '*\n');
    
    end % if isempty(stat_id_list)
    
return;

