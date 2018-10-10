% -------------------------------------------------------------------------
%
%                              vex_antenna
%
%   Writes the $ANTENNA Block to a VEX File.
%       - 
%
%   Author: 
%       2013-11-14 : Andreas Hellerschmied
%   
%   changes       :
%   - 2015-06-18, A. Hellerschmied: Updated for VieVS 2.3
%   - 2015-06-29, A. Hellerschmied: Cable wrap sector can now be written to the VEX file
%   - 2015-10-20, A. Hellerschmied: Possibility added to write combined VEX
%               files: Antenna blocks are written for all stations defined by the input
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

function [error_code, error_msg] = vex_antenna(fid_vex, vex_para, sched_data, stat_id_list)

    % ##### Init #####
    error_code = 0;
    error_msg = '';
    rad2deg = 180/pi;
    
    
    % ##### Wite block #####
    if isempty(stat_id_list) % Error
        error_code = 2;
        error_msg = 'Station ID list is empty';  
        return;
    else
        
        fprintf(fid_vex, '$ANTENNA;\n');
        fprintf(fid_vex, '*\n');
        
        
        % ##### loop over all stations in ID list ##### 
        for i_stat = 1 : length(stat_id_list)
            
            stat_id      = stat_id_list(i_stat);
            % station_name = sched_data.stat(stat_id).name; 
            
            % Define Antenna:
            try
                fprintf(fid_vex, 'def %s;\n', sched_data.stat(stat_id).antenna_label);
            catch error
                error_code = 1;
                error_msg = 'Station Name not available';  
                return;
            end

            % axis_type:
            try
                fprintf(fid_vex, '    axis_type = %s : %s;\n', sched_data.stat(stat_id).axis_type_1, sched_data.stat(stat_id).axis_type_2);
            catch error
                error_code = 1;
                error_msg = 'Axis Type not available';  
                return;
            end

            % antenna_motion, axis 1:
            try
                fprintf(fid_vex, '    antenna_motion = %s : %5.1f deg/min : %1.0f sec;\n', sched_data.stat(stat_id).axis_type_1, sched_data.stat(stat_id).axis_1_slew_rate, sched_data.stat(stat_id).axis_1_settling_time);
            catch error
                error_code = 1;
                error_msg = 'Antenna motion data for axis 1 not available';  
                return;
            end


            % antenna_motion, axis 2:
            try
                fprintf(fid_vex, '    antenna_motion = %s : %5.1f deg/min : %1.0f sec;\n', sched_data.stat(stat_id).axis_type_2, sched_data.stat(stat_id).axis_2_slew_rate, sched_data.stat(stat_id).axis_2_settling_time);
            catch error
                error_code = 1;
                error_msg = 'Antenna motion data for axis 2 not available';  
                return;
            end


            % axis_offset: 
            try
                fprintf(fid_vex, '    axis_offset = %6.5f m;\n', sched_data.stat(stat_id).axis_offset);
            catch error
                error_code = 1;
                error_msg = 'Axis offset not available';  
                return;
            end


            % Cable wrap pointing sector (only for AZEL antenna mount type and "write_cable_wrap_sector" flag in VEX parameter file = 1)
            if (strncmpi(sched_data.stat(stat_id).axis_type_1, 'az', 2) && strncmpi(sched_data.stat(stat_id).axis_type_2, 'el', 2)) && strncmpi(vex_para.write_cable_wrap_sector, '1', 1)

                % pointing_sector (counter-clockwise, &ccw):
                try
                    fprintf(fid_vex, '    pointing_sector = &ccw, : %s : %6.1f deg : %6.1f deg : %s : %6.1f deg : %6.1f deg;\n', sched_data.stat(stat_id).axis_type_1, sched_data.stat(stat_id).az_ccw1*rad2deg, sched_data.stat(stat_id).az_ccw2*rad2deg, sched_data.stat(stat_id).axis_type_2, sched_data.stat(stat_id).lim21*rad2deg, sched_data.stat(stat_id).lim22*rad2deg);
                catch error
                    error_code = 1;
                    error_msg = 'Pointing_sector (counter-clockwise) not available';  
                    return;
                end

                % pointing_sector (neutral, &n):
                try
                    fprintf(fid_vex, '    pointing_sector = &n,   : %s : %6.1f deg : %6.1f deg : %s : %6.1f deg : %6.1f deg;\n', sched_data.stat(stat_id).axis_type_1, sched_data.stat(stat_id).az_n1*rad2deg, sched_data.stat(stat_id).az_n2*rad2deg, sched_data.stat(stat_id).axis_type_2, sched_data.stat(stat_id).lim21*rad2deg, sched_data.stat(stat_id).lim22*rad2deg);
                catch error
                    error_code = 1;
                    error_msg = 'Pointing_sector (neutral) not available';  
                    return;
                end

                % pointing_sector (clockwise, &cw):
                try
                    fprintf(fid_vex, '    pointing_sector = &cw,  : %s : %6.1f deg : %6.1f deg : %s : %6.1f deg : %6.1f deg;\n', sched_data.stat(stat_id).axis_type_1, sched_data.stat(stat_id).az_cw1*rad2deg, sched_data.stat(stat_id).az_cw2*rad2deg, sched_data.stat(stat_id).axis_type_2, sched_data.stat(stat_id).lim21*rad2deg, sched_data.stat(stat_id).lim22*rad2deg);
                catch error
                    error_code = 1;
                    error_msg = 'Pointing_sector (clockwise) not available';  
                    return;
                end

            end


            % Additional notes, if available:
            try
                fprintf(fid_vex, '*   %s\n', vex_para.antenna_note);
            catch error

            end

            fprintf(fid_vex, 'enddef;\n');
            fprintf(fid_vex, '*\n');
    
    
        end % for i_stat = 1 : length(stat_id_list)

        fprintf(fid_vex, '* ---------------------------------------------------\n');
        fprintf(fid_vex, '*\n');
    
    end % if isempty(stat_id_list)
    
return;

