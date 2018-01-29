% -------------------------------------------------------------------------
%
%                              vex_station
%
%   Writes the $STATION Block to a VEX File.
%       - 
%
%   Author: 
%       2013-11-14 : Andreas Hellerschmied (heller182@gmx.at)
%   
%   changes       :
%   - 2015-06-16, A. Hellerschmied: Updated for VieVS 2.3
%   - 2015-10-19, A. Hellerschmied: Added possibility to print several station definitions within the $STATION block
%           
%
%   inputs        :
%   - fid_vex           : File ID of VEX File
%   - vex_para          : Data from VEX paramer file, structure
%   - sched_data        : scheduling data structure
%   - stat_id_list      : List of IDs of stations which should be treated
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

function [error_code, error_msg] = vex_station(fid_vex, vex_para, sched_data, stat_id_list)

    % Init
    error_code = 0;
    error_msg = '';
    
    if isempty(stat_id_list) % Error
        error_code = 2;
        error_msg = 'Station ID list is empty';  
        return;
    else
        fprintf(fid_vex, '$STATION;\n');
        fprintf(fid_vex, '*\n');

        % ##### loop over all stations in ID list ##### 
        for i_stat = 1 : length(stat_id_list)

            stat_id      = stat_id_list(i_stat);
            station_name = sched_data.stat(stat_id).name;

            % Define STATION:
            try
                fprintf(fid_vex, 'def %s;\n', sched_data.stat(stat_id).station_label);
            catch error
                error_code = 1;
                error_msg = 'Station label not available';  
                return;
            end

            % ref $SITE:
            try
                fprintf(fid_vex, '    ref $SITE = %s;\n', sched_data.stat(stat_id).site_label);
            catch error
                error_code = 1;
                error_msg = 'Site label not available';  
                return;
            end

            % ref $ANTENNA:
            try
                fprintf(fid_vex, '    ref $ANTENNA = %s;\n', sched_data.stat(stat_id).antenna_label);
            catch error
                error_code = 1;
                error_msg = 'Antenna label not available';  
                return;
            end


            % ref $DAS:
            try
                temp_str = eval(['vex_para.', station_name , '.das(1).label;']);
                fprintf(fid_vex, '    ref $DAS = %s;\n', temp_str);
            catch error
                error_code = 1;
                error_msg = 'DAS label not available';  
                return;
            end


            % Additional notes, if available:
            try
                fprintf(fid_vex, '*   %s\n', vex_para.station_note);
            catch error

            end

            fprintf(fid_vex, 'enddef;\n');
            
            if i_stat ~= length(stat_id_list)
               fprintf(fid_vex, '*\n');
            end


        end % for i_stat = 1 : length(stat_id_list)

        fprintf(fid_vex, '* ---------------------------------------------------\n');
        fprintf(fid_vex, '*\n');
        
        
    end % if isempty(stat_id_list)
    
return;

