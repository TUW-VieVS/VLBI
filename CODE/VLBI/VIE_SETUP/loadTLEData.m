% -----------------------------------------------------------------------------
%
%                              function loadTLEData.m
%
%   This function loads TLE files and stores the data as a struct.
%
%   Author: 
%   
%   changes :
%           
%
%   inputs        :
%   - tle_file_name:    Filename and path of TLE file. 
%   - mode         :    Switch between different work-modes.
%       mode = 1    :   Read Satellite name (TLE Line 0) and NORAD number.
%       mode = 2    :   Read complete TLE dataset.
%     
%
%   outputs       :
%   - tle           : Tle dats structure
%   - error_code    : Error Code.
%   - error_msg     : Error message.
%    
%
%   locals        :
% 
%
%   coupling      :
%       TLE_checksum.m  - TLE Checksum calculation and check.
%   
%
%   references    :
%
%--------------------------------------------------------------------------

function [tle, error_code, error_msg] = loadTLEData(tle_file_name, mode)

% preallocation
    error_code = 0;         % 0 = no errors occured
    error_msg = '';         % Error Message. Set, if any error occures.
    n = 0;                  % Total number of TLE datasets in file
    flag_checksum_OK = 0;   % 0 = Checksum check failed.
    
    
    % Open TLE file:
    if ~isempty(tle_file_name)
        fid = fopen(tle_file_name, 'r');
        if (fid == -1)
            %fprintf(1,'ERROR: Failed to open file: %s\n', tle_file_name); 
            error_code = 1;
        end;
    else
        error_code = 2;
    end;
    
    if ~(   (mode == 1) || ...
            (mode == 2)     )
        error_code = 3;
    end;
     
    if (error_code == 0)
        
        while (~feof(fid))

           n=n+1;

           if (feof(fid) == 0) 
               line0 = fgetl(fid);
               tle(n).sat_name = line0(1:24);
           end;

           if (feof(fid) == 0)
               line1 = fgetl(fid);
               
               if ~(length(line1) < 69)
                   
                   if (line1(1) == '1')
               
                       [flag_checksum_OK] = TLEChecksum(line1);

                       if flag_checksum_OK
                           
                           switch(mode)
                               
                               case 1
                                   tle(n).line1.sat_number = line1(3:7);
                                   
                               case 2
                                   tle(n).line1.sat_number = line1(3:7);
                                   tle(n).line1.class = line1(8);
                                   tle(n).line1.int_design_launch_year = line1(10:11);
                                   tle(n).line1.int_design_launch_num = line1(12:14);
                                   tle(n).line1.int_design_launch_piece = line1(15:17);
                                   tle(n).line1.epoch_year = str2num(line1(19:20));
                                   tle(n).line1.epoch_day = str2num(line1(21:32));
                                   tle(n).line1.mean_mot_1st = line1(34:43);
                                   tle(n).line1.mean_mot_2nd = line1(45:52);
                                   tle(n).line1.bstar_drag = line1(54:61);
                                   tle(n).line1.numb_0 = line1(63);
                                   tle(n).line1.El_set_num = line1(65:68);
                                   tle(n).line1.checksum = line1(69);
                           
                           end  % switch(mode) 
                       else
                           error_code = 5;
                       end
                   else
                       error_code = 6;
                   end
               else
                   error_code = 4;
               end
               
           end

           if (feof(fid) == 0)
               line2 = fgetl(fid);
               
               if ~(length(line2) < 69)
                   
                   if (line2(1) == '2')

                       [flag_checksum_OK] = TLEChecksum(line2);

                       if flag_checksum_OK
                           switch(mode)
                               case 1
                                   
                               case 2
                                   tle(n).line2.sat_number = line1(3:7);
                                   tle(n).line2.inc = str2num(line2(9:16));                 %! [degree]
                                   tle(n).line2.RA = str2num(line2(18:25));                 %! [degree]
                                   tle(n).line2.e = str2num(strcat('0.', line2(27:33)));    %! []
                                   tle(n).line2.omega = str2num(line2(35:42));              %! [degree]
                                   tle(n).line2.M = str2num(line2(44:51));                  %! [degree]
                                   tle(n).line2.n = str2num(line2(53:63));                  %! [rev/day]
                                   tle(n).line2.rev_num_epoch = line2(64:68);
                                   tle(n).line2.checksum = line2(69);
                                   
                           end  % switch(mode)
                       else
                           error_code = 5;
                       end
                   else
                       error_code = 7;
                   end
               else
                   error_code = 4;
               end
               
           end
           
           if (error_code ~= 0) % ERROR => Break!
               break;
           end

        end %while
        
    end;
    
    if (error_code ~= 0) % In case of ERROR => Assignme Error Message!
        
        tle = 0;    % Error initialization.

        if (error_code == 1)
            error_msg = 'Failed to open TLE file.';

        elseif (error_code == 2)
            error_msg = 'Filename is missing.';

        elseif (error_code == 3)
            error_msg = 'Unknown work-mode.';

        elseif (error_code == 4)
            error_msg = 'TLE Line is too short (< 69 characters).';
            
        elseif (error_code == 5)
            error_msg = 'TLE Ckecksum error.';
            
        elseif (error_code == 6)
            error_msg = 'TLE Line 1 does not start with "1".';
            
        elseif (error_code == 7)
            error_msg = 'TLE Line 2 does not start with "2".';

        end

    end
    
    if ( (exist('fid', 'var') ) && (fid ~= -1) )
        fclose(fid);
    end
    

return;
