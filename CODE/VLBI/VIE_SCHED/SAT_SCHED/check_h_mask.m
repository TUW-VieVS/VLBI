% -------------------------------------------------------------------------
%
%                              check_h_mask.m
%
%   This function checks, if the target (az, el coordinates) is above the min. elevation angle defined ba the horizontal mask (if available for ths station).
%   Most of the code is taken from the function "zlup.m" of VIE_SCHED V2.3 from J. Sun.
%
%   Author: 
%       Andreas Hellerschmied, 2015-07-29
%   
%   changes       : 
%   - 
%           
%   inputs        :
%       - stat_data         : station data structure
%       - i_stat            : pointer to one station in "stat_data"
%       - az_rad            : Target azimuth [rad]
%       - el_rad            : Target elevation [rad]
%
%   outputs       :
%   - flag_above_min_el     : flag: = 1, if teh target is above the min elevation defined by the horizontal mask 
%    
%
%   locals        :
% 
%
%   coupling      :
%       - 
%   
%
%   references    :
%
%-------------------------------------------------------------------------



function [flag_above_min_el] = check_h_mask(stat_data, i_stat, az_rad, el_rad) 

    hmasknum = stat_data.stat(i_stat).horizontal_mask_num;
    hmask    = stat_data.stat(i_stat).horizontal_mask;
    
    % ##### If a h-mask is available for this station: Check it! #####
    if hmasknum > 0
        % Check the h-mask type:
        if (mod(hmasknum,2) ~= 0)
            hmasktype = 1;   % step functions
        else
            hmasktype = 2;   % line segments
        end
        
        % Loop over h-mask segments and get the elevation angle, defined for the given azimuth:
        for i = 1 : floor(hmasknum/2)   
            ib = i * 2;
            if ((az_rad >= hmask(ib-1)) & (az_rad <= hmask(ib+1)))
                if (hmasktype == 1)
                    elmask = hmask(ib);
                elseif (hmasktype == 2)
                    elmask = ((hmask(ib+2) - hmask(ib)) / (hmask(ib+1) - hmask(ib-1))) * (az_rad - hmask(ib-1)) + hmask(ib);
                end
                break;
            end
        end
        
        
        flag_above_min_el = (el_rad > elmask);
    
    % ##### If there is no h-mask available for this station #####
    else 
        flag_above_min_el = 1;
    end
    
  
return
         
    

