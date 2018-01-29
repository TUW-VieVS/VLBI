% -------------------------------------------------------------------------
%
%                              function TLEChecksum.m
%
%   Calculates checksum for TLE Line 1 or 2.
%
%   Author: 
%       Andreas Hellerschmied  - 13.10.2013
%   
%   changes       :
%           
%
%   inputs        :
%       tle_line:   String, containing TLE Line 1 or 2.
%     
%
%   outputs       :
%       flag_checksum_OK:  1, if checksum is OK, otherwise 0.
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

function [flag_checksum_OK] = TLEChecksum(tle_line)

    checksum = 0;

    for i = 1 : 68
        
        if ((tle_line(i) >= '0') && (tle_line(i) <= '9'))   % Is figure?
            checksum = checksum + str2num(tle_line(i));
        elseif (tle_line(i) == '-')
            checksum = checksum + 1;
        end
        
    end
    
    if (str2num(tle_line(69)) == mod(checksum, 10) )
        flag_checksum_OK = 1;
    else
        flag_checksum_OK = 0;
    end

return;

