% -------------------------------------------------------------------------
%
%                              jdutc2ttt
%
%   Calculation of the Julan centuries of Terrestrial Time
%
%   Author: 
%       Andreas Hellerschmied (2015-02-19)
%   
%   changes       :
%           
%
%   inputs        :
%   - jdutc :   Julian date of UTC [d]
%     
%
%   outputs       :
%    
%
%   locals        :
% 
%
%   coupling      :
%   - tai_utc   : Find the number of leap seconds between TA1 and UTC
%   
%   
%   references    :
%   - D. Vallada, 2007, Fundaments of Astrodynamics and Applications,
%     Springer, p. 197, Equ. 3-47.
%
%-------------------------------------------------------------------------

function [ttt] = jdutc2ttt(jdutc)

    sec2day = 1 / (24*60*60);
    
    mjdutc = jdutc - 2.400000500000000e+006;
    [leap_seconds] = tai_utc(mjdutc);
    
    % tai = utc + leap_seconds
    % tt = tai + 32.184s
    
    jdtt = jdutc + leap_seconds * sec2day + 32.184 * sec2day;

    ttt = (jdtt - 2451545) / 36525;
    
return;

