% -------------------------------------------------------------------------
%
%                              procedure JD_2_GMST.m
%
%   This function calculates the Greenwich mean siderial time (GMST) of the
%   given Julian Date (JD).
%
%   Author: 
%       Andreas Hellerschmied, 16.9.2013
%   
%   changes       :
%           
%
%   inputs        :
%               JD - Julan Date 
%
%   outputs       :
%               GMST - Greenwich Mean Siderial Time [rad]
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


function [GMST] = JD_2_GMST(JD)

    UT = mod(JD + 0.5, 1);
    JD = JD - UT;
    TU = (JD - 2451545.0) / 36525;
    GMST = (- 6.2e-6 * TU^3+ 0.093104 * TU^2 + 8640184.812866 * TU + 24110.54841);
    GMST = mod(GMST + 86400 * 1.00273790934 * UT, 86400);
    GMST = GMST * 2 * pi() / 86400;
    
return;

