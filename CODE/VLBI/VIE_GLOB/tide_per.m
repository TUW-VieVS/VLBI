% ************************************************************************
%   Description:
%   calculates the periods of tidal terms
% 
%   Input:
%     N                    Matrix with multipliers of fundamental arguments
%                
%   Output:
%     per                  Line vector with tidal periods in hours
%
%   Coded for VieVS: 
%   23 Nov 2015 by Sigrid Boehm
%
% *************************************************************************
function [per] = tide_per(N)

%   GMST + pi
      FA(:,1) = (876600*3600 + 8640184.812866)*15;

%   Mean Anomaly of the Moon.
      FA(:,2) = 1717915923.2178;
                                            
%   Mean Anomaly of the Sun.
      FA(:,3) = 129596581.0481;
                         
%   Argument of latitude of the Moon
      FA(:,4) = 1739527262.8478;
                         
%   Mean Elongation of the Moon from the Sun.
      FA(:,5) = 1602961601.2090;
                         
%   Mean Longitude of the Ascending Node of the Moon.
      FA(:,6) = -6962890.5431;
                         
 freq = FA*N';
 per = 360./(freq/36525/3600)*24;
 