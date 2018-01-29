% ************************************************************************
%   Description:
%   gives fundamental arguments of Lunisolar Nutation (FA(1)-FA(5))
%   and of Planetary Nutation (FA(6)-FA(14)).
% 
%   Reference: 
%   IERS Conv. 2003, Chapt. 5, (40) and (41)
%
%   Input:
%     T                    time vector (t,1)
%                          [julian centuries of TDB/TT since J2000]
%     opt = 1              FA = (t, 5) --> only Lunisolar arguments
%     opt = 2              FA = (t,14) --> all 14 arguments
%     opt = 3              FA = (t,6) --> GMST+pi, lunisolar arguments
%                
%   Output:
%     FA (t,5) or (t,14) or (t,6)   fundamental arguments in [rad]
% 
%   External calls: 	
%     as2rad.m
%
%   Coded for VieVS: 
%   02 Feb 2009 by Lucia Plank
%
%   Revision: 
%   26 May 2010 by Lucia Plank: option 3 added
%
% *************************************************************************
function [FA] = fund_arg (T,opt)

T2 = T.^2;
T3 = T.^3;
T4 = T.^4;
p2 = 2*pi;

i =1;

if opt == 3
%   GMST + pi
      FA(:,1) = (67310.54841 + ...
                (876600*3600 + 8640184.812866) * T        +...
                                      0.093104 * T2       +...
                                      -6.2e-6  * T3) * 15 +...
                                      648000;
      i=2;
end

%   Mean Anomaly of the Moon.
      FA(:,i) =    485868.249036D0       +...
                (1717915923.2178D0 * T  +...
                        (31.8792D0 * T2 +...
                       (0.051635D0 * T3 +...
                     -0.00024470D0 * T4)));
                                            
%   Mean Anomaly of the Sun.
      FA(:,i+1) =  1287104.793048D0       +...
                  (129596581.0481D0 * T  +...
                         (-0.5532D0 * T2 +...
                        (0.000136D0 * T3 +...
                      -0.00001149D0 * T4)));
                         
%   Mean Longitude of the Moon minus Mean Longitude of the Ascending
%   Node of the Moon.
      FA(:,i+2) =    335779.526232D0       +...
                  (1739527262.8478D0 * T  +...
                         (-12.7512D0 * T2 +...
                        (-0.001037D0 * T3 +...
                        0.00000417D0 * T4)));
                         
%   Mean Elongation of the Moon from the Sun.
      FA(:,i+3) =   1072260.703692D0       +...
                  (1602961601.2090D0 * T  +...
                          (-6.3706D0 * T2 +...
                         (0.006593D0 * T3 +...
                       -0.00003169D0 * T4)));
                         
%   Mean Longitude of the Ascending Node of the Moon.
      FA(:,i+4) =  450160.398036D0       +...
                  (-6962890.5431D0 * T  +...
                         (7.4722D0 * T2 +...
                       (0.007702D0 * T3 +...
                     -0.00005939D0 * T4)));
                         
 FA = mod(FA,1296000);
 FA = as2rad(FA);
                        
 if opt == 2
      FA(:, 6) = ( 4.402608842D0 + 2608.7903141574D0 * T );
      FA(:, 7) = ( 3.176146697D0 + 1021.3285546211D0 * T );
      FA(:, 8) = ( 1.753470314D0 +  628.3075849991D0 * T );
      FA(:, 9) = ( 6.203480913D0 +  334.0612426700D0 * T );
      FA(:,10) = ( 0.599546497D0 +   52.9690962641D0 * T );
      FA(:,11) = ( 0.874016757D0 +   21.3299104960D0 * T );
      FA(:,12) = ( 5.481293872D0 +    7.4781598567D0 * T );
      FA(:,13) = ( 5.311886287D0 +    3.8133035638D0 * T );
      FA(:,14) =  0.024381750D0 * T +    0.00000538691D0 * T2;
 end     
      FA = mod(FA,p2);