% -------------------------------------------------------------------------
%
%                              function eci2topo.m
%
%   This function calculates Antenna pointing data in terms of azimut-,
%   elevation- angls and range, topocentric declination and right ascension,
%   local hour angle, plus all the correspondig change rates.
%   Note: The reference ellipsoid is WGS72 (hard coded!).
%   For scheduling purposes.
%
%   Author: 
%       Andreas Hellerschmied, 22.9.2013
%   
%   changes       :
%   - 2013-10-25 : Andras Hellerschmied. Calculation of rates implemented.
%   - 2014-01-18 : A. Hellerschmied: Calculation of topocentric Ra/Dev and
%       their change Rates according to D. Vallado (2007), utilizing the 
%       matlab procedure "rv2radec.m".
%   - 2014-01-18 : A. Hellerschmied: Added calculation of the local hour
%       angle. Didn't check yet, if the values and the calc.-approach is
%       correct!!!!
%   - 2014-01-21 : A. Hellerschmied: tRa - rate calculation corrected!
%   - 2014-01-21 : A. Hellerschmied: local LHa rate calc. corrected!
%   - 2014-01-21 : A. Hellerschmied: Use azel2hadec.m to calc ha/dec + rates  
%           
%
%   inputs        :
%   - jd                : Julian Dateof the considered epoch.
%   - X_sat_eci_km      : Position Vector in the ECI System [km]
%   - V_sat_eci_km      : Velocity Vector in the ECI System [km/sec]
%   - stat_lon          : Station longitude [deg]
%   - stat_lat          : Geodetic Station latitude [deg]
%   - stat_alt          : Station altitude (ellipsoid height) [m]
%     
%
%   outputs       :
%   - az                : Azimuth (Station -> Source) [deg]
%   - el                : Elevation (Station -> Source) [deg]
%   - range_s_obs       : Range (Station -> Source) [m]
%   - az_rate           : Azimuth rate (Station -> Source) [deg/sec]
%   - el_rate           : Elevation rate (Station -> Source) [deg/sec]
%   - r_rate            : Range rate (Station -> Source) [m/sec]
%   - tra               : Topocentric right ascension [deg]
%   - tdec              : Topocentric declination [deg]
%   - tra_rate          : Topocentric right ascension rate [deg/sec]
%   - tdec_rate         : Topocentric declination rate [deg/sec]
%   - local_hour_angle  : Local Hour Angle [deg]
%   - local_hour_angle_rate : Local Hour Angle change rate [deg/sec]
%    
%
%   locals        :
% 
%
%   coupling      :
%   - JD_2_GMST : Conversion from JD to Greenwich mean sidereal time.
%   - rv2radec.m : Cal. of: Ra/Dec and their Rates from satellite state
%       vectors (position and velocity vector).
%   - azel2hadec.m : Transformation Az/El (+ rates) => ha/dec (+ rates).
%   
%
%   references    :
%   - T.S. Kelso, Orbital Coordinate Systems, Part I, Satellite Times,
%     Sept/Oct, 1995
%   - T.S. Kelso, Orbital Coordinate Systems, Part II, Satellite Times,
%     Nov/Dec, 1995
%   - T.S. Kelso, Orbital Coordinate Systems, Part III, Satellite Times,
%     Jan/Feb, 1996
%   - Orbital Mechanics with Matlab, Oct. 2013,
%     Web: http://www.cdeagle.com/html/ommatlab.html
%   - Escobal, P., 1965, Methods of orbit determination, Wiley, p. 29, pp.
%       397-398.
%   - D. Vallado, 2007, Fundamentals of Astrodynamics and Applications, 3rd
%     Edition, Space Technology Library, pp. 263-264.
%
%-------------------------------------------------------------------------


function [az, el, range_s_obs, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate] = eci2topo(jd, X_sat_eci_km, V_sat_eci_km, stat_lon, stat_lat, stat_alt)

    % Earth rotaion rate [rad/sec], T.S. Kelso, 1995, Satellite Times,
    % Orbital Coordinate Systems, Part I
    omeg = 7.29211510e-005;
    
    % Init
    deg_2_rad = pi()/180;
    rad_2_deg = 180/pi();
    
    stat_lat = stat_lat * deg_2_rad;
    stat_lon = stat_lon * deg_2_rad;

    sin_stat_lat = sin(stat_lat);
    cos_stat_lat = cos(stat_lat);
    
    % Ellipsoid data:
    a_WGS72 = 6378135;
    f = 1/298.26;

    % Greenwich Mean Siderial Time [rad]
    GMST = JD_2_GMST(jd);

    C = 1/(sqrt(1 + f * (f - 2) * sin_stat_lat^2));
    S = (1 - f)^2 * C;

    % Local Mean Siderial Time [rad]
    LMST = mod(GMST + stat_lon, (2*pi));
   
    sin_LMST = sin(LMST);
    cos_LMST = cos(LMST);

    % ECI station vector [m]
%     X_obs = (a_WGS72 + stat_alt) *   [  C * cos_stat_lat * cos_LMST; ...  
%                                          C * cos_stat_lat * sin_LMST;...
%                                          S * sin_stat_lat                              ];
                                                               
    % Ref.: Escobal, P., 1965, Methods of orbit determination, p. 29
    G1 = a_WGS72 * C + stat_alt;
    G2 = a_WGS72 * S + stat_alt;
                                     
    X_obs = [  G1 * cos_stat_lat * cos_LMST; ...  
               G1 * cos_stat_lat * sin_LMST;...
               G2 * sin_stat_lat                              ];
            
                                     
    % Satellite vector (ECI) [m]
    X_sat = X_sat_eci_km' .* 10^3;

    % Range vector (ECI) [m]
    X_r = X_sat - X_obs;

    % Topocentric - Horizon Coordinate System (south, east, zenith) [m]
    
    % Transformation Matrix:  
    T = [   sin_stat_lat*cos_LMST,    sin_stat_lat*sin_LMST , -cos_stat_lat    ;
            -sin_LMST            ,    cos_LMST              , 0                ;
            cos_stat_lat*cos_LMST,    cos_stat_lat*sin_LMST , sin_stat_lat         ];
        
    X_top = T * X_r;
    
%     X_top = [   sin(stat_lat)*cos(LMST)*X_r(1) + sin(stat_lat)*sin(LMST)*X_r(2) - cos(stat_lat)*X_r(3);
%                 -sin(LMST)*X_r(1) + cos(LMST)*X_r(2);
%                 cos(stat_lat)*cos(LMST)*X_r(1) + cos(stat_lat)*sin(LMST)*X_r(2) + sin(stat_lat)*X_r(3)  ];        

    % Range
    range_s_obs = norm(X_top);

    % Elevation
    el_rad = asin(X_top(3)/range_s_obs);    % -69.214187110489718
    el = el_rad * rad_2_deg;

    % Azimuth (North-Azimuth!)
    az = (atan2(-1*X_top(2), X_top(1)) * rad_2_deg) + 180; % 3.040554500625361e+002
  

    % ECI velocity vector of the satellite [m/sec] 
    V_sat_ECI = V_sat_eci_km' .* 10^3;

    % ECI range rate vector [m/sec]
    V_ECI = V_sat_ECI - cross([0;0;omeg], X_sat);

    % Topocentric range rate vector [m/sec]
%     V_top = [   sin(stat_lat)*cos(LMST)*V_ECI(1) + sin(stat_lat)*sin(LMST)*V_ECI(2) - cos(stat_lat)*V_ECI(3);
%                 -sin(LMST)*V_ECI(1) + cos(LMST)*V_ECI(2);
%                 cos(stat_lat)*cos(LMST)*V_ECI(1) + cos(stat_lat)*sin(LMST)*V_ECI(2) + sin(stat_lat)*V_ECI(3)  ];  

    V_top = T * V_ECI;

    % range rate [m/sec] 
    r_rate = dot(X_top,V_top)/norm(X_r);
    % equal to: r_rate_1 = (X_top' * V_top)/norm(X_r);

    % Azimuth rate
    az_rate = (V_top(1)*X_top(2) - V_top(2)*X_top(1)) / (X_top(1)^2 + X_top(2)^2);
    
    % Elevation rate
    el_rate =  (V_top(3) - r_rate*sin(el_rad)) / sqrt(X_top(1)^2 + X_top(2)^2);

    az_rate = az_rate * rad_2_deg;
    el_rate = el_rate * rad_2_deg;
        

    % Topocentric Rightaszension and Declination
%     X_r = X_r./norm(X_r);
%     dec = asin(X_r(3,1));
%     cos_dec = sqrt(1 - X_r(3)^2);
%     ra = acos(X_r(1) / cos_dec);

    % Calculation of Ra/Dec [P. Escobal (1965), pp. 397-398.]:
    % atan2 => quadrant problem resolved!
%     tra = atan2(X_r(2), X_r(1));
%     tdec = atan2(X_r(3), sqrt(X_r(1)^2 + X_r(2)^2));
    
    
    % ##### Calculation of topocentric Ra/Dec and their Rates: #####
    % (Using function from D. Vallado (2007): rv2radec.m)
    X_r_km = X_r / 1000; % [km]
    V_ECI_km = V_ECI / 1000; % [km/sec]
    
    V = (V_sat_ECI - cross([0;0;omeg], X_obs)) / 1000;
    
    [rr, tra, tdec, drr, tra_rate, tdec_rate] = rv2radec(X_r_km, V);
    %[rr, tra, tdec, drr, tra_rate, tdec_rate] = rv2radec(X_r_km, V_ECI_km);
    
    % rr => range [km]
    % drr => range rate [km/sec]
    
    tra = tra * rad_2_deg;
    tdec = tdec * rad_2_deg;
    tra_rate = tra_rate * rad_2_deg;
    tdec_rate = tdec_rate * rad_2_deg;

    
    % ##### Calculation of Local Hour Angle + Rate #####
    
    % Local hour angle [deg]:
    local_hour_angle = LMST * rad_2_deg - tra;
    
    % LHA => [180 deg, +180 deg]
    while (local_hour_angle > 180)
        local_hour_angle = local_hour_angle - 360;
    end
    while (local_hour_angle < -180)
        local_hour_angle = local_hour_angle + 360;
    end
    
    % local hour angle rate = omega - Ra_rate [deg/sec]:
    local_hour_angle_rate = omeg * rad_2_deg - tra_rate;


    
% !!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!
    % Following function gives the same results for LHA + rate as the
    % values calculated above. Therefore the function invocation is not
    % necessaryly needed in this place. It is just used to cross-check the
    % results from above!!!!
    %stat_lat = stat_lat * rad_2_deg;
    %[ha, tdec1, ha_rate, tdec1_rate] = azel2hadec(az, el, az_rate, el_rate, stat_lat);
    
% !!!!!!!!! Test output: !!!!!!!!!!!!!!!
    %fprintf(1,'Ra= %f, rate= %f; LHA= %f, rate= %f, tdec= %f, rate= %f, tdec1= %f, rate= %f , ha= %f, rate= %f\n',tra, tra_rate, local_hour_angle, local_hour_angle_rate, tdec, tdec_rate, tdec1, tdec1_rate, ha, ha_rate);
    
  
    
return;





