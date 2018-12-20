% -------------------------------------------------------------------------
%
%                              procedure calc_sundist.m
%
%   Caluculation of the angular distance between a source given by azimuth
%   (src_az) and elevation (src_el) in a topocentric system at a station,
%   defined by latitude and longitude, at a certain epoch given by JD
%   (jd).
%   The calculations here follow simple models and gives results with 
%   limited accuracy!
%
%   Author: 
%       Andreas Hellerschmied (heller182@gmx.at) - 2013.10.23
%   
%   changes       :
%           
%
%   inputs        :
%   - jd        : Julain Date of the considered epoch.
%   - stat_lon  : Sation longitude [deg].
%   - stat_lat  : Station latitude [deg].
%   - src_az    : Source Azimuth at the Station (stat_lon, stat_lat) [deg].
%   - src_el    : Source Elevation at the Station (stat_lon, stat_lat) [deg].
%     
%
%   outputs       :
%    - sun_dist : Angle between unity vectors (station-source) and (station-sun) [deg];
%
%   locals        :
% 
%
%   coupling      :
%   - ssunpo.m      : Calculation of the sun position (Ra/Dec) at MJD.
%   - zazel_s.m     : Calculate [az,el] and [ha,dc] with a simple model.
%
%   references    :
%
%-------------------------------------------------------------------------

function [sun_dist] = calc_sundist(jd, stat_lon, stat_lat, src_az, src_el)

    % Init
    rad2deg = 180/pi;
    deg2rad = pi/180;

    % Unit Conversions
    stat_lon = stat_lon * deg2rad;
    stat_lat = stat_lat * deg2rad;
    src_az = src_az * deg2rad;
    src_el = src_el * deg2rad;
    mjd = jd - 2.400000500000000e+006;

    % Calculation of the sun position (Ra/Dec) at MJD
    [sun_ra, sun_de] = ssunpo(mjd);
    
    % Calculation of Az/El of the sun at the given station (stat_lon, stat_lat)
    [sun_az, sun_el] = zazel_s(mjd, stat_lon, stat_lat, sun_ra, sun_de);
    
    % Angular distance between sun (v_sun) and source (v_src)
    
    v_src = [cos(src_el)*cos(src_az);cos(src_el)*sin(src_az);sin(src_el)]; % Unity vector station-source
    v_sun = [cos(sun_el)*cos(sun_az);cos(sun_el)*sin(sun_az);sin(sun_el)]; % Unity vector station-sun

    sun_dist = acos((v_src'*v_sun)); % / (norm(v_src) * norm(v_sun))); % Point product
    sun_dist = sun_dist * rad2deg;
    
return;

