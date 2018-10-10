% -------------------------------------------------------------------------
%
%                              eci2tradec2000.m
%
%   This function calculates Antenna pointing data in terms of topocentric
%   Ra. and Dec. values for a given station-source (satellite)
%   constellation. 
%   Therefor the input ITRF station coordintes are transformed to the
%   celestioal system (ECI J2000.0), whereby an approximated transformation
%   matrix (t2c) is used.
%
%   Author: 
%       Andreas Hellerschmied, 6.11.2013
%   
%   changes       :
%   - 2014-01-14 : Andreas Hellerschmied: Resolved a quadrant problem at
%     calculating the topocentric right ascension.
%           
%   inputs        :
%   - jd                    : Julian Date of the considered epoch (UTC).
%   - X_sat_eci_j2000_km    : Position Position Vector in the ECI J2000.0 System [km]
%   - stat_x                : ITRF x-coordinate [m]
%   - stat_y                : ITRF y-coordinate [m]
%   - stat_z                : ITRF z-coordinate [m]
%   - eop_data              : EOP data structure
%     
%   outputs       :
%   - tra            : Topocentric right ascension (J2000.0) [deg]
%   - tdec           : Topocentric declination (J2000.0)[deg]
%    
%   locals        :
%
%   coupling      :
%   - ctrs2crs_sched(mjd) : Finds the approximated transformation matrix
%     from the ITRF to the ECI J2000.0 frame. 
%       -  
%   references    :
%   - D. Vallado, 2007, Fundamentals of Astrodaynamics and Applications, Third
%     Edition, Springer, NY, pp. 236-240.
%   - P. Escobal, 1965(?), Methods of Orbit Determination, John Wiley and
%     Sons, pp. 397-398.
%
%-------------------------------------------------------------------------


function [tra, tdec] = eci2tradec2000(jdutc, X_sat_eci_j2000_km, stat_x, stat_y, stat_z, eop_data)
    
    % Init
    % deg_2_rad = pi()/180;
    rad_2_deg = 180/pi();
    
    % ITRF Station vector
    X_stat = [stat_x; stat_y; stat_z];
    
    mjdutc = jdutc - 2.400000500000000e+006;
    
    % #### Approximative transformation ####
    %[t2c_approx] = ctrs2crs_sched(mjdutc);
    
    
    % #### Rigorous transformation approach, following the IERS conventions chapter 5 ####

	% A priori EOP values are determined with linear interpolation between
	% the value of midnight before and after observation time.
	% for a session from 18:00 to 18:00 this means, that there are 2 a
	% priori lines and a break at midnight

    % linear interpolation of EOP values for the required epoch:
	dut1  = interp1(eop_data.mjd, eop_data.dut1, mjdutc, 'linear', 'extrap');  
	xp    = interp1(eop_data.mjd, eop_data.xp, mjdutc, 'linear', 'extrap');  
	yp    = interp1(eop_data.mjd, eop_data.yp, mjdutc, 'linear', 'extrap');  
	dX    = interp1(eop_data.mjd, eop_data.dX, mjdutc, 'linear', 'extrap');
	dY    = interp1(eop_data.mjd, eop_data.dY, mjdutc, 'linear', 'extrap');
    
    % Define nutation model:
    nutmod = 'IAU_2006/2000A'; % 'IAU_2000A' or 'IAU_2006/2000A' 

    % Get transformation matix:
    [t2c,dQdx,dQdy,dQdut,dQdX,dQdY,X,Y,era] = trs2crs(mjdutc,xp,yp,dut1,dX,dY,nutmod);
    
    % Transformation: TRF => ECI J2000.0 [m]
    X_stat_eci_j2000 = t2c * X_stat;
                    
    
    % Satellite vector (ECI J2000.0) [km] => [m]
    X_sat_eci_j2000_m = X_sat_eci_j2000_km' .* 10^3;

    % Range vector (ECI J2000.0) [m]
    X_r_j2000 = X_sat_eci_j2000_m - X_stat_eci_j2000;
    
    % Topocentric Rightaszension and Declination (J2000.0)
     X_r_j2000 = X_r_j2000./norm(X_r_j2000);
     
%     tdec2 = asin(X_r_j2000(3,1));
%     cos_dec = sqrt(1 - X_r_j2000(3)^2);
%     tra2 = acos(X_r_j2000(1) / cos_dec);
% In the Ra/Dec calculation approach above the quadrant of the Ra is not 
% unambiguously defined! 
    
    % Calculation of Ra/Dec [P. Escobal (1965)]:
    % atan2 => quadrant problem resolved!
    tra = atan2(X_r_j2000(2), X_r_j2000(1));
    tdec = atan2(X_r_j2000(3), sqrt(X_r_j2000(1)^2 + X_r_j2000(2)^2));
    
    tdec = tdec * rad_2_deg;
    tra = tra * rad_2_deg;

return;





