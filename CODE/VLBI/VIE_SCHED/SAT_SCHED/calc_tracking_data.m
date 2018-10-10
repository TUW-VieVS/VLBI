% -------------------------------------------------------------------------
%
%                              calc_tracking_data.m
%
%   This procedure calculates the stepwise antenna pointing coordinates in
%   terms of RA/DEC positions in the J2000.0 reference frame.
%
%   Author: 
%       Andreas Hellerschmied, 19.01.2014
%   
%   changes       : 
%   - 2015-06-10 A. Hellerschmied:  - Renamed (previous name: calc_stepwise_antenna_pointing),
%                                   - Function modified to be compatible with the newest version of the sat. sched. module (V2.3)
%   - 2015-07-20: A. Hellerschmied: New TLE/SGP4 orbit propagation function used (tle_propagation.m)
%   - 2015-08-24: A. Hellerschmied: Bug-fix. Bug in errror treatment routine at the function-call of tle_propagation.m
%   - 2016-05-??, A. Hellerschmied: - Calc. tracking points (topo. ra/dec) for stepwise satellite tracking for epoch in the middle of the following repos. interval, if "PARA.CENTER_TRACK_EPOCHS" is set.
%                                   - Waitbar added.
%   - 2016-05-23, A. Hellerschmied: Flag for centering the tracjking epochs can be set in param.txt (PARA.CENTER_TRACK_EPOCHS)
%   - 2016-06-16, A. Hellerschmied: Waitbar is closed, in case there is an error in tle_propagation.m 
%   - 2016-06-20, A. Hellerschmied: - Write warning to CW, if "PARA.CENTER_TRACK_EPOCHS" is missing and give the possibility to set this parameter via CW input.
%                                   - Write un_az and el angles to sched_data.scan.stat.epoch for all antenne amount types, not only for AzEl antennas!
%   - 2016-06-29, A. Hellerschmied: Notification text for shifting the tracking point epochs changed.
%   - 2016-07-15, A. Hellerschmied: Comsetic change, without influencing the results: corrected units [km, m] of r and v vectors 
%   - 2016-10-31, A. Hellerschmied: Bug-fix: sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).un_az was not calculated for the correct epoch, if the CENTER_TRACK_EPOCHS option was used
%   - 2016-11-08, A. Hellerschmied: Option added to use the function "rv2tradc.m" to calc. topo. RaDec values (for test purposes!) => set flag "flag_calc_tradec_vallado" accordingly
%   - 2017-01-09, A. Hellerschmied: Minor changes (lines added to set up Vallado's functions for topo. RaDec calculation)
%   - 2017-01-20, A. Hellerschmied: Vallado functions not used per default for calculating tracking points
%
%   inputs        :
%   - source            : source data structure
%   - stat_data         : station data struct 
%   - sched_data        : Schedule data structure
%   - PARA              : scheduling parameter structure
%
%
%   outputs       :
%   - error_code        : Error Code (0 = no erros occured)
%   - error_msg         : Error Message (empty, if no errors occured)
%   - sched_data        : Scheduling data struct
%    
%
%   locals        :
% 
%
%   coupling      :
%   - tle_propagation
%   - jdut12ttt
%   - teme2eci
%   - eci2tradec2000
%   - invjday
%   - load_eop_finals
%   - sdfaz => Calculation of un_az
%   - eci2topo => az/el
%   
%
%   references    :
%   - David A. Vallado (2007), Fundamentals of Astrodynamics and Applications, 3rd Edition, Springer NY, pp. 236-240.
%
%-------------------------------------------------------------------------


function [sched_data, error_code, error_msg] = calc_tracking_data(stat_data, source, sched_data, PARA)

    % Init
    error_code = 0;
    error_msg = '';
    
    % Option to use function "rv2tradc.m" to calc. topo. RaDec values (for test purposes)
    % => In this case "lod" has to be defined (= input argument of "rv2tradc.m")
    % Results (in [deg]) saved in:
    %    - sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).ra_2
	%    - sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).dec_2
    flag_calc_tradec_vallado = false;
    
    % Prepare input for Vallado's functions, if required:
    %  - Set EOPs (here!)
    if flag_calc_tradec_vallado
        warning('Vallado functions are used to calculate topo. RaDec angles: Please check the LOD value defined in the code!!!"')
        % warning('Vallado functions are used to calculate topo. RaDec angles: Add the path to the directory with these functions!!"')
        terms = 0;
        lod = 0; %1.6487-6; % LOD !!!!!!!!!!!!! => https://datacenter.iers.org/eop/-/somos/5Rgv/latest/9     
        ddpsi = 0;
        ddeps = 0;
        % addpath('../COMPILE/VIE_SCHED_V30/SAT_SCHED/VALLADO') % path to rewuired m-files (by D. Vallado) % Alre4ady done in vie_batch3_0.m
    end
    
    
    disp(' ');
    disp('#### Calculation of stepwise tracking data for the scheduled satellite observations  ####.');
        
    % Center tracking point epochs:
    if isfield(PARA, 'CENTER_TRACK_EPOCHS')
        if PARA.CENTER_TRACK_EPOCHS
            t_shift = PARA.REPOS_INT/(3600*24) / 2; % [years], shifted for half the repos. interval
            fprintf(1, '  => The epochs for calculating the tracking points are shifted by half of the antenna repos. interval (= %1.1f sec)\n', PARA.REPOS_INT/2);
            fprintf(1, '\n');
        else
            fprintf(1, '\n');
        end
    else
        % WARNING
        fprintf(1, '  => WARNING: PARA.CENTER_TRACK_EPOCHS is missing!\n');
        fprintf(1, '     Do you want to shift the epochs for calculating the tracking points by half of the antenna repos. interval (= %1.0f sec)?\n', PARA.REPOS_INT/2);
        input_str = input('     Type "yes", if you want to. Typ anything else, to skip this option: ', 's');
        if strcmp(input_str, 'yes')
            t_shift = PARA.REPOS_INT/(3600*24) / 2; % [years], shifted for hafl the repos. interval
            PARA.CENTER_TRACK_EPOCHS = true;
        else
            PARA.CENTER_TRACK_EPOCHS = false;
        end
        fprintf(1, '\n');
    end
    
    
    % ##### Check, if scans are available in the schedule #####
    if sched_data.number_of_scans < 1
        error_code = 1;
        error_msg = 'No scheduled scans available.';
        return;
    end
    
    
    % ##### Get EOP data for rigorous transformation approach (TRS => CRS) #####
    mjd_start = sched_data.t_nominal_start_jd  - 2.400000500000000e+006;    % JD => MJD
    mjd_end = sched_data.t_nominal_end_jd - 2.400000500000000e+006;         % JD => MJD

    [eop_data.mjd, eop_data.xp, eop_data.yp, eop_data.dut1, eop_data.dX, eop_data.dY, error_code, error_msg] = load_eop_finals(PARA.EOP_FILEPATH, PARA.EOP_FILENAME, mjd_start, mjd_end);
    if error_code ~= 0
        error_msg = ['load_eop_finals: ', error_msg];
        return;
    end
    
    
    % ##### Open Waitbar #####
    h_wb = waitbar(0, '1' ,'Name', 'Calultaion of tracking points'); 
    
    
    % ##### Loop over all scans #####
    for i_scan = 1 : sched_data.number_of_scans
        
        % Update waitbar
        waitbar(i_scan/sched_data.number_of_scans, h_wb, sprintf('Scan %1.0f of %1.0f processed.',i_scan, sched_data.number_of_scans))
        
        % Loop init.:
        un_az_old = 0;
        
        % Satellite or quasar scan
        switch(sched_data.scan(i_scan).obs_type)
            
            % ##########################
            % ##### Satellite scan #####
            % ##########################
            case 'sat'
                
                % Loop over all involved stations:
                for i_stat = 1 : length(sched_data.scan(i_scan).stat)
                    
                    % ### Get station data ###
                    station_id = sched_data.scan(i_scan).stat(i_stat).stat_id;
                    % Get ITRF station location:
                    stat_x = stat_data.stat(station_id).location.TRF.x;
                    stat_y = stat_data.stat(station_id).location.TRF.y;
                    stat_z = stat_data.stat(station_id).location.TRF.z;
                    % long, lat, alt
                    stat_lon = stat_data.stat(station_id).location.ellipsoid.long;
                    stat_lat = stat_data.stat(station_id).location.ellipsoid.lat;
                    stat_alt = stat_data.stat(station_id).location.ellipsoid.altitude;
                    
                    % Loop over all epochs:
                    for i_epoch = 1 : length(sched_data.scan(i_scan).stat(i_stat).epoch)
                        
                        % Get epoch fo calculating the tracking points
                        if PARA.CENTER_TRACK_EPOCHS
                            t_epoch_jd = sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).jd + t_shift;
                        else
                            t_epoch_jd = sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).jd;
                        end

						if flag_calc_tradec_vallado
							jdutc  = t_epoch_jd;
							mjdutc = jdutc - 2.400000500000000e+006;
							dut1  = interp1(eop_data.mjd, eop_data.dut1, mjdutc, 'linear', 'extrap');  
							xp    = interp1(eop_data.mjd, eop_data.xp, mjdutc, 'linear', 'extrap');  
							yp    = interp1(eop_data.mjd, eop_data.yp, mjdutc, 'linear', 'extrap');
							jdut1 = jdutc + dut1 / (60*60*24);
						end
                        % ##### Calculation of satellite position (ECI TEME frame) for this epoch #####
                        [temp_sat_data, error_code, error_msg] = tle_propagation(t_epoch_jd, t_epoch_jd, PARA.REPOS_INT / 60, PARA, sched_data.scan(i_scan).sat_name);
                        % [temp_sat_data, error_code, error_msg] = tle_propagation(jdut1, jdut1, PARA.REPOS_INT / 60, PARA, sched_data.scan(i_scan).sat_name);

                        if error_code ~= 0;
                            error_msg = ['tle_propagation: ', error_msg];
                            % Close waitbar (if it is still there...)
                            try
                                close(h_wb)
                            catch
                            end
                            return;
                        end
                        
                        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Transformation of the satellite position vector ECI TEME => ECI J2000.0
                        %
                        % According to David A. Vallado (2007), Fundamentals of
                        % Astrodynamics and Applications, 3rd Ed., Springer NY, pp. 236-240.
                        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        % Setup for function "teme2eci":
                        eqeterms    = 2;        % number of terms for eqe        0, 2
                        opt         = 'a';      % option for processing          a - complete nutation
                        order       = 106;      % number of terms for nutation   4, 50, 106, ...

                        % JD of epoch => Julian Centuries of Terrestrial Time TT of epoch:
                        [ttt] = jdutc2ttt(t_epoch_jd);
                        
                        % position vector of date - true equator, mean equinox (TEME) [km]:
                        rteme = (temp_sat_data.sat(1).prop(1).r)';  
                        vteme = (temp_sat_data.sat(1).prop(1).v)'; 
                        
                        % TEME of date frame => ECI (J2000.0) frame:
                        [reci, veci] = teme2eci(rteme, vteme, ttt, order, eqeterms, opt); % Function by D.A. Vallado (2007)


                        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Calculation of topocentric Ra/Dec for
                        % the given epoch and station - satellite  constellation:
                        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

						if ~flag_calc_tradec_vallado
							jdutc  = t_epoch_jd;
						end
                        
                        X_sat_eci_J2000_km =  reci';
                        
                        % ECI J2000.0 => topo. Ra/Dec:
                        [tra, tdec] = eci2tradec2000(jdutc, X_sat_eci_J2000_km, stat_x, stat_y, stat_z, eop_data); % Function by D.A. Vallado (2007)
                        
                        % Save topo. Ra/Dec values:
                        sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).ra = tra;
                        sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).dec = tdec;
                        
                        if flag_calc_tradec_vallado
							% Calc topo Radec with function by Vallado:
                        
							[latgc,latgd,lon,hellp] = ijk2ll ([stat_x/1000; stat_y/1000; stat_z/1000]);
							alt = hellp; % alt in [km]
                            
							[rho, trtasc, tdecl, drho, dtrtasc, dtdecl] = rv2tradc(reci, veci, latgd, lon, alt, ttt, jdut1, lod, xp, yp, terms, ddpsi, ddeps);
							
							% Save topo. Ra/Dec values:
							sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).ra_2 = trtasc * 180/pi;
							sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).dec_2 = tdecl * 180/pi;
                        end
                        
                        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Calculation of unambigous azimuth (un_az) for
                        % the given epoch and station - satellite  constellation:
                        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % !!! Only for AzEl-mount antennas !!!
                        % For all other mont types the function sdfaz.m may generates an error!
                        % Not requires for other antenna types than AzEl. un_az only needed to determine the cable wrap section.
                        if strcmp(sched_data.stat(i_stat).axis_type, 'AZEL')

                            % First epoch
                            if i_epoch == 1
                                if PARA.CENTER_TRACK_EPOCHS
                                    % Determine un_az for the shifted time epoch:
                                    [az_new_deg, el_deg] = eci2topo(jdutc, temp_sat_data.sat(1).prop(1).r, temp_sat_data.sat(1).prop(1).v, stat_lon, stat_lat, stat_alt);
                                    [un_az_new_rad] = sdfaz(sched_data.stat(i_stat).lim11, sched_data.stat(i_stat).lim12, sched_data.scan(i_scan).stat(i_stat).start.un_az, az_new_deg * pi/180);
                                    un_az_old                                                   = un_az_new_rad; % [rad]
                                    sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).un_az   = un_az_new_rad; % [rad]
                                else
                                    % Take un_az from scan start epoch:
                                    un_az_old                                                   = sched_data.scan(i_scan).stat(i_stat).start.un_az; % [rad]
                                    sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).un_az   = sched_data.scan(i_scan).stat(i_stat).start.un_az; % [rad]
                                end
                                % Save el [rad]:
                                sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).el = sched_data.scan(i_scan).stat(i_stat).start.el;
                            else
                                % Calculation of topocentric coordinatee:
                                [az_new_deg, el_deg] = eci2topo(jdutc, temp_sat_data.sat(1).prop(1).r, temp_sat_data.sat(1).prop(1).v, stat_lon, stat_lat, stat_alt);

                                % Calc. un_az:
                                [un_az_new_rad] = sdfaz(sched_data.stat(i_stat).lim11, sched_data.stat(i_stat).lim12, un_az_old, az_new_deg * pi/180);

                                % Save un_az/el [rad]:
                                sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).el = el_deg * pi/180;
                                sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).un_az = un_az_new_rad;
                                un_az_old = un_az_new_rad;
                            end
                        else % All other mount types:
                            % Calculation of topocentric coordinatee:
                            [az_deg, el_deg] = eci2topo(jdutc, temp_sat_data.sat(1).prop(1).r, temp_sat_data.sat(1).prop(1).v, stat_lon, stat_lat, stat_alt);
                            
                            % Save un_az/el [rad]:
                            sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).el = el_deg * pi/180;
                            sched_data.scan(i_scan).stat(i_stat).epoch(i_epoch).un_az = az_deg * pi/180;
                        end

                    end % for i_epoch = 1 : length(sched_data.scan(i_scan).stat(i_stat).epoch)
                
                end % for i_stat = 1 : length(sched_data.scan(i_scan).stat)
                
                
                
            % #######################
            % ##### Quasar scan #####
            % #######################
            case 'quasar'
                
                % ### Get Ra/Dec values of the observed source: ###
                quasar_id = sched_data.scan(i_scan).quasar_id;
                ra = source(quasar_id).ra * 180/pi;     % [rad] => [deg]
                dec = source(quasar_id).de * 180/pi;    % [rad] => [deg]
                
                % ### Loop over all involved stations: ###
                for i_stat = 1 : length(sched_data.scan(i_scan).stat)
                    
                    % Save Ra/Dec values:
                    sched_data.scan(i_scan).stat(i_stat).epoch(1).ra = ra;
                    sched_data.scan(i_scan).stat(i_stat).epoch(1).dec = dec;
                    
                    % Save un_az/el [rad]:
                    sched_data.scan(i_scan).stat(i_stat).epoch(1).un_az = sched_data.scan(i_scan).stat(i_stat).start.un_az;
                    sched_data.scan(i_scan).stat(i_stat).epoch(1).el = sched_data.scan(i_scan).stat(i_stat).start.el;
                    
                end
                
                
        end % switch(sched_data.scan(i_scan).obs_type)
        
    end % for i_scan = 1 : sched_data.number_of_scans
    
    % ##### Close waitbar (if it is still there...) #####
    try
        close(h_wb)
    catch
    end
    
    % Save Ra/Dec epoch:
    sched_data.ra_dec_epoch = 'J2000.0';
    
return;
