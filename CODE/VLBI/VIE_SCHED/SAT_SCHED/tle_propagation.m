% #########################################################################
% #     tle_propagation
% #########################################################################
%
% DESCRIPTION
%   Prapagation of satellite orbits using the SGP4 models and TLE datasets
%
%   Options:
%   - A ASCII Textfile with time-sereis ofcalculated ECI satellite coordinates is written, 
%     if the write_file flag is set.
%   - There is a verification mode to check the correctness of the SGP4 results
%
% CREATED  
%   2013-07-20     Andreas Hellerschmied
%
% REFERENCES
%   - SGP4/SDP4 orbit propagation models (e.g. D. Vallado et al., 2006, Revisited Spacetrack Report Number 3)
%
%
% COUPLING
%   - twoline2rv
%   - rv2coe
%
% INPUT
%   - PARA                : scheduling parameter structure
%   - sat_names_to_prop   : Names of satellites which should be treated
%
%
% OUTPUT
%   - sat_data            : Structure containing all calculated data + auxilliary information
%   - error_code          : Error Code (0 = no erros occured)
%   - error_msg           : Error Message (empty, if no errors occured)
%
% CHANGES:
%   - 2016-06-16, A. Hellerschmied: Error msg corrected, in case TLE data for a specific satellite is not available in PARA
%


function [sat_data, error_code, error_msg] = tle_propagation(t_start, t_stop, delta_t, PARA, sat_names_to_prop)

% ##### preallocating #####
sat_data = struct('sat', [], 'prop_setup', []);
sat_data.sat = struct('sat_number', [], 'prop', [], 'epoch', [], 'TLE_header_line', []);
sat_data.sat.epoch = struct('r', [], 'v', [], 't_since', [], 'year', [],...
    'mon', [], 'day', [], 'h', [], 'min', [], 'sec', [], 'jd', []);
sat_data.sat.prop = struct('r', [], 'v', [], 't_since', [], 'year', [],...
    'mon', [], 'day', [], 'h', [], 'min', [], 'sec', [], 'jd', []);
sat_data.prop_setup = struct('delta_t_min', [], 'start_epoch', [], 'stop_epoch', [],...
    'tle_filepath', [], 'tle_filename', []);
sat_data.prop_setup.start_epoch = struct('year', [],...
    'mon', [], 'day', [], 'h', [], 'min', [], 'sec', [], 'jd', []);
sat_data.prop_setup.stop_epoch = struct('year', [],...
    'mon', [], 'day', [], 'h', [], 'min', [], 'sec', [], 'jd', []);

% ##### Const. #####
mjd2jd  = 2.400000500000000e+006;
rad     = 180.0 / pi;

% ##### init #####
error_code = 0;                 % 0 = No error
error_msg = '';                 % '' = No error
i_sat = 0;                      % counts the number of processed satellites/TLE datasets
i_epoch = 0;                    % counts the number of processed epochs per satelite/TLE dataset
flag_tle_dataset_found = 0;    
header_str = '';
numb_of_sats_to_prop = 0;
min_length_to_compare = 0;


% ##### SGP4 propagation settings #####
path_tle            = PARA.TLE_FILEPATH;
filename_tle        = PARA.TLE_FILENAME;
% t_start             = PARA.startmjd + mjd2jd;
% t_stop              = PARA.endmjd + mjd2jd;
% delta_t             = PARA.INIT_PROP_INTERVAL;
grav_const          = PARA.SGP4_GRAV_CONST;
write_file          = PARA.SGP4_WRITE_FILE;
path_out            = PARA.SGP4_PATH_OUT;
filename_out        = PARA.SGP4_FILENAME_OUT;
verification_mode   = PARA.SGP4_VERIFICATION_MODE;


% ##### Check input #####
numb_of_sats_to_prop = size(sat_names_to_prop, 1);
if (numb_of_sats_to_prop == 0) || verification_mode
    error_code = 1;
    error_msg = 'There is not any satellite defined to propagate the orbit!';
    return;
end


% ###############################################
% ############ Run SGP4 program code ############ 
% ###############################################


% #### Set SGP4 configuration ####

% opsmode:
% (a) calculation of GST with gstime.m, 
% (i) modern aproach of calculating GST
opsmode = 'a'; 

% whichconst:
% Set gravity constants.
% 72 = WGS72, 84 = WGS84, 721 = WGS72 low precission
whichconst = grav_const;



% ##### If an output file should be written #####
if write_file
    % If out folder doesn't exst, create it!
    if ~exist(path_out, 'dir')
        mkdir(path_out);
        fprintf(1,'Created output folder: %s\n', path_out);
    end
    % Open output file
    outfile = fopen([path_out, filename_out], 'wt');
    fprintf(1,'Created output file: %s\n\n', filename_out);
end



% ##### Save data to "sat_data" structure #####
sat_data.prop_setup.tle_filepath = path_tle;
sat_data.prop_setup.tle_filename = filename_tle;



% ##### SGP4 Propagation for defined satellites #####

if ~verification_mode
    
    % Loop over all satellites defined for propagation:
    for i_sat = 1 : numb_of_sats_to_prop

        sat_name_str = sat_names_to_prop(i_sat,:);

        % ##### Search for according TLE dataset in PARA structure #####
        flag_tle_dataset_found = 0;

        % Loop over all TLE datasets in PARA:
        for i_tle = 1 : length(PARA.tle_data.sat)
            header_str = PARA.tle_data.sat(i_tle).header_str;
            min_length_to_compare = min(length(sat_name_str), length(header_str));
            if (strncmp(sat_name_str, header_str, min_length_to_compare))
                flag_tle_dataset_found = 1;
                sat_data.sat(i_sat).TLE_header_line = header_str;    % TLE Line 0
                longstr1 = PARA.tle_data.sat(i_tle).line_1_str;
                longstr2 = PARA.tle_data.sat(i_tle).line_2_str;
                break;
            end
        end
% TEST TEST TEST  
% fprintf('tle_propagation: sat_name: %s, found: %d \n', header_str(1:24), flag_tle_dataset_found);
% TEST TEST TEST  
        
        % ##### If matching TLE dataset was found #####
        if flag_tle_dataset_found

            % Calculate JD and ydhms of start epoch
            [year,mon,day,hr,minute,sec] = invjday (t_start);

            % Save propagation setup data (start epoch) in struct
            sat_data.prop_setup.start_epoch.jd = t_start;
            sat_data.prop_setup.start_epoch.year = year;
            sat_data.prop_setup.start_epoch.mon = mon;
            sat_data.prop_setup.start_epoch.day = day;
            sat_data.prop_setup.start_epoch.h = hr;
            sat_data.prop_setup.start_epoch.min = minute;
            sat_data.prop_setup.start_epoch.sec = sec;

            % Calculate JD and ydhms of stop epoch
            [year,mon,day,hr,minute,sec] = invjday (t_stop);

            % Save propagation setup data (start epoch) in struct
            sat_data.prop_setup.stop_epoch.jd = t_stop;
            sat_data.prop_setup.stop_epoch.year = year;
            sat_data.prop_setup.stop_epoch.mon = mon;
            sat_data.prop_setup.stop_epoch.day = day;
            sat_data.prop_setup.stop_epoch.h = hr;
            sat_data.prop_setup.stop_epoch.min = minute;
            sat_data.prop_setup.stop_epoch.sec = sec;

            sat_data.prop_setup.delta_t_min = delta_t;

            % convert the char string to sgp4 elements
            % includes initialization of sgp4
            [satrec, startmfe, stopmfe, deltamin] = twoline2rv( whichconst, longstr1, longstr2, t_start, t_stop, delta_t, verification_mode);

            if write_file % If an output file should be written
                fprintf(outfile, '%d xx\n', satrec.satnum);
            end;

            % call the propagator to get the initial state vector value
            % H: Calculation of state vector for the given TLE epoch (tsince = 0.0):
            [satrec, ro ,vo] = sgp4 (satrec,  0.0);

            % Calculate JD and ydhms of given epoch
            jd = satrec.jdsatepoch; %
            [year,mon,day,hr,minute,sec] = invjday ( jd );
            tsince = 0;

            % Store data from given TLE epoch in struct
            sat_data.sat(i_sat).sat_number = satrec.satnum;
            sat_data.sat(i_sat).epoch.r = ro;
            sat_data.sat(i_sat).epoch.v = vo;
            sat_data.sat(i_sat).epoch.t_since = tsince;
            sat_data.sat(i_sat).epoch.year = year;
            sat_data.sat(i_sat).epoch.mon = mon;
            sat_data.sat(i_sat).epoch.day = day;
            sat_data.sat(i_sat).epoch.h = hr;
            sat_data.sat(i_sat).epoch.min = minute;
            sat_data.sat(i_sat).epoch.sec = sec;
            sat_data.sat(i_sat).epoch.jd = jd;

            if write_file % If an output file should be written
                fprintf(outfile, ' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f %5i%3i%3i %2i:%2i:%9.6f\n', tsince,ro(1),ro(2),ro(3),vo(1),vo(2),vo(3),year,mon,day,hr,minute,sec );
            end

            % Set time since epoch [min]
            tsince = startmfe;

            % check so the first value isn't written twice
            if ( abs(tsince) > 1.0e-8 )
                tsince = tsince - deltamin;
            end;

            i_epoch = 0;

            % #### loop to perform the propagation ####
            while ((tsince < stopmfe) && (satrec.error == 0))

                % counter for the number of current epoch
                i_epoch = i_epoch + 1;

                % Set time since epoch [min] 
                tsince = tsince + deltamin;
                if(tsince > stopmfe)
                    tsince = stopmfe;
                end;

                % Calculate propagation 
                [satrec, ro, vo] = sgp4 (satrec, tsince);

                % Error routine
                if (satrec.error > 0)
                   fprintf(1,'# *** error: t:= %f *** code = %3i\n', tsince, satrec.error);
                end 

                if (satrec.error == 0)

                    % Calculate absolute time epoch [JD and ydhms] of propagated data
                    jd = satrec.jdsatepoch + tsince/1440.0; % 1440 = 24*60
                    [year,mon,day,hr,minute,sec] = invjday ( jd );

                    % Store data from current epoch to struct
                    sat_data.sat(i_sat).prop(i_epoch).r = ro;
                    sat_data.sat(i_sat).prop(i_epoch).v = vo;
                    sat_data.sat(i_sat).prop(i_epoch).t_since = tsince;
                    sat_data.sat(i_sat).prop(i_epoch).year = year;
                    sat_data.sat(i_sat).prop(i_epoch).mon = mon;
                    sat_data.sat(i_sat).prop(i_epoch).day = day;
                    sat_data.sat(i_sat).prop(i_epoch).h = hr;
                    sat_data.sat(i_sat).prop(i_epoch).min = minute;
                    sat_data.sat(i_sat).prop(i_epoch).sec = sec;
                    sat_data.sat(i_sat).prop(i_epoch).jd = jd;

                    if write_file % If an output file should be written (or ver. mode)
                        fprintf(outfile, ' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f %5i%3i%3i %2i:%2i:%9.6f\n', tsince,ro(1),ro(2),ro(3),vo(1),vo(2),vo(3),year,mon,day,hr,minute,sec );
                    end

                end % if (satrec.error == 0)

            end; % while ((tsince < stopmfe) && (satrec.error == 0))


        % ##### If matching TLE dataset was NOT found ==>  Error! #####
        else
            error_code = 2;
            error_msg = ['TLE dataset was not found in the PARA structure for satellite: ', sat_name_str];
            return;
        end % if flag_tle_dataset_found

    end % for i_sat = 1 : numb_of_sats_to_prop
    
end % if ~verification_mode



% ##################################
% ##### Verification Mode Init #####
% ##################################
if verification_mode
    
    % Get VIE_SCHED version:
    path_str = path;
    search_str = '..\COMPILE\VIE_SCHED_V';
    path_str_pointer = strfind(path_str,search_str); 
    path_str_pointer = path_str_pointer(:,2);
    vievs_version = path_str(path_str_pointer + length(search_str): path_str_pointer + length(search_str) + 1);

    % propagation setup:
    filename_out = 'SGP4-VER.OUT';
    path_out = [pwd, '../COMPILE/VIE_SCHED_V', vievs_version ,'/SAT_SCHED/SGP4/VER_OUT/'];
    path_tle = [pwd, '../COMPILE/VIE_SCHED_V', vievs_version ,'/SAT_SCHED/SGP4/VER_TLE/'];
    filename_tle = 'SGP4-VER.TLE';
    
    % these are set in sgp4init
    global mu
    
    % Open TLE input file
    infile = fopen([path_tle, filename_tle], 'r');
    if (infile == -1)
        error_code = 1;
        error_msg = ['Failed to open TLE file for verification mode: ', filename_tle];
        return;
    end
    
    % ##### Go through TLE file #####
    while (~feof(infile))

        % number of TLE dataset currently processed
        i_sat = i_sat + 1;

        % Load 2 lines of TLE data:
        longstr1 = fgets(infile, 130);
        while ( (longstr1(1) == '#') && (feof(infile) == 0) )   % # at the beginning of a line = comment!
            longstr1 = fgets(infile, 130);
        end;
        if (feof(infile) == 0) % Not end of file
            longstr2 = fgets(infile, 130);
        end
        
        % Calculate JD and ydhms of start epoch
        [year,mon,day,hr,minute,sec] = invjday (t_start);

        % Save propagation setup data (start epoch) in struct
        sat_data.prop_setup.start_epoch.jd = t_start;
        sat_data.prop_setup.start_epoch.year = year;
        sat_data.prop_setup.start_epoch.mon = mon;
        sat_data.prop_setup.start_epoch.day = day;
        sat_data.prop_setup.start_epoch.h = hr;
        sat_data.prop_setup.start_epoch.min = minute;
        sat_data.prop_setup.start_epoch.sec = sec;

        % Calculate JD and ydhms of stop epoch
        [year,mon,day,hr,minute,sec] = invjday (t_stop);

        % Save propagation setup data (start epoch) in struct
        sat_data.prop_setup.stop_epoch.jd = t_stop;
        sat_data.prop_setup.stop_epoch.year = year;
        sat_data.prop_setup.stop_epoch.mon = mon;
        sat_data.prop_setup.stop_epoch.day = day;
        sat_data.prop_setup.stop_epoch.h = hr;
        sat_data.prop_setup.stop_epoch.min = minute;
        sat_data.prop_setup.stop_epoch.sec = sec;

        sat_data.prop_setup.delta_t_min = delta_t;

        % convert the char string to sgp4 elements
        % includes initialization of sgp4
        [satrec, startmfe, stopmfe, deltamin] = twoline2rv( whichconst, longstr1, longstr2, t_start, t_stop, delta_t, verification_mode);

        % Write to file:
        fprintf(outfile, '%d xx\n', satrec.satnum);

        % call the propagator to get the initial state vector value
        % ...Calculation of state vector for the given TLE epoch (tsince = 0.0):
        [satrec, ro ,vo] = sgp4 (satrec,  0.0);

        % Calculate JD and ydhms of given epoch
        jd = satrec.jdsatepoch; %
        [year,mon,day,hr,minute,sec] = invjday ( jd );
        tsince = 0;

        % Store data from given TLE epoch in struct
        sat_data.sat(i_sat).sat_number = satrec.satnum;
        sat_data.sat(i_sat).epoch.r = ro;
        sat_data.sat(i_sat).epoch.v = vo;
        sat_data.sat(i_sat).epoch.t_since = tsince;
        sat_data.sat(i_sat).epoch.year = year;
        sat_data.sat(i_sat).epoch.mon = mon;
        sat_data.sat(i_sat).epoch.day = day;
        sat_data.sat(i_sat).epoch.h = hr;
        sat_data.sat(i_sat).epoch.min = minute;
        sat_data.sat(i_sat).epoch.sec = sec;
        sat_data.sat(i_sat).epoch.jd = jd;

        % Write to file:
        fprintf(outfile, ' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n', tsince,ro(1),ro(2),ro(3),vo(1),vo(2),vo(3));

        % Set time since epoch [min]
        tsince = startmfe;

        % check so the first value isn't written twice
        if ( abs(tsince) > 1.0e-8 )
            tsince = tsince - deltamin;
        end

        i_epoch = 0;

        % loop to perform the propagation
        while ((tsince < stopmfe) && (satrec.error == 0))

            % counter for the number of current epoch
            i_epoch = i_epoch + 1;

            % Set time since epoch [min] 
            tsince = tsince + deltamin;
            if(tsince > stopmfe)
                tsince = stopmfe;
            end

            % Calculate propagation 
            [satrec, ro, vo] = sgp4 (satrec, tsince);

            % Error routine
            if (satrec.error > 0)
               fprintf(1,'# *** error: t:= %f *** code = %3i\n', tsince, satrec.error);
            end  

            if (satrec.error == 0)
                
                % Calculate absolute time epoch [JD and ydhms] of propagated data
                jd = satrec.jdsatepoch + tsince/1440.0; % 1440 = 24*60
                [year,mon,day,hr,minute,sec] = invjday ( jd );

                % Store data from current epoch to struct
                sat_data.sat(i_sat).prop(i_epoch).r = ro;
                sat_data.sat(i_sat).prop(i_epoch).v = vo;
                sat_data.sat(i_sat).prop(i_epoch).t_since = tsince;
                sat_data.sat(i_sat).prop(i_epoch).year = year;
                sat_data.sat(i_sat).prop(i_epoch).mon = mon;
                sat_data.sat(i_sat).prop(i_epoch).day = day;
                sat_data.sat(i_sat).prop(i_epoch).h = hr;
                sat_data.sat(i_sat).prop(i_epoch).min = minute;
                sat_data.sat(i_sat).prop(i_epoch).sec = sec;
                sat_data.sat(i_sat).prop(i_epoch).jd = jd;

                fprintf(outfile, ' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f', tsince,ro(1),ro(2),ro(3),vo(1),vo(2),vo(3));

                % Calculation classical orbital elements given the
                % geocentric equatorial position and velocity
                % vectors:
                [p,a,ecc,incl,node,argp,nu,m,arglat,truelon,lonper ] = rv2coe (ro,vo,mu);
                fprintf(outfile, ' %14.6f %8.6f %10.5f %10.5f %10.5f %10.5f %10.5f %5i%3i%3i %2i:%2i:%9.6f \n', a, ecc, incl*rad, node*rad, argp*rad, nu*rad, m*rad,year,mon,day,hr,minute,sec );

            end; % if (satrec.error == 0)
        end % while ((tsince < stopmfe) && (satrec.error == 0))
    end % while (~feof(infile))
   
    % Close files:
    fclose(infile); % Close TLE input file
    fclose(outfile); % Close output file
    
end % if verification_mode

return;

