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
%   - 2025-01-13, H. Wolf: added end to function block

function [sat_data, error_code, error_msg] = tle_propagation(t_start, t_stop, delta_t, PARA)

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

% ##### init #####
error_code = 0;                 % 0 = No error
error_msg = '';                 % '' = No error

t_start = t_start - (1/24);
t_stop = t_stop + (1/24);

% ##### SGP4 propagation settings #####
path_tle            = PARA.TLE_FILEPATH;
filename_tle        = PARA.TLE_FILENAME;
grav_const          = PARA.SGP4_GRAV_CONST;
%write_file          = PARA.SGP4_WRITE_FILE;
%path_out            = PARA.SGP4_PATH_OUT;
%filename_out        = PARA.SGP4_FILENAME_OUT;
verification_mode   = PARA.SGP4_VERIFICATION_MODE;

numb_of_sats_to_prop = size(PARA.tle_data.sat, 2);

% ###############################################
% ############ Run SGP4 program code ############ 
% ###############################################

% #### Set SGP4 configuration ####

% opsmode:
% (a) calculation of GST with gstime.m, 
% (i) modern aproach of calculating GST
% opsmode = 'a'; 

% whichconst:
% Set gravity constants.
% 72 = WGS72, 84 = WGS84, 721 = WGS72 low precission
whichconst = grav_const;

% ##### Save data to "sat_data" structure #####
sat_data.prop_setup.tle_filepath = path_tle;
sat_data.prop_setup.tle_filename = filename_tle;


% ##### SGP4 Propagation for defined satellites #####

   % Loop over all satellites defined for propagation:
     for i_sat = 1 : numb_of_sats_to_prop

        header_str = PARA.tle_data.sat(i_sat).header_str;
        sat_data.sat(i_sat).TLE_header_line = header_str;  
        longstr1 = PARA.tle_data.sat(i_sat).line_1_str;
        longstr2 = PARA.tle_data.sat(i_sat).line_2_str;
      
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

        % Set time since epoch [min]
        tsince = startmfe;

        % check so the first value isn't written twice
        if ( abs(tsince) > 1.0e-8 )
            tsince = tsince - deltamin;
        end

        i_epoch = 0;

        % #### loop to perform the propagation ####
        while ((tsince < stopmfe) && (satrec.error == 0))
            i_epoch = i_epoch + 1;

            % Set time since epoch [min] 
            tsince = tsince + deltamin;
            if(tsince > stopmfe)
                tsince = stopmfe;
            end

            [satrec, ro, vo] = sgp4 (satrec, tsince);

            if (satrec.error == 0)
                % Calculate absolute time epoch [JD and ydhms] of propagated data
                jd = satrec.jdsatepoch + tsince/1440.0; % 1440 = 24*60
                [year,mon,day,hr,minute,sec] = invjday ( jd );

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
            end 
        end 
     end 
end