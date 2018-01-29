% #########################################################################
% #     station_based_quasar_scheduling
% #########################################################################
%
% DESCRIPTION
%   This function is requried to use standard VIE_SCHED scheduling functions within the satellite scheduling environment.
%
%   NOTE: The following VIE_SCHED features are not available here:
%    - Tagalong mode
%    - Source bases strategy
%    - Twin telescope features
%
%   Sky coverage treatment:
%   - The sky coverage coefficeints (obs.staobs.ncat, obs.staobs.cattn, obs.staobs.catsn) are reseted, whenever this function is called, if "flag_reset_sky_coverage" is set!
%   - Sky coverage data determined in this function is written to output structures
%
% CREATED  
%   2015-07-06     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING
%   - vie_sched_s2_sat_sched
%   - vie_sched_src.m
%   - sdfaz.m
%   - checkSkyCoverage.m
%   - sched_analyser.m
%   - checkAndStat.m
%
%
% INPUT
%   - 
%
%   - stat_data                     : station data structure
%   - station_id_list               : List of IDs of the quasar observing stations (referring to the required entry of "stat_data.stat" and "obs_data.stat")
%   - source                        : Source structure
%   - obs_data                      : Observation data structure (preallocated in "find_observation_times.m")
%   - t_start_jd                    : Calculation epoch
%   - t_end_jd                      : Calculation epoch
%   - PARA                          : Global scheduling parameter strucutre
%
%
% OUTPUT
%   - error_code                    : Error Code (0 = no erros occured)
%   - error_msg                     : Error Message (empty, if no errors occured)
%
% CHANGES:
% - 2016-04-21, A. Hellerschmied: Bug fixes; Sub-netting added.
% - 2016-10-27, A. Hellerschmied: Call "vie_sched_src.m" with full functionality instead of "vie_sched_s2_sat_sched.m"
% - 2016-11-12, A. Hellerschmied: Bug-fix at using station indices
% - 2016-11-25, A. Hellerschmied: - PARA.fid_body set to "1" if the field does not exist (no file ID defined).
%                                 - Option to open sched_analyzer added
%                                 - sched structures are saved at ../DATA/LEVEL5/ and ../DATA/SCHED/<sub_dir>/LEVEL5/
%                                 - Added "checkSkyCoverage.m"
%                                 - Added "checkAndStat"
% - 2016-12-22, A. Hellerschmied: Added options to write CW text to log file (PARA.fid_body)
%

function [sched_data, obs_data, error_code, error_msg] = station_based_quasar_scheduling(source, station, PARA, sched_data, obs_data, station_id_list, t_start_jd, t_end_jd)

    % ##### Init #####
    error_code = 0;
    error_msg = '';
       
    number_of_quasars  = length(source);
    
    PARA_temp = PARA; % copy PARA and apply modifications to the copied version (keep orig. struct)
    
    if PARA.SAVEOUTPUT
        if sched_data.number_of_vie_sched_calls == 0
            fprintf(PARA.fid_body, '########## Quasar block %d ##########\n', sched_data.number_of_vie_sched_calls + 1);
        else
            fprintf(PARA.fid_body, '\n\n########## Quasar block %d ##########\n', sched_data.number_of_vie_sched_calls + 1);
        end
    end

    
    % ##### Options #####
    flag_reset_sky_coverage              = false; % true: Sky coverage data from "obs_data" is considered when starting station based scheduling with "vie_sched_src.m"
%     flag_write_sky_cov_to_out_struct    = true; % true: Sky coverage data determined in this function is written zu the output structure ("sched_data" and "obs_data")
    
    % ##### Prepare station structure #####
    % - Bring stations into the same order as in "sched_data" and "obs_data"
    % - Only consider the stations defined in "station_id_list"
    i_stat = 1;
    for i_stat_1 = 1 : length(sched_data.stat(station_id_list))
        for i_stat_2 = 1 : length(station)
            if strcmpi(sched_data.stat(station_id_list(i_stat_1)).name, station(i_stat_2).name)
                station_temp(i_stat) = station(i_stat_2);
                i_stat = i_stat + 1;
                break;
            end
        end
    end
    station = station_temp;
    
    
    % ##### Special settings used for quasar observation blocks #####
    flag_use_special_q_param = false;
    
    % ##### Prepare twin structure #####
    % Twin antenna functions are disable here
    twin.num    = 0;
    twin.sn     = [0, 0];

    
    % ##### Adjust PARA structure #####
    % Set start and end time:
	PARA_temp.startmjd  = t_start_jd - 2.400000500000000e+006;
	PARA_temp.endmjd    = t_end_jd - 2.400000500000000e+006;
    % Choose station based scheduling approach:
    PARA_temp.OPTIMIZATION = 2;
    
    % Debug output:
    if ~isfield(PARA_temp, 'DEBUG_FLAG')
        PARA_temp.DEBUG_FLAG = false;
    end
    if ~isfield(PARA_temp, 'fid_body')
        PARA_temp.fid_body = 1;
    end
    if ~isfield(sched_data, 'number_of_vie_sched_calls')
        sched_data.number_of_vie_sched_calls = 0;
    end
    if ~isfield(PARA_temp, 'openAnalyser')
        PARA_temp.openAnalyser = 0;
    end
    

    % Apply spezial parameterization:
    if flag_use_special_q_param
        fprintf('\n');
        fprintf(' ============================================================================\n');
        fprintf(' ==> NOTE: special parameterisation in PARA used for quasar scheduling!   <==\n');
        fprintf(' ==> Defined in the code of function "station_based_quasar_scheduling.m"  <==\n');
        fprintf(' ==> ....by setting "flag_use_special_q_param = true"                     <==\n');
        fprintf(' ============================================================================\n');
        PARA_temp.CALIBRATION    = 10;
        PARA_temp.SOURCE         = 5;
        PARA_temp.TAPETM         = 1;
        PARA_temp.CORSYNCH       = 3;
        disp('PARA_temp.CALIBRATION    = 10;');
        disp('PARA_temp.SOURCE         = 5');
        disp('PARA_temp.TAPETM         = 1');
        disp('PARA_temp.CORSYNCH       = 3;');
        fprintf(' ============================================================================\n');
        fprintf('\n');
    end


    
    % ##### Prepare "obs" struct #####
    % ### Initialize "obs.staobs" ###
    number_of_stat              = length(station);
    obs_data_stat_temp  = obs_data.stat(station_id_list);
    if obs_data.number_of_scans == 0 % If there wasn't any scan before...
        for i_stat = 1 : number_of_stat
            obs.staobs(i_stat).ifirst   = 1; % Has to be = 1, if there weren't any scans before...
            obs.staobs(i_stat).endmjd   = PARA_temp.startmjd;
            obs.staobs(i_stat).az       = obs_data_stat_temp(i_stat).end_of_last_obs.un_az; % [rad]
            obs.staobs(i_stat).el       = obs_data_stat_temp(i_stat).end_of_last_obs.el;    % [rad]
            obs.staobs(i_stat).ha       = obs_data_stat_temp(i_stat).end_of_last_obs.ha;    % [rad]
            obs.staobs(i_stat).dc       = obs_data_stat_temp(i_stat).end_of_last_obs.dc;    % [rad]
            % Sky coverage data:
            % - is resetted
            obs.staobs(i_stat).ncat     = 0;                      % Number of observed fields within PARA.SKYDT; = length of fields "cattn" and "catsn"
            obs.staobs(i_stat).cattn(1) = PARA_temp.startmjd;     % Vector with epochs of last observations in the sky cov. fields indicated in "catsn"
            obs.staobs(i_stat).catsn(1) = 0;                      % Vector with sky cov. field numbers
        end
    else % If there was at least one scan before => initialize with correct pointing angles!
        for i_stat = 1 : number_of_stat
            obs.staobs(i_stat).ifirst   = 0; % Has to be = 1, if there weren't any scans before...
            obs.staobs(i_stat).endmjd   = PARA_temp.startmjd;
            obs.staobs(i_stat).az       = obs_data_stat_temp(i_stat).end_of_last_obs.un_az; % [rad]
            obs.staobs(i_stat).el       = obs_data_stat_temp(i_stat).end_of_last_obs.el;    % [rad]
            obs.staobs(i_stat).ha       = obs_data_stat_temp(i_stat).end_of_last_obs.ha;    % [rad]
            obs.staobs(i_stat).dc       = obs_data_stat_temp(i_stat).end_of_last_obs.dc;    % [rad]
            % Sky coverage data:
            % Init. sky coverage:
            obs.staobs(i_stat).ncat     = 0;                      % Number of observed fields within PARA.SKYDT; = length of fields "cattn" and "catsn"
            obs.staobs(i_stat).cattn(1) = PARA_temp.startmjd;     % Vector with epochs of last observations in the sky cov. fields indicated in "catsn"
            obs.staobs(i_stat).catsn(1) = 0;                      % Vector with sky cov. field numbers
            if ~flag_reset_sky_coverage
                % => Here, the sky coverage data is initialized properly, based on the data in "obs_data.stat(i_stat).sky_cov"
                % Loop over all entries of sky coverage vector of the current station:
                for i_sc = 1 : length(obs_data_stat_temp(i_stat).sky_cov)
                    % Check, if there is an entry in the field with the index "i_sc":
                    if ~isempty(obs_data_stat_temp(i_stat).sky_cov(i_sc))
                        % Check, if the entry is still valid (ref. epoch = end of last scan)
                        if (((obs.staobs(i_stat).endmjd - (obs_data_stat_temp(i_stat).sky_cov(i_sc) - 2400000.5) ) * 1440) <= PARA_temp.SKYDT)
                            % Sky cov. is reseted:
                            obs.staobs(i_stat).ncat                             = obs.staobs(i_stat).ncat + 1;                          % Number of observed fields within PARA.SKYDT; = length of fields "cattn" and "catsn"
                            obs.staobs(i_stat).cattn(obs.staobs(i_stat).ncat)   = obs_data_stat_temp(i_stat).sky_cov(i_sc) - 2400000.5; % Vector with epochs of last observations in the sky cov. fields indicated in "catsn"
                            obs.staobs(i_stat).catsn(obs.staobs(i_stat).ncat)   = i_sc;                                                 % Vector with sky cov. field numbers
                        end
                    end
                end
            end % if ~flag_reset_sky_coverage
        end
    end
    
    % ### Initialize "obs.srcobs" ###
    % - initialize "obs.srcobs" with the values from "obs_data.quasars", if "obs_data.quasars" was already initialized:
    if (length(obs_data.quasars) == 1) && isempty(obs_data.quasars(1).number_of_scans) % Was NOT initialized
        for i_q = 1 : number_of_quasars
            obs.srcobs(i_q).obsmjd                  = 0.0;
            obs.srcobs(i_q).nscan                   = 0.0;
            obs_data.quasars(i_q).last_obs_jd       = 0;
            obs_data.quasars(i_q).number_of_scans   = 0;
        end
    else % Was already initialized
        for i_q = 1 : number_of_quasars
            obs.srcobs(i_q).obsmjd = obs_data.quasars(i_q).last_obs_jd - 2.400000500000000e+006;
            obs.srcobs(i_q).nscan  = obs_data.quasars(i_q).number_of_scans;
        end
    end

    
    
    % ### source based obs parameters ###
    obs.srcat       = 0;
    obs.srcn        = 0;
    obs.icatpair    = 0;
    catpair         = [];
    srcat           = 0;
    
    
    
    % ##### Load ephem. data #####
    jpl = load_jpl_Earth('jpl_421');
    
    
    % ##### Start station based scheduling #####
%      [sched] = vie_sched_s2_sat_sched(station, twin, source, PARA_temp.obsmode, PARA_temp, obs_data, station_id_list);
    try
        [sched_new, obs_new] = vie_sched_src(station, twin, source, PARA_temp.obsmode, srcat, catpair, PARA_temp, obs, jpl);
    catch
        fprintf(' ===>> Error in "vie_sched_src.m"!\n')
        keyboard 
    end
    % Assign output structures, if there was not error in "vie_sched_src":
    obs = obs_new;
    sched = sched_new;
    
    
    % 
    
    % #####################################################
    % ##### Save scan: Update obs_data and sched_data #####
    % #####################################################
    
    
    % loop init.:
    scan_id = sched_data.number_of_scans; % Continuous scan number within this function to adress a scan in "sched_data"
    t_end_last_scan_jd = 0;
    
    
    % #### Loop over all scans in "sched" structure ####
    for i_scan = 1 : length(sched)
        
        % #### Loop over all sub-scans in "sched.scan" structure ####
        for i_sub_scan = 1 : sched(i_scan).nscan
        
            % Scan ID in "sched_data" structure
            scan_id = scan_id + 1;
            sched_data.number_of_scans = scan_id;
            
            % Number of scans:
            obs_data.number_of_scans = obs_data.number_of_scans + 1;

            % ##### Get scan start [JD] #####
            t_start_jd = sched(i_scan).scan(i_sub_scan).startmjd + 2.400000500000000e+006;

            % ##### Get scan end [JD] #####
            t_end_mjd_list = sortrows([sched(i_scan).scan(i_sub_scan).sta.endmjd]);
            t_scan_end_jd   = t_end_mjd_list(end) + 2.400000500000000e+006;
                       
            % Nominal session starte/end time:
            if isempty(sched_data.t_nominal_start_jd)
                sched_data.t_nominal_start_jd = t_start_jd;
            elseif sched_data.t_nominal_start_jd > t_start_jd
                sched_data.t_nominal_start_jd = t_start_jd;
            end
            if isempty(sched_data.t_nominal_end_jd)
                sched_data.t_nominal_end_jd = t_scan_end_jd;
            elseif sched_data.t_nominal_end_jd < t_scan_end_jd
                sched_data.t_nominal_end_jd = t_scan_end_jd;
            end

            sched_data.scan(scan_id).number_of_stat     = sched(i_scan).scan(i_sub_scan).nsta;
            sched_data.scan(scan_id).repos_int_sec      = PARA_temp.REPOS_INT;
            sched_data.scan(scan_id).flag_record_scan   = 1;                       % !!!! Hard coded here! If it is required, change that later on.... 
            sched_data.scan(scan_id).t_start_jd         = t_start_jd;
            sched_data.scan(scan_id).t_end_jd           = t_scan_end_jd;

            % Get quasar ID:
            quasar_id = sched(i_scan).scan(i_sub_scan).srcid;
            
            if t_scan_end_jd > obs_data.quasars(quasar_id).last_obs_jd
                obs_data.quasars(quasar_id).last_obs_jd = t_scan_end_jd;
            end

            % Preallocate "epoch" data:
            sched_data.scan(scan_id).number_of_epochs = 1;

            % Save some data:
            sched_data.scan(scan_id).sat_id         = [];
            sched_data.scan(scan_id).sat_name       = [];
            sched_data.scan(scan_id).sat_number     = [];
            sched_data.scan(scan_id).quasar_id      = quasar_id;
            sched_data.scan(scan_id).quasar_name    = source(quasar_id).name;
            sched_data.scan(scan_id).obs_type       = 'quasar';

            obs_data.quasars(quasar_id).number_of_scans = obs_data.quasars(quasar_id).number_of_scans + 1;

            % Topocentric pointing angles at scan start and end:
            for i_stat = 1 : sched(i_scan).scan(i_sub_scan).nsta

                station_id  = station_id_list(sched(i_scan).scan(i_sub_scan).sta(i_stat).staid);    % ref. to "obs_data" and is stored in "sched_data.stat.stat_id"
                staid       = sched(i_scan).scan(i_sub_scan).sta(i_stat).staid;                     % re. to "station" 
                % i_stat                                                                            % ref. to "sched_data.scan(scan_id).stat(i_stat)", "sched(i_scan).scan(i_sub_scan).sta(i_stat)", "obs.staobs(i_stat)"

                % Preallocate "epoch" data for stepwise tracking
                sched_data.scan(scan_id).stat(i_stat).epoch(1).jd   = t_start_jd;
                sched_data.scan(scan_id).stat(i_stat).stat_id       = station_id;
                
                sched_data.scan(scan_id).stat(i_stat).duration_sec  = sched(i_scan).scan(i_sub_scan).sta(i_stat).duration;
                sched_data.scan(scan_id).stat(i_stat).slew_time_sec = sched(i_scan).scan(i_sub_scan).sta(i_stat).slewtime;
                
                

                
                % #### Scan end time: ####

                % ## topocentric poitng data at scan end time (individual for each station!): ##
%                 [az, el, ha, dc] = zazel_s(sched(i_scan).scan(i_sub_scan).sta(i_stat).endmjd, station(station_id).llh(1), station(station_id).llh(2), source(quasar_id).ra, source(quasar_id).de);
                [az, el, ha, dc] = zazel_s(sched(i_scan).scan(i_sub_scan).sta(i_stat).endmjd, station(staid).llh(1), station(staid).llh(2), source(quasar_id).ra, source(quasar_id).de);
%                 [az, el, ha, dc] = zazel_s(sched(i_scan).scan(i_sub_scan).startmjd, station(station_id).llh(1), station(station_id).llh(2), source(quasar_id).ra, source(quasar_id).de);

                % NOTE: un_az is not available here for scan start time!
                % => empty field!
                sched_data.scan(scan_id).stat(i_stat).end.az        = az;
                sched_data.scan(scan_id).stat(i_stat).end.el        = el;
                sched_data.scan(scan_id).stat(i_stat).end.ha        = ha;
                sched_data.scan(scan_id).stat(i_stat).end.dc        = dc;
                sched_data.scan(scan_id).stat(i_stat).end.un_az     = []; % Not available! If this value is required, it has to be added here!
                sched_data.scan(scan_id).stat(i_stat).end.jd        = sched(i_scan).scan(i_sub_scan).sta(i_stat).endmjd + 2.400000500000000e+006;
                
                 if PARA_temp.DEBUG_FLAG   
                    % Test ++++
                    fprintf('sched.sta.endmjd  .: %2.15f   %2.15f   %2.15f   %2.15f\n', az, el, ha, dc);
                    fprintf('sched.scan.stat   .: %2.15f   %2.15f   %2.15f   %2.15f\n', sched(i_scan).scan(i_sub_scan).sta(i_stat).az, sched(i_scan).scan(i_sub_scan).sta(i_stat).el, sched(i_scan).scan(i_sub_scan).sta(i_stat).ha, sched(i_scan).scan(i_sub_scan).sta(i_stat).dc);
                    % [az, el, ha, dc] = zazel_s(sched(i_scan).scan(i_sub_scan).startmjd, station(station_id).llh(1), station(station_id).llh(2), source(quasar_id).ra, source(quasar_id).de); 
%                     [az, el, ha, dc] = zazel_s(sched(i_scan).scan(i_sub_scan).sta(i_stat).startmjd, station(station_id).llh(1), station(station_id).llh(2), source(quasar_id).ra, source(quasar_id).de); 
                    [az, el, ha, dc] = zazel_s(sched(i_scan).scan(i_sub_scan).sta(i_stat).startmjd, station(staid).llh(1), station(staid).llh(2), source(quasar_id).ra, source(quasar_id).de); 
                    fprintf('sched.sta.startmjd.: %2.15f   %2.15f   %2.15f   %2.15f\n', az, el, ha, dc);
                    % Test ----
                 end

% Calc. pointing angels at scan start time (t_start_jd):
% [az, el, ha, dc] = zazel_s(sched(i_scan).scan(i_sub_scan).startmjd, station(station_id).llh(1), station(station_id).llh(2), source(quasar_id).ra, source(quasar_id).de); 
[az, el, ha, dc] = zazel_s(sched(i_scan).scan(i_sub_scan).startmjd, station(staid).llh(1), station(staid).llh(2), source(quasar_id).ra, source(quasar_id).de); 
if PARA_temp.DEBUG_FLAG  
    fprintf('sched.startmjd    .: %2.15f   %2.15f   %2.15f   %2.15f\n', az, el, ha, dc);
end
% Calc. un_az at scan start time (t_start_jd):
[un_az_new_rad] = sdfaz(sched_data.stat(station_id).lim11, sched_data.stat(station_id).lim12, sched(i_scan).scan(i_sub_scan).sta(i_stat).az, az);
% [un_az_new_rad] = sdfaz(sched_data.stat(i_stat).lim11, sched_data.stat(i_stat).lim12, sched(i_scan).scan(i_sub_scan).sta(i_stat).az, az);

                if PARA_temp.DEBUG_FLAG   
                    % Test ++++
                    fprintf('factor pi         .: %2.15f\n', abs(un_az_new_rad - az)/pi);
                    fprintf('un_az_new_rad     .: %2.15f\n', un_az_new_rad);
                    fprintf('\n');
                    % Test ----
                end

                

                % #### Scan start time: ####
%                 sched_data.scan(scan_id).stat(i_stat).start.az    = sched(i_scan).scan(i_sub_scan).sta(i_stat).az; % = unaz
%                 sched_data.scan(scan_id).stat(i_stat).start.el    = sched(i_scan).scan(i_sub_scan).sta(i_stat).el;
%                 sched_data.scan(scan_id).stat(i_stat).start.ha    = sched(i_scan).scan(i_sub_scan).sta(i_stat).ha;
%                 sched_data.scan(scan_id).stat(i_stat).start.dc    = sched(i_scan).scan(i_sub_scan).sta(i_stat).dc;
%                 sched_data.scan(scan_id).stat(i_stat).start.un_az = sched(i_scan).scan(i_sub_scan).sta(i_stat).az; % = unaz
%                 sched_data.scan(scan_id).stat(i_stat).start.jd    = t_start_jd;


% Mit neuem Update von Matthias => Verwende auskommentierten Block von oben!!!!
% => Dann: "sched.scan.sta.az/el/ha/dc" auf die epoche "sched.scan.startmjd bezogen"
% => Cable wrap lösen mit sdfaz.m unnötig für scan start time!!

                sched_data.scan(scan_id).stat(i_stat).start.az    = az;
                sched_data.scan(scan_id).stat(i_stat).start.el    = el;
                sched_data.scan(scan_id).stat(i_stat).start.ha    = ha;
                sched_data.scan(scan_id).stat(i_stat).start.dc    = dc;
                sched_data.scan(scan_id).stat(i_stat).start.un_az = un_az_new_rad; % = unaz
                sched_data.scan(scan_id).stat(i_stat).start.jd    = t_start_jd;
                
                % Calc. sky coverage number at scan start for this station and source:
                [sky_cov_num] = zcoverage(az, el);
                sched_data.scan(scan_id).stat(i_stat).sky_cov_num  = sky_cov_num;
                
                obs_data.stat(station_id).sky_cov(sky_cov_num)          = t_start_jd; % Update sky-coverage
                obs_data.stat(station_id).end_of_last_obs.quasar_id     = quasar_id;

                if t_scan_end_jd > t_end_last_scan_jd
                    t_end_last_scan_jd = t_scan_end_jd;
                end
                
            end % for i_stat = 1 : length(station_id_list)
        end % for i_sub_scan = 1 : sched(10).nscan
    end % for i_scan = 1 : length(sched)
    

    
    % ##### Save scan data of last scan to "obs_data" and reset temp. data #####
     
    obs_data.end_of_last_scan_jd = t_end_last_scan_jd;
    
    % Loop over all stations in station_id_list
    for i_stat = 1 : length(station_id_list)
    
        station_id = station_id_list(i_stat);
        % station_id ... index for "obs_data"
        % i_stat     ... index for staobs, station, sched, etc...

        % #### Save antenna data to "obs_data" structure ####
        obs_data.stat(station_id).end_of_last_obs.jd            = obs.staobs(i_stat).endmjd + 2.400000500000000e+006;    % [JD]
        obs_data.stat(station_id).end_of_last_obs.az            = obs.staobs(i_stat).az;            % [rad]
        obs_data.stat(station_id).end_of_last_obs.el            = obs.staobs(i_stat).el;            % [rad]
        obs_data.stat(station_id).end_of_last_obs.ha            = obs.staobs(i_stat).ha;            % [rad]
        obs_data.stat(station_id).end_of_last_obs.dc            = obs.staobs(i_stat).dc;            % [rad]
        obs_data.stat(station_id).end_of_last_obs.un_az         = obs.staobs(i_stat).az;            % [rad]
        obs_data.stat(station_id).end_of_last_obs.sat_id        = [];

        % ##### Reset temp. "obs_data" #####
        obs_data.stat(station_id).begin_of_new_obs_temp.un_az   = [];
        obs_data.stat(station_id).begin_of_new_obs_temp.el      = [];
        obs_data.stat(station_id).begin_of_new_obs_temp.jd      = [];
        obs_data.stat(station_id).end_of_new_obs_temp.un_az     = [];
        obs_data.stat(station_id).end_of_new_obs_temp.el        = [];
        obs_data.stat(station_id).end_of_new_obs_temp.jd        = [];
        
    end % for i_stat = 1 : length(station_id_list)
    
    
    % saves structures in /SCHED/<sub_dir>/LEVEL5/
    % => Can be used later on e.g. as input for the sched analyzer!
    poutfile_sub = ['../DATA/SCHED/', PARA.OUTPUT_SUBDIR, '/'];
    if ~isdir([poutfile_sub, 'LEVEL5'])
        mkdir([poutfile_sub, 'LEVEL5']);
    end
    
    sched_data.number_of_vie_sched_calls = sched_data.number_of_vie_sched_calls + 1;
    
    obsmode = PARA_temp.obsmode;
    
    level5_subdir_str = [poutfile_sub, 'LEVEL5/', sprintf('q_block_%0d', sched_data.number_of_vie_sched_calls), '/'];
    if ~isdir(level5_subdir_str)
        mkdir(level5_subdir_str);
    end
    
    save([level5_subdir_str ,'schedparam.mat'], 'PARA_temp'); 
    save([level5_subdir_str, 'station.mat'], 'station'); 
    save([level5_subdir_str, 'sched.mat'], 'sched'); 
    save([level5_subdir_str, 'source.mat'], 'source'); 
    save([level5_subdir_str, 'obsmode.mat'], 'obsmode'); 
    
    
    save([PARA.pfolder 'schedparam.mat'], 'PARA');
    save([PARA.pfolder 'station.mat'], 'station');
    save([PARA.pfolder 'sched.mat'], 'sched');
    save([PARA.pfolder 'source.mat'], 'source');
    save([PARA.pfolder 'obsmode.mat'], 'obsmode');
    
    [meanSky,h] = checkSkyCoverage( sched,PARA_temp,station );
    checkAndStat( sched,station,source,PARA,meanSky );
    h.Name = sprintf('Sky Coverage of quasar block %d',sched_data.number_of_vie_sched_calls);
    print(h ,sprintf('%sSky Coverage of quasar block %d',level5_subdir_str , sched_data.number_of_vie_sched_calls),'-dpng','-r150') % PNG
    savefig(h, [sprintf('%sSky Coverage of quasar block %d',level5_subdir_str , sched_data.number_of_vie_sched_calls), '.fig']);
    
    % Open sched analyser:
    if PARA_temp.openAnalyser
        sched_analyser();
    end
    
    
return;