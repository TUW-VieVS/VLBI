% #########################################################################
% #     sort_scan_cons
% #########################################################################
%
% DESCRIPTION
%   This function sorts a number scan configurations (contained in structure "scan_cons"). 
%   Satellite and quasar scans are supported.
%
%   The following criteria are considered:
%   - Quasar scans:
%       - 
%       - 
%       - 
%       - 
%
%   - Satellite scans:
%       - 
%       - 
%       - 
%
%
% CREATED  
%   2015-11-02     Andreas Hellerschmied
%
% REFERENCES
%
%
% COUPLING (external function calls)
%   - zcoverage
%
%
% SUB-ROUTINES defined within this file
%   - 
%
%
% INPUT
%   - scan_cons                     : Structure containing all possible scan-configurations for the defined time and observation constellation
%   - stat_data                     : station data structure
%   - source                        : Source structure
%   - PARA                          : Global scheduling parameter strucutre
%   - t_epoch_in_jd                 : Scan start time
%   - scan_type                     : Scan type (sat. or quasar)
%   - obs_data                      : Observation data structure (preallocated in "find_observation_times.m")
%   - stat_id_list                  : IDs of all stations of the current station network (statellites/quasars)
%
%
% OUTPUT
%   - sort_index                    : Sort index for the structure "scan_cons"
%   - scan_cons                     : Structure containing all possible scan-configurations for the defined time and observation constellation (weight factors added)
%   - obs_data                      : Observation data structure (preallocated in "find_observation_times.m"), sky-coverage number updated
%   - error_code                    : Error Code (0 = no erros occured)
%   - error_msg                     : Error Message (empty, if no errors occured)
%
% CHANGES:
% - 2016-05-09: A. Hellerschmied: Separate sky coverage treatment for satellites added.  
%

function [sort_index, scan_cons, obs_data, error_code, error_msg] = sort_scan_cons(scan_cons, stat_data, source, PARA, t_epoch_in_jd, scan_type, obs_data, stat_id_list)

    % ##### Init #####
    error_code = 0;
    error_msg = '';
    
    number_of_scan_cons = length(scan_cons);
    number_of_all_stations = length(stat_id_list);
    
    % ##### Constants #####
    mjd2jd = 2.400000500000000e+006;

    % ##### Preallocation #####
    scon_weights            = zeros(number_of_scan_cons, 1);
    scon_weights_obs_num    = zeros(number_of_scan_cons, 1);
    scon_weights_sky_cov    = zeros(number_of_scan_cons, 1);
    scon_weights_scan_end   = zeros(number_of_scan_cons, 1);
    sky_coverage            = zeros(number_of_scan_cons, 1);
    number_of_obs           = zeros(number_of_scan_cons, 1);
    t_scon_end_jd           = zeros(number_of_scan_cons, 1);

    
    % ############################################################################
    % ##### 1.) Calc all required factors for the weighting of the scan cons #####
    % ############################################################################
    
    % ##### Loop over all possible scan-cons #####
    for i_scon = 1 : number_of_scan_cons
    
        
        % #### Loop over all scans (1 or 2) within one scan-con ####
        for i_scan = 1 : scan_cons(i_scon).number_of_scans
            
            % ##### Define the reference epoch for the cal. of the sky coverage #####
            if isempty(t_epoch_in_jd)
                t_epoch_jd = scan_cons(i_scon).scan(i_scan).scan_start_jd; % Start of the curren scan (used for the final sorting!)
            else
                t_epoch_jd = t_epoch_in_jd;
            end
            
            
            number_of_stat  = length(scan_cons(i_scon).scan(i_scan).list_stat_id_obs);
            
            
            % ##### Number of observations #####
            number_of_obs(i_scon) = number_of_obs(i_scon) + number_of_stat * (number_of_stat - 1)/2;
            
            
            % ##### Loop over all stations within this scan #####
            for i_stat = 1 : number_of_stat

                stat_id = scan_cons(i_scon).scan(i_scan).list_stat_id_obs(i_stat);

                % #### Get the latest "observation end time" of all stations within all scans of this scan-con ####
                if scan_cons(i_scon).scan(i_scan).stat(i_stat).end_of_new_obs_temp.jd > t_scon_end_jd(i_scon)
                    t_scon_end_jd(i_scon) = scan_cons(i_scon).scan(i_scan).stat(i_stat).end_of_new_obs_temp.jd;
                end

                
                % #### Calc. sky coverage ####
                switch(scan_type)

                    % #############################
                    % ####     Quasar scan     ####
                    % #############################
                    case 'q'

                        % AzEl:
                        source_quasar = source(scan_cons(i_scon).scan(i_scan).quasar_id);
                        [az_rad, el_rad, ha, dc] = zazel_s(t_epoch_jd - mjd2jd, stat_data.stat(i_stat).location.ellipsoid.long * pi/180, stat_data.stat(i_stat).location.ellipsoid.lat * pi/180, source_quasar.ra, source_quasar.de);

                        % Calc. sky coverage number for this station and ssource:
                        [cat] = zcoverage(az_rad, el_rad);
                        % Check, if there was an observation with the sky coverage number "cat" within the time "PARA.SKYDT" => No => increase sky coverage counter:
                        if obs_data.stat(stat_id).sky_cov(cat) < (t_epoch_jd - PARA.SKYDT/1440)
                            sky_coverage(i_scon) = sky_coverage(i_scon) + 1;
                        end
                        
                    % ################################
                    % ####     Satellite scan     ####
                    % ################################
                    case 's'
                        
                        % AzEl:
                        [az_deg, el_deg, range, az_rate, el_rate, r_rate, tra, tdec, tra_rate, tdec_rate, local_hour_angle, local_hour_angle_rate, error_code, error_msg] = tle_2_topo(stat_data, stat_id, PARA, scan_cons(i_scon).scan(i_scan).sat_id, t_epoch_jd);
                        if error_code > 0
                            error_msg = ['tle_2_topo: ', error_msg];
                            return; 
                        end
                        az_rad = az_deg * pi/180;
                        el_rad = el_deg * pi/180;
                        
                        % Calc. sky coverage number for this station and ssource:
                        [cat] = zcoverage(az_rad, el_rad);
                        % Check, if there was an observation with the sky coverage number "cat" within the time "PARA.SKYDT(_SAT)" => No => increase sky coverage counter:
                        if PARA.CALC_SKYCOV_SAT
                            if obs_data.stat(stat_id).sky_cov_sat(cat) < (t_epoch_jd - PARA.SKYDT_SAT/1440)
                                sky_coverage(i_scon) = sky_coverage(i_scon) + 1;
                            end
                        else
                            if obs_data.stat(stat_id).sky_cov(cat) < (t_epoch_jd - PARA.SKYDT/1440)
                                sky_coverage(i_scon) = sky_coverage(i_scon) + 1;
                            end
                        end

                end % switch(scan_type)

%                 % Calc. sky coverage number for this station and ssource:
%                 [cat] = zcoverage(az_rad, el_rad);
%                 % Check, if there was an observation with the sky coverage number "cat" within the time "PARA.SKYDT" => No => increase sky coverage counter:
%                 if obs_data.stat(stat_id).sky_cov(cat) < (t_epoch_jd - PARA.SKYDT/1440)
%                     sky_coverage(i_scon) = sky_coverage(i_scon) + 1;
%                 end
                scan_cons(i_scon).scan(i_scan).stat(i_stat).sky_cov_num = cat;
% TEST TEST TEST      
% switch(scan_type)
%    case 'q'
%        src_name = source_quasar.name;
%     case 's'
%         src_name = obs_data.sat(scan_cons(i_scon).scan(i_scan).sat_id).name(1:24);
% end
% fprintf('- src: %10s, %8s, cat: %2d, scon_id: %3d, sky_cov_w: %1d\n', src_name, stat_data.stat(stat_id).name, cat, i_scon, sky_coverage(i_scon));
% TEST TEST TEST      

            end % for i_stat = 1 : number_of_stat

        end % for i_scan = 1 : scan_cons(i_scon).number_of_scans 
 
    end % for i_scon = 1 : length(scan_cons)
    
    
    
    % ############################################################################
    % ##### 2.) Calc the weights                                             #####
    % ############################################################################
    
    max_scon_end_time = max(t_scon_end_jd);
    min_scon_end_time = min(t_scon_end_jd);
     
    % ##### Loop over all possible scan-cons #####
    for i_scon = 1 : number_of_scan_cons
        
        % ### Numbe of observation ###
        % Nomralized by the maximal number of observations with all stations partizipating
        scon_weights_obs_num(i_scon) = scon_weights_obs_num(i_scon) + number_of_obs(i_scon) / (number_of_all_stations * (number_of_all_stations - 1)/2);
        
        
        % ### Sky coverage ###
        % Nomralized by the maximal number of all partizipating stations 
        scon_weights_sky_cov(i_scon) = scon_weights_sky_cov(i_scon) + sky_coverage(i_scon) / number_of_all_stations;
        
        
        % ### Scan end time ###
        % Only apply this property for weighting, if the end times differ between scan configurations!
        if (max_scon_end_time ~= min_scon_end_time)
            scon_weights_scan_end(i_scon) = scon_weights_scan_end(i_scon) +((max_scon_end_time - t_scon_end_jd(i_scon)) / (max_scon_end_time - min_scon_end_time))/10; % max = 0.1
        end
        
        % Sum up the individual weights:
        scon_weights(i_scon) = scon_weights_obs_num(i_scon) + scon_weights_sky_cov(i_scon) +  scon_weights_scan_end(i_scon);
        
        % Save wheight for each scan con.:
        scan_cons(i_scon).weight = scon_weights(i_scon);
        
% TEST TEST TEST      
% fprintf('- scon_id: %3d, sky_cov_w: %2d, obs_num_w: %2d, scan_end_w: %2d => sum: %2d\n', i_scon, scon_weights_sky_cov(i_scon), scon_weights_obs_num(i_scon), scon_weights_scan_end(i_scon), scon_weights(i_scon));
% TEST TEST TEST
        
    end % for i_scon = 1 : number_of_scan_cons
    
    
    
    % ############################################################################
    % ##### 3.) Sort all scan cons by applying the weigths                   #####
    % ############################################################################
    [xx, sort_index] = sort(scon_weights);

    
end % function


%% #################### Sub routines #########################


