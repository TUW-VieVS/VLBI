% #########################################################################
% #     read_vso
% #########################################################################
%
% DESCRITPION
% This function reads files in the VieVS simple observation file format (.vso)
% and creates the antenna, sources and scan structures.
%
%
% Possible VSO formats:
%
% 1.) ##### Baseline delay only: #####
% | Delay reference epoch            | Stat. 1 name  | Stat. 2 name  | source name | obs. type | baseline delay [nsec]/[sec] |
% | YYYY MM DD hh mm ss.ssssssssssss | <max 8 char.> | <max 8 char.> | <string>    |   s/q     |      float                  |
%
% 2.) ##### Baseline delay + formal error: #####
% | Delay reference epoch            | Stat. 1 name  | Stat. 2 name  | source name | obs. type | baseline delay [nsec]/[sec] | delay formal error [nsec]/[sec] |
% | YYYY MM DD hh mm ss.ssssssssssss | <max 8 char.> | <max 8 char.> | <string>    |   s/q     |      float                  |      float                      |
%
% 3.) ##### Baseline delay + formal error  and Ion. delay + formal error: #####
% | Delay reference epoch            | Stat. 1 name  | Stat. 2 name  | source name | obs. type | baseline delay [nsec]/[sec] | delay formal error [nsec]/[sec] | Ion delay [nsec] | Ion formal error [nsec] |
% | YYYY MM DD hh mm ss.ssssssssssss | <max 8 char.> | <max 8 char.> | <string>    |   s/q     |      float                  |      float                      |      float       |      float              |
%
% 4.) ##### No observationas (used e.g. to calculate a priori delays in vie_mod for the correlator input model) #####
% - If the input vso file has only 10 colums, this formati is selected automatically
% - If the second vso header line (starting with "#") contains the phrase "noDelay", all entries after colums 10 are skipped (= vso format 4 is selected)
% | Delay reference epoch            | Stat. 1 name  | Stat. 2 name  | source name | obs. type |
% | YYYY MM DD hh mm ss.ssssssssssss | <max 8 char.> | <max 8 char.> | <string>    |   s/q     |
%
% 5.) ##### Baseline delay + formal error, Ion. delay + formal error, cable corr: #####
% | Delay reference epoch            | Stat. 1 name  | Stat. 2 name  | source name | obs. type | baseline delay [nsec]/[sec] | delay formal error [nsec]/[sec] | Ion delay [nsec] | Ion formal error [nsec] | cable corr. stat. 1 [ns] | cable corr. stat. 2 [ns] |
% | YYYY MM DD hh mm ss.ssssssssssss | <max 8 char.> | <max 8 char.> | <string>    |   s/q     |      float                  |      float                      |      float       |      float              |  float                   |  float                   |
%
% 6.) ##### Baseline delay + formal error, Ion. delay + formal error, cable corr, met. data (T, p, e): #####
% | Delay reference epoch            | Stat. 1 name  | Stat. 2 name  | source name | obs. type | baseline delay [nsec]/[sec] | delay formal error [nsec]/[sec] | Ion delay [nsec] | Ion formal error [nsec] | cable corr. stat. 1 [ns] | cable corr. stat. 2 [ns] | T stat. 1 [°C] | p stat. 1 [mbar] | Rel. humidity e stat. 1 |T stat. 2 [°C] | p stat. 2 [mbar] | Rel. humidity e stat. 2 |
% | YYYY MM DD hh mm ss.ssssssssssss | <max 8 char.> | <max 8 char.> | <string>    |   s/q     |      float                  |      float                      |      float       |      float              |  float                   |  float                   | float          | float            | float                   |float          | float            | float                   |    
%
% Note 1:
% ### Units in the scan structure: ###
% delay                         : [sec]
% delay formal error            : [sec]
% ion. delay                    : [nsec]
% ion. delay formal error       : [nsec]
% cable delays                  : [nsec]
%
% Note 2:
% Baseline delays + their formal errors can have the units [sec] or [nsec] in the input vso files!
% => Please select the correct option with the flag "input_delay_format".
%
% Note 3:
% Different VSO formats are recognized automatically. Notification text is written to the CW.
% => Comments start with "#"!
% There is also the possibility to define tags in the header lines of the input vso file to declare a specific format. 
% The following tags are availabel:
% "# noDelay"                       => VSO format 4 is used (all further columns, if available, are ignored)
% "# geocentricAprioriDelay"        => VSO format 4 is used (all further columns, if available, are ignored) + scans in the scan structure will not be sorted (flag_sort_scans_by_epoch = false)
%
% Note 4:
% Observations with the same ref. epoch and the same observation target (in the input vso file) are assumed to be part of the same scan!!
% => sub-netting is supported: Observations with the same epoch, but a different observation target belong to different scans!
%
%
%
% AUTHOR 
%   A. Hellerschmied, 2016-07-12
%
% INPUT
%  - vso_file_path          - filepath of VSO file (string)
%  - vso_file_name          - filename of VSO file (string)
%  - trf                    - VieVS trf data strucure
%  - trf_name_str           - Name of selected TRF realiation within the trf structure (string)
%  - crf                    - VieVS crf data strucure
%  - crf_name_str           - Name of selected CRF realiation within the crf structure (string)
%  - ini_opt                - input options (structure)
%  - sat_orbit_file_path    - 
%  - sat_orbit_file_name    - 
%  - sat_orbit_file_type    - 
%  - parameter              - VieVS parameter structure
%
% OUTPUT
%  - antenna                - VieVS antenna structure
%  - sources                - VieVS sources structure
%  - scan                   - VieVS scan structure
%  - parameter              - VieVS parameter structure
%
% COUPLING
%  - read_sp3.m
%  - orbit_data2sources.m
%  - read_sat_ephem_trf.m 
%  - modjuldat.m
%
% CHANGES
%  - 2016-08-01, A. Hellerschmied:  - Orbit data can be read from SP3 files (read_sp3.m) and written to the sources.s structure (orbit_data2sources.m) 
%                                   - possibility added to have VSO files with different columns (VSO versions 1, 2 and 3; see header for description) 
%  - 2016-08-04, A. Hellerschmied:  CW notifications changed. 
%  - 2016-08-08, A. Hellerschmied:  Fields added to sources.s structure: x_crf, y_crf, z_crf (filled in vie_mod)
%  - 2016-09-01, A. Hellerschmied:  Unit of delay values in scan.obs.obs cahnged from [ns] to [sec].
%                                   - Added option to select wheter the input delays (+ formal errors, + ion. delays, + Ion. formal error) are given in [sec] or in [nsec]; var.: delay_format
%  - 2016-09-21, A. Hellerschmied:  - Possibility added to read space craft ephemerids (TRF posit vions and velocities) from an external file usign the function read_sat_ephem_trf.m
%                                   - Possibility added to load vso files without observation data (= vso_format 4) => see info in the header!
%  - 2016-09-26, H. Krasna: Changes related to updated supersource file: CommonNames do not exist anymore
%  - 2016-10-24, A. Hellerschmied: - Added VSO format 5 and 6, see description in header
%                                  - Corrected inout units (see header)
%  - 2016-11-28, A. Hellerschmied  - modjuldat.m used instead of juliandate.m (more significant decimal figures)
%  - 2016-12-05, A. Hellerschmied: Added fields 'mjd', 'sec_of_day', 'year', 'month', 'day', 'hour', 'minu', 'sec' to sources structure preallocation
%  - 2017-01-24, D. Landskron: preallocation of GPT2 changed to GPT3, GPT removed
%  - 2017-02-06, A. Hellerschmied: - Support of sub-netting added => Observations belonging to two scans distinguish in "scan epoch" AND/OR "observed source"!
%                                  - Option added to sort scans (read from the input vso file) by the scan epoch ("flag_sort_scans_by_epoch"). Default: "flag_sort_scans_by_epoch = false"
%                                  - parameter struct added as in/output argument
%  - 2017-02-09, D. Landskron: Preallocation extended
%  - 2017-02-22, A. Hellerschmied: antenna.psd initialized
%
%
function [antenna, sources, scan, parameter] = read_vso(vso_file_path, vso_file_name, trf, trf_name_str, crf, crf_name_str, ini_opt, sat_orbit_file_path, sat_orbit_file_name, sat_orbit_file_type, parameter)

% ##### Load constants #####


% ##### Options #####
error_code_invalid_met_data = -999; % Error corde for missing met. data in NGS file (numerical)
flag_sort_scans_by_epoch    = true; % Sort scans (read from the input vso file) by the scan epoch? [true=t/false]

% => Define parameters (for all scans) which are not available in the VSO fils:
scan_obs_q_code_ion     = 0;
scan_obs_q_code         = 0;

% => Define formal error of baseline delay and ion. delay [sec], if these values are not available in the input vso file:
% - For satellite observations (formal error of baseline delay):
formal_error_sat_obs_sec = 30e-12;
% - For quasar observations (formal error of baseline delay):
formal_error_quasar_obs_sec = 30e-12;


% => Select wheter the input delays (+ formal errors) are given in [sec] or in [nsec]
% - sec     : [sec]
% - nsec    : [nsec] (e.g. as provided in NGS files)
input_delay_format = 'nsec';
switch(input_delay_format)
    case 'sec'
        sec_unit_conversion = 1;        
    case 'nsec'
        sec_unit_conversion = 1e-9;
end


% Init.:
small = 1e-11; % for matching mjd epochs (e.g. outliers)
vso_format = 0;

% ##### Check input options (ini_opt) and initialise them, if required #####
if ~isfield(ini_opt,'minel') % Min. elevation
    ini_opt.minel=0;
end
if ~isfield(ini_opt,'Qlim') % Quality code limit (NOT CONDSIDERED HERE!)
    ini_opt.Qlim=0;
end
%if sta_excl, sour_excl, or scan_excl not specified, initialize them as empty matrices
if ~isfield(ini_opt,'sta_excl')     % Excluded stations from OPT file
    ini_opt.sta_excl=[];
end
if ~isfield(ini_opt,'sour_excl')    % Excluded sources from OPT file
    ini_opt.sour_excl=[];
end
if ~isfield(ini_opt,'scan_excl')    % From outlier file: | stat. 1 name | stat.2 name | MJD |
    ini_opt.scan_excl=[];
end
if ~isfield(ini_opt,'no_cab')       % There is no cable cal. information in the vso format so far
    ini_opt.no_cab=[];
end
if ~isfield(ini_opt,'bas_excl')     % Excluded baselines from OPT file
    ini_opt.bas_excl=[];
end


% ##### preallocate structures #####
% scan        = struct('mjd', []);
% scansource  = [];
% sources     = struct('name', [], 'IERSname', [], 'ICRFdes', [], 'ra2000', [], 'de2000', [], 'ra_sigma', [], 'de_sigma', [], 'corr', [], 'in_crf', [], 'flag_defining', []);


% ###############################
% ##### Load all data files #####
% ###############################

% #### Load VSO file ####
fid_vso = fopen([vso_file_path, vso_file_name], 'r');
if fid_vso == -1
    error(['Cannot open VSO file: ', vso_file_path, vso_file_name]);
end 

% Distinguish between VSO formats and recognize flags in the file header:
in_str = fgetl(fid_vso);
while strcmp(in_str(1), '#')
    in_str = fgetl(fid_vso);
    if logical(strfind(in_str, 'noDelay'))
        vso_format = 4;                     % no delay
    end
    if logical(strfind(in_str, 'geocentricAprioriDelay'))
        flag_sort_scans_by_epoch = false;    % Do not sort scans by epoch
        vso_format = 4;                     % no delay
        
        % Set special handling tag in parameter struct.:
        parameter.vie_mod.special_handling_tag = 'geocentricAprioriDelay';
    end
end
vso_data = textscan(in_str, '%s');
num_of_col_in_vso_file = size(vso_data{1}, 1);
if vso_format == 0
    switch(num_of_col_in_vso_file)
        case 11
            vso_format = 1;
            fprintf(1, '\n');
            fprintf(1, 'VSO format 1: Only baseline delay.\n');
            fprintf(1, '   - Formal errors are set to predefined values (sat.: %e sec, quasars: %e sec)\n', formal_error_sat_obs_sec, formal_error_quasar_obs_sec);
            fprintf(1, '   - Ion. delay is set to zero.\n');
        case 12
            vso_format = 2;
            fprintf(1, '\n');
            fprintf(1, 'VSO format 2: Baseline delay + formal error.\n');
            fprintf(1, '   - Ion. delay is set to zero.\n');
        case 14
            vso_format = 3;
            fprintf(1, '\n');
            fprintf(1, 'VSO format 3: Baseline delay + formal error and Ion. delay + formal error.\n');
        case 10
            vso_format = 4;
            fprintf(1, '\n');
            fprintf(1, 'VSO format 4:\n');
            fprintf(1, '   - NO OBSERVATIONS INCLUDED!\n');
            fprintf(1, '   - Formal errors are set to predefined values (sat.: %e sec, quasars: %e sec)\n', formal_error_sat_obs_sec, formal_error_quasar_obs_sec);
            fprintf(1, '   - Ion. delay is set to zero.\n');
        case 16
            vso_format = 5;
            fprintf(1, '\n');
            fprintf(1, 'VSO format 5: Baseline delay + formal error; Ion. delay + formal error; cable corrections\n');
        case 22
            vso_format = 6;
            fprintf(1, '\n');
            fprintf(1, 'VSO format 6: Baseline delay + formal error; Ion. delay + formal error; cable corrections; met. data (T, p, e)\n');
        otherwise
            error('Unknown vso format.');
    end
    
    frewind(fid_vso);

    % vso_data = textscan(fid_vso, '%d %d %d %d %d %f %s %s %s %s %f', 'CommentStyle', '#');
    % vso_data = textscan(fid_vso, '%f %f %f %f %f %f %s %s %s %s %f', 'CommentStyle', '#');
    % vso_data = textscan(fid_vso, '%f %f %f %f %f %f %s %s %s %s %f %f %f %f', 'CommentStyle', '#');
    vso_data = textscan(fid_vso, '%f %f %f %f %f %f %s %s %s %s %f %f %f %f %f %f %f %f %f %f %f %f', 'CommentStyle', '#');
    
    
elseif vso_format == 4
    fprintf(1, '\n');
    fprintf(1, '   - NO OBSERVATIONS INCLUDED!\n');
    fprintf(1, '   - Formal errors are set to predefined values (sat.: %e sec, quasars: %e sec)\n', formal_error_sat_obs_sec, formal_error_quasar_obs_sec);
    fprintf(1, '   - Ion. delay is set to zero.\n');
    
    frewind(fid_vso);
    
    % Skip all columns atfter column 10
    vso_data = textscan(fid_vso, '%f %f %f %f %f %f %s %s %s %s %*[^\n]', 'CommentStyle', '#');
end

if ~flag_sort_scans_by_epoch
    fprintf(1, '\n');
    fprintf(1, 'PLEASE NOTE: Observations will not be sorted by reference time in the source structure (flag_sort_scans_by_epoch = false)!\n');
    fprintf(1, '\n');
end


fclose(fid_vso);

% Distinguish between VSO formats:


% Prepare data from vso file:
% | Delay reference epoch            | Stat. 1 name  | Stat. 2 name  | source name | obs. type | baseline delay [nsec] | delay formal error [nsec] | Ion delay [nsec] | Ion formal error [nsec] |
% | YYYY MM DD hh mm ss.ssssssssssss | <max 8 char.> | <max 8 char.> | <string>    |   s/q     |      float            |      float                |      float       |      float              |

year 			= vso_data{1};
mon 			= vso_data{2};
day 			= vso_data{3};
hour 			= vso_data{4};
minu 			= vso_data{5};
sec 			= vso_data{6};
stat_1_name_str = vso_data{7};
stat_2_name_str = vso_data{8};
source_name_str = vso_data{9};
obs_type_str 	= vso_data{10};

if vso_format == 1
    delay           = vso_data{11};	
elseif vso_format == 2
    delay           = vso_data{11};	
    delay_sig       = vso_data{12};
elseif vso_format == 3
    delay           = vso_data{11};	
    delay_sig       = vso_data{12};
    ion_del         = vso_data{13};
    ion_sig         = vso_data{14};
elseif vso_format == 5
    delay           = vso_data{11};	
    delay_sig       = vso_data{12};
    ion_del         = vso_data{13};
    ion_sig         = vso_data{14};
    cab_stat1       = vso_data{13}; % cable correction of station 1 [ns]
    cab_stat2       = vso_data{14}; % cable correction of station 2 [ns]
elseif vso_format == 6
    delay           = vso_data{11};	
    delay_sig       = vso_data{12};
    ion_del         = vso_data{13};
    ion_sig         = vso_data{14};
    cab_stat1       = vso_data{15}; % cable correction of station 1 [ns]
    cab_stat2       = vso_data{16}; % cable correction of station 2 [ns]
    t_stat1         = vso_data{17};	% Ambient atmospheric temperature at site 1 [deg. C]
    p_stat1         = vso_data{18}; % Ambient atmospheric barometric pressure at site 1 [mbar]
    e_stat1         = vso_data{19}; % Realative humitity at site 1 
    t_stat2         = vso_data{20}; % Ambient atmospheric temperature at site 2 [deg. C]
    p_stat2         = vso_data{21}; % Ambient atmospheric barometric pressure at site 2 [mbar]
    e_stat2         = vso_data{22}; % Realative humitityat site 2
    
end

number_of_all_obs   = length(year);

% Check obs_type_str  
if sum(strcmp(obs_type_str, 'q') | strcmp(obs_type_str, 's')) ~= number_of_all_obs
    error('Invalid entries for "obs_type_str" in input VSO file!');
end


% Convert epochs to MJD:
% mjd = juliandate([year, mon, day, hour, minu, sec]) -  2400000.5;
mjd = modjuldat(year, mon, day, hour, minu, sec);


% DOY:
[itim,doy] = dday(year, mon, day, hour, minu);


[~, ~, local_soure_id] = unique(source_name_str);
number_of_all_scans = size(unique([mjd, local_soure_id], 'rows'), 1); % Individual scans differ in epoch ("mjd") AND observed source ("local_soure_id")!

% ###############################
% ##### Check scans         #####
% ###############################

% Init. index of valid scans:
index_valid_scans = ones(number_of_all_obs, 1);


% Throw out all invalid scans
% Checks:
% - Outliers
% - Excluded stations
% - Excluded baselines
% - Excluded sources
% - 

% => Just get the indices of invalid scans and exclude them! 

% #### Check if outliers should be excluded ####
excl_outliers_ind = false(number_of_all_obs, 1);
if ~isempty(ini_opt.scan_excl)
    for i_tmp = 1 : length(ini_opt.scan_excl)
        
        % Check outlier epoch:
        excl_out_mjd_tmp_ind = abs(mjd - ini_opt.scan_excl(i_tmp).mjd) < small;
        
        % Check stations of baseline:
        excl_out_bl_ind = (strcmp(stat_1_name_str, deblank(ini_opt.scan_excl(i_tmp).sta1)) & strcmp(stat_2_name_str, deblank(ini_opt.scan_excl(i_tmp).sta2))) | (strcmp(stat_1_name_str, deblank(ini_opt.scan_excl(i_tmp).sta2)) & strcmp(stat_2_name_str, deblank(ini_opt.scan_excl(i_tmp).sta1)));        
   
        excl_outliers_tmp_ind = excl_out_mjd_tmp_ind & excl_out_bl_ind;
        excl_outliers_ind = excl_outliers_ind | excl_outliers_tmp_ind;
    end
end

% #### Check if any of the stations are in the list of stations to be excluded ####
excl_stat_ind = false(number_of_all_obs, 1);
if ~isempty(ini_opt.sta_excl)
    for i_tmp = 1 : size(ini_opt.sta_excl, 1)
        
        % Get obs. indizes for station to be excluded:
        excl_stat_1_ind = strcmp(stat_1_name_str, deblank(ini_opt.sta_excl(i_tmp,:)));
        excl_stat_2_ind = strcmp(stat_2_name_str, deblank(ini_opt.sta_excl(i_tmp,:)));
        excl_stat_ind_tmp = excl_stat_1_ind | excl_stat_2_ind;
        
        % Check, if only a certain time span should be excluded:
        if ini_opt.sta_excl_start(i_tmp) ~= 0
           excl_time_ind = (mjd > ini_opt.sta_excl_start(i_tmp)) & (mjd < ini_opt.sta_excl_end(i_tmp));
           excl_stat_ind_tmp = excl_stat_ind_tmp & excl_time_ind;
        end
        excl_stat_ind = excl_stat_ind | excl_stat_ind_tmp;
    end
end

% #### Check if any of the sources are in the list of sources to be excluded ####
excl_src_ind = false(number_of_all_obs, 1);
if ~isempty(ini_opt.sour_excl)
    for i_tmp = 1 : size(ini_opt.sour_excl, 1)
        
        % Get obs. indizes for sources to be excluded:
        excl_src_ind_tmp = strcmp(source_name_str, deblank(ini_opt.sour_excl(i_tmp,:)));
        
        % Check, if only a certain time span should be excluded:
        if ini_opt.sour_excl_start(i_tmp) ~= 0
           excl_time_ind = (mjd > ini_opt.sour_excl_start(i_tmp)) & (mjd < ini_opt.sour_excl_end(i_tmp));
           excl_src_ind_tmp = excl_src_ind_tmp & excl_time_ind;
        end
        excl_src_ind = excl_src_ind | excl_src_ind_tmp;
    end
end

% #### Check if baselines should be excluded ####
excl_bl_ind = false(number_of_all_obs, 1);
if ~isempty(ini_opt.bas_excl)
    for i_tmp = 1 : length(ini_opt.bas_excl)
        excl_bl_ind = excl_bl_ind | ((strcmp(stat_1_name_str, deblank(ini_opt.bas_excl(i_tmp).sta1)) & strcmp(stat_2_name_str, deblank(ini_opt.bas_excl(i_tmp).sta2))) | (strcmp(stat_1_name_str, deblank(ini_opt.bas_excl(i_tmp).sta2)) & strcmp(stat_2_name_str, deblank(ini_opt.bas_excl(i_tmp).sta1))));
    end
end

% Combine all indices to exclude observations:
index_valid_scans = ~(excl_outliers_ind | excl_stat_ind | excl_src_ind | excl_bl_ind);

% Apply indices to exclude observations:
year 			= year(index_valid_scans);
mon 			= mon(index_valid_scans);
day 			= day(index_valid_scans);
hour 			= hour(index_valid_scans);
minu 			= minu(index_valid_scans);
sec 			= sec(index_valid_scans);
stat_1_name_str = stat_1_name_str(index_valid_scans);
stat_2_name_str = stat_2_name_str(index_valid_scans);
source_name_str = source_name_str(index_valid_scans);
obs_type_str 	= obs_type_str(index_valid_scans);
% delay           = delay(index_valid_scans);
mjd             = mjd(index_valid_scans);  
doy             = doy(index_valid_scans); 

% if vso_format == 2
%     delay_sig       = delay_sig(index_valid_scans);
% elseif vso_format == 3
%     delay_sig       = delay_sig(index_valid_scans);
%     ion_del         = ion_del(index_valid_scans);
%     ion_sig         = ion_sig(index_valid_scans);
% end

if vso_format == 1
    delay           = delay(index_valid_scans);	
elseif vso_format == 2
    delay           = delay(index_valid_scans);	
    delay_sig       = delay_sig(index_valid_scans);
elseif vso_format == 3
    delay           = delay(index_valid_scans);	
    delay_sig       = delay_sig(index_valid_scans);
    ion_del         = ion_del(index_valid_scans);
    ion_sig         = ion_sig(index_valid_scans);
elseif vso_format == 5
    delay           = delay(index_valid_scans);	
    delay_sig       = delay_sig(index_valid_scans);
    ion_del         = ion_del(index_valid_scans);
    ion_sig         = ion_sig(index_valid_scans);
    cab_stat1       = cab_stat1(index_valid_scans);
    cab_stat2       = cab_stat2(index_valid_scans);
elseif vso_format == 6
    delay           = delay(index_valid_scans);	
    delay_sig       = delay_sig(index_valid_scans);
    ion_del         = ion_del(index_valid_scans);
    ion_sig         = ion_sig(index_valid_scans);
    cab_stat1       = cab_stat1(index_valid_scans);
    cab_stat2       = cab_stat2(index_valid_scans);    
    t_stat1         = t_stat1(index_valid_scans);	% Ambient atmospheric temperature at site 1 [deg. C]
    p_stat1         = p_stat1(index_valid_scans);   % Ambient atmospheric barometric pressure at site 1 [mbar]
    e_stat1         = e_stat1(index_valid_scans);   % Realative humitity at site 1 
    t_stat2         = t_stat2(index_valid_scans);   % Ambient atmospheric temperature at site 2 [deg. C]
    p_stat2         = p_stat2(index_valid_scans);   % Ambient atmospheric barometric pressure at site 2 [mbar]
    e_stat2         = e_stat2(index_valid_scans);   % Realative humitityat site 2
end

% Get names of all stations with observations left:
stat_names_remaining_str            = unique([stat_1_name_str; stat_2_name_str]);
[source_names_remaining_str, ia, local_soure_id]    = unique(source_name_str);

% Get the type of the sources:
source_type_str = obs_type_str(ia);

number_of_remaining_obs = sum(index_valid_scans);



% Get number of remaining scans and the indices to assign observations to scans (obs2scan_ind)
% => Observations with the same epoch and the same observation target are assumed to be part of the same scan (=> sub-netting is supported!)!!
if(flag_sort_scans_by_epoch) % Sort input scans by epoch?
    [mjd_scan, scan2obs_ind, obs2scan_ind] = unique([mjd, local_soure_id], 'rows', 'sorted'); % Individual scans differ in epoch ("mjd") AND observed source ("local_soure_id")!
else
    [mjd_scan, scan2obs_ind, obs2scan_ind] = unique([mjd, local_soure_id], 'rows', 'stable'); % Individual scans differ in epoch ("mjd") AND observed source ("local_soure_id")!
end
mjd_scan = mjd_scan(:, 1);
number_of_remaining_scans = length(mjd_scan);



% #### Load sat orbit file ####
% Only, if at least one satellite was observed!


% ##### Write some information on excluded observations to CW #####
fprintf(1, '\n');
fprintf(1, 'Number of excluded observations due to...\n');
fprintf(1, '  ...outliers:          %d\n', sum(excl_outliers_ind));
fprintf(1, '  ...excluded stations: %d\n', sum(excl_stat_ind));
fprintf(1, '  ...excluded sources:  %d\n', sum(excl_src_ind));
fprintf(1, '  ...excluded baslines: %d\n', sum(excl_bl_ind));
fprintf(1, '%d observations of %d excluded. %d observations left.\n', sum(~index_valid_scans), number_of_all_obs, number_of_remaining_obs);
fprintf(1, '%d scans of %d excluded. %d scans left.\n', number_of_all_scans - number_of_remaining_scans, number_of_all_scans, number_of_remaining_scans);
fprintf(1, '\n');




% ###############################
% #####  Create structures  #####
% ###############################

% ###############################################################
% ##### ANTENNA                                             #####
% ###############################################################

number_of_remaining_stat = length(stat_names_remaining_str);

% preallocate the structure:
antenna(number_of_remaining_stat) = struct('IDsuper', [], 'in_trf', [], 'name', [], 'x', [], 'y', [], 'z', [], 'firstObsMjd', [], 'vx', [], 'vy', [], 'vz', [], 'epoch', [], 'start', [], 'end', [],...
                                            'x_sigma', [], 'y_sigma', [], 'z_sigma', [], 'vx_sigma', [], 'vy_sigma', [], 'vz_sigma', [], 'thermal', [], 'comments', [], 'domes', [], 'code', [],... 
                                            'ecc', [], 'ecctype', [], 'axtyp', [], 'offs', [], 'gpt3pres', [], 'gpt3temp', [], 'gpt3e', [], 'gpt3', [], 'noGrad', [], ...
                                            'cto', [], 'cta', [], 'cnta_dx', [], 'vm1', [], 'opl', [], 'numobs', [], 'lastObsMjd', [], 'psd', []);
[antenna.gpt3] = deal(struct('pres', [], 'temp', [],'lapserate', [], 'e', [], 'ah', [], 'aw', [], 'undulation', []));



for i_stat = 1 : number_of_remaining_stat
    
    stat_name_str       = stat_names_remaining_str{i_stat};
    stat_name_8char_str = sprintf('%s%s', stat_names_remaining_str{i_stat}, blanks(8 - length(stat_names_remaining_str{i_stat})));
    trf_id              = find(strcmpi({trf.name}, stat_name_8char_str));
    
    % #### Check, if there is an entry for the current station in the superstation file. If not => Error Msg. and abort! ####
    if isempty(trf_id)
        error('ERROR (read_vso.m): No entry for station %s in the superstation file! Add it!\n', stat_name_8char_str);
    end
    
    antenna(i_stat).IDsuper     = trf_id;
    antenna(i_stat).name        = trf(trf_id).name;
    
    % Get all valid scans, where the current station participated:
    valid_obs_with_cur_stat_ind = strcmp(stat_1_name_str, stat_name_str) | strcmp(stat_2_name_str, stat_name_str);
    
    % #### Get the number of observations of the current station ####
     antenna(i_stat).numobs     = sum(valid_obs_with_cur_stat_ind);
    
    % #### Get mjd of first and last observation ####
    mjd_tmp = mjd(valid_obs_with_cur_stat_ind);
    antenna(i_stat).firstObsMjd = mjd_tmp(1);
    antenna(i_stat).lastObsMjd  = mjd_tmp(end);
    
    
    % #### Get break (coordinate epoch) ####
    
    % ### Check, if there are coordinates for the chosen TRF available and get the "break_id" ###
    if ~isempty(trf(trf_id).(trf_name_str))
        
        % if start/end fields are missing => add them!
        if ~isfield(trf(trf_id).(trf_name_str).break(1),'start') % 25/06/2014
            trf(trf_id).(trf_name_str).break.start=[];
            trf(trf_id).(trf_name_str).break.end=[];
        end
        % if epoch information is missing, put 99999 to start and 0 to end
        emptyStartLog = cellfun(@isempty,{trf(trf_id).(trf_name_str).break.start});
        emptyEndLog   = cellfun(@isempty,{trf(trf_id).(trf_name_str).break.end});
        if sum(emptyStartLog) > 0
            try
                trf(trf_id).(trf_name_str).break(emptyStartLog).start = zeros(sum(emptyStartLog), 1); % 25/06/2014
            catch
                keyboard;
            end
        end
        if sum(emptyEndLog) > 0
            trf(trf_id).(trf_name_str).break(emptyEndLog).end = repmat(99999, sum(emptyEndLog), 1); % 25/06/2014
        end

        break_id = find(antenna(i_stat).firstObsMjd >= [trf(trf_id).(trf_name_str).break.start] & antenna(i_stat).firstObsMjd <= [trf(trf_id).(trf_name_str).break.end]);

    else % Not found...
        fprintf('Station not found in %s - write coordinates to ASCII file and renew superstation file!\n', trf_name_str);
        break_id = [];
    end
    
    
    % #### If the "break_id" was not found in selected TRF ####
    % => Get coordinates from VieVS TRF (= backup TRF!)
    if isempty(break_id)

        fprintf('No %s coordinates for in %s ... get vievsTrf coordinates (backup).\n', trf_name_str, stat_name_str);

        % if there is no start break - only approx coords (e.g. station "VLA     ")
        if ~isfield(trf(trf_id).vievsTrf.break, 'start')
            break_id = 1;
        else
            break_id = find(antenna(i_stat).firstObsMjd >= [trf(trf_id).vievsTrf.break.start] & antenna(i_stat).firstObsMjd <= [trf(trf_id).vievsTrf.break.end]);
        end

        % ### Check, if VieVS TRF coordinates were found: ###
        if isempty(break_id)
            warning('Station not found in vievs TRF - write coordinates to ASCII file and renew superstation file!\n');
            keyboard;
        end

        % ### Get the break sub-structure from the trf structure where coordinatess should be taken from ###
        curBreak_substruct = trf(trf_id).vievsTrf.break(break_id);

        % ### Define, that the station (coordiantes from VieVS TRF, as BACKUP!) is NO datum station!
        antenna(i_stat).in_trf = 0;

    else % Station coordinates were taken from selected TRF file

        % ### Get the break sub-structure from the trf structure where coords should be taken ###
        curBreak_substruct = trf(trf_id).(trf_name_str).break(break_id);

        % ### Define if the current station should be a datum station ###
        % => since station is found in chosen trf, it would be 1 by default, but it could be 0, if it is set to 0 in "manual" trf file
        % => If the "indatum" flag has a "NaN" value => in_trf = 1
        if isfield(curBreak_substruct, 'indatum') && ~isnan(curBreak_substruct.indatum)
            antenna(i_stat).in_trf = curBreak_substruct.indatum;
        else
            antenna(i_stat).in_trf = 1;
        end

    end % if isempty(break_id)
    
    % ##### Add data to antenna structure #####
    antenna(i_stat).x              = curBreak_substruct.x;
    antenna(i_stat).y              = curBreak_substruct.y;
    antenna(i_stat).z              = curBreak_substruct.z;
    
    % Check, if velocity information is available:
    if isfield(curBreak_substruct, 'vx')
        antenna(i_stat).vx         = curBreak_substruct.vx;
        antenna(i_stat).vy         = curBreak_substruct.vy;
        antenna(i_stat).vz         = curBreak_substruct.vz;
        antenna(i_stat).epoch      = curBreak_substruct.epoch;
        antenna(i_stat).start      = curBreak_substruct.start;
        antenna(i_stat).end        = curBreak_substruct.end;
    else % if coords are taken from blokq.dat: not even that info is given...
        antenna(i_stat).vx         = 0;
        antenna(i_stat).vy         = 0;
        antenna(i_stat).vz         = 0;
        antenna(i_stat).epoch      = 0;
        antenna(i_stat).start      = 0;
        antenna(i_stat).end        = 99999;
    end
    % check if error measures are available:
    if isfield(curBreak_substruct,'x_sigma')
        antenna(i_stat).x_sigma    = curBreak_substruct.x_sigma;
        antenna(i_stat).y_sigma    = curBreak_substruct.y_sigma;
        antenna(i_stat).z_sigma    = curBreak_substruct.z_sigma;
        antenna(i_stat).vx_sigma   = curBreak_substruct.vx_sigma;
        antenna(i_stat).vy_sigma   = curBreak_substruct.vy_sigma;
        antenna(i_stat).vz_sigma   = curBreak_substruct.vz_sigma;
    else
        antenna(i_stat).x_sigma    = [];
        antenna(i_stat).y_sigma    = [];
        antenna(i_stat).z_sigma    = [];
        antenna(i_stat).vx_sigma   = [];
        antenna(i_stat).vy_sigma   = [];
        antenna(i_stat).vz_sigma   = [];
    end

    antenna(i_stat).thermal        = trf(trf_id).antenna_info;
    antenna(i_stat).comments       = trf(trf_id).comments;
    antenna(i_stat).domes          = trf(trf_id).domes;
    antenna(i_stat).code           = trf(trf_id).code;

    % Eccentricity from superstation file
    if isempty(trf(trf_id).ecc)
        antenna(i_stat).ecc        = [0 0 0];
        antenna(i_stat).ecctype    = 'NEU';
    else
        for i_tmp = 1 : length(trf(trf_id).ecc.break)
            tecc        = trf(trf_id).ecc.break(i_tmp).starting;
            y           = str2num(tecc(1:4));
            m           = str2num(tecc(6:7));
            d           = str2num(tecc(9:10));
            h           = str2num(tecc(12:13));
            mi          = str2num(tecc(15:16));
            s           = 0;
            mjd_start   = modjuldat(y,m,d,h,mi,s);
            tecc        = [];
            tecc = trf(trf_id).ecc.break(i_tmp).ending;
            y           = str2num(tecc(1:4));
            m           = str2num(tecc(6:7));
            d           = str2num(tecc(9:10));
            h           = str2num(tecc(12:13));
            mi          = str2num(tecc(15:16));
            s           = 0;
            mjd_end     = modjuldat(y,m,d,h,mi,s);
            tecc        = [];
            if (floor(antenna(i_stat).firstObsMjd) >= mjd_start) && (floor(antenna(i_stat).firstObsMjd) <= mjd_end)
                break;
            end
        end
        antenna(i_stat).ecc        = [trf(trf_id).ecc.break(i_tmp).FCE, trf(trf_id).ecc.break(i_tmp).SCE, trf(trf_id).ecc.break(i_tmp).TCE];
        antenna(i_stat).ecctype    = trf(trf_id).ecc.break(i_tmp).type_e;
    end

    % mounting type from superstation file
    antenna(i_stat).axtyp      = '';
    if ~isempty(trf(trf_id).antenna_info) % if there was information in antenna-info.txt
        antenna(i_stat).axtyp  = trf(trf_id).antenna_info.mount(4:end);
    end
    if strcmp(antenna(i_stat).axtyp, 'XYNO')
        antenna(i_stat).axtyp  = 'X-Y1';
    end

    % axis offset from superstation file
    antenna(i_stat).offs = 0;
    if ~isempty(trf(trf_id).antenna_info) % if there was information in antenna-info.txt
        antenna(i_stat).offs = trf(trf_id).antenna_info.axis_offset; %m
    end

    % axis offset from superstation file
    antenna(i_stat).offs = 0;
    if ~isempty(trf(trf_id).antenna_info) % if there was information in antenna-info.txt
        antenna(i_stat).offs = trf(trf_id).antenna_info.axis_offset; % m
    end

    % Init. flags:
    antenna(i_stat).gpt3pres   = 0;
    antenna(i_stat).gpt3temp   = 0;
    antenna(i_stat).gpt3e      = 0;
                        
end % for i_stat = 1 : length(stat_names_remaining_str)



% #################################################################
% ##### SOURCES                                               #####
% #################################################################

% Init.:
number_of_remaining_sources = length(source_names_remaining_str);
i_quasar    = 0;
i_sat       = 0;

number_of_remainging_sat_sources        = sum(strcmp(source_type_str, 's'));
number_of_remainging_quasar_sources    = sum(strcmp(source_type_str, 'q'));

% Preallocate structure:
% Also nit. lists of source names for quasars and satellites (same names as in input file):
% => Same order as in sources struct to be match the names to get the right source ID for the scan structure

% For quasars:
if number_of_remainging_quasar_sources > 0
    sources.q(number_of_remainging_quasar_sources) = struct('name', [], 'IERSname', [], 'IVSname', [], 'ICRFdes', [], 'ra2000', [], 'de2000', [], 'ra_sigma', [], 'de_sigma', [], ...
                                              'corr', [], 'in_crf', [], 'flag_defining', [], 'firstObsMjd', [], 'numobs', [], 'lastObsMjd', []);
    quasar_name_list{number_of_remainging_quasar_sources}   = []; 
else
    sources.q = [];
end
% For satellites:
if number_of_remainging_sat_sources > 0
    sources.s(number_of_remainging_sat_sources) = struct('name', [], 'numobs', [], 'firstObsMjd', [], 'lastObsMjd', [], 'orbit_file_type', [], 'x_trf', [], 'y_trf', [], 'z_trf', [], 'x_crf', [], 'y_crf', [], 'z_crf', [], 'mjd', [], 'sec_of_day', [], 'year', [], 'month', [], 'day', [], 'hour', [], 'minu', [], 'sec', []);
    sat_name_list{number_of_remainging_sat_sources}         = []; 
else
    sources.s = [];
end


% Loop over all sources:
for i_src = 1 : number_of_remaining_sources

    source_name = source_names_remaining_str{i_src};
    
    % Distinguish between quasar and satellite scan:
    switch(source_type_str{i_src})

        case 'q' % quasar
            
            i_quasar = i_quasar + 1;
            quasar_name_list{i_quasar} = source_name;
            
            % ##### Find index of current source in crf (crf = superstation file) #####
            
            % search at first within IERS names
            crf_ind = strcmp(source_name, {crf.IERSname});
            
            % if not found within IERS names, search within IVS names
            actIVSName = 0;
            if sum(crf_ind) == 0
                crf_ind = strcmp(strtrim(cellstr(char(crf.IVSname))), strtrim(source_name));
                if sum(crf_ind) > 0
                    actIVSName = 1; % 1 if the name of the actual source is the IVS one
                end
            end
            
            % if not existing in supersource file (should not happen)
            if sum(crf_ind) == 0
                error('ERROR: Source %s does not exist in supersource file - add it!\n', source_name);
            end

            % check if chosen-catalog coordinates exist
            if isempty(crf(crf_ind).(crf_name_str))
                fprintf('No %s coordinates for %s ... get vievsCrf coordinates\n', crf_name_str, source_name);
                crfCatalogToTake = 'vievsCrf';
            else
                crfCatalogToTake = crf_name_str;
            end
            
            % add information
            sources.q(i_quasar).name = sprintf('%s%s', source_name, blanks(8 - length(source_name)));
            if actIVSName == 1
                sources.q(i_quasar).IERSname = crf(crf_ind).IERSname;
            else
                sources.q(i_quasar).IERSname = source_name;
            end
            sources.q(i_quasar).IVSname = crf(crf_ind).IVSname;
            if ~isempty(crf(crf_ind).designation)
                sources.q(i_quasar).ICRFdes = crf(crf_ind).designation(6:end);
            else
                sources.q(i_quasar).ICRFdes = 'J               ';
            end
            sources.q(i_quasar).ra2000 = crf(crf_ind).(crfCatalogToTake).ra;
            sources.q(i_quasar).de2000 = crf(crf_ind).(crfCatalogToTake).de;
            
            % add sigmas if exist
            if isfield(crf(crf_ind).(crfCatalogToTake), 'ra_sigma')
                sources.q(i_quasar).ra_sigma = crf(crf_ind).(crfCatalogToTake).ra_sigma;
                sources.q(i_quasar).de_sigma = crf(crf_ind).(crfCatalogToTake).de_sigma;
            else
                sources.q(i_quasar).ra_sigma = [];
                sources.q(i_quasar).de_sigma = [];
            end
            
            % add correlation if exist
            if isfield(crf(crf_ind).(crfCatalogToTake), 'corr')
                sources.q(i_quasar).corr = crf(crf_ind).(crfCatalogToTake).corr;
            else
                sources.q(i_quasar).corr = [];
            end
            
            % in_crf (important for estimating sources)
            if strcmp(crfCatalogToTake, 'vievsCrf') % if the user has chosen vievsCrf
                sources.q(i_quasar).in_crf = crf(crf_ind).(crfCatalogToTake).in_crf; % take the in_crf from the "source.cat" "textfile"
            else
                % if we take backup coords for the current source
                if strcmp(crfCatalogToTake, 'vievsCrf')
                    sources.q(i_quasar).in_crf = 0;
                else
                    sources.q(i_quasar).in_crf = 1;
                end
            end
            
            % David - always add information about defining sources from the ICRF2
            if ~isempty(crf(crf_ind).icrf2)
                sources.q(i_quasar).flag_defining=crf(crf_ind).icrf2.defining;
            else
                sources.q(i_quasar).flag_defining=0;
            end

            % Index of observations the current source:
            src_obs_ind = strcmp(source_name_str, source_name);
            
            % Number of observations
            sources.q(i_quasar).numobs = sum(src_obs_ind);
            
            % MJD of first and last observation:
            mjd_tmp = mjd(src_obs_ind);
            sources.q(i_quasar).firstObsMjd     = mjd_tmp(1);
            sources.q(i_quasar).lastObsMjd      = mjd_tmp(end);
            
        case 's' % satellites
            i_sat = i_sat + 1;
            sat_name_list{i_sat}            = source_name;
            sources.s(i_sat).name           = source_name;
            
            % Index of observations the current source:
            src_obs_ind = strcmp(source_name_str, source_name);
            
            % Number of observations
            sources.s(i_sat).numobs = sum(src_obs_ind);
            
            % MJD of first and last observation:
            mjd_tmp                         = mjd(src_obs_ind);
            sources.s(i_sat).firstObsMjd    = mjd_tmp(1);
            sources.s(i_sat).lastObsMjd     = mjd_tmp(end);
            
    end % switch(obs_type_str)

end % for i_src = 1 : number_of_remaining_sources

% #### Add orbit data to sources.s ####
if ~isempty(sources.s)
    if ~isempty(sat_orbit_file_type)
        switch(sat_orbit_file_type)
            case 'sp3'
                % Read SP3 file and writ data to orbiot_data strucutre
                % - GPS time epochs in SP3 files are converted to UTC 
                [orbit_data] = read_sp3(sat_orbit_file_path, sat_orbit_file_name,{sources.s.name});
                [sources.s.orbit_file_type] = deal('sp3');
                
            case 'sat_ephem_trf'
                % Read sate emphemeris file with ITRF positions (and velocities) and writ data to orbiot_data strucutre
                [orbit_data] = read_sat_ephem_trf(sat_orbit_file_path, sat_orbit_file_name);
                [sources.s.orbit_file_type] = deal('sat_ephem_trf');

            otherwise
                error('Unknown orbit file type!');
        end
    else
        error('No orbit data file specified (VieVS GUI menu: Models/Space Crafts)!');
    end
    
    % Check, if all observations are coverd by the orbit data time series:
    if (max([sources.s.lastObsMjd] > orbit_data.mjd_end) || (min([sources.s.firstObsMjd]) < orbit_data.mjd_start))
        error('Satellite observation eopochs are not covered by the orbit data time series!');
    end
    
    % Write orbit data to the source structure:
    [sources] = orbit_data2sources(orbit_data, sources);
    
end



% #################################################################
% ##### SCAN                                                  #####
% #################################################################

% Preallocate structure:
scan(number_of_remaining_scans) = struct('mjd', [], 'stat', [], 'tim', [], 'nobs', [], 'space', [], 'obs', [], 'iso', []);
[scan.space] = deal(struct('source', zeros(3,3), 'xp', 0,'yp', 0, 'era', 0, 'xnut', 0, 'ynut', 0, 't2c', zeros(3,3)));
[scan.stat] = deal(struct('x', [], 'temp', [], 'pres', [], 'e', [], 'az', [], 'zd', [], 'zdry', [], 'cab', [], 'axkt', [], 'therm', [], 'pantd', [], 'trop', [], 'psd', []));
[scan.obs] = deal(struct('i1', [], 'i2', [], 'obs', [], 'sig', [], 'com', 0, 'q_code', [], 'q_code_ion', [], 'delion', [], 'sgdion', []));

year_scan   = year(scan2obs_ind);
mon_scan    = mon(scan2obs_ind);
day_scan    = day(scan2obs_ind);
hour_scan   = hour(scan2obs_ind);
min_scan    = minu(scan2obs_ind);
sec_scan    = sec(scan2obs_ind);
doy_scan    = doy(scan2obs_ind);

% loop over all scans:

for i_scan = 1 : number_of_remaining_scans
    
    % Obs. in scan:
    obs_in_scan_ind         = obs2scan_ind == i_scan;
    number_of_obs_in_scan   = sum(obs_in_scan_ind);
    obs_in_scan_id_list     = find(obs_in_scan_ind);
    % Stations in scan:
    stations_in_scan        = unique([stat_1_name_str(obs_in_scan_ind); stat_2_name_str(obs_in_scan_ind)]);
    number_of_stat_in_scan  = length(stations_in_scan);
    
    % #### Add data to scan structure: ####
    scan(i_scan).nobs       = number_of_obs_in_scan;
    scan(i_scan).obs_type   = obs_type_str{scan2obs_ind(i_scan)};
    scan(i_scan).src_name   = source_name_str{scan2obs_ind(i_scan)};
    % Get source ID:
    switch(scan(i_scan).obs_type)
        case 'q'
            scan(i_scan).iso        = find(strcmp(quasar_name_list, scan(i_scan).src_name));
        case 's'
            scan(i_scan).iso        = find(strcmp(sat_name_list, scan(i_scan).src_name));
    end
    scan(i_scan).mjd        = mjd_scan(i_scan);
    scan(i_scan).tim        = [year_scan(i_scan); mon_scan(i_scan); day_scan(i_scan); hour_scan(i_scan); min_scan(i_scan); sec_scan(i_scan); doy_scan(i_scan)];
    
    % #### Add stat sub-structure: ####
    for i_stat = 1 : number_of_stat_in_scan
        stat_name_8char_str = sprintf('%s%s', stations_in_scan{i_stat}, blanks(8 - length(stations_in_scan{i_stat})));
        % Antenna ID (refering to the antenna structure)
        antenna_ind = strcmp({antenna.name}, stat_name_8char_str);
        scan(i_scan).stat(antenna_ind) = struct('x', [0 0 0], 'temp', error_code_invalid_met_data, 'pres', error_code_invalid_met_data, 'e', error_code_invalid_met_data, 'az', 0, 'zd', [], 'zdry', 0, 'cab', 0, 'axkt', 0, 'therm', 0, 'pantd', [0 0 0], 'trop', [], 'psd', 0);
    end
    
    % #### Add obs sub-structure ####
    % Loop over observations in scan:
    for i_obs = 1 : number_of_obs_in_scan
        % Get observation ID:
        obs_id = obs_in_scan_id_list(i_obs);
        
        % Get antenna IDs:
        stat_name_8char_str = sprintf('%s%s', stat_1_name_str{obs_id}, blanks(8 - length(stat_1_name_str{obs_id})));
        antenna_1_id = find(strcmp({antenna.name}, stat_name_8char_str));
        stat_name_8char_str = sprintf('%s%s', stat_2_name_str{obs_id}, blanks(8 - length(stat_2_name_str{obs_id})));
        antenna_2_id = find(strcmp({antenna.name}, stat_name_8char_str));
        
        switch(vso_format)
            case 1 % Baseline delay only
                obs     = delay(obs_id) * sec_unit_conversion;
                switch(scan(i_scan).obs_type)
                    case 'q'
                        sig = formal_error_quasar_obs_sec;
                    case 's'
                        sig = formal_error_sat_obs_sec;
                end
                delion  = [];
                sgdion  = [];
                
            case 2 % Baseline delay + formal error
                obs     = delay(obs_id) * sec_unit_conversion;
                sig     = delay_sig(obs_id) * sec_unit_conversion;
                delion  = [];
                sgdion  = [];
                
            case 3 % Baseline delay + formal error  and Ion. delay + formal error
                sig     = sqrt((ion_sig(obs_id) * sec_unit_conversion)^2 + (delay_sig(obs_id)*1e-9)^2);
                obs     = delay(obs_id) * sec_unit_conversion - ion_del(obs_id) * 1e-9; 
                delion  = ion_del(obs_id);
                sgdion  = ion_sig(obs_id);
            
            case 4 % No observations at all
                obs     = [];
                
                switch(scan(i_scan).obs_type)
                    case 'q'
                        sig = formal_error_quasar_obs_sec;
                    case 's'
                        sig = formal_error_sat_obs_sec;
                end
                delion  = [];
                sgdion  = [];
                
            case 5 % Baseline delay + formal error; Ion. delay + formal error, cable corrections
                cor_cable_cal   = cab_stat2(obs_id) - cab_stat1(obs_id);    
                sig     = sqrt((ion_sig(obs_id) * sec_unit_conversion)^2 + (delay_sig(obs_id)*1e-9)^2);
                obs             = delay(obs_id) * sec_unit_conversion - ion_del(obs_id) * 1e-9 + cor_cable_cal * 1e-9; 
                delion          = ion_del(obs_id);
                sgdion          = ion_sig(obs_id);  % Has to be stored in [ns] in the scan struct!!!
                % Add data to "scan.stat":
                scan(i_scan).stat(antenna_1_id).cab = cab_stat1(obs_id);
                scan(i_scan).stat(antenna_2_id).cab = cab_stat2(obs_id);

            case 6 % Baseline delay + formal error; Ion. delay + formal error, cable corrections, met. data (T, p, e)
                cor_cable_cal   = cab_stat2(obs_id) - cab_stat1(obs_id);    
                sig     = sqrt((ion_sig(obs_id) * sec_unit_conversion)^2 + (delay_sig(obs_id)*1e-9)^2);
                obs             = delay(obs_id) * sec_unit_conversion - ion_del(obs_id) * 1e-9 + cor_cable_cal * 1e-9; 
                delion          = ion_del(obs_id);
                sgdion          = ion_sig(obs_id);  % Has to be stored in [ns] in the scan struct!!!
                % Add data to "scan.stat":
                scan(i_scan).stat(antenna_1_id).cab     = cab_stat1(obs_id);
                scan(i_scan).stat(antenna_1_id).temp    = t_stat1(obs_id);
                scan(i_scan).stat(antenna_1_id).pres    = p_stat1(obs_id);
                scan(i_scan).stat(antenna_1_id).e       = e_stat1(obs_id);
                scan(i_scan).stat(antenna_2_id).cab     = cab_stat2(obs_id);
                scan(i_scan).stat(antenna_2_id).temp    = t_stat2(obs_id);
                scan(i_scan).stat(antenna_2_id).pres    = p_stat2(obs_id);
                scan(i_scan).stat(antenna_2_id).e       = e_stat2(obs_id);
        end
        
        scan(i_scan).obs(i_obs) = struct('i1', antenna_1_id, 'i2', antenna_2_id, 'obs', obs, 'sig', sig, 'com', 0, 'q_code', scan_obs_q_code, 'q_code_ion', scan_obs_q_code_ion, 'delion', delion, 'sgdion', sgdion);
    end

end % for i_scan = 1 : number_of_remaining_scans

disp('finished');




