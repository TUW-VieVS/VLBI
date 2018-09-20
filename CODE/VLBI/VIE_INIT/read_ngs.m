% ************************************************************************
%   Description:
%   File to read the VLBI observations from a NGS-file. A structure array
%   scan containing all scans found in the NGS file will be created, along
%   with antenna and sources structure arrays for the stations and the
%   sources respectively. The positions etc. of the stations/sources are
%   read from the TRF/STA_INFO/CRF files. If a station/source is not found
%   in these files the information in the header of the NGS-file will be used.
%
%   Reference: Description of the NGS file format:
%       http://lacerta.gsfc.nasa.gov/mk5/help/dbngs_format.txt
%
%
%   Input:
%      ngsfil   (string)  The name of the NGS file, including the path
%      trffil   (string)  The name of the mat-file containing the station
%                          positions and velocities, including the path
%      crffil   (string)  The name of the mat-file containing the source
%                          positions, including the path
%      sta_excl (n*8 string) List of stations to be excluded
%      sour_excl (n*8 string) List of sources to be excluded
%      scan_excl (length n struct array) List of outliers. The structure
%                         array should contain the fields sta1, sta2 (name
%                         of the stations in observation) and mjd (time of
%                         observation in modified julian day.
%      minel    (double)  Minimum elevation angle.
%
%   Output:
%      antenna    Structure array containing information about the stations
%      sources    Structure array containing information about the sources
%      scan       Structure array containing the VLBI observations and
%                 other scan-specific parameters.
%
%   External calls:
%   	read_ngs_hdr.m,	modjuldat.m (in VIE_MOD directory),
%       cart2llh.m, elev.m, constants.m
%
%   Coded for VieVS:
%   May 2009 by Tobias Nilsson
%
%   Revision:
%   30 Jun 2009 by Tobias Nilsson: Now source and station coordinates are
%        taken from the NGS file header if not found in the CRF/TRF file.
%   13 Oct 2009 by Tobias Nilsson: Problem when there are asterisks in the
%        NGS file fixed. In these cases the values are set to zero and a
%        warning message is displayed.
%   05 Nov 2009 by Tobias Nilsson: Added the possibility to exclude
%        stations.
%   18 Nov 2009 by Tobias Nilsson: Possibility to exclude sources. Added
%        removal of outliers and bad quality code observations.
%   03 Dec 2009 by Tobias Nilsson: Possible to specify an elevation angle
%        cutoff.
%   17 Dec 2009 by Lucia Plank: replacement of cart2llh.m with xyz2ell;
%   12 Jan 2010 by Tobias Nilsson: Added possibility to include velocities
%        and breaks in a "special" TRF-file (.txt)
%   04 Feb 2010 by Tobias Nilsson: Some NGS files contains
%        #INF...-entries. Fixed the error caused by this when reading these
%        files.
%   05 Feb 2010 by Tobias Nilsson: Fixed a bug occuring sometimes when
%        there are *'s in the NGS-file.
%   01 Mar 2010 by Tobias Nilsson: Minor modifications, the program now
%     seems to work with Octave.
%   01 Mar 2010 by Lucia Plank: Minor changes concerning the scan array.
%   30 Mar 2010 by Lucia Plank: scan.tim --> civil date
%   22 Jun 2010 by Hana Spicakova: save start and end of the epoch of
%        station coordinates for global solution
%   27 Aug 2010 by Tobias Nilsson: Added the possibility to use a different
%        limit for the Q-codes than zero
%   31 Aug 2010 by Sigrid Böhm and Tobias Nilsson: Check for blanks in the
%        middle of station names and consider them as underscore in the
%        station search in the TRF
%   07 Sep 2010 by Tobias Nilsson: Fixed a bug with removal of
%        outliers
%   10 Nov 2010 by Tobias Nilsson: Fixed a bug occuring when there
%        are TRF coordinates but not for the epoch of the session.
%   19 Nov 2010 by Tobias Nilsson: (Re-)added the possibility to exclude
%       baselines again!!!!
%   25 Nov 2010 by Tobias Nilsson: added the fields firstObsMjd and
%       lastObsMjd to the antenna structure
%   04 May 2011 by Matthias Madzak: Iono delay is now only applied when in
%       vie_mod-Gui "from NGS" is ticked. If "from external file" is chosen,
%       the iono delay is applied in vie_mod.
%   05 May 2011 by Matthias Madzak: Speeding up the search if the scan
%       already exists in the scan struct, and preallocating.
%   18 Aug 2011 by Lucia Plank: sta0.trop = [];
%   10 Jan 2012 by Jing Sun: added ionospheric delay to scan
%   10 Apr 2012 by Matthias Madzak: Changed TRF a priori coords to be taken
%       from the superstations file instead of the different mat file (e.g.
%       VTRF2008.mat)
%   01 Jun 2012 by Lucia Plank: antenna(st(stnr)).opl=[] + in_trf for
%        manual trf file
%   06 Nov 2012 by Tobias Nilsson: include more information from superstation
%        in antenna, added gpt2
%   08 Nov 2012 by Hana Krasna: id of the antenna in superstation file saved in antenna struct
%   08 Nov 2012 by Hana Krasna: eccentricity read in from superstation file
%       (not '../TRF/ECC.mat anymore);
%   08 Nov 2012 by Hana Krasna: GPT2 function added for pressure and
%       temperature (constant value saved scan-wise)
%   13 Nov 2012 by Hana Krasna: antenna offset and mounting type taken from
%       superstatin file (antenna_info). No more from NGS header.
%   ----------  by Matthias Madzak: supersource file
%   27 March 2013 by Hana Krasna: changes according to supersource file
%   13 Jun 2013 by Hana Krasna: bug corrected for sources which have more
%       than 1 common name
%   18 Dec 2013 by Hana Krasna: ICRF designation added into source.mat
%   05 Feb 2014 by Lucia Plank: exclude scans in JET file
%   25 Jun 2014 by Hana Krasna: bug found and corrected for breaks in
%       station coordinates
%   03 Sep 2014 by Hana Krasna: bug found and corrected for excluding cable cal in
%       OPT file
%   15 Jan 2015 by Nastya&Matthias: Now more than one downtime per station
%       is possible (at stations to be excluded)
%   30 Jan 2015 by Lucia: only defining sources in NNR condition. THIS IS A
%       PRELIMINARY IMPLEMENTATION ONLY. Not included in the GUI yet. Search
%       for "LuciaNNR" in this file.
%   03 Mar 2015 by Caroline Schönberger: exclude observations between twin telescopes
%   20 Apr 2015 by Daniel Landskron: water vapor pressure e is now also
%       read out from the NGS file and saved in the "scan" structure
%   21 Jul 2015 by A.Girdiuk: finding of outliers and ngs-card coincidences
%   3  Nov 2015 by A.Girdiuk: duplicated fields (antenna_info and thermal) are in antenna-structure
%							  replaced as: thermal. To get access you need to use: antenna(:).thermal.
%							  Note: the field is generated with the trf variable
%   04 Dec 2015 by A. Hellerschmied: Minor bug-fix.
%   19 Jan 2016 by M. Madzak: Smarter antenna names ('_' instead of ' ' in
%                             beginning of name) are written to antenna
%                             struct
%   29 Apr 2016 by M. Madzak: Added post-seismic deformation to
%                             preallocation of superstation file
%   28 Jun 2016 by A. Hellerschmied: General revision
%      - Calculation of met. data (with GPT model) is moved to vie_init
%      - Waning appears, if the hunidity is not given as relative humidity as expected
%      - Lots of comments added to improve redability
%      - Minort bug-fixes and changes
%   06 Jul 2016 by A. Hellerschmied: - Corrected minor bugs
%                                    - Rigorous checks, if e is valid
%   07 Jul 2016 by D. Mayer: included NNR for only defining sources
%   15 Jul 2016 by D. Mayer: Added IVS source name
%   18 Jul 2016 by D. Mayer: include fistObsMjd and lastObsMjd into sources
%   22 Jul 2016 by D. Mayer: bug fix; excluded baselines were not read
%        correctly
%   09 Aug 2016 by A. Hellerschmied: Field scan(i_scan).obs_type added.
%   10 Sep 2016 by A. Girdiuk: bug-fix: sources can be excluded in selected time interval
%   26 Sep 2016 by H. Krasna: changes related to updates in supersource file
%   12 Oct 2016 by A. Girdiuk: bug-fix: sources can be excluded when
%                               NO time span is written down to OPT-file
%   25 Oct 2016 by A. Hellerschmied: - Changed warning msg. when den obs. delay is decreased/increased here by 1 sec
%                                    - Added additional status msg. to CW (total number of obs., number of valid/excluded obs.)
%                                    - Added field "scan.src_name"
%   23 Jan 2017 by M. Schartner: bug-fix: if you find multiple breaks in your TRF the first one is taken and a warning is displayed. 
%   24 Jan 2017 by D. Landskron: preallocation of GPT2 changed to GPT3
%   09 Feb 2017 by D. Landskron: Preallocation extended
%   14 Feb 2017 by M. Schartner: changes to improve speed
%   22 Feb 2017 by A. Hellerschmied: antenna.psd initialized
%   05 Jul 2018 by D. Landskron: vm1 renamed to vmf1 and VMF3 added to the troposphere models 
% ************************************************************************


% READ_NGS version 1.00827
function [antenna,sources,scan]=read_ngs(ngsfil,trffil,crffil,ini_opt, trf, crf)

% ##### Load constants #####
% constants; % ???

% ##### Options #####
url_vievswiki_create_superstation   = 'http://vievswiki.geo.tuwien.ac.at/doku.php?id=public:vievs_manual:data#create_a_superstation_file';
error_code_invalid_met_data         = -999; % Error corde for missing met. data in NGS file (numerical)

% ##### Check inoput options (ini_opt) #####
if ~isfield(ini_opt,'minel')
    ini_opt.minel=0;
end
if ~isfield(ini_opt,'Qlim')
    ini_opt.Qlim=0;
end
%if sta_excl, sour_excl, or scan_excl not specified, initialize them as empty matrices
if ~isfield(ini_opt,'sta_excl')
    ini_opt.sta_excl=[];
end
if ~isfield(ini_opt,'sour_excl')
    ini_opt.sour_excl=[];
end
if ~isfield(ini_opt,'scan_excl')
    ini_opt.scan_excl=[];
end
if ~isfield(ini_opt,'no_cab')
    ini_opt.no_cab=[];
end
if ~isfield(ini_opt,'bas_excl')
    ini_opt.bas_excl=[];
end

% ##### preallocate structures #####
scan        = struct('mjd', []);
scansource  = [];
sources     = struct('name', [], 'IERSname', [], 'IVSname', [], 'ICRFdes', [], 'ra2000', [], 'de2000', [], 'ra_sigma', [], 'de_sigma', [], 'corr', [], 'in_crf', [], 'flag_defining', []);



% ##### Initialisations #####
num_of_stat             =0;
num_of_sources          =0;
num_of_scans            =0;
num_of_obs              =0;

sta0.x          =[0 0 0];
sta0.temp       =0;
sta0.pres       =0;
sta0.e          =0;
sta0.az         =0;
sta0.zd         =[];
sta0.zdry       =0;
sta0.cab        =0;
sta0.axkt       =0;
sta0.therm      =0;
sta0.pantd      =[0 0 0];
sta0.trop       =[];
sta0.psd        =0;

space0.source   = zeros(3,3);
space0.xp       =0;
space0.yp       =0;
space0.era      =0;
space0.xnut     =0;
space0.ynut     =0;
space0.t2c      =zeros(3,3);

flag_ast_found  =0;

% Counters for warning msg. when the obs. delas are in/decreased by i sec:
count_obs_delay_decreased = 0;
count_obs_delay_increased = 0;
num_of_obs_in_ngs_file    = 0;

% ##### Open NGS file #####
fid_ngs = fopen(ngsfil,'r');

% #### Loop Init.: ####
% Get the first line:
% input_str = fgetl(fid_ngs);

% ##### Read the NGS file #####
% Loop over all lines of NGS file
if fid_ngs ~= -1
    wholeFile = textscan(fid_ngs,'%s','delimiter', '\n', 'whitespace', '');
    wholeFile = wholeFile{1};
    idx_line = 1;
    nlines = length(wholeFile);
else
    error(' couln''t find NGS-file %s',ngsfil)
end


while 1
    if length(wholeFile{idx_line})>= 80
        break;
    end
    idx_line = idx_line+1;
end

while (idx_line <= nlines)
    
    input_str = wholeFile{idx_line};
    idx_line = idx_line+1;
    
    % Get sequence and card number
    sequ_num        = input_str(71:78); % Sequence number (= observation number)
    ngs_card_num    = sscanf(input_str(79:80),'%d'); % NGS card number

    % #### Read station and source names of current observation/sequence ####
    sta_names(1,:)  = input_str(1:8);
    sta_names(2,:)  = input_str(11:18);
    source_name     = input_str(21:28);

    % check if there is blanks in the station name and replace them with "_";
    sta_names(1, sta_names(1, 1:max(find(sta_names(1,:) ~= ' '))) == ' ') = '_';
    sta_names(2, sta_names(2, 1:max(find(sta_names(2,:) ~= ' '))) == ' ') = '_';


    % #### Read epoch and convert to MJD ####
    tim = sscanf(input_str(30:60),'%f');
    [~,doy] = dday(tim(1), tim(2), tim(3), tim(4), tim(5));
    mjd = modjuldat(tim(1), tim(2), tim(3), tim(4), tim(5), tim(6));
    mjd_floor = floor(mjd);

    % Look init.:
    flag_next_sequ = 0;

    % ##############################################################################
    % ##### Loop over all other NGS cards of current sequence and read in data #####
    % ##############################################################################
    while (flag_next_sequ == 0) && (idx_line <= nlines)

        % Get next line:
        input_str = wholeFile{idx_line};
        idx_line = idx_line+1;
        
        new_sequ_num    = input_str(71:78);
        ngs_card_num    = input_str(79:80);

        if all(new_sequ_num == sequ_num) % Still the same sequence?

            % #### Check the current line ####
            % Check, if there are asterisks in the input line:
            id_ast = find(input_str=='*');
            if ~isempty(id_ast)
                if flag_ast_found == 0
                    flag_ast_found = 1;
                    fprintf('Asterisk(s) found in NGS file!!! Value(s) treated as zero!\n');
                end
                numinrow = 0;
                for i_tmp = 1:length(id_ast) - 1
                    if id_ast(i_tmp+1)==id_ast(i_tmp)+1;
                        input_str(id_ast(i_tmp))=' ';
                        numinrow=numinrow+1;
                    else
                        input_str(id_ast(i_tmp))='0';
                        numinrow=0;
                    end
                    if (numinrow==21)
                        input_str(id_ast(i_tmp-1))='0';
                    end
                end
                input_str(id_ast(i_tmp+1))='0';
            end

            % Check if there is a #INF.... in the line (occurs in some NGS
            % files). If so, replace with inf to avoid error later.
            id_inf = strfind(input_str,'#INF');
            if ~isempty(id_inf)
                for i_tmp=1:length(id_inf)
                    input_str(id_inf(i_tmp))=' ';
                    input_str(id_inf(i_tmp)+1:id_inf(i_tmp)+3)='inf';
                    b=id_inf(i_tmp)+4;
                    while input_str(b)~=' '
                        input_str(b)=' ';
                        b=b+1;
                    end
                end
            end

            % ##### Card 2: observed values #####
            if all(ngs_card_num == '02')
                trmp_input  = sscanf(input_str(1:62),'%f');
                delay       =trmp_input(1); % Observed delay (ns)
                sigdel      =trmp_input(2); % Formal error for the observed delay (ns)
                rate        =trmp_input(3); % Observed delay rate (ps/sec)
                sigrat      =trmp_input(4); % Formal error for the observed delay rate (ps/sec)
                q_code      =trmp_input(5); % Data quality flag (blank or 0 indicates good data)

                % Lucia:
                if delay < -4e8
%                             warning('PLUSSSSSSS: delay = delay + 1e9');
                    delay = delay + 1e9;
                    count_obs_delay_increased = count_obs_delay_increased + 1;
                elseif delay > 4e8
%                             warning('MINUSSSSSSSS: delay = delay - 1e9');
                    delay = delay - 1e9;
                    count_obs_delay_decreased = count_obs_delay_decreased + 1;
                end
                num_of_obs_in_ngs_file = num_of_obs_in_ngs_file + 1; 

                % ##### Card 5: cable cal & wvr corrections #####
            elseif all(ngs_card_num == '05')
                trmp_input  = sscanf(input_str(1:20),'%f');
                cab1        = trmp_input(1); % Cable calibration correction (one-way) for site 1 (ns)
                cab2        = trmp_input(2); % Cable calibration correction (one-way) for site 2 (ns)

                % Apply OPT file options:
                % Set cable cal. value =0, if defined in OPT file!
                for i_tmp = 1 : size(ini_opt.no_cab,1)
                    if strcmp(ini_opt.no_cab(i_tmp, :), sta_names(1,:))
                        cab1 = 0;
                    end
                    if strcmp(ini_opt.no_cab(i_tmp,:), sta_names(2,:))
                        cab2 = 0;
                    end
                end
                cor_cabel_cal = cab2 - cab1;

                % ##### Card 6: met data #####
            elseif all(ngs_card_num == '06')
                trmp_input  = sscanf(input_str(1:64),'%f');
                tdry1       = trmp_input(1); % Ambient atmospheric temperature at site 1 (deg. C)
                tdry2       = trmp_input(2); % Ambient atmospheric temperature at site 2 (deg. C)
                press1      = trmp_input(3); % Ambient atmospheric barometric pressure at site 1 (mb)
                press2      = trmp_input(4); % Ambient atmospheric barometric pressure at site 2 (mb)
                hum1        = trmp_input(5); % Ambient atmospheric humidity at site 1
                hum2        = trmp_input(6); % Ambient atmospheric humidity at site 2
                hum_par_1   = trmp_input(7); % Humidity parameter definition code for site 1
                hum_par_2   = trmp_input(8); % Humidity parameter definition code for site 2
                % 0=humidity parameter is relative humidity (%)
                % 1=humidity parameter is dew point (deg. C)
                % 2=humidity parameter is wet bulb temperature (deg. C)
                % Calculate relative humidity with the formula by Magnus:
                e1 = 6.1078 * exp((17.1 * tdry1) / (235 + tdry1)) * hum1/100;   % formula by Magnus * relative humidity
                e2 = 6.1078 * exp((17.1 * tdry2) / (235 + tdry2)) * hum2/100;   % formula by Magnus * relative humidity
                if (hum_par_1 ~= 0) || (hum_par_2 ~= 0)
                    warning('The unit of the hunidity parameter in the NGS file is NOT relative humidity (%) as expected!\n');
                end

                % ##### Card 8: ionospheric corrections #####
            elseif all(ngs_card_num == '08')
                % Check, if there is an "*":
                if input_str(51) == '*'
                    input_str(51:60) = '         0';
                end
                trmp_input  = sscanf(input_str(1:63),'%f');
                if strcmp(ini_opt.iono, 'ngs')
                    delion  = trmp_input(1); % Delay ionosphere correction (ns)
                else
                    delion  = 0; % E.g. use ion. correction from external source
                end
                sgdion      = trmp_input(2); % Delay ionosphere correction formal error (ns)
                ration      = trmp_input(3); % Delay rate ionosphere correction (ps/s)
                sgrion      = trmp_input(4); % Delay rate ionosphere correction formal error (ps/s)
                q_iono      = trmp_input(5); % Ionosphere error flag (0=ionosphere correction OK)
                coride      = -1 * delion; % change sign
                corira      = -1 * ration; % change sign
                corsgd      = sgdion;
                corsgr      = sgrion;
            end
        else
            flag_next_sequ = 1;
            idx_line = idx_line-1;
        end % if new_sequ_num == sequ_num
    end % while (flag_next_sequ == 0) && (~feof(fid_ngs))



    % #####################################################################
    % ##### Check if the observation is OK and apply OPT file options #####
    % #####################################################################

    % Init.:
    flag_obs_ok = 1;

    % #### Check the quality code ####
    if (q_code > ini_opt.Qlim) || (q_iono > ini_opt.Qlim)
        continue;
    end
    % ???
    % #### Exclude observations between twin telescopes ####
    if strcmp(sta_names(1,1:7), sta_names(2,1:7))
        fprintf(1, '=> Exclude observation between twin telescopes\n');
        continue;
    end
    % ???
    % #### Check if any of the stations are in the list of stations to be excluded
    if ~isempty(ini_opt.sta_excl)
        for i_tmp = 1 : size(ini_opt.sta_excl, 1)
            if strcmp(sta_names(1,:), ini_opt.sta_excl(i_tmp,:)) || strcmp(sta_names(2,:), ini_opt.sta_excl(i_tmp,:))
                if ini_opt.sta_excl_start(i_tmp) ~= 0
                    if (mjd > ini_opt.sta_excl_start(i_tmp)) && (mjd < ini_opt.sta_excl_end(i_tmp))
                        flag_obs_ok = 0;
                        break;
                    end
                else
                    flag_obs_ok = 0;
                    break;
                end
            end
        end
        if ~flag_obs_ok
            continue;
        end
    end

    % #### Check if the sources should be excluded ####
    if ~isempty(ini_opt.sour_excl)
        for i_tmp = 1 : size(ini_opt.sour_excl,1)
            if strcmp(source_name, ini_opt.sour_excl(i_tmp,:))
                if (mjd > ini_opt.sour_excl_start(i_tmp)) && (mjd < ini_opt.sour_excl_end(i_tmp))
                    flag_obs_ok = 0;
                    break;
                end
                if ~logical(ini_opt.sour_excl_start(i_tmp))
                    flag_obs_ok = 0;
                    break;
                end
            end
        end
        if ~flag_obs_ok
            continue;
        end
    end

    % #### Check if scan is in the list of outliers ####
    if ~isempty(ini_opt.scan_excl)
        %                idout=find(abs([ini_opt.scan_excl.mjd]-mjd)<10^-6);
        idout = find([ini_opt.scan_excl.mjd] == mjd);            %%%>> A. Girdiuk 2015-07-21
        for i_tmp = 1 : length(idout)
            if (strcmp(ini_opt.scan_excl(idout(i_tmp)).sta1,sta_names(1,:)) && ...
                    strcmp(ini_opt.scan_excl(idout(i_tmp)).sta2,sta_names(2,:))) || ...
                    (strcmp(ini_opt.scan_excl(idout(i_tmp)).sta1,sta_names(2,:)) && ...
                    strcmp(ini_opt.scan_excl(idout(i_tmp)).sta2,sta_names(1,:)))
                flag_obs_ok = 0;
                break;
            end
        end
        if ~flag_obs_ok
            continue;
        end
    end

    % #### Check if scan is in the list of removed baselines ####
    if ~isempty(ini_opt.bas_excl)
        for i_tmp = 1 : length(ini_opt.bas_excl)
            if (strcmp(sta_names(1,:), ini_opt.bas_excl(i_tmp).sta1) &&...
                    strcmp(sta_names(2,:), ini_opt.bas_excl(i_tmp).sta2)) ||...
                    (strcmp(sta_names(1,:), ini_opt.bas_excl(i_tmp).sta2) &&...
                    strcmp(sta_names(2,:), ini_opt.bas_excl(i_tmp).sta1))
                flag_obs_ok = 0;
                break;
            end
        end
        if ~flag_obs_ok
            continue;
        end
    end

    % #### Check for bad jet angles ####
    if ~isempty(ini_opt.scan_jet)
        idout=find(abs([ini_opt.scan_jet.mjd]-mjd)<10^-6);
        for i_tmp=1:length(idout)
            if (strcmp(ini_opt.scan_jet(idout(i_tmp)).sta1,sta_names(1,:)) && ...
                    strcmp(ini_opt.scan_jet(idout(i_tmp)).sta2,sta_names(2,:))) || ...
                    (strcmp(ini_opt.scan_jet(idout(i_tmp)).sta1,sta_names(2,:)) && ...
                    strcmp(ini_opt.scan_jet(idout(i_tmp)).sta2,sta_names(1,:)))
                flag_obs_ok = 0;
                break;
            end
        end
        if ~flag_obs_ok
            continue;
        end
    end



    % ########################################################################################
    % ##### Update the antenna-, scan-, and sources-structure, if the curretn scan is OK #####
    % ########################################################################################

    % Look for the stations in the antenna struct. array. If not found, add it.
    if flag_obs_ok == 1


        % ##### Get station IDs #####
        % Loop init.:
        stat_id_vec = zeros(2,1); % Station IDF vector of curretn baseline

        % Loop over two station of baseline
        for i_stat = 1 : 2
            %                     stat_id_vec(i_stat) = 0;
            if num_of_stat>0
                foundID = find(strcmp({sta_names(i_stat,:)}, {antenna.name}));
                if ~isempty(foundID)
                    stat_id_vec(i_stat) = foundID;
                end
            end
            % #################################################################
            % ##### ANTENNA                                               #####
            % #################################################################

            % If the station ID is not already in the antenna structure ==> add it!
            if stat_id_vec(i_stat) == 0 % If station not in antenna, add it

                % Increase of the station index:
                num_of_stat = num_of_stat + 1;
                stat_id_vec(i_stat)=num_of_stat; % New station ID




                % ##### Get data from superstation file #####
                % The data from "Manual TRF files" is also available in the trf structure (as the data from the superstation file.), TRF-name: "manualTrf":
                % To disinguish between different TRF realisations take the trffil variable:
                %  - trffil{1} ... path + name of TRF file (superstaion.mat or "manual" TRF file (.txt)).
                %  - trffil{2} ... Name of TRF realisation (e.g. vieTrf, etc.); "manualTrf", if a "manual" TRF file (.txt) was selected in the GUI.


                trf_id = find(strcmpi({trf.name}, sta_names(i_stat,:)));

                % #### Check, if there is an entry for the current station in the superstation file. If not => Error Msg. and abort! ####
                if isempty(trf_id)
                    error('Station %s not found in the superstation file. Add this station to the superstation file by following the steps described at %s\n', sta_names(i_stat,:), url_vievswiki_create_superstation);
                end

                % Save the ID of the antenna in the superstation file (needed for finding of the corrections in Vie_MOD)
                antenna(stat_id_vec(i_stat)).IDsuper = trf_id;

                % #### Get break (coordinate epoch) ####

                % ### Check, if there are coordinates for the chosen TRF available and get the "break_id" ###
                if ~isempty(trf(trf_id).(trffil{2}))

                    % if a break is empty - put 99999 to start and 0 to end
                    if ~isfield(trf(trf_id).(trffil{2}).break(1),'start') % 25/06/2014
                        trf(trf_id).(trffil{2}).break.start=[];
                        trf(trf_id).(trffil{2}).break.end=[];
                    end
                    emptyStartLog = cellfun(@isempty,{trf(trf_id).(trffil{2}).break.start});
                    emptyEndLog   = cellfun(@isempty,{trf(trf_id).(trffil{2}).break.end});
                    if sum(emptyStartLog) > 0
                        try
                            trf(trf_id).(trffil{2}).break(emptyStartLog).start = zeros(sum(emptyStartLog), 1); % 25/06/2014
                        catch
                            keyboard;
                        end
                    end
                    if sum(emptyEndLog) > 0
                        trf(trf_id).(trffil{2}).break(emptyEndLog).end = repmat(99999, sum(emptyEndLog), 1); % 25/06/2014
                    end

                    break_id = find(mjd >= [trf(trf_id).(trffil{2}).break.start] & mjd <= [trf(trf_id).(trffil{2}).break.end]);
                    if length(break_id) > 1
                        break_id = break_id(1);
                        warning(['Multiple breaks found in TRF data for station: ' trf(trf_id).name]);
                    end

                else % Not found...
                    fprintf('Station not found in %s!\n', trffil{2});
                    break_id = [];
                end


                % #### If the "break_id" was not found in selected TRF ####
                % => Get coordinates from VieVS TRF (= backup TRF!)
                if isempty(break_id)

                    fprintf('No %s coordinates for %s in %s ... get vievsTrf coordinates\n', trffil{2},sta_names(i_stat,:),trffil{1});

                    % if there is no start break - only approx coords (e.g. station "VLA     ")
                    if ~isempty(trf(trf_id).vievsTrf)
                        if ~isfield(trf(trf_id).vievsTrf.break, 'start')
                            break_id = 1;
                        else
                            break_id = find(mjd >= [trf(trf_id).vievsTrf.break.start] & mjd <= [trf(trf_id).vievsTrf.break.end]);
                            if length(break_id) > 1
                                break_id = break_id(1);
                                warning(['Multiple breaks found in TRF data for station: ' trf(trf_id).name]);
                            end
                        end
                    else
                        error('Station %s has no vievsTRF coordinates in the superstation file. Add this station to the superstation file (to vievsTRF.txt) by following the steps described at %s\n', sta_names(i_stat,:), url_vievswiki_create_superstation);
                    end

                    % ### Check, if VieVS TRF coordinates were found: ###
                    if isempty(break_id)
                        error('Station %s not found in the superstation file (vievsTRF). Add this station to the superstation file by following the steps described at %s\n', trf(trf_id).name, url_vievswiki_create_superstation);
                    end

                    % ### Get the break sub-structure from the trf structure where coords should be taken ###
                    curBreak_substruct = trf(trf_id).vievsTrf.break(break_id);

                    % ### Define, that the station (coordiantes from VieVS TRF, as BACKUP!) is NO datum station!
                    antenna(stat_id_vec(i_stat)).in_trf = 0;

                else % Station coordinates were taken from selected TRF file


                    % ### Get the break sub-structure from the trf structure where coords should be taken ###
                    curBreak_substruct = trf(trf_id).(trffil{2}).break(break_id);




                    % ### Define if the current station should be a datum station ###
                    % => since station is found in chosen trf, it would be 1 by default, but it could be 0, if it is set to 0 in "manual" trf file
                    % => If the "indatum" flag has a "NaN" value => in_trf = 1
                    if isfield(curBreak_substruct, 'indatum') && ~isnan(curBreak_substruct.indatum)
                        antenna(stat_id_vec(i_stat)).in_trf = curBreak_substruct.indatum;
                    else
                        antenna(stat_id_vec(i_stat)).in_trf = 1;
                    end

                end % if isempty(break_id)
                %                         end



                % ##### Add data to antenna structure #####

                antenna(stat_id_vec(i_stat)).name           = trf(trf_id).name;
                antenna(stat_id_vec(i_stat)).x              = curBreak_substruct.x;
                antenna(stat_id_vec(i_stat)).y              = curBreak_substruct.y;
                antenna(stat_id_vec(i_stat)).z              = curBreak_substruct.z;
                
                % ### Transform Cartesian coordinates X,Y,Z to ellipsoidal to improve speed
                [phi(stat_id_vec(i_stat)),lam(stat_id_vec(i_stat))] =xyz2ell([antenna(stat_id_vec(i_stat)).x, antenna(stat_id_vec(i_stat)).y, antenna(stat_id_vec(i_stat)).z]);

                antenna(stat_id_vec(i_stat)).firstObsMjd    = mjd;
                % Check, if velocity information is available:
                if isfield(curBreak_substruct, 'vx')
                    antenna(stat_id_vec(i_stat)).vx         = curBreak_substruct.vx;
                    antenna(stat_id_vec(i_stat)).vy         = curBreak_substruct.vy;
                    antenna(stat_id_vec(i_stat)).vz         = curBreak_substruct.vz;
                    antenna(stat_id_vec(i_stat)).epoch      = curBreak_substruct.epoch;
                    antenna(stat_id_vec(i_stat)).start      = curBreak_substruct.start;
                    antenna(stat_id_vec(i_stat)).end        = curBreak_substruct.end;
                else % if coords are taken from blokq.dat: not even that info is given...
                    antenna(stat_id_vec(i_stat)).vx         = 0;
                    antenna(stat_id_vec(i_stat)).vy         = 0;
                    antenna(stat_id_vec(i_stat)).vz         = 0;
                    antenna(stat_id_vec(i_stat)).epoch      = 0;
                    antenna(stat_id_vec(i_stat)).start      = 0;
                    antenna(stat_id_vec(i_stat)).end        = 99999;
                end
                % check if error measures are available:
                if isfield(curBreak_substruct,'x_sigma')
                    antenna(stat_id_vec(i_stat)).x_sigma    = curBreak_substruct.x_sigma;
                    antenna(stat_id_vec(i_stat)).y_sigma    = curBreak_substruct.y_sigma;
                    antenna(stat_id_vec(i_stat)).z_sigma    = curBreak_substruct.z_sigma;
                    antenna(stat_id_vec(i_stat)).vx_sigma   = curBreak_substruct.vx_sigma;
                    antenna(stat_id_vec(i_stat)).vy_sigma   = curBreak_substruct.vy_sigma;
                    antenna(stat_id_vec(i_stat)).vz_sigma   = curBreak_substruct.vz_sigma;
                else
                    antenna(stat_id_vec(i_stat)).x_sigma    = [];
                    antenna(stat_id_vec(i_stat)).y_sigma    = [];
                    antenna(stat_id_vec(i_stat)).z_sigma    = [];
                    antenna(stat_id_vec(i_stat)).vx_sigma   = [];
                    antenna(stat_id_vec(i_stat)).vy_sigma   = [];
                    antenna(stat_id_vec(i_stat)).vz_sigma   = [];
                end

                antenna(stat_id_vec(i_stat)).thermal        = trf(trf_id).antenna_info;	%%% Girdiuk 3 Nov 2015
                antenna(stat_id_vec(i_stat)).comments       = trf(trf_id).comments;
                antenna(stat_id_vec(i_stat)).domes          = trf(trf_id).domes;
                antenna(stat_id_vec(i_stat)).code           = trf(trf_id).code;


                % Eccentricity from superstation file
                if isempty(trf(trf_id).ecc)
                    antenna(stat_id_vec(i_stat)).ecc        = [0 0 0];
                    antenna(stat_id_vec(i_stat)).ecctype    = 'NEU';
                else
                    for i_tmp = 1 : length(trf(trf_id).ecc.break)
                        tecc        = trf(trf_id).ecc.break(i_tmp).starting;
                        y           = str2double(tecc(1:4));
                        m           = str2double(tecc(6:7));
                        d           = str2double(tecc(9:10));
                        h           = str2double(tecc(12:13));
                        mi          = str2double(tecc(15:16));
                        s           = 0;
                        mjd_start   = modjuldat(y,m,d,h,mi,s);
                        tecc = trf(trf_id).ecc.break(i_tmp).ending;
                        y           = str2double(tecc(1:4));
                        m           = str2double(tecc(6:7));
                        d           = str2double(tecc(9:10));
                        h           = str2double(tecc(12:13));
                        mi          = str2double(tecc(15:16));
                        s           = 0;
                        mjd_end     = modjuldat(y,m,d,h,mi,s);
                        if (mjd_floor >= mjd_start) && (mjd_floor <= mjd_end)
                            break;
                        end
                    end
                    antenna(stat_id_vec(i_stat)).ecc        = [trf(trf_id).ecc.break(i_tmp).FCE, trf(trf_id).ecc.break(i_tmp).SCE, trf(trf_id).ecc.break(i_tmp).TCE];
                    antenna(stat_id_vec(i_stat)).ecctype    = trf(trf_id).ecc.break(i_tmp).type_e;
                end

                % mounting type from superstation file
                antenna(stat_id_vec(i_stat)).axtyp      = '';
                if ~isempty(trf(trf_id).antenna_info) % if there was information in antenna-info.txt
                    antenna(stat_id_vec(i_stat)).axtyp  = trf(trf_id).antenna_info.mount(4:end);
                end
                if strcmp(antenna(stat_id_vec(i_stat)).axtyp, 'XYNO')
                    antenna(stat_id_vec(i_stat)).axtyp  = 'X-Y1';
                end

                % axis offset from superstation file
                antenna(stat_id_vec(i_stat)).offs = 0;
                if ~isempty(trf(trf_id).antenna_info) % if there was information in antenna-info.txt
                    antenna(stat_id_vec(i_stat)).offs = trf(trf_id).antenna_info.axis_offset; %m
                end


                % Init. flags:
                antenna(stat_id_vec(i_stat)).gpt3pres   = 0;
                antenna(stat_id_vec(i_stat)).gpt3temp   = 0;
                antenna(stat_id_vec(i_stat)).gpt3e      = 0;
                antenna(stat_id_vec(i_stat)).noGrad     = 0;
                
                % Init. emptx fields => They will be filled in vie_mod (calc_met_data.m)
                antenna(stat_id_vec(i_stat)).gpt3.p      = [];
                antenna(stat_id_vec(i_stat)).gpt3.T      = [];
                antenna(stat_id_vec(i_stat)).gpt3.dT     = [];
                antenna(stat_id_vec(i_stat)).gpt3.Tm     = [];
                antenna(stat_id_vec(i_stat)).gpt3.e      = [];
                antenna(stat_id_vec(i_stat)).gpt3.ah     = [];
                antenna(stat_id_vec(i_stat)).gpt3.aw     = [];
                antenna(stat_id_vec(i_stat)).gpt3.lambda = [];
                antenna(stat_id_vec(i_stat)).gpt3.undu   = [];
                antenna(stat_id_vec(i_stat)).gpt3.Gn_h   = [];
                antenna(stat_id_vec(i_stat)).gpt3.Ge_h   = [];
                antenna(stat_id_vec(i_stat)).gpt3.Gn_w   = [];
                antenna(stat_id_vec(i_stat)).gpt3.Ge_w   = [];

                antenna(stat_id_vec(i_stat)).cto        = [];
                antenna(stat_id_vec(i_stat)).cta        = [];
                antenna(stat_id_vec(i_stat)).cnta_dx    = [];
                antenna(stat_id_vec(i_stat)).vmf3       = [];
                antenna(stat_id_vec(i_stat)).vmf1       = [];
                antenna(stat_id_vec(i_stat)).opl        = [];

                antenna(stat_id_vec(i_stat)).numobs     = 0;
                
                antenna(stat_id_vec(i_stat)).lastObsMjd = [];
                
                antenna(stat_id_vec(i_stat)).psd        = [];

            end %
        end


        % #################################################################
        % ##### SOURCES                                               #####
        % #################################################################

        % Find source in sources structure array. If not found, add it!
        foundSource = strcmp(source_name, {sources.name});
        if sum(foundSource) == 0

            % find index of current source in crf (crf ~=superstation file)
            % search at first within IERS names
            curSourceInCrf = strcmp(source_name, {crf.IERSname});

            % if not found within IERS names, search within 
            % IVS names
            actIVSName = 0;
            if sum(curSourceInCrf) == 0
%                 curSourceInCrf = strcmp(strtrim(cellstr(char(crf.IVSname))), strtrim(source_name));
                curSourceInCrf = strcmp(strtrim({crf.IVSname}),strtrim({source_name}));
                if sum(curSourceInCrf) > 0
                    actIVSName = 1; % 1 if the name of the actual source is the IVS one
                end
            end

            % number of sources we already have put to sources
            if length(sources) > 0 % in case the first source has zero observations,
                % the field structure is deleted and the source
                % field has dimention of 1x0!!!!
                if isempty(sources(1).name)
                    nOldSources = 0;
                else
                    nOldSources = length(sources);
                end
            else
                nOldSources = 0;
            end

            % get index of new source in sources struct
            indOfNewSourceInSources = nOldSources + 1;

            % if not existing in supersource file (should not
            % happen)
            if sum(curSourceInCrf) == 0
                fprintf('ERROR: Source %s does not exist in supersource file (%s) - add there!\n', source_name,crffil{1});
            end

            % check if chosen-catalog coordinates exist
            if isempty(crf(curSourceInCrf).(crffil{2}))
                fprintf('No %s coordinates for %s in %s ... get vievsCrf coordinates\n',crffil{2},source_name,crffil{1});
                crfCatalogToTake = 'vievsCrf';
            else
                crfCatalogToTake = crffil{2};
            end

            % add information
            sources(indOfNewSourceInSources).name = source_name;
            if actIVSName == 1
                sources(indOfNewSourceInSources).IERSname = crf(curSourceInCrf).IERSname;
            else
                sources(indOfNewSourceInSources).IERSname = source_name;
            end
            sources(indOfNewSourceInSources).IVSname = crf(curSourceInCrf).IVSname;

            % Check if source.name from NGS file is equal to
            % source.IVSname from translation table
            if ~strcmp(sources(indOfNewSourceInSources).name,sources(indOfNewSourceInSources).IVSname)
                fprintf('WARNING: Source name from NGS file %s does not correspond to the IVS name %s from translation table!\n', sources(indOfNewSourceInSources).name,sources(indOfNewSourceInSources).IVSname);
            end


            if ~isempty(crf(curSourceInCrf).designation)
                sources(indOfNewSourceInSources).ICRFdes = crf(curSourceInCrf).designation(6:end);
            else
                sources(indOfNewSourceInSources).ICRFdes = 'J               ';
            end
            sources(indOfNewSourceInSources).ra2000 = crf(curSourceInCrf).(crfCatalogToTake).ra;
            sources(indOfNewSourceInSources).de2000 = crf(curSourceInCrf).(crfCatalogToTake).de;

            % add sigmas if exist
            if isfield(crf(curSourceInCrf).(crffil{2}), 'ra_sigma')
                sources(indOfNewSourceInSources).ra_sigma = crf(curSourceInCrf).(crfCatalogToTake).ra_sigma;
                sources(indOfNewSourceInSources).de_sigma = crf(curSourceInCrf).(crfCatalogToTake).de_sigma;
            else
                sources(indOfNewSourceInSources).ra_sigma = [];
                sources(indOfNewSourceInSources).de_sigma = [];
            end

            sources(indOfNewSourceInSources).firstObsMjd    = mjd;

            % add correlation if exist
            if isfield(crf(curSourceInCrf).(crfCatalogToTake), 'corr')
                sources(indOfNewSourceInSources).corr = crf(curSourceInCrf).(crfCatalogToTake).corr;
            else
                sources(indOfNewSourceInSources).corr = [];
            end

            % in_crf (important for estimating sources)
            if strcmp(crffil{2}, 'vievsCrf') % if the user has chosen vievsCrf
                sources(indOfNewSourceInSources).in_crf = crf(curSourceInCrf).(crffil{2}).in_crf; % take the in_crf from the "source.cat" "textfile"
            else
                % if we take backup coords for the current source
                if strcmp(crfCatalogToTake, 'vievsCrf')
                    sources(indOfNewSourceInSources).in_crf = 0;
                else
                    sources(indOfNewSourceInSources).in_crf = 1;
                end
            end

            %LuciaNNR + (uncomment this and also change
            %LSMopt_INTERNAL)
            % add information about defining
%                     if isfield(crf(curSourceInCrf).(crfCatalogToTake),'defining')
%                         sources(indOfNewSourceInSources).defining=crf(curSourceInCrf).(crfCatalogToTake).defining;
%                     else
%                         sources(indOfNewSourceInSources).defining=[];
%                     end
            %LuciaNNR -

            %David - always add information about defining sources
            %from the ICRF2
            if ~isempty(crf(curSourceInCrf).icrf2)
                sources(indOfNewSourceInSources).flag_defining=crf(curSourceInCrf).icrf2.defining;
            else
                sources(indOfNewSourceInSources).flag_defining=0;
            end



            % set number of observations to zero
            sources(indOfNewSourceInSources).numobs = 0;

        else
            % the source of this obs was already put to sources struct
            indOfNewSourceInSources = find(foundSource);
        end




        % #################################################################
        % ##### SCAN                                                  #####
        % #################################################################
        %Look for the scan in the scan structure array. If not found, add it

        % Calc. elevation of observed source:
        el1 = elev(mjd, [antenna(stat_id_vec(1)).x, antenna(stat_id_vec(1)).y, antenna(stat_id_vec(1)).z], sources(indOfNewSourceInSources).de2000, sources(indOfNewSourceInSources).ra2000,phi(stat_id_vec(1)),lam(stat_id_vec(1)));
        el2 = elev(mjd, [antenna(stat_id_vec(2)).x, antenna(stat_id_vec(2)).y, antenna(stat_id_vec(2)).z], sources(indOfNewSourceInSources).de2000, sources(indOfNewSourceInSources).ra2000,phi(stat_id_vec(2)),lam(stat_id_vec(2)));

        % Check min elevation:
        if (el1 >= ini_opt.minel) && (el2 >= ini_opt.minel)

            % Init.:
            i_scan = 0; % Scan index

            % Check, if we have that scan date already and at same epoch the current source was observed:
            id_sc = find((scansource == indOfNewSourceInSources) & ([scan.mjd] == mjd));
            if ~isempty(id_sc)
                i_scan = id_sc;
            end

            % if we have a new scan:
            if i_scan == 0
                num_of_scans = num_of_scans + 1;
                scansource(num_of_scans) = indOfNewSourceInSources;
                i_scan = num_of_scans;

                % Init. scan structure:
                scan(i_scan).stat(stat_id_vec(1))   = sta0;
                scan(i_scan).stat(stat_id_vec(2))   = sta0;
                scan(num_of_scans).mjd              = mjd;
                scan(num_of_scans).tim              = [tim; doy];
                scan(num_of_scans).nobs             = 0;
                scan(num_of_scans).space            = space0;
                scan(num_of_scans).src_name         = sources(indOfNewSourceInSources).name;
%                         fprintf('%3d : %s\n', num_of_scans, scan(num_of_scans).src_name);
            end

            % Write observation data to scan:
            num_of_obs          = num_of_obs + 1;
            delay               = (delay + cor_cabel_cal + coride) / 1000;
            rate                = rate + corira;
            sigdel              = sqrt(sigdel^2 + corsgd^2);
            sigrat              = sqrt(sigrat^2+corsgr^2);
            scan(i_scan).nobs   = scan(i_scan).nobs + 1;

            scan(i_scan).obs(scan(i_scan).nobs).i1          = stat_id_vec(1);
            scan(i_scan).obs(scan(i_scan).nobs).i2          = stat_id_vec(2);
            scan(i_scan).obs(scan(i_scan).nobs).obs         = delay * 1e-6;
            scan(i_scan).obs(scan(i_scan).nobs).sig         = sigdel * 1e-9;
            scan(i_scan).obs(scan(i_scan).nobs).com         = 0;
            scan(i_scan).iso                                = indOfNewSourceInSources;
            scan(i_scan).obs(scan(i_scan).nobs).q_code      = q_code;
            scan(i_scan).obs(scan(i_scan).nobs).q_code_ion  = q_iono;

            scan(i_scan).obs(scan(i_scan).nobs).delion      = delion;   % Jing SUN, Jan 10, 2012
            scan(i_scan).obs(scan(i_scan).nobs).sgdion      = sgdion;   % Jing SUN, Jan 10, 2012

            scan(i_scan).stat(stat_id_vec(1))               = sta0;
            scan(i_scan).stat(stat_id_vec(2))               = sta0;

            scan(i_scan).stat(stat_id_vec(1)).cab           = cab1;
            scan(i_scan).stat(stat_id_vec(2)).cab           = cab2;

            scan(i_scan).space                              = space0;

            scan(i_scan).obs_type                           = 'q'; % Specify observation type (q = quasar)


            % Met. data:
            %  => Save valid measured data
            %  => Init. empty fields, if no valid data is available

            % Check, if the temp. value at station 1 from NGS file is valid:
            if tdry1 > -99
                scan(i_scan).stat(stat_id_vec(1)).temp = tdry1;
            else
                scan(i_scan).stat(stat_id_vec(1)).temp = error_code_invalid_met_data;
            end

            % Check, if the temp. value at station 2 from NGS file is valid:
            if tdry2 > -99
                scan(i_scan).stat(stat_id_vec(2)).temp = tdry2;
            else
                scan(i_scan).stat(stat_id_vec(2)).temp = error_code_invalid_met_data;
            end

            % Check, if the pres. value at station 1 from NGS file is valid:
            if press1 > 0
                scan(i_scan).stat(stat_id_vec(1)).pres = press1;
            else
                scan(i_scan).stat(stat_id_vec(1)).pres = error_code_invalid_met_data;
            end

            % Check, if the pres. value at station 2 from NGS file is valid:
            if press2 > 0
                scan(i_scan).stat(stat_id_vec(2)).pres = press2;
            else
                scan(i_scan).stat(stat_id_vec(2)).pres = error_code_invalid_met_data;
            end

            % Check, if the water vapor pressure at station 1 from NGS file is valid:
            % => Also check temp. here, because "e" is calculated from "tdry" and "hum".
            if (e1 > 0) && (tdry1 > -99) && (hum1 > 0)
                scan(i_scan).stat(stat_id_vec(1)).e = e1;
            else
                scan(i_scan).stat(stat_id_vec(1)).e = error_code_invalid_met_data;
            end

            % Check, if the water vapor pressure at station 2 from NGS file is valid:
            % => Also check temp. here, because "e" is calculated from "tdry" and "hum".
            if (e2 > 0) && (tdry2 > -99) && (hum2 > 0)
                scan(i_scan).stat(stat_id_vec(2)).e = e2;
            else
                scan(i_scan).stat(stat_id_vec(2)).e = error_code_invalid_met_data;
            end


            % #### Update antenna and sources structures: ####
            antenna(stat_id_vec(1)).numobs          = antenna(stat_id_vec(1)).numobs + 1;
            antenna(stat_id_vec(2)).numobs          = antenna(stat_id_vec(2)).numobs + 1;
            antenna(stat_id_vec(1)).lastObsMjd      = mjd;
            antenna(stat_id_vec(2)).lastObsMjd      = mjd;

            sources(indOfNewSourceInSources).numobs = sources(indOfNewSourceInSources).numobs + 1;
            sources(indOfNewSourceInSources).lastObsMjd = mjd;
        else % Invalid scan => Remove it!

            if antenna(stat_id_vec(2)).numobs == 0
                antenna(stat_id_vec(2))     = [];
                num_of_stat                 = num_of_stat - 1;
            end
            if antenna(stat_id_vec(1)).numobs == 0
                antenna(stat_id_vec(1))     = [];
                num_of_stat                 = num_of_stat - 1;
            end
            if sources(indOfNewSourceInSources).numobs == 0
                sources(indOfNewSourceInSources)    = [];
                num_of_sources                      = num_of_sources - 1;
            end

        end % if (el1 >= ini_opt.minel) && (el2 >= ini_opt.minel)
    end % if ngs_card_num == 1
end


% ##### Write warnings to CW: #####
if count_obs_delay_decreased > 0
    fprintf('WARNING: %d observed delay values from the NGS file were decreased by 1 sec!\n', count_obs_delay_decreased);
end

if count_obs_delay_increased > 0
    fprintf('WARNING: %d observed delay values from the NGS file were increased by 1 sec!\n', count_obs_delay_increased);
end

if count_obs_delay_increased || count_obs_delay_decreased
    
end

fprintf('Total number of observations in NGS file: %d\n', num_of_obs_in_ngs_file);
fprintf('Number of valid observations:             %d\n', num_of_obs);
fprintf('Number of excluded observations:          %d\n', num_of_obs_in_ngs_file - num_of_obs);


fclose(fid_ngs);

