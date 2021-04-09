% ************************************************************************
%   Description:
%   Function to estimate geodetic parameters of VLBI with least squares adjustment
%
%   Reference:
%   [1] J. Boehm, H. Spicakova, L. Plank, K. Teke, A. Pany, J. Wresnik, S.
%       Englich, H. Schuh, T. Hobiger, R. Ichikawa, Y. Koyama, T. Gotoh,
%       T. Otsubo, T. Kubooka (2009), ?Plans for the Vienna VLBI Software VieVS?,
%       Proceedings of the 19th European VLBI for Geodesy and Astrometry
%       Working Meeting, edited by G. Bourda, P. Charlot, A. Collioud,
%       Universite Bordeaux1-CNRS, 23-28 March 2009, Bordeaux, France,
%       pp. 161-164 (Oral presentation and paper).
%
%   [2]K. Teke, J. Boehm, H. Spicakova, A. Pany, and H. Schuh, (2009),
%      ?Sub-daily Parameter Estimation for the Continuous VLBI Campaign,
%      CONT05?, European Geosciences Union General Assembly, 19-24 April 2009,
%      Vienna, Austria (Poster).
%
%   [3]K. Teke, J. Boehm, E. Tanir, H. Schuh, (2009), ?Piecewise Linear
%      Offsets for VLBI Parameter Estimation?, Proceedings of the 19th
%      European VLBI for Geodesy and Astrometry Working Meeting,
%      edited by G. Bourda, P. Charlot, A. Collioud, Universite Bordeaux1-CNRS,
%      23-28 March 2009, Bordeaux, France, pp. 63-67 (Poster and paper).
%
%
%   Input:
%      'VieVS/DATA/LEVEL1/parameter.session_name_antenna.mat'   structure array   (for info. /DOC/structures.xls)
%      'VieVS/DATA/LEVEL1/parameter.session_name_parameter.mat'   structure array   (for info. /DOC/structures.xls)
%      'VieVS/DATA/LEVEL1/parameter.session_name_scan.mat'   structure array   (for info. /DOC/structures.xls)
%      'VieVS/DATA/LEVEL1/parameter.session_name_sources.mat'   structure array   (for info. /DOC/structures.xls)
%
%   Output:
%      'VieVS/DATA/LEVEL3/opt_',parameter.session_name,'.mat'    structure array  (for info. ../DOC/structures.xls)
%      'VieVS/DATA/LEVEL3/x_',parameter.session_name,'.mat'      structure array  (for info. ../DOC/structures.xls)
%      'VieVS/DATA/LEVEL3/atpa_',parameter.session_name,'.mat'   structure array  (for info. ../DOC/structures.xls)
%      'VieVS/DATA/LEVEL3/atpl_',parameter.session_name,'.mat'   structure array  (for info. ../DOC/structures.xls)
%      'VieVS/DATA/LEVEL2/' parameter.session_name '_glob.mat'   structure array  (for info. ../DOC/structures.xls)
%
%   External calls:
%   a_love.m, a_shida.m, a_source.m, a_xyz.m, ahp_dut1.m, ahp_nutdx.m,
%   ahp_nutdy.m, ahp_xpol.m, ahp_ypol.m, apw_egr.m, apw_ngr.m, apw_source.m,
%   apw_xyz.m, apw_zwd.m, apwq_clk.m, delmodel.m, delparam.m, delref.m,
%   delsource.m, hpoc_sources.m,
%   hpoc_xyz.m, lsmopt.m, helmert.m, reduce_oc.m,
%   sourcewisepar.m, splitx.m, stwisepar.m, vie_lsm_gui_clock.fig,
%   vie_lsm_gui_clock.m, vie_lsm_gui_eop.fig, vie_lsm_gui_eop.m,
%   vie_lsm_gui_first.fig, vie_lsm_gui_first.m, vie_lsm_gui_global.fig,
%   vie_lsm_gui_global.m, vie_lsm_gui_sourcoor.fig, vie_lsm_gui_sourcoor.m,
%   vie_lsm_gui_statcoor.fig, vie_lsm_gui_statcoor.m, vie_lsm_gui_tropo.fig,
%   vie_lsm_gui_tropo.m, a_tidpm, a_tidut
%   write_outlier_file.m
%   a_love.m, a_shida.m, a_stsespos.m, a_rg.m,
%   satellitewisepar.m, apw_satellite.m, hpoc_satellites.m
%
%   Coded for VieVS:
%   12 May 2009 by Johannes Boehm and Kamil Teke
%
%   Revision:
%   04 Jun 2009 by Johannes Boehm: output for global solutions and
%   the paths of output changed
%   04 Aug 2009 by Kamil Teke: the parametrization of lsmopt.m changed
%   24 Aug 2009 by Johannes Boehm: a function out of it created
%   12 Sep 2009 by Kamil Teke: single and global solutions for sources
%   added
%   28 Sep 2009 by Kamil Teke: GUIs added
%   06 Oct 2009 by Kamil Teke: header added
%   18 Nov 2009 by Hana Spicakova: partials for station velocities for
%   global solution added
%   05 Dec 2009 by Kamil Teke: outlier reduction (elemination) and iteration related
%   codes deleted, outliers are book kept to an ASCII file but not
%   eleminated anymore in this subroutine of VieVS
%   16 Dec 2009 by Hana Spicakova: changes related to global solution
%   28 Jan 2010 by Tobias Nilsson: outlier now written to files in yearly
%   directories (../DATA/OUTLIERS/1994, ../DATA/OUTLIERS/2008, etc).
%   05 Feb 2010 by Tobias Nilsson: Removed a bug in the outlier test.
%   18 Feb 2010 by Kamil Teke: Absolute constraints for troposphere north
%   and east gradients added.
%   23 Feb 2010 by Kamil Teke: file paths for input and output are modified
%   user defined sub-directory for input and output is available
%   21 Jun 2010 by Hana Spicakova: changes related to global solution
%   02 Sep 2010 by Sigrid Boehm: added condition number to opt_
%   15 Sep 2010 by Johannes Boehm: the zenith distance is now provided to
%   reduce_oc instead of the hydrostatic mapping function
%   05 Oct 2010 by Hana Spicakova: final N matrix (datum free) and b vector
%   can be reduced and stored for sinex output (option in the last GUI)
%   10 Nov 2010 by Hana Spicakova: it is possible only to build N and b
%   for global solution (LEVEL2) or for sinex output (LEVEL3/SINEX) without
%   estimating parameters with LSM from the single session (option in the last GUI)
%   12 Nov 2010 by Hana Spicakova: reference epoch for estimation of
%   station velocities is taken from GUI
%   17 Jan 2010 by Hana Spicakova: parameters for SOLUTION/STATISTIC block
%   in sinex file are stored in LEVEL3/SINEX/col_sinex_*.mat
%   18 Apr 2011 by Tobias Nilsson: Changed input variables
%   19 May 2011 by Hana Spicakova: Sinex output for zwd and tgr added
%   06 Jun 2012 by Matthias Madzak: Residual parameters are now saved in
%   DATA/LEVEL3/ as res_*.mat. Those are needed to plot residuals later
%   (e.g. in the interface).
%   19 Jun 2012 by Hana Krasna:
%       Love and Shida numbers
%       FCN period from solid Earth tides
%       acceleration of SSB
%       velocity of sources
%       amplitudes of seasonal variation in station positions
%       pole tide Love and Shida number
%   25 Sep 2012 by Hana Krasna: partials for global solution w.r.t. relativistic parameter
%          gamma added (only Sun's contribution)
%   24 Mar 2013 by Hana Krasna: IERS name of source is added to x_ (glob2)
%   26 Mar 2013 by Hana Krasna: info added to opt_.source even if only N
%           matrix is produced - needed for Vie_Glob
%   24 Apr 2013 by Matthias Madzak: screen output of m0 changed to chi^2
%   24 Apr 2013 by Hana Krasna: variance factor for Sinex output corrected
%   26 Apr 2013 by Sigrid Boehm: output of condition number of N matrix removed
%   04 Oct 2013 by Hana Krasna: estimation of antenna axis offset added -
%          only for global solution
%   05 Dec 2013 by Hana Krasna: APL regression coefficients added
%   07 Jan 2014 by Hana Krasna: bug corrected: if requested LEVEL2 data can be
%          stored in a different subdirectory
%   16 May 2014 by Monika Tercjak: change the way of preparing weight
%           matrix - possibility to down-weight stations
%   28 May 2014 by Younghee Kwak: WRMS of post-fit residuals added
%                                 unit of residual v changed to [cm] at the comment
%   18 Jun 2014 by Hana Krasna: sources can be estimated in vie_lsm with NNR condition
%   29 Jul 2014 by Hana Krasna: baseline dependent weighting
%           (Minttu Uunila coded the function calc_bas_weights.m)
%   13 Jan 2015 by Matthias Madzak: Contact details now given as input for
%           writing SINEX files
%   14 Jan 2015 by Daniel Landskron: shape of output in Command Window
%           slightly changed
%   03 Feb 2015 by Andreas Hellerschmied: LEVEL1 subfilder ("opt.level1OutDir") is handed over
%       to the function "write_sinex_vievs".
%   19 Feb 2015 by Andreas Hellerschmied: Bug fixing (MATLAB2014b compatibility issues)
%   10 Apr 2015 by David Mayer: Bugs in first solution fixed
%   28 May 2015 by Lucia Plank: hardcoded option of constraining HbHo; to
%       apply, search for 'sibling' and set the parameter =1
%   19 Jun 2015 by David Mayer:bug fix the condition number of the N matrix is now calculated correctly
%   21 Jul 2015 by A. Girdiuk: explanation report about outliers
%   24 Aug 2015 by S. Boehm: estimation of tidal ERP variations
%                           - only for global solution
%   21 Oct 2015 by Hana Krasna: seasonal variations in station positions (annual, semi-annual), pole
%               tide Love and Shida numbers, APL regression coefficients added
%   04 Dec 2015 by A. Hellerschmied: OPT file is now loaded in VIE_INIT only!
%   04 Feb 2016 by Y. Kwak: common parameter constraints (ngr&egr) at the co-located sites, currently only for Hobart
%   19 Feb 2016 by Y. Kwak: common parameter constraints (clk&zwd) and local ties at the co-located sites,
%                           currently only for Hobart; to apply, set 'common = 1'
%   11 Jun 2016 by H. Krasna: implementation of partials for the
%          FCN period and amplitudes
%   14 Jul 2016 by A. Hellerschmied: - "opt" isn't a global variable any more
%                                    - Content of "sources" structure (natural sources, quasars) is now stored in the sub-structure "sources.q"
%   23 Sep 2016 by D. Mayer: added significance test of estimates
%   10 Oct 2016 by A. Girdiuk: code-optimized concerning loop over stations and A matrix of partial derivatives
%   10 Oct 2016 by A. Hellerschmied: Enabled support of observations to satellites.
%                                    - In case the required field "obs_type" in the scan structure is missing and no satellites were observed, obs. type "q" is assigned here!
%   12 Oct 2016 by A. Hellerschmied: Separate function to write outlier files (write_outlier_file)
%   19 Oct 2016 by A. Girdiuk: initialization of parameters for sinex part; in part for glob: order of parameters is now consistent with concatenation in matrix
%   04 Nov 2016 by H. Krasna: bug in computing the number of constraints for the global solution corrected (w.o. abs. constr. on source coordinates)
%   30 Nov 2016 by A. Girdiuk: helmert.m is introduced for station coordinates constraints instead of nnt.m and
%                               opt.fixed_station message is modified accordingly
%   15 Dec 2016 by H. Krasna: introducing of elevation dependent noise in the P matrix
%   17 Jan 2017 by H. Krasna: bug by creating the A matrix for station seasonal
%                     harmonic variation corrected (probably Matlab version dependent)
%   17 Feb 2017 by Y. Kwak: hardcoded option of common parameter constraints (ngr,egr,zwd)
%                           three sites added: Hartebeesthoek, Wettzell, Yebes
%   23 Feb 2017 by A. Hellerschmied: constant noise added to sig. for P matrix is now taken from GUI
%   01 Mar 2017 by A. Hellerschmied: call of "write_opt.m" removed along with the "manually find clock break" functions
%   04 May 2017 by A. Hellerschmied: - Added some functions to estimate pwl offsets to a priori satellite positions (preliminary implementation!)
%                                    - Some minor general revisions were done
%   21 Sep 2017 by H. Krasna: bugs corrected which appeared after the implemenation of satellite positions update w.r.t. global solution
%   06 Mar 2019 by D. Landskron: suffix checkbox added to the sinex files
% ************************************************************************

function vie_lsm(antenna, sources, scan, parameter, dirpth, dirpthL2)

tic

fprintf('---------------------------------------------------------------\n')
fprintf('|                   Welcome to VIE_LSM!!!!!                   |\n')
fprintf('---------------------------------------------------------------\n\n')

% ##### Init.: #####
format compact;

% Constants:
c = 299792458; % velocity in m/s
rad2mas = (180/pi)*3600*1000; % radian to milli arc second

% Preallocate variables:
x   = [];
mi  = [];
tso = []; % +hana 10Nov10

% Check sub-directories:
if ~exist('dirpth')
  dirpth='';
end
if ~exist('dirpthL2')
  dirpthL2='';
end


% #######################################################################
% # Clean VieVS structures                                              #
% #######################################################################
% Based on
% - OPT files
% - Outliers
% - Observation restrictions
fprintf('0. APPLY OUTLIERS AND ANALYSIS OPTIONS\n');

% ########################
% ##### Load Options #####
% ########################

% Init.:
clean_opt.sta_excl          = ''; % Exclude stations
clean_opt.sta_excl_start    = 0;  % Exclude stations start epoch
clean_opt.sta_excl_end      = 0;  % Exclude stations end epoch

clean_opt.sour_excl         = ''; % Exclude sources
clean_opt.bas_excl          = []; % Exclude baselines

clean_opt.num_clk_breaks        = 0;  % Number of clock breaks
clean_opt.refclock              = ''; % Reference clock
clean_opt.clk_break.stat_name   = []; % Clock break stations
clean_opt.clk_break.mjd         = []; % Clock break epochs

clean_opt.stat_dw       = ''; % Downweight stations
clean_opt.no_cab        = ''; % Cable cal

clean_opt.bdco_est.sta1 = '';
clean_opt.bdco_est.sta2 = '';

clean_opt.scan_excl     = []; % Outliers

parameter.opt.options = clean_opt; % Init


% ini_opt.bas_excl(nex).sta1 = station_name_1_str;
% ini_opt.bas_excl(nex).sta2 = station_name_2_str;
% bas_excl(nex,:) = [station_name_1_str, ' ', station_name_2_str];

remove_sprecial_stations = false;
stations_to_be_removed = {''; ''; ''; ''}; % Can be set here in cell array!


% ##### OPT files #####

% Init.:
parameter.vie_init.stat_dw = [];

% read OPT-file
if length(parameter.session_name) == 14
    opt_file_path_name = ['../../VLBI_OPT/', parameter.vie_init.diropt, '/', parameter.year, '/', parameter.session_name(1:end-5), '.OPT'];
elseif length(parameter.session_name) == 9
    opt_file_path_name = ['../../VLBI_OPT/', parameter.vie_init.diropt, '/', parameter.year, '/', parameter.session_name, '.OPT'];
elseif length(parameter.session_name) ~= 19
    warning('Session name does not follow convention');
end
if parameter.opt.use_opt_files
    if exist(opt_file_path_name, 'file')
        [clean_opt, ~] = readOPT(opt_file_path_name,remove_sprecial_stations,stations_to_be_removed);
        parameter.opt.options = clean_opt;
    else
        fprintf(' - No OPT file applied: %s\n', opt_file_path_name);
        if remove_sprecial_stations
            clean_opt.sta_excl = char(stations_to_be_removed);
            clean_opt.sta_excl_start = zeros(1,size(stations_to_be_removed,1));
        end
    end
end

% ### write info about excluded baselines, stations, sources to CW ###
if parameter.opt.use_opt_files
    fprintf('Stations to be excluded: %1.0f\n', size(clean_opt.sta_excl,1))
    for k=1:size(clean_opt.sta_excl,1)
        if clean_opt.sta_excl_start(k)==0
            fprintf('%s \n', clean_opt.sta_excl(k,:));
        else
            fprintf('%s %f %f\n', clean_opt.sta_excl(k,:), clean_opt.sta_excl_start(k),clean_opt.sta_excl_end(k));
        end
    end

    fprintf('Stations to be down-weighted: %1.0f\n', size(clean_opt.stat_dw,1))
    if size(clean_opt.stat_dw,1) == 0
        parameter.vie_init.stat_dw = [];
    else
        parameter.vie_init.stat_dw = {};
        for k=1:size(clean_opt.stat_dw,1)
            fprintf('%s', clean_opt.stat_dw(k,:),' ', clean_opt.stat_co(k,:))
            fprintf('\n')
            parameter.vie_init.stat_dw(k,:) = {clean_opt.stat_dw(k,:)};
            parameter.vie_init.stat_co(k,:) = str2num(clean_opt.stat_co(k,:));
        end
    end

    fprintf('Sources to be excluded: %1.0f\n', size(clean_opt.sour_excl,1))
    for k=1:size(clean_opt.sour_excl,1)
        %fprintf('%s\n', ini_opt.sour_excl(k,:))
        if clean_opt.sour_excl_start(k)==0
            fprintf('%s \n', clean_opt.sour_excl(k,:));
        else
            fprintf('%s %f %f\n', clean_opt.sour_excl(k,:), clean_opt.sour_excl_start(k),clean_opt.sour_excl_end(k));
        end
    end

    remove_sources_from_list = false;
    if remove_sources_from_list
        path2sourcelist = '';     %add the path of your .txt file here. Format is the same as glob input .txt files
        fid = fopen(path2sourcelist);
        if fid == -1
            warning('File with list of removed sources can not be found\n');
        else
            remove_sources = textscan(fid, '%8s','Delimiter','\n');
            remove_sources = remove_sources{1};
            clean_opt.sour_excl = [clean_opt.sour_excl;char(remove_sources)];
            disp('+ sources from external file removed');
            if isfield(clean_opt, 'sour_excl_start')
                clean_opt.sour_excl_start = [clean_opt.sour_excl_start, zeros(1,length(remove_sources))];
                clean_opt.sour_excl_end = [clean_opt.sour_excl_end, zeros(1,length(remove_sources))];
            else
                clean_opt.sour_excl_start = zeros(1,length(remove_sources));
                clean_opt.sour_excl_end = zeros(1,length(remove_sources));
            end
        end
        fclose(fid);
    end

    fprintf('Baselines to be excluded: %1.0f\n', size(clean_opt.bas_excl, 2))
    for k=1:size(clean_opt.bas_excl, 2)
        fprintf('%8s - %8s \n', clean_opt.bas_excl(k).sta1, clean_opt.bas_excl(k).sta2)
    end
    fprintf('No cable calibration: %1.0f\n', size(clean_opt.no_cab,1))
    for k=1:size(clean_opt.no_cab,1)
        fprintf('%s\n', clean_opt.no_cab(k,:))
    end
    fprintf('\n')
end

% ##### Outlier files #####
parameter.outlier.obs2remove = []; % init empty field for outliers

outlier_filename_path = ['../DATA/OUTLIER/', parameter.outlier.out_file_dir, '/', parameter.year, '/', parameter.session_name, '.OUT'];
if parameter.outlier.flag_remove_outlier
    if exist(outlier_filename_path, 'file')
        [parameter.outlier.obs2remove] = readOUT(outlier_filename_path);
        fprintf('%d outliers will be removed:\n',size(parameter.outlier.obs2remove,2));        %%%=> A. Girdiuk 2015-07-21
        for k=1:size(parameter.outlier.obs2remove,2)
            fprintf(' - %8s %8s %5.2f %8s\n', parameter.outlier.obs2remove(k).sta1, parameter.outlier.obs2remove(k).sta2, parameter.outlier.obs2remove(k).mjd, parameter.outlier.obs2remove(k).sou);
        end
    else
        fprintf('Outlier list not available: %s\n', outlier_filename_path);
    end
else
    fprintf('Outliers will not be removed\n');
end
fprintf('\n')


% ############################
% ##### Clean structures #####
% ############################
% => Remove scans, observations, antennas and soures based on analysis options (obs. restrictions, OPT and OUTLIER files)
[scan, sources, antenna] = cleanScan(scan, sources, antenna, parameter);

% ##### Check minimum number of observations per estimated parameter #####

% The following options have to be set in the GUI in futur (preliminary here):

% Required number of observations per source:
% - If this threshold is not met => Remove observation(s) from vievs data structures (scan, antenna, sources)
% - Only remove observations, if source coordinates are estimated!Otherwise they are fixed anyway.

% Check, if source coordinates are estimated and if a threchold for the num of observations per source is set:
if (parameter.lsmopt.min_num_obs_per_est_source > 1) && (parameter.lsmopt.pw_sou || parameter.lsmopt.est_sourceNNR)
    fprintf(1, 'Remove observations to sources with less then %d observations in total.\n', parameter.lsmopt.min_num_obs_per_est_source);
    % Check threchold and remove observations to sources, where the criteria is not met:
    [scan, sources, antenna] = check_num_of_obs_per_parameter(scan, sources, antenna, parameter);
end


fprintf('1. LOAD AND PREPARE DATA\n');

opt = parameter.lsmopt;

% to ensure compatibility with parameter files before the vievs update (10/2020)
% where the estimation of scale offset was added. Can be removed in the future!
if ~isfield(parameter.lsmopt, 'est_scale') 
    opt.est_scale=0;
end

% to ensure compatibility with parameter files before the vievs update (10/2020)
% where the bas-dep clock offsets were added. Can be removed in the future!
if ~isfield(parameter.lsmopt, 'est_bdco') 
   opt.est_bdco=0; 
end

if opt.addSnxSource==1
    opt.est_source=1;
end

% ##### Initial checks: #####
% Check, if the sources structure has the correct format (with fields "q" and "s" for quasars and satellites respectively).
% ==> Older format (VieVS 2.3): copy content from "sources" to "sources.q"
% !!! Remove this check after the release of VieVS V3.0!!
if ~isfield(sources, 'q') && ~isfield(sources, 's')
%     error('Please run VIE_MOD and VIE_INIT again for session %s! The format of the "sources" structure is not supported any more (sub-strucutre sources.q missing).', parameter.session_name);
    sources_tmp = sources;
    clear sources;
    sources.q = sources_tmp;
    clear sources_tmp;
    % In case an older scan structure (from LEVEL1) is loaded without the field scan.obs_type:
    % ==> Assume, that only quasars were observed and add the field!
    [scan.obs_type] = deal('q');
    warning('The format of the "sources" structure of session %s is not supported any more => "sources" copied to "sources.q".', parameter.session_name);
end

% ##### Get session info #####
n_scan      = length(scan);     % number of scans
na          = length(antenna);  % number of antennas
n_observ    = sum([scan.nobs]); % number of observations
% number of sources
% - Satellites
if isfield(sources, 's')
    ns_s = length(sources.s);
else
    ns_s = 0;
end
% - Quasars
if isfield(sources, 'q')
    ns_q = length(sources.q);
else
    ns_q = 0;
end


% Write info to CW:
fprintf('number of scans            : %d\n',n_scan);
fprintf('number of antennas         : %d\n',na);
if ns_q > 0
    fprintf('number of sources (quasars): %d\n',ns_q);
end
if ns_s > 0
    fprintf('number of sources (sc)     : %d\n',ns_s);
end
fprintf('number of obs.             : %d\n',n_observ);

mjd1            = min([scan.mjd]);     % The time of the first scan in mjd, UTC
mjd0            = floor(mjd1);         % The midnight (the beginning of the day of the session), UTC
mjd2            = max([scan.mjd]);     % The time of the last scan in mjd, UTC
opt.first_scan  = mjd1;
opt.last_scan   = mjd2;

% ASSIGNING VARIABLES FROM 'scan' TO SOURCE 'obs_per_source'
% - Quasars

if ns_q > 0
%     % preallocation:
    obs_per_source(1,ns_q) = struct('mjd', [], 'iso', [], 'i1', [], 'i2', []);
    for isou = 1 : ns_q
        n_obs_per_src = 0;                  % Observation of the specific source
        for i_scan = 1 : n_scan               % Number of scans per session
            if strcmp(scan(i_scan).obs_type, 'q') && (scan(i_scan).iso == isou)
                for iobs = 1 : scan(i_scan).nobs  % Number of observations per scan
                    n_obs_per_src                           = n_obs_per_src + 1;
                    obs_per_source(isou).mjd(n_obs_per_src) = scan(i_scan).mjd;
                    obs_per_source(isou).iso(n_obs_per_src) = scan(i_scan).iso;
                    obs_per_source(isou).i1(n_obs_per_src)  = scan(i_scan).obs(iobs).i1;
                    obs_per_source(isou).i2(n_obs_per_src)  = scan(i_scan).obs(iobs).i2;
                end
            end
        end
    end
else
    obs_per_source = [];
end
% - Satellites
if ns_s > 0
    % preallocation:
    obs_per_satellite(1, ns_s) = struct('mjd', [], 'iso', [], 'i1', [], 'i2', []);
    for isou = 1 : ns_s % Loop over all satellites (isou = sat. ID)
        n_obs_per_src = 0;                  % Observation of the specific source
        for i_scan = 1 : n_scan               % Number of scans per session
            if strcmp(scan(i_scan).obs_type, 's') && (scan(i_scan).iso == isou)
                for iobs = 1 : scan(i_scan).nobs  % Number of observations per scan
                    n_obs_per_src                           = n_obs_per_src + 1;
                    obs_per_satellite(isou).mjd(n_obs_per_src) = scan(i_scan).mjd;
                    obs_per_satellite(isou).iso(n_obs_per_src) = scan(i_scan).iso;
                    obs_per_satellite(isou).i1(n_obs_per_src)  = scan(i_scan).obs(iobs).i1;
                    obs_per_satellite(isou).i2(n_obs_per_src)  = scan(i_scan).obs(iobs).i2;
                end
            end
        end
    end
else
    obs_per_satellite = [];
end

% fprintf('3. READING OPTIONS FILE "lsm_opt.dat" FOR SPECIFIYING the LSM PARAMETERS\n');
% -------------------------------------------------------------------------
% READING THE FILE #lsm_opt.dat# specifiy the LSM parameters
% SUBROUTINES THAT ARE USED "read_lsmopt"
fprintf('2. CREATING DEFAULT OPTIONS\n');

opt = lsmopt(antenna, sources, na, ns_q, ns_s, obs_per_source, obs_per_satellite, parameter, opt);
opt.scans_total = n_scan;
opt.midnight = ceil(mjd1);

fprintf('3. FORMING THE WEIGHT MATRIX OF THE OBSERVATIONS "Pobserv"\n');
% -------------------------------------------------------------------------
% FORMING THE WEIGHT MATRIX OF THE OBSERVATIONS (Pobserv)
temp = [scan.obs];
mi_observ = [temp.sig]; % [seconds]
stations =  parameter.vie_init.stat_dw; %list of down-weighted stations


% Add noise
if opt.eldep_noise==1 % Elevation dependent noise
    fprintf('INTRODUCING ELEVATION DEPENDENT NOISE \n');
    k = 0;
    for j = 1 : length(scan)
        for i = 1 : scan(j).nobs

            e1 = pi/2-scan(j).stat(scan(j).obs(i).i1).zd;
            e2 = pi/2-scan(j).stat(scan(j).obs(i).i2).zd;

            addnoise_ps2= (6/sin(e1))^2 + (6/sin(e2))^2; % ps^2

            k = k + 1;
            addnoise(k)=sqrt(addnoise_ps2*(10^-24)) * c; % m
        end
    end
else
    % add constand noise (noise level defined in GUI => VieVS estimation settings'; default: 1 cm)
    addnoise = opt.add_const_noise_cm * 1e-2;
    fprintf('Constant noise added to input sig.: %4.2f cm\n', opt.add_const_noise_cm);
end

%---------------------------
if ~isempty(stations)
    nmi_observ = (addnoise/c).^2+(mi_observ).^2;  % [seconds2]
    an_weight = parameter.vie_init.stat_co;
    numbers = find(ismember ({antenna.name},stations)==1);
    i12 =  [[temp.i1]; [temp.i2] ];
    for j = 1:size(numbers,2)
        aa = find(ismember(i12(1,:),numbers(j))==1);  bb = find(ismember(i12(2,:),numbers(j))==1); AB = [aa,bb]; AB = sort(AB);
        for i = 1:size(AB,2)
            nmi_observ(AB(i)) = nmi_observ(AB(i)) + (an_weight(j)/c)^2;
        end
    end
    nmi_observ = sqrt(nmi_observ); % [seconds]
else
    nmi_observ = sqrt((addnoise/c).^2+(mi_observ).^2); % [seconds]
end
%-------------------
so = sqrt(((nmi_observ.*c*100)*(nmi_observ.*c*100)')/n_observ); % [cm] apriori std. dev. of unit weight
opt.so = so;
Pobserv = diag(sparse(1./((nmi_observ.^2).*c^2*100^2))); % [1/cm^2]
%Qll = diag(sparse((nmi_observ.^2).*c^2*100^2)); % [cm^2]
fprintf('apriori std. dev. of unit weight. : %3.4f\n',so);

% preallocation
if ~opt.est_stsespos
    per_stat(1,na)=struct('mjd', [], 'oc_nob', [], 'zd', [], 'first', [], 'other', [], 'mf', [], 'az', [],...
                 'dx', [],     'dy', [], 'dz', [],    'xo', [],    'yo', [], 'zo', [],...
                 'dAO', [], 'drg', []);
else
    per_stat(1,na)=struct('mjd', [], 'oc_nob', [], 'zd', [], 'first', [], 'other', [], 'mf', [], 'az', [],...
                 'dx', [],     'dy', [], 'dz', [],    'xo', [],    'yo', [], 'zo', [],...
                 'dAO', [], 'drg', [], 'pAcr', [], 'pAce', [], 'pAcn', [], 'pAsr', [], 'pAse', [], 'pAsn', []);
end
% -------------------------------------------------------------------------
% ASSIGNING THE INPUT PARAMETERS TO EACH STATION
for istat = 1:na % number of stations
    i = 0; n_obs_per_src = 0;
    for itim = 1:n_scan % number of scans per session
        for iobs = 1:scan(itim).nobs % number of observations per scan
            i = i + 1; % i : i. observation in the session
            i1 = scan(itim).obs(iobs).i1;
            i2 = scan(itim).obs(iobs).i2;
            if i1 == istat || i2 == istat
                n_obs_per_src = n_obs_per_src + 1; % k : k. observation of the specific station
                per_stat(istat).mjd(n_obs_per_src) = scan(itim).mjd; % The times of scans [day]
                per_stat(istat).oc_nob(n_obs_per_src) = i;
                per_stat(istat).zd(n_obs_per_src) = scan(itim).stat(istat).zd;  % Boehm 21 Aug 2009, 15 Sep 2010
                obs_mjd(i) = (scan(itim).mjd - mjd0)*24*60; % [minute]
            end
            if i1 == istat
                per_stat(istat).first(n_obs_per_src) = -1;
                per_stat(istat).other(n_obs_per_src) = i2;  % Boehm 21 Aug 2009

                %
                %   partials prepared
                %

                per_stat(istat).mf(n_obs_per_src) = scan(itim).stat(i1).mfw; % The mapping function value
                per_stat(istat).az(n_obs_per_src) = scan(itim).stat(i1).az; % Azimuth [radians]
%                    per_stat(istat).zd(k) = scan(itim).stat(i1).zd; % Zenith distance [radians]

                per_stat(istat).dx(n_obs_per_src) = -scan(itim).obs(iobs).pstat1(1); % The partial derivatives of delay wrt antenna coordinates
                per_stat(istat).dy(n_obs_per_src) = -scan(itim).obs(iobs).pstat1(2);
                per_stat(istat).dz(n_obs_per_src) = -scan(itim).obs(iobs).pstat1(3);

                %obs_per_stat.dx(k) = scan(itim).stat(i1).pantd(1); % The partial derivatives of delay wrt antenna coordinates
                %obs_per_stat.dy(k) = scan(itim).stat(i1).pantd(2);
                %obs_per_stat.dz(k) = scan(itim).stat(i1).pantd(3);

                per_stat(istat).xo(n_obs_per_src) = scan(itim).stat(i1).x(1); % Apriori coordinates of the antennas
                per_stat(istat).yo(n_obs_per_src) = scan(itim).stat(i1).x(2); % y
                per_stat(istat).zo(n_obs_per_src) = scan(itim).stat(i1).x(3); % z

                per_stat(istat).dAO(n_obs_per_src) = -scan(itim).obs(iobs).pAO_st1;

                per_stat(istat).drg(n_obs_per_src)= -scan(itim).obs(iobs).prg_st1;

                 if opt.est_stsespos ==1
    		        per_stat(istat).pAcr(n_obs_per_src,:) = -scan(itim).obs(iobs).pAcr_st1;
    		        per_stat(istat).pAce(n_obs_per_src,:) = -scan(itim).obs(iobs).pAce_st1;
    		        per_stat(istat).pAcn(n_obs_per_src,:) = -scan(itim).obs(iobs).pAcn_st1;
    		        per_stat(istat).pAsr(n_obs_per_src,:) = -scan(itim).obs(iobs).pAsr_st1;
    		        per_stat(istat).pAse(n_obs_per_src,:) = -scan(itim).obs(iobs).pAse_st1;
		            per_stat(istat).pAsn(n_obs_per_src,:) = -scan(itim).obs(iobs).pAsn_st1;
                 end
            end
            if i2 == istat
                per_stat(istat).first(n_obs_per_src) = +1;
                per_stat(istat).other(n_obs_per_src) = i1;  % Boehm 21 Aug 2009

                per_stat(istat).mf(n_obs_per_src) = scan(itim).stat(i2).mfw; % The mapping function value
                per_stat(istat).az(n_obs_per_src) = scan(itim).stat(i2).az; % Azimuth [radians]
%                    per_stat(istat).zd(k) = scan(itim).stat(i2).zd; % Zenith distance [radians]

                per_stat(istat).dx(n_obs_per_src) = scan(itim).obs(iobs).pstat2(1); % The partial derivatives of delay wrt antenna coordinates
                per_stat(istat).dy(n_obs_per_src) = scan(itim).obs(iobs).pstat2(2);
                per_stat(istat).dz(n_obs_per_src) = scan(itim).obs(iobs).pstat2(3);

                %obs_per_stat.dx(k) = scan(itim).stat(i2).pantd(1); % The partial derivatives of delay wrt antenna coordinates
                %obs_per_stat.dy(k) = scan(itim).stat(i2).pantd(2);
                %obs_per_stat.dz(k) = scan(itim).stat(i2).pantd(3);

                per_stat(istat).xo(n_obs_per_src) = scan(itim).stat(i2).x(1); % Apriori coordinates of the antennas
                per_stat(istat).yo(n_obs_per_src) = scan(itim).stat(i2).x(2);
                per_stat(istat).zo(n_obs_per_src) = scan(itim).stat(i2).x(3);

                per_stat(istat).dAO(n_obs_per_src) = scan(itim).obs(iobs).pAO_st2;

                per_stat(istat).drg(n_obs_per_src) = scan(itim).obs(iobs).prg_st2;

                 if opt.est_stsespos ==1
        		    per_stat(istat).pAcr(n_obs_per_src,:) = scan(itim).obs(iobs).pAcr_st2;
        		    per_stat(istat).pAce(n_obs_per_src,:) = scan(itim).obs(iobs).pAce_st2;
        		    per_stat(istat).pAcn(n_obs_per_src,:) = scan(itim).obs(iobs).pAcn_st2;
        		    per_stat(istat).pAsr(n_obs_per_src,:) = scan(itim).obs(iobs).pAsr_st2;
        		    per_stat(istat).pAse(n_obs_per_src,:) = scan(itim).obs(iobs).pAse_st2;
		            per_stat(istat).pAsn(n_obs_per_src,:) = scan(itim).obs(iobs).pAsn_st2;
                 end
            end
        end
    end
    fprintf('obs. of the antenna %s : %4d\n',antenna(istat).name,n_obs_per_src);
end

fprintf('4. FORMING THE REDUCED OBSERVATION VECTOR "oc_observ"\n');
% -------------------------------------------------------------------------
% FORMING THE REDUCED OBSERVATION VECTOR (oc_observ)
% oc_observ = [];
numberOfLSMs = length(temp(1).obs);
oc_observ = ([temp.obs]-repmat([temp.com],numberOfLSMs,1))'.*c*100; % [cm]

if opt.first ~= 1
    first_solution.ref_st = [];
elseif opt.first == 1
    [oc_observ,first_solution,opt] = reduce_oc(n_observ,na,n_scan,scan,Pobserv,oc_observ,opt,per_stat,parameter.session_name,dirpth);
end

if opt.second == 0
    return
end

% clear temp

fprintf('5. FORMING THE DESIGN MATRICES "A(i).sm" ...\n' );
% -------------------------------------------------------------------------
% FORMING THE DESIGN MATRICES OF THE REAL OBSERVATIONS ([A1clk A2clk A3zwd])
% SUBROUTINES THAT ARE USED "stwisepar","apwq_clk","apw_zwd"
number_of_estimated_parameters = 20;
sum_.clk(1)     = 0;
sum_.qclk(1)    = 0;
sum_.zwd(1)     = 0;
sum_.ngr(1)     = 0;
sum_.egr(1)     = 0;
sum_.xyz(1)     = 0;
A(1, number_of_estimated_parameters)    = struct('sm', []);
H(1, number_of_estimated_parameters)    = struct('sm', []);
Ph(1, number_of_estimated_parameters)   = struct('sm', []);
och(1, number_of_estimated_parameters)  = struct('sv', []);
dj                                      = zeros(1, number_of_estimated_parameters);
num_psob                                = zeros(1, number_of_estimated_parameters);
sum_dj                                  = zeros(1, number_of_estimated_parameters + 1);

if opt.global_solve == 1 || opt.ascii_snx ==1
    A_Acr=[]; A_Ace=[]; A_Acn=[]; A_Asr=[]; A_Ase=[]; A_Asn=[]; A_AO=[]; A_rg=[];
end

t(1,na) = struct('clk',[],'zwd',[],'ngr',[],'egr',[],'xyz',[]);             % Estimation intervals
xo(1,na)=0; yo=xo; zo=xo;                                                   % A priori antenna coordinates
n_(1,na) = struct('clk',[],'qclk',[],'zwd',[],'ngr',[],'egr',[],'xyz',[]);  % Number of estimated offsets per estimation interval

for istat = 1 : na
    a = struct('pwclk',[],'rqclk',[],'zwd',[],'ngr',[],'egr',[],'x',[],'y',[],'z',[]);
    %clkbreak = opt.stat(istat).clkbreak; % The vector of clock break per clock [minutes]
    mjdstat = per_stat(istat).mjd; % time of the observations per station [day]
    int.zwd = opt.stat(istat).int_zwd; % zwd estimation interval per station [minute]
    int.clk = opt.stat(istat).int_clk; % clock estimation interval per station [minute]
    int.egr = opt.stat(istat).int_egr; % east gradients estimation interval per station [minute]
    int.ngr = opt.stat(istat).int_ngr; % north gradients estimation interval per station [minute]
    int.xyz = opt.stat(istat).int_xyz; % coordinates estimation interval per station [minute]

    % ASSIGN ALL THE PARAMETERS FROM SCANWISE TO STATIONWISE - FORM1
    [obs_per_stat,n_unk,T_] = stwisepar(per_stat,mjdstat,istat,mjd0,int);

    t(istat).clk = T_.clk; % estimation intervals of clock function offsets
    t(istat).zwd = T_.zwd; % estimation intervals of zenith wet delay offsets
    t(istat).ngr = T_.ngr; % estimation intervals of tropospheric north gradient offsets
    t(istat).egr = T_.egr; % estimation intervals of tropospheric east gradient offsets
    t(istat).xyz = T_.xyz; % estimation intervals of antenna coordinate offsets

    xo(istat) = per_stat(istat).xo(1);
    yo(istat) = per_stat(istat).yo(1);
    zo(istat) = per_stat(istat).zo(1);

    % FORMING THE DESIGN MATRIX FOR THE CLOCKS AS CONTINUOUS PIECEWISE
    % LINEAR FUNCTIONS
    [Apwclk, Arqclk] = apwq_clk(obs_per_stat,n_observ,per_stat(istat).oc_nob,n_unk,T_,opt);

    % Concatenating
    a(istat).pwclk = Apwclk;   % Piecewise linear clock offsets (nobs of stat. x numb. of est. intervals)
    A(1).sm = horzcat(A(1).sm, a(istat).pwclk);
    % Concatenating
    a(istat).rqclk = Arqclk;   % Rate and quadratic terms of the clock function(nobs of stat. x 2)
    A(2).sm = horzcat(A(2).sm, a(istat).rqclk);

    % FORMING THE DESIGN MATRIX FOR THE WET ZENITH DELAYS AS CONTINUOUS
    % PIECEWISE LINEAR FUNCTIONS
     [Azwd] = apw_zwd(obs_per_stat,n_observ,per_stat(istat).oc_nob,n_unk,T_);
    % Concatenating
    a(istat).zwd = Azwd;        % (nobs of stat. x numb. of est. intervals)
    A(3).sm = horzcat(A(3).sm, a(istat).zwd);

    % FORMING THE DESIGN MATRIX FOR THE TROPOSPHERIC NORTH  GRADIENTS AS
    % CONTINUOUS PIECEWISE LINEAR FUNCTIONS
    [Apwngr] = apw_ngr(obs_per_stat,n_observ,per_stat(istat).oc_nob,n_unk,T_);
    % Concatenating
    a(istat).ngr = Apwngr;      % (nobs of stat. x numb. of est. intervals)
    A(4).sm = horzcat(A(4).sm, a(istat).ngr); % Concatenating

    % FORMING THE DESIGN MATRIX FOR THE TROPOSPHERIC EAST GRADIENTS AS
    % CONTINUOUS PIECEWISE LINEAR FUNCTIONS
    [Apwegr] = apw_egr(obs_per_stat,n_observ,per_stat(istat).oc_nob,n_unk,T_);
    % Concatenating
    a(istat).egr = Apwegr;      % (nobs of stat. x numb. of est. intervals)
    A(5).sm = horzcat(A(5).sm, a(istat).egr); % Concatenating

    if opt.pw_stc == 0 % if NO station coor. is estimated
        % OR if ONE OFFSET per session is estimated
        [Ax,Ay,Az] = a_xyz(obs_per_stat,n_observ,per_stat(istat).oc_nob);
        % Concatenating
        a(istat).x = Ax; a(istat).y = Ay; a(istat).z = Az;  % (nobs of stat. x numb. of est. intervals (=1))
    elseif opt.pw_stc == 1 % if station coor. is estimated as PWL OFFSETS
        [Apwx,Apwy,Apwz] = apw_xyz(obs_per_stat,n_observ,per_stat(istat).oc_nob,n_unk,T_);
        % Concatenating
        a(istat).x = Apwx; a(istat).y = Apwy; a(istat).z = Apwz;
    end
    A(13).sm = horzcat(A(13).sm, a(istat).x);
    A(14).sm = horzcat(A(14).sm, a(istat).y);
    A(15).sm = horzcat(A(15).sm, a(istat).z);

    if opt.est_vel == 1
        % +hana 18Nov09
        refvel_mjd=modjuldat(opt.refvel);
        opt.refvel_mjd=refvel_mjd;
        fctvel=(mjd0-refvel_mjd)/36525;  % [day/century] +hana 12Nov2010
        A_vx = A(13).sm; A_vx = A_vx*fctvel;
        A_vy = A(14).sm; A_vy = A_vy*fctvel;
        A_vz = A(15).sm; A_vz = A_vz*fctvel;
        % -hana
    end


    % A matrix for Axis Offset
    if opt.est_AO == 1
        [A_ao] = a_AO(obs_per_stat,n_observ,per_stat(istat).oc_nob);
        A_AO = horzcat(A_AO,A_ao);
    end

    % a_stsespos.m
    % Amplitudes of seasonal variations in the station positions
    nsewa=0;
    if opt.est_stsespos ==1
        [AAcr,AAce,AAcn,AAsr,AAse,AAsn] = a_stsespos(obs_per_stat,n_observ,per_stat(istat).oc_nob);
        % Concatenating - for global adjustment
        A_Acr = horzcat(A_Acr,AAcr);
        A_Ace = horzcat(A_Ace,AAce);
        A_Acn = horzcat(A_Acn,AAcn);
        A_Asr = horzcat(A_Asr,AAsr);
        A_Ase = horzcat(A_Ase,AAse);
        A_Asn = horzcat(A_Asn,AAsn);
        nsewa=size(obs_per_stat.pAcr,2);
    end

    if opt.est_rg == 1
        [AA_rg] = a_rg(obs_per_stat,n_observ,per_stat(istat).oc_nob);
        A_rg = horzcat(A_rg,AA_rg);
    end

    % the number of clock function offsets estimated for each estimation interval
    n_(istat).clk = size(a(istat).pwclk,2);
    sum_.clk(istat+1) = sum_.clk(istat) + n_(istat).clk;
    % the pair of rate and quadratic coefficients of the clock function
    n_(istat).qclk = size(a(istat).rqclk,2);
    sum_.qclk(istat+1) = sum_.qclk(istat) + n_(istat).qclk;
    % the number of zenith wet delay offsets estimated for each estimation interval
    n_(istat).zwd = size(a(istat).zwd,2);
    sum_.zwd(istat+1) = sum_.zwd(istat) + n_(istat).zwd;
    % the number of tropospheric north gradient offsets estimated for each estimation interval
    n_(istat).ngr = size(a(istat).ngr,2);
    sum_.ngr(istat+1) = sum_.ngr(istat) + n_(istat).ngr;
    % the number of tropospheric east gradient offsets estimated for each estimation interval
    n_(istat).egr = size(a(istat).egr,2);
    sum_.egr(istat+1) = sum_.egr(istat) + n_(istat).egr;
     % the number of site coordinate offsets estimated for each estimation interval
    n_(istat).xyz = size(a(istat).x,2);
    sum_.xyz(istat+1) = sum_.xyz(istat) + n_(istat).xyz;
end

clear Apwclock Arqclock Azwd Apwegr Apwngr Apwx Apwy Apwz Ax Ay Ay a A_ao

% -------------------------------------------------------------------------
% FORMING DESIGN MATRICES FOR SOURCE COORDINATES (QUASARS)
A_ra_glob       = [];               % For global parameter estimation in vie_glob
A_de_glob       = [];               % For global parameter estimation in vie_glob
A_vra_glob      = [];               % For global parameter estimation in vie_glob
A_vde_glob      = [];               % For global parameter estimation in vie_glob
sumso.sources   = zeros(1,ns_q);    % total vector of source coor. estimates after eliminating non-observed sources
nso.sources     = 0;                % the number of source coordinate offsets for each source

if (opt.global_solve==1 && opt.est_source==1) || (opt.ascii_snx==1 && opt.est_source==1) || opt.est_sourceNNR==1 || opt.pw_sou == 1
	for isou = 1 : ns_q  % loop over all sources (quasars)
    	mjdsource   = obs_per_source(isou).mjd;
    	int_sou     = opt.source(isou).int_rade;

    	% ASSIGN ALL THE PARAMETERS FROM SCANWISE TO SOURCEWISE
    	[per_source, n_unk,T_]   = sourcewisepar(opt,scan,n_scan,mjdsource,isou,mjd0,int_sou);
    	tso(isou).sources       = T_.source;
    	a(isou).ra              = [];
        a(isou).de              = [];

    	% FORMING DESIGN MATIRECES FOR SOURCE COORDINATES
    	% if one offset per session is estimated (for global solution)
    	if  (opt.global_solve==1 && opt.est_source==1) || (opt.ascii_snx==1 && opt.est_source==1) || opt.est_sourceNNR==1 % +hana 18Jun14
    	    [Ara,Ade] = a_source(per_source,n_observ);
    	    % Concatenating
    	    a(isou).ra = Ara; a(isou).de = Ade;
    	elseif opt.pw_sou == 1                          % if pwl offsets are estimated
    	    [Apw_ra,Apw_de] = apw_source(per_source,n_observ,n_unk,T_);
    	    % Concatenating
    	    a(isou).ra = Apw_ra; a(isou).de = Apw_de;
        end

    	nso(isou).sources = size(a(isou).ra, 2);         % the number of source coordinate offsets for each source

        % Preparations for global parameter estimation
    	if (opt.global_solve==1 && opt.est_source==1) || (opt.ascii_snx==1 && opt.est_source==1)
    	    A_ra_glob = horzcat(A_ra_glob,a(isou).ra);  % design matrix for right ascensions of sources
    	    A_de_glob = horzcat(A_de_glob,a(isou).de);  % design matrix for declination of sources
    	    % velocities of sources added
    	    if opt.est_source_velo ==1
    	        refvelsou_mjd= 51544; %2000
    	        opt.refvelsou_mjd=refvelsou_mjd;
    	        fctvelsou=(mjd0-refvelsou_mjd)/36525;   % [day/century]
    	        A_vra_glob = A_ra_glob*fctvelsou;
    	        A_vde_glob = A_de_glob*fctvelsou;
    	    end
        end

    	if opt.pw_sou == 1 || opt.est_sourceNNR==1
    	    A(11).sm = horzcat(A(11).sm, a(isou).ra);   % design matrix for right ascensions of sources
    	    A(12).sm = horzcat(A(12).sm, a(isou).de);   % design matrix for declination of sources
    	    ra(isou) = opt.source(isou).ra;
    	    de(isou) = opt.source(isou).de;
    	end
        sumso.sources(isou+1) = sumso.sources(isou) + nso(isou).sources; % total vector of source coor. estimates after eliminating non-observed sources
	end
end

ts.sources = sum([nso.sources]); %

clear Apw_ra Apw_de Ara Ade a

% -------------------------------------------------------------------------
% FORMING DESIGN MATRICES FOR SATELLITE COORDINATES

number_pwlo_per_sat = 0;

% If estimation of PWL offsets for satellite positions was selected
if opt.pw_sat
    % loop over all satellites
    for i_sat = 1 : ns_s

    	% ASSIGN ALL THE PARAMETERS FROM SCAN-WISE TO SATELLITE-WISE
    	[per_satellite, n_unk, T_] = satellitewisepar(opt, scan, n_scan, obs_per_satellite(i_sat).mjd, i_sat, mjd0, opt.satellite(i_sat).sat_pos_int);
    	tso(i_sat).sat = T_.sat;

    	% FORMING DESIGN MATRICES FOR SATELLITE COORDINATES
        [Apw_pos1, Apw_pos2, Apw_pos3] = apw_satellite(per_satellite, n_observ, n_unk, T_);
        % Concatenating
        a(i_sat).pos1 = Apw_pos1;
        a(i_sat).pos2 = Apw_pos2;
        a(i_sat).pos3 = Apw_pos3;

        number_pwlo_per_sat(i_sat) = size(a(i_sat).pos1, 2);          % number of coordinate offsets for each satellite

        A(16).sm = horzcat(A(16).sm, a(i_sat).pos1);   % design matrix for satellote coord. 1
    	A(17).sm = horzcat(A(17).sm, a(i_sat).pos2);   % design matrix for satellote coord. 2
        A(18).sm = horzcat(A(18).sm, a(i_sat).pos3);   % design matrix for satellote coord. 2

    end % for i_sat = 1 : ns_s
end % if opt.pw_sat

clear Apw_pos1 Apw_pos2 Apw_pos3 a




% -------------------------------------------------------------------------

% XPOL (piecewise) - ahp_xpol
[Apwxpol,T,Hxpol,Phxpol,oc_hxpol] = ahp_xpol(scan,mjd0,opt,c,rad2mas,obs_mjd);
A(6).sm = Apwxpol; H(6).sm = Hxpol; Ph(6).sm = Phxpol; och(6).sv = oc_hxpol;

% YPOL (piecewise) - ahp_ypol
[Apwypol,T,Hypol,Phypol,oc_hypol] = ahp_ypol(scan,mjd0,opt,T,c,rad2mas,obs_mjd);
A(7).sm = Apwypol; H(7).sm = Hypol; Ph(7).sm = Phypol; och(7).sv = oc_hypol;

% dUT1 (piecewise) - ahp_dut1
[Apwdut1,T,Hdut1,Phdut1,oc_hdut1] = ahp_dut1(scan,mjd0,opt,T,c,rad2mas,obs_mjd);
A(8).sm = Apwdut1;  H(8).sm = Hdut1; Ph(8).sm = Phdut1; och(8).sv = oc_hdut1;

% Celestial Pole Offsets dx (DEPS) (piecewise) - ahp_nutdx
[Apwnutdx,T,Hnutdx,Phnutdx,oc_hnutdx] = ahp_nutdx(scan,mjd0,opt,T,c,rad2mas,obs_mjd);
A(9).sm = Apwnutdx; H(9).sm = Hnutdx; Ph(9).sm = Phnutdx; och(9).sv = oc_hnutdx;

% Celestial Pole Offsets dy (DPSI) (piecewise) - ahp_nutdy
[Apwnutdy,T,Hnutdy,Phnutdy,oc_hnutdy] = ahp_nutdy(scan,mjd0,opt,T,c,rad2mas,obs_mjd);
A(10).sm = Apwnutdy; H(10).sm = Hnutdy; Ph(10).sm = Phnutdy; och(10).sv = oc_hnutdy;

clear Apwxpol Apwypol Apwdut1 Apwnutdx Apwnutdy


%  Correction to the scale factor
if opt.est_scale==1
    for n_obs_per_src = 1 : n_observ
        % assigning observationwise partial derivatives
        A(19).sm(n_obs_per_src,1) = temp(n_obs_per_src).pscale;  % [s] 
        H(19).sm(1) = 0;
        Ph(19).sm(1) = 0;
    end
end

% Baseline dependent clock offset
ebsl_bdco=[];
if opt.est_bdco==1
    [A_bdco,H_bdco,Ph_bdco,ebsl_bdco] = abdo_clk(scan,antenna,opt,n_observ,parameter);
    A(20).sm = A_bdco;
    H(20).sm = H_bdco;
    Ph(20).sm = Ph_bdco;
end
parameter.lsmopt.bdco_nrlist = ebsl_bdco;

%%
fprintf('6. FORMING THE CONSTRAIN MATRIX and WEIGHT MATRIX OF CONSTRAINTS\n');

% Constraints as pseudo observations - H, Ph, och
% Constraints for clocks, wet zenith delays, north gradients, east gradients:
[H,Ph,och]=hpoc(H,Ph,och,n_,opt,na,size(A(2).sm,2));

sibling=0;

if sibling==1
    sites = {antenna.name};
    HbHo(1)=strmatch('HOBART12',sites);
    HbHo(2)=strmatch('HOBART26',sites);
    [H,Ph,och,t_] = hpoc_xyz_hbho(H,Ph,och,n_,opt,na,t_,HbHo); % Constraints for antenna coordinates
else
    [H,Ph,och] = hpoc_xyz(H,Ph,och,n_,opt,na); % Constraints for antenna coordinates
end

% common parameters at the co-located sites
common=0;
if common==1
    % ! caution ! current version only supports the case that the numbers of common parameters are same at the co-located sites.
    %             therefore, if there're empty intervals, you cannot apply this option currently.
    % ! caution ! current version only supports the Hobart site.
    % [H,Ph,och,t_] = hpoc_comclk(H,Ph,och,n_,opt,na,c,antenna);    % Constraints for common clock rates at co-located sites % check the estimation setting in the subroutine file
    [H,Ph,och,t_] = hpoc_comzwd(H,Ph,och,n_,opt,na,c,antenna);    % Constraints for common zenith wet delays at co-located sites
    [H,Ph,och,t_] = hpoc_comngr(H,Ph,och,n_,opt,na,c,antenna);    % Constraints for common north gradients at co-located sites
    [H,Ph,och,t_] = hpoc_comegr(H,Ph,och,n_,opt,na,c,antenna);    % Constraints for common east gradients at co-located sites
    % [H,Ph,och,t_,nlt] = hpoc_lt(H,Ph,och,n_,opt,na,t_,antenna);   % local ties as fictitious observations at co-located sites
    nlt=0;
end

if opt.pw_sou == 1  || opt.est_sourceNNR==1
    [H,Ph,och,ts] = hpoc_sources(H,Ph,och,nso,ns_q,opt,ts); % Constraints for source coordinates
end
if opt.pw_sat
	[H, Ph, och] = hpoc_satellites(H, Ph, och, number_pwlo_per_sat, ns_s, opt); % Constraints for satellite coordinates
end


% EXCLUDING THE OFFSET PARAMETERS OF THE SPECIFIED STATIONS' ZWDs, NORTH and EAST GRADIENTS and COORDINATES
[n_,t,A,H,Ph,och] = delparam(na,t,opt,n_,sum_,A,H,Ph,och);

if opt.pw_sou == 1
    % Excluding the specified - fixed sources from the design matrix
    [nso,tso,A,H,Ph,och] = delsource(ns_q,tso,opt,nso,sumso,A,H,Ph,och);
end
% Add the following code if exclusion of single satellites from estimation is needed:
% if opt.pw_sat
% 	[H, Ph, och] = hpoc_satellites(H, Ph, och, number_pwlo_per_sat, ns_s, opt); % Constraints for satellite coordinates
% end


% DELETING THE REFERENCE CLOCK
[A,H,Ph,och,nistat,n_,t] = delref(na,opt,sum_,A,H,Ph,och,n_,t);
opt.fixed_clock = antenna(nistat).name;

% FORMING "A" and "H" ACCORDING TO THE LSM OPTIONS
[A,H,Ph,T,och,n_] = delmodel(opt,A,H,Ph,T,och,n_);


Ablk = [];
Hblk = [];
for i = 1 : length(A)
    dj(i)       = size(A(i).sm,2);
    Ablk        = sparse(horzcat(Ablk,A(i).sm));
    Hblk        = sparse(blkdiag(Hblk,H(i).sm));
    Pobserv     = sparse(blkdiag(Pobserv,Ph(i).sm));
    oc_observ   = vertcat(oc_observ,repmat(och(i).sv,1,numberOfLSMs));
	sum_dj(i+1) = sum_dj(i) + dj(i);
	num_psob(i) = length(find((any(H(i).sm,1) == 1)));      % number of pseudo-observations for each parameter (hana 17May11)
end

if sibling == 1
    if opt.pw_stc == 0
        oc_observ(size(oc_observ,1)-length(och(13).sv)-length(och(14).sv)-length(och(15).sv)+1:size(oc_observ,1)-length(och(14).sv)-length(och(15).sv)-1,:) = [];
        oc_observ(size(oc_observ,1)-length(och(14).sv)-length(och(15).sv)+1:size(oc_observ,1)-length(och(15).sv)-1,:) = [];
        oc_observ(size(oc_observ,1)-length(och(15).sv)+1:size(oc_observ,1)-1,:) = [];
    end
elseif common==1 && (nlt > 0)
    if opt.pw_stc == 0
        oc_observ(size(oc_observ,1)-length(och(13).sv)-length(och(14).sv)-length(och(15).sv)+1:size(oc_observ,1)-length(och(14).sv)-length(och(15).sv)-nlt,:) = [];
        oc_observ(size(oc_observ,1)-length(och(14).sv)-length(och(15).sv)+1:size(oc_observ,1)-length(och(15).sv)-nlt,:) = [];
        oc_observ(size(oc_observ,1)-length(och(15).sv)+1:size(oc_observ,1)-nlt,:) = [];
    end
else
    if opt.pw_stc == 0 % 1 estimate all selected station coordinates as pwl offsets (NNT/NNR is NOT available - fix at least 3 stations)
        oc_observ(size(oc_observ,1)-length(och(13).sv)-length(och(14).sv)-length(och(15).sv)+1:size(oc_observ,1),:) = []; % omc values for rel. constraints are removed again here (from the END of the vector)!
    end
end

% Eliminate elements which are zero:
Hblk(~any(Hblk,2),:)            = [];
Pobserv(~any(Pobserv, 2), :)    = [];
Pobserv(:, ~any(Pobserv, 1))    = [];

pobserv     = Pobserv(1:n_observ,1:n_observ);
mo          = [];
opt_.mo     = mo;
wrms        = [];
opt_.wrms   = wrms;
vTPv        = [];

%% Estimation of the parameters
ess = opt.est_singleses;
if ess == 1 % +hana 10Nov10

    fprintf('7. ESTIMATING THE PARAMETERS WITH LEAST SQUARES\n');
    % -------------------------------------------------------------------------
    % LEAST SQUARES PARAMETER ESTIMATION
    % Normal equation matrix
    A = vertcat(Ablk,Hblk);
    N = A'*Pobserv*A;
    N = full(N);

    % Introducing NNT/NNR for station coordinates
    if opt.stc == 1 % if station coordinates will be estimated
        if opt.nnt_stc == 0
            opt.fixed_station = 'all of them';
            fprintf('!!! Not all station coordinates are estimated (selected are fixed)!!!\n');
        elseif opt.nnt_stc == 1
            if sum([opt.stat.nnt_inc])==1 && sum([opt.stat.nnr_inc])==1
                opt.fixed_station = 'none';
            else
                opt.fixed_station = 'some of them';
            end
            [N] = helmert(n_, na, xo, yo, zo, opt, sum_dj, N);
        end
    end

    % Introducing NNR for source coordinates
    if opt.est_sourceNNR==1
        fprintf('!!! NNR condition is introduced to matrix N for source coordinates!!!\n');
        if sum([opt.source.nnr_inc]) < 3
            fprintf(1,'Not enough sources for NNR condition. All sources are used for NNR instead\n');
            [opt.source.nnr_inc] = deal(1);
        end
%         
%         
%         
%         for iSrc = 1:length(sources.q)
%             if sources.q(iSrc).numobs < 20
%                 opt.source(iSrc).nnr_inc = 0;
%                 fprintf('Removing source %s from datum\n',opt.source(iSrc).name)
%             end
%         end
%         
        [N] = nnr_cond(ns_q,ra,de,opt,sum_dj,N);
    end

    fprintf('clock %s is selected as the ref. clock for the main solution\n',antenna(nistat).name);

    % change by D. Mayer 19 Jun 2015
    % Condition number of the Normal equation matrix
    % [unitless], close to 1 indicates a well-conditioned matrix
    condn = cond(N);
    [opt.condn] = condn;

    % Applying LS
    Qxx = inv(N); % [cm2] & [mas2]
    Qxx = Qxx(1:sum_dj(length(sum_dj)),1:sum_dj(length(sum_dj)));


    n = A'*Pobserv*oc_observ; % [1/cm] & [1/mas]
    x = Qxx*n; % (16) x_c.......... [cm] & [mas]

    % -------------------------------------------------------------------------
    % ACCURACY CRITERIA
    v = A*x - oc_observ; % RESIDUALS OF THE OBSERVATIONS -------01 [cm]
    v_real = v(1:n_observ,:);

    if opt.bsldep ==1 % Do you want to apply baseline dependent weights?
        fprintf('VERSION WITH BASELINE DEPENDENT WEIGHTS\n')

        %20140708 Baseline dependent weighting coded by Minttu Uunila
        [var_bas,Pobserv_bas]=calc_bas_weights(scan,v_real,antenna);%
        Pobserv = sparse(blkdiag(Pobserv_bas,Pobserv(n_observ+1:end,n_observ+1:end)));

        fprintf('7.1 ESTIMATING THE PARAMETERS WITH LEAST SQUARES - 2nd run\n');
        % -------------------------------------------------------------------------
        % LEAST SQUARES PARAMETER ESTIMATION
        % Normal equation matrix
        %A = vertcat(Ablk,Hblk);
        N = A'*Pobserv*A;
        N = full(N);
        % Condition number of the Normal equation matrix
        % [unitless], close to 1 indicates a well-conditioned matrix


        % Introducing NNT/NNR for station coordinates
        if opt.stc == 1 % if station coordinates will be estimated
            if opt.nnt_stc == 0
                opt.fixed_station = 'all of them';
                fprintf('!!! all station coordinates are not estimated (selected as fixed)!!!\n');
            elseif opt.nnt_stc == 1
                if sum([opt.stat.nnt_inc])==1 && sum([opt.stat.nnr_inc])==1
                    opt.fixed_station = 'none';
                else
                    opt.fixed_station = 'some of them';
                end
                [N] = helmert(n_,na,xo,yo,zo,opt,sum_dj,N);
            end
        end

        % Introducing NNR for source coordinates
        if opt.est_sourceNNR==1
            fprintf('!!! NNR condition is introduced to matrix N for source coordinates!!!\n');
            [N] = nnr_cond(ns_q,ra,de,opt,sum_dj,N);
        end

        fprintf('clock %s is selected as the ref. clock for the main solution\n',antenna(nistat).name);

        condn = cond(N);
        [opt.condn] = condn;

        % Applying LS
        Qxx = inv(N); % [cm2] & [mas2]
        Qxx = Qxx(1:sum_dj(length(sum_dj)),1:sum_dj(length(sum_dj)));


        n = A'*Pobserv*oc_observ; % [1/cm] & [1/mas]
        x = Qxx*n; % (16) x_c.......... [cm] & [mas]

        % -------------------------------------------------------------------------
        % ACCURACY CRITERIA
        v = A*x - oc_observ; % RESIDUALS OF THE OBSERVATIONS -------01 [cm]
        v_real = v(1:n_observ,1);
    end


    %mo = sqrt((v_real'*pobserv*v_real)/(nobserv-length(x)));  % RMS [1]
    vTPv=diag(v'*Pobserv*v)';
    mo = sqrt(vTPv/(n_observ+size(Hblk,1)-size(x,1)));  % RMS [2] SINEX V2.01 Appendix(2) Eq.(20)
    opt.mo = mo;

    vTPv_real=diag(v_real'*Pobserv(1:n_observ,1:n_observ)*v_real);     %% numerator of wrms
    weightsum=sum(sum(Pobserv(1:n_observ,1:n_observ)));          %% denominator of wrms
    wrms=sqrt(vTPv_real/weightsum);                            %% wrms of post-fit residual
    opt.wrms = wrms;

    mi = repmat(mo,length(Qxx),1).*repmat(sqrt(diag(Qxx)),1,numberOfLSMs); % std. dev. of estimated parameters [cm,mas]
    
    % Covariance matrix a posteriori, Cvv = sigma_0^2 * Qvv, Qvv=Qll-Q~ll, Qll= inv(P), Q~ll=A*Qxx*A'
    Qll = inv(Pobserv(1:n_observ,1:n_observ));
    Qlldach = A(1:n_observ,:)*Qxx*A(1:n_observ,:)';
    Qvv = Qll-Qlldach;
    Cvv=mo^2*Qvv; 
    sigma_residuals_aposteriori = sqrt(diag(Cvv));

    % Outlier test begins here
    % DETECTING OUTLIER
    if opt.basic_outlier == 1
       qll = diag(inv(pobserv));
       AQh=Ablk*sqrtm(Qxx);
    end
    if opt.simple_outlier == 1
        qvv = 1;
    end

    count_out = 0; out_v = [];
    if opt.basic_outlier == 1 || opt.simple_outlier == 1
       for v_i = 1 : size(v_real,1)
          if opt.basic_outlier == 1
            qvv = qll(v_i) - AQh(v_i,:)*AQh(v_i,:)';
          end
          if abs(v_real(v_i)) >= opt.par_outlier*mo*sqrt(qvv)
             count_out = count_out + 1;
             out_v(count_out) = v_i; % out_v is the vector that contains all outliers in oc_observ
          end
       end
       fprintf('----------\n');
       fprintf('detecting outliers:\n');
       fprintf('num. of outliers: %2d \n',count_out);
    else
       fprintf('----------\n');
       fprintf('outlier detection test was not applied!\n')
       fprintf('----------\n');
    end

    % ##### Write outliers to ASCII file: #####
    if ~isempty(out_v)
        out = write_outlier_file(out_v, antenna, scan, parameter, sources);
    end

    % ##### write residuals/outlier info to new variable #####
    % if the res file already exist, ie if it was written in first solution
    if numberOfLSMs == 1
        resFilename = ['../DATA/LEVEL3/', dirpth, '/res_', parameter.session_name, '.mat'];
        
        if exist(resFilename, 'file') && opt.first
            load(resFilename);
        else  % else: first solution was not run - 'res' variable does not exist

            % Station names:
            res.allStatNames = {antenna.name};
            % Quasar names:
            if isfield(sources, 'q')
                if isfield(sources.q, 'name')
                    res.allSourceNames={sources.q.name};
                else
                    res.allSourceNames = {};
                end
            else
                res.allSourceNames = {};
            end
            % Satellite names:
            if isfield(sources, 's')
                if isfield(sources.s, 'name')
                    res.allSatelliteNames={sources.s.name};
                else
                    res.allSatelliteNames = {};
                end
            else
                res.allSatelliteNames = {};
            end

            lengthOfScans       = cellfun(@length, {scan.obs});
            mjdOfObs            = zeros(sum(lengthOfScans), 1);
            res.baselineOfObs   = zeros(sum(lengthOfScans),2);
            res.source          = zeros(sum(lengthOfScans),1);

            runningInd = 1;
            for iScan = 1 : size(scan, 2)

                % source index of current scan
                res.source(runningInd:runningInd+lengthOfScans(iScan)-1) = repmat(scan(iScan).iso, lengthOfScans(iScan),1);

                % obs type
                res.obs_type(runningInd:runningInd+lengthOfScans(iScan)-1) = repmat(scan(iScan).obs_type, lengthOfScans(iScan),1);

                % get station names for each observation of this scan
                res.baselineOfObs(runningInd:runningInd+lengthOfScans(iScan)-1, :) = [[scan(iScan).obs.i1]', [scan(iScan).obs.i2]'];
                for iObs=1:lengthOfScans(iScan)
                    mjdOfObs(runningInd,1)=scan(iScan).mjd;
                    runningInd=runningInd+1;
                end
            end
            res.obs_type = res.obs_type';

            % save mjd/baselines if this was not done in first solution
            res.mjd = mjdOfObs;
        end

        % main solution residuals and outliers should be saved in any case
        res.mainVal = v_real;
        res.outlier = out_v;
        res.sigma_from_fringe_fitting = mi_observ;
        res.sigma_residuals_aposteriori = sigma_residuals_aposteriori;


        if ~isempty(dirpth) && ~exist(['../DATA/LEVEL3/',dirpth])
            mkdir(['../DATA/LEVEL3/',dirpth])
        end
        % save file (again)
        save(resFilename, 'res');
    end
    



    % -------------------------------------------------------------------------
    % Dividing the x vector in to sub-vectors
    CQ = mat2cell(Qxx,dj,dj);
    cx = mat2cell(x,dj);
    opt.total_est = sum_dj;
    fprintf('---------------------------------------------------------\n');
    fprintf('total clock offsets:                          %4d\n',dj(1));
    fprintf('total rate and quad. terms of clock funct.:   %4d\n',dj(2));
    fprintf('total zenith wet delay offsets:               %4d\n',dj(3));
    fprintf('total tropo. north gradients:                 %4d\n',dj(4));
    fprintf('total tropo. east gradients:                  %4d\n',dj(5));
    fprintf('total pole coor. (x-pol) offsets:             %4d\n',dj(6));
    fprintf('total pole coor. (y-pol) offsets:             %4d\n',dj(7));
    fprintf('total dUT1 offsets:                           %4d\n',dj(8));
    fprintf('total celestial pole (nutation dx) offsets:   %4d\n',dj(9));
    fprintf('total celestial pole (nutation dy) offsets:   %4d\n',dj(10));
    fprintf('total right ascension offsets of sources :    %4d\n',dj(11));
    fprintf('total declination offsets of sources :        %4d\n',dj(12));
    fprintf('antenna coor. dx offsets:                     %4d\n',dj(13));
    fprintf('antenna coor. dy offsets:                     %4d\n',dj(14));
    fprintf('antenna coor. dz offsets:                     %4d\n',dj(15));
    if logical(opt.pw_sat)
        fprintf('satellite pos. 1 offsets:                     %4d\n',dj(16));
        fprintf('satellite pos. 2 offsets:                     %4d\n',dj(17));
        fprintf('satellite pos. 3 offsets:                     %4d\n',dj(18));
    end
    if logical(opt.est_scale)
        fprintf('scale parameter:                              %4d\n',dj(19));
    end
    if logical(opt.est_bdco)
        fprintf('total baseline dependent clock offsets:       %4d\n',dj(20));
    end

    fprintf('---------------------------------------------------------\n');
    fprintf('total number of estimated parameters:         %4d\n',sum_dj(end));
    fprintf('---------------------------------------------------------\n');

    
    [x_] = splitx(x,first_solution,mi,na,sum_dj,n_,mjd0,mjd1,t,T,opt,antenna,ns_q,nso,tso,ess, ns_s, number_pwlo_per_sat, ebsl_bdco);
    x_.mo = mo;
    x_.mo_first = first_solution.mo;
    x_.units.mo = 'chi of main solution vTPv/degOfFreedom [] (NOT SQUARED!)';
    x_.units.mo_first = 'chi of main solution vTPv/degOfFreedom [] (NOT SQUARED!)';
    x_.wrms = wrms;
    x_.units.m02 = 'WRMS of post-fit residuals sqrt(v_realTPv_real/sumOfWeights) [cm]';
    x_.nobs = n_observ;
    x_.nscans = n_scan;

    res.mo = mo;
    res.mo_first = first_solution.mo;
    res.units.mo = 'chi of main solution vTPv/degOfFreedom [] (NOT SQUARED!)';
    res.units.mo_first = 'chi of main solution vTPv/degOfFreedom [] (NOT SQUARED!)';
    res.wrms = wrms;
    res.units.m02 = 'WRMS of post-fit residuals sqrt(v_realTPv_real/sumOfWeights) [cm]';



    [atpa_.mat] = N;
    % [atpl_.vec] = n;
    [opt_] = opt;
    % save files

    fprintf('----------\n');
    test_significance(x_,opt_,5);

    % Save the "cleaned" VieVS structures (consistent with the results in x_ and res_):
    fprintf('Estimated parameters are saved as ../DATA/LEVEL3/%s/x_%s.mat\n',dirpth,parameter.session_name);
    save(['../DATA/LEVEL3/',dirpth,'/x_',parameter.session_name,'.mat'],'x_');

    fprintf('Estimation options are saved as ../DATA/LEVEL3/%s/opt_%s.mat\n',dirpth,parameter.session_name);
    save(['../DATA/LEVEL3/',dirpth,'/opt_',parameter.session_name,'.mat'],'opt_');

    fprintf('normal equation matrix is saved as ../DATA/LEVEL3/%s/atpa_%s.mat\n',dirpth,parameter.session_name);
    save(['../DATA/LEVEL3/',dirpth,'/atpa_',parameter.session_name,'.mat'],'atpa_');

%     fprintf('right hand side vector is saved as ../VieVS/DATA/LEVEL3/%s/atpl_%s.mat\n',dirpth,parameter.session_name);
%     save(['../DATA/LEVEL3/',dirpth,'/atpl_',parameter.session_name,'.mat'],'atpl_');

    fprintf('Residuals are saved as ../DATA/LEVEL3/%s/res_%s.mat\n',dirpth,parameter.session_name);
    save(['../DATA/LEVEL3/',dirpth,'/res_',parameter.session_name,'.mat'],'res');

end

if numberOfLSMs == 1
    % Save the "cleaned" VieVS structures (consistent withe the results in x_ and res_):
    fprintf('Save VieVS data structures (outliers and OPT file options applied!)\n');
    fprintf('  ../DATA/LEVEL3/%s/%s_antenna.mat\n',dirpth,parameter.session_name);
    save(['../DATA/LEVEL3/',dirpth,'/',parameter.session_name,'_antenna.mat'],'antenna');
    fprintf('  ../DATA/LEVEL3/%s/%s_source.mat\n',dirpth,parameter.session_name);
    save(['../DATA/LEVEL3/',dirpth,'/',parameter.session_name,'_sources.mat'],'sources');
    fprintf('  ../DATA/LEVEL3/%s/%s_parameter.mat\n',dirpth,parameter.session_name);
    save(['../DATA/LEVEL3/',dirpth,'/',parameter.session_name,'_parameter.mat'],'parameter');
    fprintf('  ../DATA/LEVEL3/%s/%s_scan.mat\n',dirpth,parameter.session_name);
    save(['../DATA/LEVEL3/',dirpth,'/',parameter.session_name,'_scan.mat'],'scan');
end
[opt_] = opt;

%%

% -------------------------------------------------------------------------
% +Boehm 2 July 2009
% FOR GLOBAL SOLUTION (CLOCK IS FIXED (or NNT) and WITH CONSTRAINTS BETWEEN OFFSETS)

% +hana 05Oct10
% and for N,b for SINEX output (because of source coordinates)

% +hana 10Nov10
% In case, that the single-session solution wasn't applied,
% in the x_ variable will be written only information about
% columns in the N and b and the time information

if opt.global_solve == 1 || opt.ascii_snx ==1 % +hana 05Oct10
    if opt.est_scale == 1 % estimate scale
        fprintf('\nWe are sorry, but currently it is not possible to create the sinex file or glob data if the scale is estimated!\n\n')
        return
    end


    [x_] = splitx(x,first_solution,mi,na,sum_dj,n_,mjd0,mjd1,t,T,opt,antenna,ns_q,nso,tso,ess, ns_s, number_pwlo_per_sat, ebsl_bdco);

    glob_dj = dj;

    x_.source=[''];
    if opt.est_source == 1
        % +hana 21 Jun 2010
        % Assign the common and IERS name of the sources
        for isou = 1 : length(opt.source)
            x_.source(isou).name = opt.source(isou).name;
            x_.source(isou).IERSname = opt.source(isou).IERSname; % hana 24Mar13
        end
        % -hana

        glob_dj(length(glob_dj)+1) = size(A_ra_glob,2);
        glob_dj(length(glob_dj)+1) = size(A_de_glob,2);
    elseif opt.est_source == 0
        glob_dj(length(glob_dj)+1) = 0;
        glob_dj(length(glob_dj)+1) = 0;
    end

    if opt.est_vel == 1
        glob_dj(length(glob_dj)+1) = size(A_vx,2);
        glob_dj(length(glob_dj)+1) = size(A_vy,2);
        glob_dj(length(glob_dj)+1) = size(A_vz,2);
    elseif opt.est_vel == 0
        A_vx = []; A_vy = []; A_vz = [];
        glob_dj(length(glob_dj)+1) = 0;
        glob_dj(length(glob_dj)+1) = 0;
        glob_dj(length(glob_dj)+1) = 0;
    end

    if opt.est_AO == 1
        glob_dj(length(glob_dj)+1) = size(A_AO,2);
    else
        glob_dj(length(glob_dj)+1) = 0;
    end



    % amplitudes of seasonal variations in station position
    glob_dj(length(glob_dj)+1) = size(A_Acr,2);
    glob_dj(length(glob_dj)+1) = size(A_Ace,2);
    glob_dj(length(glob_dj)+1) = size(A_Acn,2);
    glob_dj(length(glob_dj)+1) = size(A_Asr,2);
    glob_dj(length(glob_dj)+1) = size(A_Ase,2);
    glob_dj(length(glob_dj)+1) = size(A_Asn,2);

    % Love number for pole tide
    A_hpole=[];
    if opt.est_hpole==1
        for n_obs_per_src = 1 : n_observ
            A_hpole(n_obs_per_src,1) = temp(n_obs_per_src).phpole(1); % assigning observationwise partial derivatives
        end
    end
    glob_dj(length(glob_dj)+1) = size(A_hpole,2);

    % Shida number for pole tide
    A_lpole=[];
    if opt.est_lpole==1
        for n_obs_per_src = 1 : n_observ
            A_lpole(n_obs_per_src,1) = temp(n_obs_per_src).plpole(1); % assigning observationwise partial derivatives
        end
    end
    glob_dj(length(glob_dj)+1) = size(A_lpole,2);

     % APL RG
    if opt.est_rg == 1
        glob_dj(length(glob_dj)+1) = size(A_rg,2);
    else
        glob_dj(length(glob_dj)+1) = 0;
    end


    % FCN period from nutation
    A_FCNnut=[]; % period
    if opt.est_FCNnut==1
        for n_obs_per_src = 1 : n_observ
            A_FCNnut(n_obs_per_src,1) = temp(n_obs_per_src).pFCNnut(1)*c*100; %[cm*day]
        end
    end
    glob_dj(length(glob_dj)+1) = size(A_FCNnut,2);

    % Amplitudes for FCN
    A_FCNnutAc=[];
    if opt.est_FCNnutAmpl==1
        for n_obs_per_src = 1 : n_observ
            A_FCNnutAc(n_obs_per_src,1) = temp(n_obs_per_src).pFCNnutAc(1)*c*100/rad2mas; % [cm/mas]
        end
    end
    glob_dj(length(glob_dj)+1) = size(A_FCNnutAc,2);

    A_FCNnutAs=[];
    if opt.est_FCNnutAmpl==1
        for n_obs_per_src = 1 : n_observ
            A_FCNnutAs(n_obs_per_src,1) = temp(n_obs_per_src).pFCNnutAs(1)*c*100/rad2mas; % [cm/mas]
        end
    end
    glob_dj(length(glob_dj)+1) = size(A_FCNnutAs,2);



    % Tidal ERP variations
    A_tidap = []; A_tidbp = []; A_tidam = []; A_tidbm =[];
    A_tidutc = []; A_tiduts = [];
    if opt.est_erpsp == 1
        [A_tidap,A_tidbp,A_tidam,A_tidbm] = a_tidpm(temp,n_observ,c,rad2mas);
        glob_dj(length(glob_dj)+1) = size(A_tidap,2);
        glob_dj(length(glob_dj)+1) = size(A_tidbp,2);
        glob_dj(length(glob_dj)+1) = size(A_tidam,2);
        glob_dj(length(glob_dj)+1) = size(A_tidbm,2);
        [A_tidutc,A_tiduts] = a_tidut(temp,n_observ,c,rad2mas);
        glob_dj(length(glob_dj)+1) = size(A_tidutc,2);
        glob_dj(length(glob_dj)+1) = size(A_tiduts,2);
    end


    % external functions: a_love.m, a_shida.m
    A_love = []; A_shida = [];
    % a_love & a_shida
    if opt.est_love == 1
        [A_love] = a_love(temp,n_observ);
        glob_dj(length(glob_dj)+1) = size(A_love,2);
    elseif opt.est_love == 0
        glob_dj(length(glob_dj)+1) = 0;
    end

    if opt.est_shida == 1
    	[A_shida] = a_shida(temp,n_observ);
        glob_dj(length(glob_dj)+1) = size(A_shida,2);
    elseif opt.est_shida == 0
        glob_dj(length(glob_dj)+1) = 0;
    end

    % FCN period
    A_FCN=[];
    if opt.est_FCNset==1
        for n_obs_per_src = 1 : n_observ
            A_FCN(n_obs_per_src,1) = temp(n_obs_per_src).pFCN(1); % assigning observationwise partial derivatives
        end
    end
    glob_dj(length(glob_dj)+1) = size(A_FCN,2);
    % acceleration of SSB
    A_accSSB=[];
    if opt.est_accSSB==1
        for n_obs_per_src = 1 : n_observ
            % assigning observationwise partial derivatives
            A_accSSB(n_obs_per_src,1:3) = temp(n_obs_per_src).pacc;    % [sec^3/cm]     (don't use: .*(100*c)^3; %sec^3/cm --> cm^2)
        end
    end
    glob_dj(length(glob_dj)+1) = size(A_accSSB,2);

    % velocities of sources added
    glob_dj(length(glob_dj)+1) = size(A_vra_glob,2);
    glob_dj(length(glob_dj)+1) = size(A_vde_glob,2);


    % Gamma parameter
    A_gamma=[];
    if opt.est_gamma==1
        for n_obs_per_src = 1 : n_observ
            % assigning observationwise partial derivatives
            A_gamma(n_obs_per_src,1) = temp(n_obs_per_src).pGamma *c*100;  % [cm]
        end
    end
    glob_dj(length(glob_dj)+1) = size(A_gamma,2);


    oc_observ_real = oc_observ(1:n_observ);


    A_add = horzcat(A_ra_glob,A_de_glob,A_vx,A_vy,A_vz,A_AO,...
            A_Acr,A_Ace,A_Acn,A_Asr,A_Ase,A_Asn,...
            A_hpole,A_lpole,A_rg,A_FCNnut,A_FCNnutAc,A_FCNnutAs,...
            A_tidap,A_tidbp,A_tidam,A_tidbm,A_tidutc,A_tiduts,...
            A_love,A_shida,A_FCN,A_accSSB,A_vra_glob,A_vde_glob,A_gamma);
        %
        % N.B. x_.col_ order should be consistent with concatenation in the matrix
        %
    H_add = zeros(size(Hblk,1),size(A_add,2));

    A_global = horzcat(Ablk,A_add);
    H_global = horzcat(Hblk,H_add);
    A_glob = vertcat(A_global,H_global);
    oc_glob = vertcat(oc_observ_real,zeros(size(Hblk,1),1));



    N_global = A_glob' * Pobserv * A_glob;
    b_global = A_glob' * Pobserv * oc_glob;

    N_global = full(N_global);
    b_global = full(b_global);

    sum_glob_dj(1) = 0;
    for i = 1 : length(glob_dj)
        sum_glob_dj(i+1) = sum_glob_dj(i) + glob_dj(i);
    end

    IDglobdj=number_of_estimated_parameters+1;

    x_.col_soura = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_soude = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_vx    = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_vy    = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_vz    = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_AO    = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;

    nAcr=sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    nAce=sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    nAcn=sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    nAsr=sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    nAse=sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    nAsn=sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;

    % number of seasonal waves (station position variation) : nsewa
    for i=1:na % number of antennas
        x_.col_ssp_Acr(i).col=nAcr(nsewa*i-nsewa+1:nsewa*i);
        x_.col_ssp_Ace(i).col=nAce(nsewa*i-nsewa+1:nsewa*i);
        x_.col_ssp_Acn(i).col=nAcn(nsewa*i-nsewa+1:nsewa*i);
        x_.col_ssp_Asr(i).col=nAsr(nsewa*i-nsewa+1:nsewa*i);
        x_.col_ssp_Ase(i).col=nAse(nsewa*i-nsewa+1:nsewa*i);
        x_.col_ssp_Asn(i).col=nAsn(nsewa*i-nsewa+1:nsewa*i);
    end

    x_.col_hpole   = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_lpole   = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;

    x_.col_rg = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;

    % FCN from nutation
    x_.col_FCNnut     = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_FCNnutAc   = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_FCNnutAs   = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;


    % tidal ERP variations
    if opt.est_erpsp == 1
    x_.col_tidap = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_tidbp = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_tidam = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_tidbm = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_tidutc = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_tiduts = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    end

    x_.col_love  = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_shida = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;

    x_.col_FCNset= sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;

    x_.col_accSSB= sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_souvra= sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_souvde= sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_gamma = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;


    if opt.global_solve == 1
        allant= [antenna.name];
        for i=1:size(allant,2)/8
             ant(i,1:8)= allant((i*8)-7:i*8);
        end
        an_glob.name = ant;
        an_glob.in_trf = [antenna.in_trf];
        an_glob.x=[antenna.x];
        an_glob.y=[antenna.y];
        an_glob.z=[antenna.z];
        an_glob.vx=[antenna.vx];
        an_glob.vy=[antenna.vy];
        an_glob.vz=[antenna.vz];
        an_glob.epoch=[antenna.epoch];
        an_glob.start=[antenna.start];
        an_glob.end=[antenna.end];
        an_glob.firstscan_mjd=mjd1;
        % a priori APL Reg Coeff
        for i=1:length(antenna)
            an_glob.aplrg(i,:) = antenna(i).crg; % [ref.pres, a priori RgC]
        end


%         if ~isempty(dirpthL2) & ~exist(['../DATA/LEVEL2/',dirpthL2])
%             mkdir(['../DATA/LEVEL2/',dirpthL2])
%         end

        fprintf('Data for GLOBAL SOLUTION is saved as ../VieVS/DATA/LEVEL2/%s/%s_an_glob.mat\n',dirpthL2,parameter.session_name);
        fprintf('Data for GLOBAL SOLUTION is saved as ../VieVS/DATA/LEVEL2/%s/%s_par_glob.mat\n',dirpthL2,parameter.session_name);
        fprintf('Data for GLOBAL SOLUTION is saved as ../VieVS/DATA/LEVEL2/%s/%s_Nb_glob.mat\n',dirpthL2,parameter.session_name);

        % hana 26March13 - needed for Vie_GLOB
        if ess == 0
            opt_.source = opt.source;
        end

        opt_.trf=parameter.vie_init.trf;
        opt_.crf=parameter.vie_init.crf;

        opt_.total_obs = n_observ;
        % Sources: currently the sources are stored twice in the N_global matrix!
        % The H(11) and H(12) are included in N_global. - Needs further
        % improvement of the code!!!! This is just a preliminary version.
        % (In Vie_GLOB the respective columns/rows are deleted anyway.)
        % Hana 11/2016
        opt_.nconstr = size(Hblk,1)-size(H(11).sm,1)-size(H(12).sm,1); % minus abs. constr for RA and De if sources estimated session-wise
        opt_.lTPl = oc_observ_real'*pobserv*oc_observ_real; %l'Pl for global solution

        glob1.an = an_glob;
        glob2.x = x_;
        glob2.opt = opt_;
        glob3.N = N_global;
        glob3.b = b_global;

%         if ~exist(['../DATA/LEVEL2/',dirpthL2])
%             mkdir(['../DATA/LEVEL2/',dirpthL2])
%         end

        if ~isempty(dirpthL2) && ~exist(['../DATA/LEVEL2/',dirpthL2])
            mkdir(['../DATA/LEVEL2/',dirpthL2])
        end

        save(['../DATA/LEVEL2/',dirpthL2,'/',parameter.session_name,'_an_glob.mat'],'glob1');
        save(['../DATA/LEVEL2/',dirpthL2,'/',parameter.session_name,'_par_glob.mat'],'glob2');
        save(['../DATA/LEVEL2/',dirpthL2,'/',parameter.session_name,'_Nb_glob.mat'],'glob3');
    end
end

% -------------------------------------------------------------------------
% +hana 16 May 2011
% N and b FOR SINEX OUTPUT from the global matrix (CLOCK IS FIXED (or NNT) and WITH CONSTRAINTS BETWEEN OFFSETS)
% Reduce: always clock
%         optional zwd, trop. gradients
% Keep in sinex file: station coordinates, EOP
% Source coordinates can be writen into sinex file in the N-matrix, but
% cannot be estimated from the single-session solution

% Constraints of parameters which will be written into SINEX file are
% removed!

if opt.ascii_snx == 1
    clear N11 N12 N21 N22 b1 b2 col_red col_est


    outsnx=opt.outsnx;
    [col_red, col_est] = snx_split(x_,outsnx);

    % remove constraints of parameters which will be written into SINEX file
    sA=[];
    sA=size(A_global,1);
    A_glob(sA+1:end,col_est)=0;

    Nsnx=A_glob' * Pobserv * A_glob;
    bsnx=A_glob' * Pobserv * oc_glob;

    Nsnx = full(Nsnx);
    bsnx = full(bsnx);

    % REDUCTION
    N11=Nsnx(col_est,col_est);
    N12=Nsnx(col_est,col_red);
    N21=Nsnx(col_red,col_est);
    N22=Nsnx(col_red,col_red);

    b1=bsnx(col_est);
    b2=bsnx(col_red);

    N_sinex = N11 - N12*inv(N22)*N21;
    b_sinex = b1 - N12*inv(N22)*b2;


    % Reduction of lTPl:
    lTPl = oc_observ_real'*pobserv*oc_observ_real;
    lTPlreduc = lTPl - (b2'*inv(N22)*b2);

    % number of pseudo-observations which were reduced
    % they are added to the number of real observation and the sum is
    % written into the sinex file
    nconstr_red=0;
    if outsnx.clk==0; nconstr_red = nconstr_red + num_psob(1)+num_psob(2); end
    if outsnx.zwd==0; nconstr_red = nconstr_red + num_psob(3); end
    if outsnx.tgr==0; nconstr_red = nconstr_red + sum(num_psob(4:5)); end
    if outsnx.eop==0; nconstr_red = nconstr_red + sum(num_psob(6:10)); end
    if outsnx.xyz==0; nconstr_red = nconstr_red + sum(num_psob(13:15)); end


    total_est = sum_dj;
    real_obs = n_observ;
    all_obs = real_obs + nconstr_red;


    % Save info about columns
    col_sinex = snx_newcol(col_est, x_, antenna, outsnx);
    % Save info about statistic
    col_sinex.lTPlreduc = lTPlreduc;
    col_sinex.nr_unknowns = total_est(end);
    col_sinex.nr_obs = real_obs; % write only the real observations into sinex
    col_sinex.vTPv = vTPv;
    col_sinex.varfac = mo^2; % hana 24 Apr 2013
    col_sinex.outsnx = outsnx;

    % Change units of the N matrix and b vector !!!
    [N_sinex, b_sinex]=snx_changeunits(N_sinex,b_sinex,col_sinex,outsnx);


%%% Saving of variables commented by S. Boehm (slows the processing and uses a lot of disk space)
%     if exist(['../DATA/LEVEL3/',dirpth,'/SINEX/'])~=7
%         mkdir(['../DATA/LEVEL3/',dirpth,'/SINEX/'])
%     end
% 
%     fprintf('\nReduced N and b for SINEX output are saved in ../VieVS/DATA/LEVEL3/%s/SINEX/ \n',dirpth);
%     save(['../DATA/LEVEL3/',dirpth,'/SINEX/N_sinex_',parameter.session_name,'.mat'],'N_sinex');
%     save(['../DATA/LEVEL3/',dirpth,'/SINEX/b_sinex_',parameter.session_name,'.mat'],'b_sinex');
%     save(['../DATA/LEVEL3/',dirpth,'/SINEX/col_sinex_',parameter.session_name,'.mat'],'col_sinex');

    % create an ascii sinex file in DATA/SNX
    fprintf('\nWriting SINEX file ... \n');
    write_sinex_vievs(parameter.session_name, [dirpth '/'], outsnx.firstname, outsnx.lastname, outsnx.email, outsnx.suffix,...
                      N_sinex,b_sinex,col_sinex); % vars are loaded to the function directly -> saving uses disk space needlessly
end


% Command window output of resulting statistics
if numberOfLSMs == 1
    fprintf('\n');
    fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')
    fprintf('chi-squared of main solution vTPv/degOfFreedom: %.3f\n',mo.^2);
    fprintf('WRMS of post-fit residuals sqrt(v_realTPv_real/sumOfWeights): %.3f cm (%.3f ps)\n',wrms,wrms*100/2.99792458);  % for post-fit residual check
    fprintf('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n')
    fprintf('\n')
end

% -------------------------------------------------------------------------
fprintf('vie_lsm is successfully finished after %.2f seconds!\n',toc);
fprintf('\n\n')
