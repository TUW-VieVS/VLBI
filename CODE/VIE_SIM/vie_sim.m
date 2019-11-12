% ************************************************************************
%   Description:
%   In vie_sim simulated group delay observables are generated which
%   consist of three important stochastic error sources in VLBI: the wet
%   troposphere delay, station clock and thermal noise. vie_sim simulates
%   the (o-c) values and writes (o-c) + c to NGS files which can then be
%   used in the analysis as real NGS files.
%
%   References:
%    - Petrachenko et al. 2010 Design Aspects of the VLBI2010 System
%
%   Input:
%      session   -    string with session name
%      dirpt     -    name of the subdirectory
%
%   Output:
%      ---
%
%   External calls:
%      read_turbparam.m
%      sim_clk.m
%      sim_swd.m
%      sim_wn.m
%      vievs2ngs.m
%
%   Coded for VieVS:
%   July 2010 by Andrea Pany
%
%   Revision
%   18 Apr 2011 by Tobias Nilsson: Changed input variables
%   26 Apr 2011 by Andrea Pany: antenna array added to the function call in line #357
%   06 Jun 2011 by Andrea Pany: changed input variables
%   21 Jun 2011 by Andrea Pany: added antenna name to the sim output
%   24 Jan 2014 by Lucia Plank: wn from input/NGS-file corrected
%   05 Feb 2014 by Lucia Plank: source structure simulation added
%   11 Apr 2014 by Lucia Plank: source structure delay added to sim.mat
%                  output. Bug in source structure calculation corrected.
%   12 Oct 2015 by David Mayer: saving in subdirectory is now possible
%   17 Jul 2016 by A. Hellerschmied: Content of "sources" structure (natural sources, quasars) is now stored in the sub-structure "sources.q"
%   19 Oct 2016 by A. Hellerschmied: - Possibility to initialize the random number generator with a seed defined in the GUI, or with the "shuffle" option.
%                                    - Function changed to support the simulation of satellite observations
%                                        - Separate white noise parameter (simparam.wn_sat) added for sat. obs.
%   21 Oct 2016 by A. Hellerschmied: Simulation of observatons of satellites added
%   21 Nov 2016 by M. Schartner: - now saves structures in level1 directory 
% 								 - new return parameters 	
% 								 - new input parameter "parameter" 	
% 								 - standard welcome message
%   30 Nov 2016 by M. Schartner: small bugfix if you only simulate 1 day
%   27 Jan 2017 by Daniel Landskron: GPT2w changed to GPT3
%   29 Aug 2017 by A. Hellerschmied: bug-fix 
%
% ************************************************************************

function [antenna,scan,sources,session,parameter] = vie_sim(antenna,scan,sources,session,dirpt,parameter)

% Options:
% flag_write_vso = true;

dirpt0 = dirpt;

fprintf('---------------------------------------------------------------\n')
fprintf('|                   Welcome to VIE_SIM!!!!!                   |\n')
fprintf('---------------------------------------------------------------\n\n')
if ~exist('dirpt','var')
    dirpt='../DATA/LEVEL1';
    dirpt1='../DATA/LEVEL4';
else
    dirpt1=['../DATA/LEVEL4/' dirpt];
    dirpt=['../DATA/LEVEL1/' dirpt];
end

if ~exist('../DATA/LEVEL4','dir')
    mkdir('../DATA/LEVEL4');
end

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%  SET PATHS/NAMES OF .MAT-FILES and load files
% disp ' loading files...'
% pfile=[dirpt '/' session] %#ok<NOPRT>
pfile1='../DATA/LEVEL4/';
load ([pfile1 'simparam.mat'],'simparam');
% disp ' data successfully loaded'
% disp '    '
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
fprintf('##### Simulating session %s for %d days #####\n', session, simparam.idays)


% Initialize random number generator:
if simparam.rng_use_seed % Use the seed defined in the GUI to initialize the random number generator
    % Check the seed value from the GUI
    if (simparam.rng_seed <= 0) || (mod(simparam.rng_seed, 1) ~= 0)
        error('The seed for the random number generator defined in the VIE_SIM GUI has to be a nonnegative integer! Current seed: %f\n', simparam.rng_seed);
    end
    rng(simparam.rng_seed);
else % Do not use a seed to initialize the random number generator
    rng('shuffle');
end

% ##### Output options #####
flag_write_ngs_file = simparam.ngs;
flag_write_vso_file = simparam.write_vso;


% number of scans
nscan = length(scan);
% number of antennas participating
nant = length(antenna);

% preallocate structure arrays
azel = struct;
mf = struct;
st = zeros(nant,1);
% extract az, el, and mjd (stationwise)
% loop over all scans
for i = 1:nscan
    % loop over all stations in the scan
    for j = 1:length(scan(i).stat)
        if ~isempty(scan(i).stat(j).az)
            st(j) = st(j) + 1;
            azel(j).az(st(j)) = scan(i).stat(j).az;
            azel(j).el(st(j)) = pi/2 - scan(i).stat(j).zd;
            azel(j).mjd(st(j)) = scan(i).mjd;
            mf(j).mfw(st(j)) = scan(i).stat(j).mfw;
        end
    end
end

% *************************************************************************
% check if simulation parameters are read from file or entered directly
if ~isempty(simparam.turbfile)
    % load turbulence parameterfile
    disp (['loading turbulence parameter file: ../DATA/TURB/',simparam.turbfile])
    disp '    '
    simparam1 = read_turbparam(['../DATA/TURB/',strtrim(simparam.turbfile)],antenna);
    % number of days to simulate
    sim_idays = simparam.idays;
    % first index for the running number of simulated NGS files
    sind = simparam.sind;
    % set ref clock to 0?
    ref0 = simparam.ref0;
    % simulate swd?
    sswd = simparam.sswd;
    % simulate clk?
    sclk = simparam.sclk;
    % simulate wn?
    swn  = simparam.swn;
    % is a turbulence parameter file used?
    tfil_used_falg = 1;
    % simulate soustruc?
    ssou=simparam.ssou;
    if ssou==1
        sou_cat=simparam.sou_cat;
    end
    clear simparam
else
    simparam1 = simparam;
    sim_idays = simparam.idays;
    % first index for the running number of simulated NGS files
    sind = simparam.sind;
    % set reference clock to 0?
    ref0 = simparam.ref0;
    % simulate swd?
    sswd = simparam.sswd;
    % simulate clk?
    sclk = simparam.sclk;
    % simulate wn?
    swn  = simparam.swn;
    % is a turbulence parameterfile used?
    tfil_used_falg = 0;
    % simulate soustruc?
    ssou=simparam.ssou;
    if ssou==1
        sou_cat=simparam.sou_cat;
    end
    clear simparam
end
% *************************************************************************

% pre-allocate
sim = struct;
cov = struct;
for i = 1:nant
    % + 21 Jun 2011, A.Pany
    sim(i).name = antenna(i).name;
    % - A. Pany
    
    
    % ##### Simulate station clock #####
    if (tfil_used_falg == 1) && (sclk == 1)
        disp (['    simulating clock for station ',antenna(i).name])
        sim(i).clk = sim_clk(azel(i).mjd,simparam1(i).sy1,simparam1(i).sy2,sim_idays);
    elseif (tfil_used_falg == 0) && (sclk == 1)
        disp (['    simulating clock for station ',antenna(i).name])
        sim(i).clk = sim_clk(azel(i).mjd,simparam1(1).sy1,simparam1(1).sy2,sim_idays);
    else
        sim(i).clk = zeros(length(azel(i).mjd),sim_idays);
    end
    
    
    % ##### Simulate troposphere slant wet delay #####
    if (tfil_used_falg == 1) && (sswd == 1)
        disp (['    simulating wet delay for station ',antenna(i).name])
        [swd,cv] = sim_swd(azel(i).mjd,azel(i).az,azel(i).el,mf(i).mfw,sim_idays,...
            simparam1(i).dhseg,simparam1(i).dh,simparam1(i).Cn, ...
            simparam1(i).H,simparam1(i).vn,simparam1(i).ve,simparam1(i).wzd0);
    elseif (tfil_used_falg == 0) && (sswd == 1)
        disp (['    simulating wet delay for station ',antenna(i).name])
        [swd,cv] = sim_swd(azel(i).mjd,azel(i).az,azel(i).el,mf(i).mfw,sim_idays,...
            simparam1(1).dhseg,simparam1(1).dh,simparam1(1).Cn, ...
            simparam1(1).H,simparam1(1).vn,simparam1(1).ve,simparam1(1).wzd0);
    else
        swd = zeros(length(azel(i).mjd),sim_idays);
        cv = 0;
    end
    sim(i).swd = swd;
    cov(i).cov = cv;
end


% ##### Simulate White Noise #####
% simulate white noise to be added to the baseline observations
obs = [scan.obs];               % Extract obs. sub-structures
nbslobs = length([obs.obs]);    % Number of baseline observations

% shall white noise be simulated?
if swn == 1
    
    % Check, wich observations types should be simulated:
    if isempty(sources.q)
        flag_sim_quasar_obs = false;
    else
        flag_sim_quasar_obs = true;
    end
    if isempty(sources.s)
        flag_sim_satellite_obs = false;
    else
        flag_sim_satellite_obs = true;
    end
    
    % Get obs inidices for obs types:
    if flag_sim_quasar_obs && ~flag_sim_satellite_obs % Only quasar obs
        ind_q_obs = true(nbslobs, 1);
        ind_s_obs = false(nbslobs, 1);
    elseif ~flag_sim_quasar_obs && flag_sim_satellite_obs % Only sat obs
        ind_s_obs = true(nbslobs, 1);
        ind_q_obs = false(nbslobs, 1);
    else % quasar ans sat obs
        ind_s_obs = false(nbslobs, 1);
        ind_q_obs = false(nbslobs, 1);
        % Loop over scans
        last_obs_of_scan_ind = 0;
        for i_scan = 1 : length(scan) 
            first_obs_of_scan_ind   = last_obs_of_scan_ind + 1;
            last_obs_of_scan_ind    = last_obs_of_scan_ind + length(scan(i_scan).obs);
            switch(scan(i_scan).obs_type) 
                case 'q'
                    ind_q_obs(first_obs_of_scan_ind : last_obs_of_scan_ind) = true;
                case 's'
                    ind_s_obs(first_obs_of_scan_ind : last_obs_of_scan_ind) = true;
            end
        end
    end

    
    % yes
    wnoise = simparam1(1).wn;
    if isfield(simparam1(1), 'wn_sat')
        wnoise_sat  = simparam1(1).wn_sat;
    else
        wnoise_sat  = 0;
        % Check, if wnoise_sat is available in case satellites were
        % observed:
        if flag_sim_satellite_obs
           error('White noise for satellite observations should be simulated, but the parameter "wn_sat" was not defined!'); 
        end
    end
else
    % no
    wnoise = 0;
end

% set reference clock to zero if specified in the gui
if ref0 == 1
    [rr,cc] = size(sim(1).clk);
    rclk = zeros(rr,cc);
    sim(1).clk = rclk;
end

% Init. wn:
wn = zeros(nbslobs,sim_idays);

% if wnoise = 0     -> zero white noise will be added
% if wnoise = 999   -> the actual uncertainties from the NGS file will be used

% Observations of quasars:
if flag_sim_quasar_obs
    switch wnoise
        case 0
            disp '    '
            disp '    not simulating white noise per baseline observation (set to zero for all observations)(quasars)'
            % wn(ind_q_obs, :) = zeros(wnoise, sum(ind_q_obs),sim_idays);
        case 999
            disp '    '
            disp '    take the uncertaintiers from the input file (NGS, etc) as basis for the simulation of white noise (quasars)'
            sig(:,1) = [obs(ind_q_obs).sig]';
            for i_d = 1 : sim_idays
                wn(ind_q_obs, i_d) = randn(sum(ind_q_obs), 1) .* sig;
            end  
        otherwise
            disp '    '
            disp '    simulating white noise per baseline observation (quasars)'
            wn(ind_q_obs, :) = sim_wn(wnoise, sum(ind_q_obs), sim_idays);
    end
end

% Observations of satellites:
if flag_sim_satellite_obs
    switch wnoise_sat
        case 0
            disp '    '
            disp '    not simulating white noise per baseline observation (set to zero for all observations)(satellites)'
            % wn(ind_s_obs, :) = zeros(sum(ind_s_obs), sim_idays);
        case 999
            disp '    '
            disp '    take the uncertaintiers from the input file (NGS, etc) as basis for the simulation of white noise (satellites)'
            sig(:,1) = [obs(ind_s_obs).sig]';
            for i_d = 1 : sim_idays
                wn(ind_s_obs, i_d) = randn(sum(ind_s_obs), 1) .* sig;
            end                        
        otherwise
            disp '    '
            disp '    simulating white noise per baseline observation (satellites)'
            wn(ind_s_obs, :) = sim_wn(wnoise_sat, sum(ind_s_obs), sim_idays);
    end
end

% Save wn values:
sim(1).wn = wn;


% ##### Simulate source structure #####
if ssou == 1
    disp 'sim sou'
    % read in catalogue
    filename_sou_cat=strcat(['../CRF/SOURCE_STRUCTURE_CAT/',sou_cat]);
    [cat_comp.name,cat_comp.flux,cat_comp.maj,cat_comp.min,cat_comp.angle,cat_comp.dRA,cat_comp.dDec]=textread(filename_sou_cat,'%s%f%f%f%f%f%f','delimiter',',');
    sim(1).soustruccat= filename_sou_cat;
end



% ##### Compute simulated baseline delay observables #####

% compute baseline delay observables and replace the real observed and
% computed delay values with the simulated ones: the simulated o-c vector
% for a baseline is computed as [swd(2)+clk(2)] - [swd(1)+clk(1)] + wn(bsl).
% Since with real data the o-c vector is formed in the lsm part, the
% computed delay is set to zero for the simulated scenario.

st = zeros(nant,1);
k = 0;
% simulate or zero input?
if sswd == 0 && sclk == 0 && swn == 0
    % if zero input -> o-c is zero and o = c -> when the NGS files are
    % created simply write the c value to the o value
    zinp = 1;
    disp '    '
else
    % if not zero input -> compute the simulated o-c
    zinp = 0;
    disp '    '
    disp 'computing the simulated o-c values...'
    disp '    '
    % loop over all scans
    for i = 1:nscan
        % display the progress
        if mod(i,50) == 0
            disp (['    progress: scan ',num2str(i),' of ',num2str(nscan)])
        end
        % number of observations at each station (not all stations participate in all scans)
        l1 = zeros(nant,1);
        l2 = zeros(nant,1);
        l1([scan(i).obs.i1]) = 1;
        l2([scan(i).obs.i2]) = 1;
        % a station can be listed once as 1st and once as 2nd station but
        % only observes once in a scan
        st = st + or(l1,l2);
        % loop over the observations in a scan
        for j = 1:length(scan(i).obs)
            k = k + 1;
            % index of first station
            st1 = scan(i).obs(j).i1;
            % index of second station
            st2 = scan(i).obs(j).i2;
            
            % source structure; same for each simulated day
            soucorr=0;
            if ssou==1
                ind=strmatch(sources.q(scan(i).iso).name,cat_comp.name,'exact');
                sources.q(scan(i).iso).sim_model=[cat_comp.flux(ind),cat_comp.maj(ind),cat_comp.min(ind),cat_comp.angle(ind),cat_comp.dRA(ind),cat_comp.dDec(ind)];
                soucorr=modDelay(sources.q(scan(i).iso).sim_model,scan(i).stat(st1).x,scan(i).stat(st2).x,([8213 8252 8353 8513 8733 8853 8913 8933]+4), 8217, sources.q(scan(i).iso).ra2000, sources.q(scan(i).iso).de2000, scan(i).mjd) * 1E-12;
                sim(1).soucorr(k)=soucorr;
                % the CRS position would be scan(isc).stat(ist).xcrs
            end
            
            % generating the (o-c) delay observables
            sim_st1 = sim(st1);
            sim_st2 = sim(st2);
            st_st1 = st(st1);
            st_st2 = st(st2);
            for iday = 1:sim_idays
                scan(i).obs(j).obs(iday) = sim_st2.swd(st_st2,iday) + sim_st2.clk(st_st2,iday) - sim_st1.swd(st_st1,iday) - sim_st1.clk(st_st1,iday) + wn(k,iday) + soucorr; % [s]
            end
            
            % if a white noise was specified in the turbulence parameter file or in the GUI 
            % => replace the uncertainty obtained from the input data file (e.g. NGS)
            if wnoise ~= 999
                switch(scan(i).obs_type) 
                    case 'q'
                        scan(i).obs(j).sig = wnoise * 1e-12; % [s]
                    case 's'
                        scan(i).obs(j).sig = wnoise_sat * 1e-12; % [s]
                end 
            end
            
        end
    end
end

% save simulated data, covariance matrices, and the new scan structure
% array to ../DATA/LEVEL4
disp '    '
disp 'saving simulated data...'

% if ~exist(dirpt1,'dir')
%     mkdir(dirpt1)
% end
% 
% try
%     if sind < 10
%         save([dirpt1 '/' session '_S00' num2str(sind) '_sim.mat'],'sim');
%     elseif sind < 100
%         save([dirpt1 '/' session '_S0' num2str(sind) '_sim.mat'],'sim');
%     else
%         save([dirpt1 '/' session '_S' num2str(sind) '_sim.mat'],'sim');
%     end
% catch
%     try
%         if sind < 10
%             save([dirpt1 '/' session '_S00' num2str(sind) '_sim.mat'],'sim');
%         elseif sind < 100
%             save([dirpt1 '/' session '_S0' num2str(sind) '_sim.mat'],'sim');
%         else
%             save([dirpt1 '/' session '_S' num2str(sind) '_sim.mat'],'sim');
%         end
%     catch
%         if sind < 10
%             save([dirpt1 '/' session '_S00' num2str(sind) '_sim.mat'],'sim');
%         elseif sind < 100
%             save([dirpt1 '/' session '_S0' num2str(sind) '_sim.mat'],'sim');
%         else
%             save([dirpt1 '/' session '_S' num2str(sind) '_sim.mat'],'sim');
%         end
%     end
% end
% try
%     if sind < 10
%         save([dirpt1 '/' session '_S00' num2str(sind) '_cov.mat'],'cov');
%     elseif sind < 100
%         save([dirpt1 '/' session '_S0' num2str(sind) '_cov.mat'],'cov');
%     else
%         save([dirpt1 '/' session '_S' num2str(sind) '_cov.mat'],'cov');
%     end
% catch
%     try
%         if sind < 10
%             save([dirpt1 '/' session '_S00' num2str(sind) '_cov.mat'],'cov');
%         elseif sind < 100
%             save([dirpt1 '/' session '_S0' num2str(sind) '_cov.mat'],'cov');
%         else
%             save([dirpt1 '/' session '_S' num2str(sind) '_cov.mat'],'cov');
%         end
%     catch
%         if sind < 10
%             save([dirpt1 '/' session '_S00' num2str(sind) '_cov.mat'],'cov');
%         elseif sind < 100
%             save([dirpt1 '/' session '_S0' num2str(sind) '_cov.mat'],'cov');
%         else
%             save([dirpt1 '/' session '_S' num2str(sind) '_cov.mat'],'cov');
%         end
%     end
% end
% save([dirpt1 '/' session '_scan.mat'],'scan');
% try
%     save([dirpt1 '/' session '_mf.mat'],'mf');
% catch
%     try
%         save([dirpt1 '/' session '_mf.mat'],'mf');
%     catch
%         save([dirpt1 '/' session '_mf.mat'],'mf');
%     end
% end
% try
%     save([dirpt1 '/' session '_azel.mat'],'azel');
% catch
%     try
%         save([dirpt1 '/' session '_azel.mat'],'azel');
%     catch
%         save([dirpt1 '/' session '_azel.mat'],'azel');
%     end
% end

% ##### Write output #####

% ### write NGS files ###
if (flag_write_ngs_file == 1)     
    % Write warning to CW, if a NGS file is written with simulated satellite observations 
    if flag_sim_satellite_obs
       fprintf('WARNING: Session %s contains observations to satellites. The NGS format is not suitable for this observation type!\n', session) 
    end
    fprintf('Writing NGS files...\n') 
    scan2ngs(session,antenna,scan,sources.q,scan(1).tim(1),zinp,sind,sim_idays,dirpt0); 
    fprintf('...finished\n') 
end

% ### write VSO files ###
if (flag_write_vso_file)     
    fprintf('Writing VSO files...\n') 
    write_vso_sim(session, antenna, scan, sources, scan(1).tim(1), sim_idays, sind, zinp, dirpt0)
    fprintf('...finished\n') 
end

% ###### just updateing structures if you use them in vie_lsm afterwards ######
% ### calc observation and sigma
% scan structure
for iscan = 1:length(scan)
    for iobs = 1:scan(iscan).nobs        
        % zero input or simulated input?
        if zinp == 1
            % if zero input -> o = c -> o-c = 0
            scan(iscan).obs(iobs).obs = scan(iscan).obs(iobs).com;
        else
            scan(iscan).obs(iobs).obs = scan(iscan).obs(iobs).com + [scan(iscan).obs(iobs).obs]';
        end
        % the sigma of the simulated delay observable is set to
        % the value of the simulated thermal noise, ionospheric formal
        % error has to be taken into account.
        scan(iscan).obs(iobs).sig = sqrt(sqrt((1.0d9 * scan(iscan).obs(iobs).sig)^2 - ...
            scan(iscan).obs(iobs).sgdion^2)^2+scan(iscan).obs(iobs).sgdion^2)*1e-9;
    end
end

% update antenna structure
if length(antenna(1).session) ~=9 
    if strcmpi(antenna(1).session(11),'V')
        fname = antenna(1).session;
    else
        fname = antenna(1).session(1:9);
    end
else
    fname = antenna(1).session;
end

if isempty(dirpt0)
    sdir = ['../SIM/', num2str(parameter.year)];
else
    sdir = ['../SIM/', dirpt0, '/',num2str(parameter.year)];
end


antenna(1).ngsfile = [];
antenna(1).session = [];
antenna(1).ngsfile = sdir;

for isim = 1:sim_idays
    antenna(1).session{isim} = sprintf('%s_S%03d',fname,sind+isim-1);
end
for iant = 1:length(antenna)
    antenna(iant).gpt3e = 1;
end

% update parameter structure
% parameter.session_name = [];
% jetfilnam = parameter.vie_init.jetfilnam;
% parameter.vie_init.jetfilnam=[];
% jetfilnamuv = parameter.vie_init.jetfilnamuv;
% parameter.vie_init.jetfilnamuv=[];
% jetfilnamjb = parameter.vie_init.jetfilnamjb;
% parameter.vie_init.jetfilnamjb = [];
% 
% parameter.filepath = ['../DATA/' sdir(4:end)];
% 
% if sim_idays == 1
%     parameter.session_name = sprintf('%s_S%03d',fname,sind+isim-1);
%     parameter.vie_init.jetfilnam = sprintf('%s_S%03d.JET',jetfilnam(1:end-4-5),sind+isim-1);
%     parameter.vie_init.jetfilnamuv = sprintf('%s_S%03d.JETUV',jetfilnamuv(1:end-6-5),sind+isim-1);
%     parameter.vie_init.jetfilnamjb = sprintf('%s_S%03d.JETJB',jetfilnamjb(1:end-6-5),sind+isim-1);
% else
%     for isim = 1:sim_idays
%         parameter.session_name{isim} = sprintf('%s_S%03d',fname,sind+isim-1);
%         parameter.vie_init.jetfilnam{isim} = sprintf('%s_S%03d.JET',jetfilnam(1:end-4-5),sind+isim-1);
%         parameter.vie_init.jetfilnamuv{isim} = sprintf('%s_S%03d.JETUV',jetfilnamuv(1:end-6-5),sind+isim-1);
%         parameter.vie_init.jetfilnamjb{isim} = sprintf('%s_S%03d.JETJB',jetfilnamjb(1:end-6-5),sind+isim-1);
%     end
% end
% ### write output in LEVEL1 folder ###
% allScans = scan;
% allSessionName = parameter.session_name;
% parameter.session_name = [];
% allInitJetfilnam = parameter.vie_init.jetfilnam;
% parameter.vie_init.jetfilnam=[];
% allInitJetfilnamuv = parameter.vie_init.jetfilnamuv;
% parameter.vie_init.jetfilnamuv=[];
% allInitJetfilnamjb = parameter.vie_init.jetfilnamjb;
% parameter.vie_init.jetfilnamjb = [];
% for isim = 1:sim_idays
%     if sim_idays == 1
%         parameter.session_name = allSessionName;
%         parameter.vie_init.jetfilnam   = allInitJetfilnam;
%         parameter.vie_init.jetfilnamuv = allInitJetfilnamuv;
%         parameter.vie_init.jetfilnamjb = allInitJetfilnamjb;
%         antenna(1).session = allSessionName;
%     else
%         parameter.session_name = allSessionName{isim};
%         parameter.vie_init.jetfilnam   = allInitJetfilnam{isim};
%         parameter.vie_init.jetfilnamuv = allInitJetfilnamuv{isim};
%         parameter.vie_init.jetfilnamjb = allInitJetfilnamjb{isim};
%         antenna(1).session = allSessionName{isim};
%     end
%     
%     for iscan = 1:length(scan)
%         for iobs = 1:length(scan(iscan).obs)
%             scan(iscan).obs(iobs).obs = allScans(iscan).obs(iobs).obs(isim);
%         end
%     end
%     save(sprintf('%s/%s_parameter.mat',dirpt,parameter.session_name),'parameter');
%     save(sprintf('%s/%s_antenna.mat',dirpt,parameter.session_name),'antenna');
%     save(sprintf('%s/%s_scan.mat',dirpt,parameter.session_name),'scan');
%     save(sprintf('%s/%s_sources.mat',dirpt,parameter.session_name),'sources');
% end
% 
% parameter.session_name = allSessionName;
% parameter.vie_init.jetfilnam = allInitJetfilnam;
% parameter.vie_init.jetfilnamuv = allInitJetfilnamuv;
% parameter.vie_init.jetfilnamjb = allInitJetfilnamjb;
% scan = allScans;
% antenna(1).session = allSessionName;
% 
% disp '    '
disp 'vie_sim successfully finished!'