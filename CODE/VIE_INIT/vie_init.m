% ************************************************************************
%   Description:
%   The first part of VieVS. Reads the data from the NGS files and some
%   additional inforamtion.
%
%   Input:	
%      Uses the parameter file stored in currentf_pa.mat. The string
%      ngsfile, containing the name of NGS-file to be processed must exist
%      in the workspace.
%   Input arguments:
%    - obs_file_name            - directory and name of the input file (e.g. "2005/05APR04XA_N004")
%    - parameter                - VieVS parameter structure (contains GUI parameters)
%    - out_vie_init_subdir      - sub-directory for VIE_INIT (LEVEL0)
%    - obs_file_dir             - Input file directiory (if it is empty the default path is used for different input file formats: "obs_fle_dir_tmp")
%    - trf                      - TRF data structure (optional)
%    - crf                      - CRF data structure (optional)
%
%   Output:
%      The antenna, scan, and sources structure arrays, saved in the
%      in the LEVEL0 directory.
% 
%   External calls: 	
%   	read_ngs.m	  
%       read_vso.m
%       get_trf_and_crf.m
%       constants.m
%       read_vgosdb_input_settings.m
%       
%   Coded for VieVS: 
%   May 2009 by Tobias Nilsson
%
%   Revision: 
%   18 Nov 2009 by Tobias Nilsson: the program now reads the outlier-file
%     for the session if it exist and specified in the option. The
%     information from the outlier-file is used in read_ngs.
%   28 Jan 2010 by Tobias Nilsson: Outliers now read from yearly
%     directories (../DATA/OUTLIERS/YYYY/).
%   24 Feb 2010 by Tobias Nilsson: file paths for input and output are modified 
%     user defined sub-directory for input and output is available
%   26 Feb 2010 by Tobias Nilsson: Possible to use different subdirectories
%     for outliers.
%   24 Jan 2011 by Hana Spicakova: choice between GPT function
%      (Global Pressure and Temperature) for all observations and measured 
%      meteorological data from NGS file with GPT only as backup
%   18 Apr 2011 by TObias Nilsson: Changed input/output variables
%   04 May 2011 by Matthias Madzak: New parameter variables
%       (parameter.vie_init.iono, _.ionoFolder) to ini_opt (they are needed
%       in read_ngs).
%   09 May 2011 by Tobias Nilsson: Now also the baselines excluded
%               are printed on the screen.
%   13 Nov 2012 by Hana Kr???sn???: antenna offset and mounting type taken from
%               superstatin file (antenna_info). No more from NGS header.
%   14 Mar 2013 by Matthias Madzak: implementation of supersource file.
%   05 Feb 2014 by Lucia Plank: read JET files
%   16 May 2014 by Monika Tercjak: reading and saving names of stations
%               which should be down-weighted
%   30 May 2014 by David Mayer: a message when no OPT file was found 
%   26 Sep 2014 by Hana Krasna: bug fixed, by simulations the Outlier
%               directory is now recognised
%   11 Nov 2014 by A. Hellerschmied: added option to decide wether to use
%               OPT files or not.
%   12 Jan 2015 by Caroline/Matthias: NO CABLE CAL Info print at command 
%               window
%   14 Jan 2015 by Daniel Landskron: output in Command Window slightly
%               modified
%   20 Jul 2015 by A. Girdiuk: output start and end time moments for
%       excluded station
%   21 Jul 2015 by A.Girdiuk: output for outliers list
%   04 Dec 2015 by A. Hellerschmied: OPT file is now loaded in VIE_INIT only!
%   10 Dec 2015 by A. Hellerschmied: Bug fix: Init. of clock break and ref. 
%       station data added in case no OPT file is available.
%   08 Jan 2016 by M. Madzak: Added read possibility for reading vgos-db
%       (netCDF) files
%   28 Jan 2016 by C. Sch???nberger: sources can be excluded for a certain time span.
%   28 Jun 2016 by A. Hellerschmied: Changes in call of read_ngs.m
%   12 Jul 2016 by A. Hellerschmied: Revised vie_init.m for VieVS 3.0:
%                   - enhanced support of vgosDB files as data input
%                   - support of VSO files as data input
%                   - Many small changes and fixes
%   17 Jul 2016 by A. Hellerschmied: - Content of "sources" structure (natural sources, quasars) is now stored in the sub-structure "sources.q"
%                                    - The sub-structure "sources.sc" can be used to define space-crafts
%   27 Jul 2016 by D. Mayer: Bug fix: Sessions with 2 underscores are now read correctly 
%   03 Aug 2016 by A. Hellerschmied: Enhanced support of .vso files as input data.
%   08 Aug 2016 by A. Girdiuk: bug-fix: parameter initialization is added for stations to be down-weighted
%   09 Aug 2016 by A. Hellerschmied: "sources.sc" changed to "sources.s"
%   21 Sep 2016 by A. Hellerschmied: Possibility added to read space craft ephemerids (TRF positions and velocities) from an external file
%   26 Oct 2016 by A. Hellerschmied: When loading vso files the fileapth is now taken from "parameter.filepath" (defined in vie_batchX_X.m)
%   08 Nov 2016 by D. Mayer: Sources can be deleted from .txt file
%   10 Nov 2016 by H. Krasna: exclude sources with <3 observations from the
%               NNR constraints, or fixed them if the sources are estimated as 
%               pwl offsets to avoid singularity
%   07 Feb 2017 by A. Hellerschmied: "parameter" structure added as in/output arguments for read_vso.m
%   14.Feb 2017 by M. Schartner: new optional parameters: trf and crf.
%   22 Feb 2017 by A. Hellerschmied: Manual TRF file support established
%   23 Feb 2017 by A. Hellerschmied: get_trf_and_crf.m used for loading CRF and TRF data
%   14 Mar 2017 by A. Hellerschmied: call of function "constants.m"
%   31 Mar 2017 by D. Mayer: added the possibility to remove list of
%   			station from every session in the code
%   13 Nov 2017 by J. Gruber: new functions to specify insitute provider,
%   			version name and frequency band (read_vgosdb_input_settings.m)
%   06 Jan 2018 by J. Gruber: Changed call of cleanScan.m
%   28 Aug 2018 by D. Landskron: shape of output slightly changed
%   28 Nov 2018 by D. Landskron: workaround concerning OPT files changed: now NO observations are excluded, because everything will be done in vie_lsm later
%   05 Dec 2018 by D. Landskron: command window output slightly extended

% ************************************************************************
%
function [antenna,sources,scan,parameter]=vie_init(obs_file_name, parameter, out_vie_init_subdir, obs_file_dir, varargin)
disp('---------------------------------------------------------------')
disp('|                     Welcome to VIE_INIT                     |')
disp('---------------------------------------------------------------')
disp(' ')

% Call the constants function to get the global variables for constants:
constants;
% Get the session name:
switch(parameter.data_type)
    case 'ngs'
        obs_fle_dir_tmp              = 'NGS';
    case 'vso'
        obs_fle_dir_tmp              = 'VSO';
    case 'vgosdb'
        obs_fle_dir_tmp              = 'vgosDB';
end % switch(parameter.data_type)

session = parameter.session_name;

% check if ngsdir exists:
if exist('ngsdir', 'var')
    if isempty(obs_file_dir)
        obs_file_dir = obs_fle_dir_tmp;
    end
else
    obs_file_dir = obs_fle_dir_tmp;
end


%% ##### CRF + TRF #####
trffile = parameter.vie_init.trf;
crffile = parameter.vie_init.crf;

% Load TRF and CRF data, if the data are not provided as input arguments
if isempty(varargin)
    [trf, crf, parameter] = get_trf_and_crf(parameter);
else
    trf = varargin{1};
    crf = varargin{2};
end



%% ##### Load observation data #####

% Distinguish between different input data types:
switch(parameter.data_type)

    % #############################
    % #####     vgosDB        #####
    % #############################
    case 'vgosdb'
    
        % folder of e.g. Head.nc file
        curNcFolder = ['../DATA/vgosDB/',parameter.year, '/',parameter.session_name,'/'];

        % uncompress vgosDB *.tar.gz or *.tgz file
        wasCompressed = false;
        if ~isfolder(curNcFolder)
            wasCompressed = true;
            
            % uncompress vgosDB *.tar.gz or *.tgz file
            vgosTargz = [curNcFolder(1:end-1),'.tar.gz'];
            vgosTgz = [curNcFolder(1:end-1),'.tgz'];
            curSlash = sort([strfind(curNcFolder,'/'), strfind(curNcFolder,'\')]);
            vgosTgzFolder = curNcFolder(1:curSlash(end-1));

            if exist(vgosTgz,'file')
                untar(vgosTgz,vgosTgzFolder);
            elseif exist(vgosTargz,'file')
                untar(vgosTargz,vgosTgzFolder);
            else
                fprintf('ERROR: %s does not exist!\n',vgosTgz);
            end       
        end
        
        % The folder within the .tgz file might have an name extension e.g.:
        % 19JUN11VG.tgz contains 19JUN11VG_Wf_noOw
        % check if the folder with the expected name exists and if not try
        % to update the curNcFolder var.
        if ~isfolder(curNcFolder)
            fprintf('No folder with the name %s was found after uncompressing the .tgz file.\nChecking for name extensions.\n', parameter.session_name)
            
            potentialFolders = dir(['../DATA/vgosDB/',parameter.year, '/',parameter.session_name,'*']);
            dirFlags = [potentialFolders.isdir];
            potentialFolders = potentialFolders(dirFlags);

            if numel(potentialFolders) == 0
                fprintf('No folders found.\n')
            elseif numel(potentialFolders) > 1
                fprintf('Multiple folders with the name %s* found:\n', parameter.session_name)
                fprintf('%s\n', potentialFolders.name)
                fprintf('Aborting.\n')
            else
                fprintf('Found one matching folder: %s. Using this folder.\n', potentialFolders.name)
                curNcFolder = ['../DATA/vgosDB/',parameter.year, '/', potentialFolders.name,'/'];
            end
        end
        
        try
        % read netCDF data
        [out_struct, nc_info]=read_nc(curNcFolder);
                  
        % convert gui values to vie_init values
        [fb,in,ambcorr,ioncorr,wrapper_v,wrapper_k] = vgosdbinGUI2vieinit(parameter.vie_init.vgosDb_observation_parameter,parameter.vie_init.vgosDb_institute,parameter.vie_init.ambiguity_correction,parameter.vie_init.iono_correction,parameter.vie_init.vgosDb_wrapper_version);
                        
        % Standard settings, which are used if not defined differently in settings file:
        if isempty(in{1}) % institute
            in = {'IVS'};
            fprintf('Set institute to default: %s\n', in{1})
        end
        if isempty(fb) % frequency band
            fb = 'GroupDelayFull_bX';
            fprintf('Set frequency band to default: %s\n', fb)
        end
        if isempty(wrapper_k) % wrapper tag
            wrapper_k = 'all';
            fprintf('Set wrapper tag to default: %s\n', wrapper_k)
        end
        if isempty(wrapper_v) % wrapper version
            wrapper_v = 'highest_version';
            fprintf('Set wrapper version to default: %s\n', wrapper_v)
        end
        if isempty(ioncorr) % ionosphere correction
            ioncorr = 'on';
            fprintf('Set ionosphere correction: %s (default)\n', ioncorr)
        end
        if isempty(ambcorr) % ambiguity correction
            ambcorr = 'on';
            fprintf('Set ambiguity correction: %s (default)\n', ambcorr)
        end
        
        % Read wrapper
        wrapper_data = read_vgosdb_wrapper(curNcFolder, parameter.session_name, in, wrapper_k, wrapper_v);
        
        % get scan, antenna, source struct from netCDF files
        scan        = nc2scan(out_struct, nc_info, fb, ioncorr, ambcorr, wrapper_data, parameter);
        antenna     = nc2antenna(out_struct, trf, trffile{2}, wrapper_data);
        sources     = nc2sources(out_struct, crf, crffile{2}, wrapper_data);
        
        % "clean" scan struct (because of exclusions)
        % [scan, sources, antenna] = cleanScan(scan, sources, antenna, out_struct.head.StationList.val', out_struct.head.SourceList.val', ini_opt, bas_excl, parameter.obs_restrictions.Qlim, parameter.obs_restrictions.cut_off_elev);

    
        % Create a sub-structure in "sources" for quasars sources:
        % => If observations of near field sources will be supported by vgosDB in future, the nc2sources.m script has to be adopted!
        q = sources;
        clear sources,
        sources.q   = q;
        sources.s 	= [];
        fprintf('...reading the vgosDB file finished!\n');
        
        catch ME
           rmdir(curNcFolder, 's');
           rethrow(ME)
        end
        
        % remove the unpacked vgosDB folder
        if wasCompressed
            rmdir(curNcFolder, 's');
        end
    

    % #############################
    % #####        NGS        #####
    % #############################
    case 'ngs'
        
         
        % +++++++++++++++++++ read jet ang file ++++++++++++++++++++++++++++++++++
        %  => Only works for NGS files!
        % parameter.vie_init.jetfilnam=['../DATA/JETANG/',session(1:min([length(session),14])),'.JETUV'];
        parameter.vie_init.jetfilnam    = ['../DATA/JETANG/',session(1:min([length(session),14])),'.JET'];
        parameter.vie_init.jetfilnamuv  = ['../DATA/JETANG/',session(1:min([length(session),14])),'.JETUV'];
        parameter.vie_init.jetfilnamjb  = ['../DATA/JETANG/',session(1:min([length(session),14])),'.JETJB'];
        if parameter.vie_init.ex_jet == 1
            jamax = parameter.vie_init.jetang;
            [ini_opt.scan_jet] = readJET(parameter.vie_init.jetfilnam,jamax);
            fprintf('read JET file: %s\n',parameter.vie_init.jetfilnam)
            fprintf('Number of obs to be excluded due to jet angle: %2.1f %3.0f\n', jamax,length(ini_opt.scan_jet))
        else
            ini_opt.scan_jet = [];  
        end
        % ------------------- read jet ang file -------------------------------------
        
        
        % this is necessary in read_NGS, for whatever reason
        ini_opt.iono = parameter.vie_init.iono; % Use Iono. corrections from NGS file or from .ion files?
        
        
        fprintf(' => Start reading %s\n',obs_file_name);
        if isnan(str2double(obs_file_name(1))) % if first element in ngsfile is character - absolute path of NGS file is given
            [antenna,sources,scan] = read_ngs(obs_file_name, trffile, crffile, ini_opt, trf, crf);
        else
            [antenna,sources,scan] = read_ngs(['../DATA/' obs_file_dir '/' obs_file_name], trffile, crffile, ini_opt, trf, crf);
        end
        fprintf('...reading the NGS file finished!\n');
        
        % Create a sub-structure in "sources" for quasars sources:
        q = sources;
        clear sources
        sources.q   = q;
        sources.s 	= [];
        


    % #############################
    % #####     VSO           #####
    % #############################
    case 'vso'
        vso_file_path = parameter.filepath; % ['../DATA/VSO/',parameter.year, '/'];
        vso_file_name = parameter.session_name;
        fprintf(' => Start reading %s\n', [vso_file_path, vso_file_name]);
        
        % Get satellite ephem. file path and name
        str_ind = max([strfind(parameter.vie_init.sc_orbit_file_path_name, '\'), strfind(parameter.vie_init.sc_orbit_file_path_name, '/')]);
        if ~isempty(str_ind)
            sat_orbit_file_path = parameter.vie_init.sc_orbit_file_path_name(1 : str_ind); 
            sat_orbit_file_name = parameter.vie_init.sc_orbit_file_path_name(str_ind+1 : end);
% +++ Workaround solution! Has to be implemented properly in the GUI!
%   => Currently the file type is derived from the file extension 
            % Get file type:
            if strcmp(parameter.vie_init.sc_orbit_file_path_name(end-3 : end), '.sp3')
                parameter.vie_init.sc_orbit_file_type = 'sp3';
            else
                parameter.vie_init.sc_orbit_file_type = 'sat_ephem_trf';
            end
% --- Workaround solution!
            sat_orbit_file_type = parameter.vie_init.sc_orbit_file_type;
        else
            sat_orbit_file_name = '';
            sat_orbit_file_path = '';
            sat_orbit_file_type = '';
        end
        
        [antenna, sources, scan, parameter] = read_vso(vso_file_path, vso_file_name, trf, trffile{2}, crf, crffile{2}, sat_orbit_file_path, sat_orbit_file_name, sat_orbit_file_type, parameter);
        fprintf('...reading the VSO file finished!\n');
        
end % switch(parameter.data_type)


% ##### Write info to CW #####
fprintf('\n');
fprintf('A total of %d stations, %d sources (quasars), %d scans and %d observations were found.\n', length(antenna), length(sources.q), length(scan), length([scan.obs]));
disp('The following stations were found:')
for i_ant = 1:length(antenna)
  fprintf('%2.0f%s%s\n',i_ant,'. ',antenna(i_ant).name)
end

% ????????? Needed: % Yes, currently (2018-10-04) it is still needed for VIE_SIM!!! This is stupid => Theses parameters are stored in the "parameters" struct!
% => Problem: If antenna(1) is removed via "cleanScan.m" in VIE_LSM, these parameters are not accessible any more!
antenna(1).ngsfile=obs_file_name;
antenna(1).session=session;
% ???????????????????????????????????????

fprintf('\nvie_init successfully finished!\n');

