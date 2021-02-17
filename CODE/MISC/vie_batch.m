% #########################################################################
% #     vie_batch3_0
% #########################################################################
%
% DESCRIPTION
%   This batch file is used to run the different VieVS modules of version 3.0.
%   It also manages parallel processing features and the according error handling.
%
% CREATED
%   ????
%
% COUPLING
% - vie_sched
% - vie_lsm
% - vie_mod
% - vie_sim
% - vie_lsm_scanwise
% - vie_glob
% - get_trf_and_crf
% - writeLEVEL3toXLSX
%
%
% CHANGES
% - May 2013 by Lucia & Matthias, parameter file
% - Jan 2014 by Matthias & Hana, bug corrected for LEVEL2 subdirectory
% - 2015-03-26: Andreas Hellerschmied: Parallel computing compatibility issue: use of "parpool" instead of "matlabpool" for MATLAB version R2013b and later.
% - Apr 2015 by David, added exception handling
% - May 11 2015 by David, added progress output in percentage for parfor loop
% - 2015-03-26: Andreas Hellerschmied: added filepath to ../COMPILE/VIE_SCHED_V23/SAT_SCHED
% - 2015-07-20: Andreas Hellerschmied: added filepath to ../COMPILE/VIE_SCHED_V23/SAT_SCHED/SGP4
% - 2015-08-20: D. Mayer: Bug fix
% - 2015-12-14: S. Boehm: Fixed the unnecessary going through all sessions of the process_list, if only vie_glob is run.
% - 2015-12-18, A. Hellerschmied: - Adapted for VieVS 3.0
%                                 - Renamed: vie_batch2_3.m => vie_batch3_0.m
% - 2016-01-05, A. Hellerschmied: Minor changes (code for pre 2.0 VieVS versions removed) 
% - 2016-07-12, A. Hellerschmied: Revised this function
%                   - Support of different types of input data (vgosDB, NGS, VSO)
%                   - New field in the parameter strucutre: parameter.year, parameter.data_type, parameter.session_name
% - 2016-08-30, A. Hellerschmied: filepath of input file is now saved in parameter.filepath
% - 2016-09-29, H. Krasna: save('PROCESSLIST/failed_sessions.mat','process_list'); \ change to / to be compatible with Linux
% - 2016-10-13, A. Hellerschmied: Field "parameter.fielpath" was not written for vgosDB files
% - 2016-10-19, A. Hellerschmied: Warning messages corrected
% - 2016-10-26, A. Hellerschmied: Correct "parameter.filepath" and "parameter.year" for VSO files
% - 2016-11-28, M. Schartner: rearranging modules
%                             - deeper integration of vie_sched an vie_sim
%                             - parallel pool is created at the beginning of vie_batch3_0
% - 2016-11-30, A. Hellerschmied: Added exception for empty process lists (as may returned from vie_sched)
% - 2017-01-09, A. Hellerschmied: Added dir. to search path: ../COMPILE/VIE_SCHED_V30/SAT_SCHED/VALLADO/
% - 2017-02-14, M. Schartner: Now loads CRF and TRF here, if you run vie_init
% - 2017-02-23, A. Hellerschmied:   - get_trf_and_crf.m used for loading CRF and TRF data
%                                   - parfor_progress removed 
% - 2017-02-24, A. Hellerschmied: Fixed some bugs and removed unnecessary code parts
% - 2017-02-27, M. Schartner: Now writes LEVEL3 statistics for multisched sessions to xlsx file
% - 2017-03-10, M. Schartner: small bugfix
% - 2017-03-21, A. Hellerschmied: Minor bug-fix 
% - 2017-05-31, A. Hellerschmied: - Adapted for VieVS 3.1
%                                 - Renamed: vie_batch3_0.m => vie_batch3_1.m
% - 2017-12-14, A. Hellerschmied: Path to /COMPILE/MISC/ added
% - 2018-01-18, A. Hellerschmied    - Revised for handling of VieVS with GIT

%*************************************************************************

function vie_batch(varargin)

%% load('currentf_pa.mat','parameter');
load('process_list','process_list');
load('runp','runp')
load('guiparameter.mat');
cleanVgosDB();

pthDALE = '../'; % path to DAta LEvel (only LEVEL0 and LEVEL1 at the moment)
% pthDALE = '/data/VIEVS/' 

% parameter.pthDALE = pthDALE; % can be put in GUI


startTime = tic;
%% if you passed a argument interpret it as session name
if(nargin > 0)
    
    % get session
    session = varargin{1};
    found = false;
    for i=1:size(process_list,1)
        if(contains(process_list(i,:),session))
            found = true;
            break
        end
    end
    if(~found)
        error('Session not found in process_list');
    end
    process_list = process_list(i,:);
    fprintf('Only running VIE_LSM for session: %s',process_list)
    
    % outlier test
    outlierTest = varargin{2};
    if strcmp(outlierTest,'simple')
        parameter.lsmopt.simple_outlier = 1;
        parameter.lsmopt.basic_outlier = 0;
    elseif strcmp(outlierTest,'normal')
        parameter.lsmopt.simple_outlier = 0;
        parameter.lsmopt.basic_outlier = 1;
    else
        parameter.lsmopt.simple_outlier = 0;
        parameter.lsmopt.basic_outlier = 0;
    end
        
    % eliminate outlier
    outlierDirectory = varargin{3};
    if strcmp(outlierDirectory,'none')
        parameter.outlier.flag_remove_outlier = 0;
    else
        parameter.outlier.flag_remove_outlier = 1;
        parameter.outlier.out_file_dir = outlierDirectory;
    end
    
    runp.init = 0;
    runp.mod = 0;
    runp.lsm = 1;
    runp.parallel = 0;
end

fid = fopen('failed_sess_temp.txt', 'w'); %delete content of file - this file is only temporary. It's purpose is to save the failed sessions even when Matlab is aborted by hand
fclose(fid);
% guiparameter=parameter;

if ~isfield(runp,'sim')
    runp.sim=0;
end

%% get parallel
if isfield(runp, 'parallel')
    VieVS_parallel=runp.parallel;
else
    VieVS_parallel=0;
end

if runp.parallel
    % Get MATLAB version (release name)
    % Use "parpool" instead of "matlabpool" for version "2013b" and following releases
    matlab_version = version('-release');
    release_year = str2double(matlab_version(1:4));
    if (release_year > 2013) || ((release_year == 2013) && (strncmp(matlab_version(5), 'b', 1)))
        flag_release_r2013b_or_later = 1;
    else
        flag_release_r2013b_or_later = 0;
    end
    % get number of cores
    nCores=runp.nCores{1};
    % write user info
    fprintf('Starting parallel VieVS...\n');
    
    % ##### Parallel computing before release R2013b #####
    if flag_release_r2013b_or_later == 0
        % close matlabpool if already open
        matlabpool size;
        if ans~=0
            matlabpool close;
        end
        % open matlabpool
        if strcmp(nCores, 'auto')
            matlabpool open
        else
            eval(['matlabpool open ', nCores]);
        end
        
    % ##### Parallel computing for release R2013b and later #####
    elseif flag_release_r2013b_or_later == 1
        % close matlabpool if it is already open
        poolobj = gcp('nocreate'); % Get current pool (object)
        if isempty(poolobj)
%             delete(poolobj)
            % Create a parallel pool
            % ...with the default profile:
            if strcmp(nCores, 'auto')
                poolobj = parpool;
                %...with the specified number of workers:
            else
                nCores = str2double(nCores);
                poolobj = parpool(nCores);
            end
        end
    end
end



%% VIE_SCHED

if isfield(runp, 'sched')
    if runp.sched
        process_list = vie_sched();
    end
end

%% if you run vie_init get TRF and CRF
    [trf, crf, parameter] = get_trf_and_crf(parameter);


% setup for parallel counter in case MATLAB version > 9.1 (2017a). 
number_of_sessions = size(process_list, 1);
if ~verLessThan('Matlab','9.2') && runp.parallel
    D = parallel.pool.DataQueue;
    h = waitbar(0, 'Please wait ...');
    afterEach(D, @nUpdateWaitbar);
    counter = 0;
else
    h = [];
    D = [];
end

% define function for increment the parallel counter. Only used if
% Matlab version is > 9.1 (2017a)
function nUpdateWaitbar(~)
    counter = counter + 1;
    time = toc(startTime);
    frac = counter/number_of_sessions;
    dtime = time/frac*(1-frac)/60;

    if isvalid(h)
        waitbar(counter/number_of_sessions, h, sprintf('finished: %d of %d (%.0f min left)',counter,number_of_sessions, dtime));
    end
    dispSessionNumber(counter, number_of_sessions)
end

%% number of sessions to process:
% Only start VIE_INT, VIE_MOD, VIE_LSM, VIE_SIM and VIE_GLOB, if the process list is not empty!
if ~isempty(process_list)
    
    % Get the number of sessions in the current process list:
    % - Remove invalid entries (only one blank)
    while isempty(find(process_list(number_of_sessions,:)~=' ',1))
        number_of_sessions  = number_of_sessions - 1;
    end

    sess_err = cell(number_of_sessions,1); % initialize cell array for failed sessions
    
    % Write parameters to a new variable for par. proc.
    if runp.parallel 
        parameter_tmp = parameter;
    end


    %% main 
    if runp.init || runp.mod || runp.lsm || runp.lsm_scanwise || runp.sim
        %% default - no parallel:
        if VieVS_parallel==0
            delete(h);
            for isess=1:number_of_sessions
                %% parameters for each session 
                fprintf('%s %4.0f %s %4.0f\n','session ',isess,' of ',number_of_sessions)
                session_name = process_list(isess,:);

                while session_name(length(session_name))==' '
                    session_name=session_name(1:length(session_name)-1);
                end

                % Get file format of input data file:
                parameter.data_type = 'ngs';
                if strfind(session_name, ' [vgosDB]')
                    parameter.data_type = 'vgosdb'; 
                elseif strfind(session_name, ' [VSO]')
                    parameter.data_type = 'vso'; 
                end

                %% Get the session name, the fielpath and the year (string)
                switch(parameter.data_type)
                    case 'ngs'
                        flag_absolut_path = false;
                        ind_tmp = [strfind(session_name, '\'), strfind(session_name, '/')];
                        if isempty(ind_tmp)
                            % Format: <session_name>
                            parameter.session_name  = session_name;
                            parameter.filepath      =  '../DATA/NGS/';
                        else
                            parameter.session_name  = session_name(max(ind_tmp) + 1 : end);
                            if (length(ind_tmp) == 1) && (length(session_name(1:min(ind_tmp-1))) == 4)
                                % Format: yyy/<session_name>
                                if ~isnan(str2double(session_name(1:4)))
                                    if (str2double(session_name(1:4)) < 2025) && (str2double(session_name(1:4)) > 1979)
                                        parameter.filepath =  ['../DATA/NGS/', session_name(1:4), '/'];
                                    else
                                        flag_absolut_path = true;
                                    end
                                else
                                    flag_absolut_path = true;
                                end
                            else
                                flag_absolut_path = true;
                            end
                        end

                        % Format: <absolut path>/<session_name>
                        if flag_absolut_path
                            parameter.filepath =  session_name(1 : max(ind_tmp));
                        end

                        year_tmp = str2double(parameter.session_name(1:2));
                        % Check if converversion was sucessfull:
                        if isnan(year_tmp)
                            error('Invalid NGS file name: The first two letters have to represent the year of the session! Please keept the standard naming convention, e.g. "16AUG26XU_N004"!');
                        else
                            % Convert two digit year to four digits!
                            if year_tmp < 79
                                year_tmp = year_tmp + 2000;
                            elseif year_tmp >= 79
                                year_tmp = year_tmp + 1900;
                            end
                        end
                        parameter.year = num2str(year_tmp); % year has to be saved as string!
                        fprintf(' Input file format: NGS\n');

                    case 'vso'
                        flag_absolut_path = false;
                        ind_tmp = [strfind(session_name, '\'), strfind(session_name, '/')];
                        if isempty(ind_tmp)
                            % Format: <session_name>
                            parameter.session_name  = session_name(1 : (strfind(session_name, ' [VSO]')-1));
                            parameter.filepath      =  '../DATA/VSO/';
                            parameter.year          = [];       % year is not included in the session name
                            % NOTE: The name of the folder containing the vso file has to be the session year!
                            % => parameter.year is required to look get the paths to OPT and OUITLIER files in vie_init!
                            error('Invalid filepath of input vso file(%s)! NOTE: The name of the folder containing the vso file has to be the session year!', parameter.filepath)
                        else
                            parameter.session_name  = session_name(max(ind_tmp) + 1 : (strfind(session_name, ' [VSO]')-1));
                            if (length(ind_tmp) == 1) && (length(session_name(1:min(ind_tmp-1))) == 4)
                                % Format: yyy/<session_name>
                                if ~isnan(str2double(session_name(1:4)))
                                    if (str2double(session_name(1:4)) < 2025) && (str2double(session_name(1:4)) > 1979)
                                        parameter.filepath  =  ['../DATA/VSO/', session_name(1:4), '/'];
                                        parameter.year      = session_name(1:4);
                                    else
                                        flag_absolut_path = true;
                                    end
                                else
                                    flag_absolut_path = true;
                                end
                            else
                                flag_absolut_path = true;
                            end
                        end

                        % Format: <absolut path>/<year>/<session_name>
                        % NOTE: The name of the folder containing the vso file has to be the session year!
                        % => parameter.year is required to look get the paths to OPT and OUITLIER files in vie_init!
                        if flag_absolut_path
                            parameter.filepath  = session_name(1 : max(ind_tmp));
                            % Get year:
                            if (length(ind_tmp) >= 2) && (length(session_name(ind_tmp(end-1)+1 : ind_tmp(end)-1)) == 4)
                                if ~isnan(str2double(session_name(ind_tmp(end-1)+1 : ind_tmp(end)-1)))
                                    if (str2double(session_name(ind_tmp(end-1)+1 : ind_tmp(end)-1)) < 2025) && (str2double(session_name(ind_tmp(end-1)+1 : ind_tmp(end)-1)) > 1979)
                                        parameter.year = session_name(ind_tmp(end-1)+1 : ind_tmp(end)-1);
                                    else
                                        error('Invalid filepath of input vso file(%s)! NOTE: The name of the folder containing the vso file has to be the session year!', parameter.filepath)
                                    end
                                else
                                    error('Invalid filepath of input vso file(%s)! NOTE: The name of the folder containing the vso file has to be the session year!', parameter.filepath)
                                end
                            else
                                error('Invalid filepath of input vso file(%s)! NOTE: The name of the folder containing the vso file has to be the session year!', parameter.filepath)
                            end
                        end
                        fprintf(' Input file format: VSO\n');

                    case 'vgosdb'
                        parameter.session_name  = session_name(6 : (strfind(session_name, ' [vgosDB]')-1));
                        parameter.year          = session_name(1:4);
                        parameter.filepath      = ['../DATA/vgosDB/', parameter.year, '/'];
                        fprintf(' Input file format: vgosDB\n');
                end % switch(parameter.data_type)
                session = parameter.session_name;

                try
                    %% VIE_INIT
                    if runp.init
                        fil=[pthDALE 'DATA/LEVEL0/' runp.init_path '/' session];
                        fprintf('Current file: %s\n', fil);
                        if isfield(runp,'ngsdir')
                            [antenna,sources,scan,parameter]=vie_init(session_name,parameter,runp.init_path,runp.ngsdir,trf,crf);
                        else
                            [antenna,sources,scan,parameter]=vie_init(session_name,parameter,runp.init_path,[],trf,crf);                            
                        end
                        if ~exist([pthDALE 'DATA/LEVEL0/' runp.init_path], 'dir')
                            mkdir([pthDALE 'DATA/LEVEL0/' runp.init_path])
                        end
                        savestruct(fil,parameter,antenna,scan,sources);
                    end

                    %% VIE_MOD

                    if runp.mod
                        fil=[pthDALE 'DATA/LEVEL0/' runp.init_path '/' session];
                        if (exist([fil '_parameter.mat'], 'file')&&...
                                exist([fil '_scan.mat'], 'file')&&...
                                exist([fil '_antenna.mat'], 'file')&&...
                                exist([fil '_sources.mat'], 'file'))

                            if ~runp.init
                                tmp=load([fil '_parameter.mat']);
                                parameter.vie_init=tmp.parameter.vie_init;
                                tmp=load([fil '_scan.mat']);scan=tmp.scan;
                                tmp=load([fil '_antenna.mat']);antenna=tmp.antenna;
                                tmp=load([fil '_sources.mat']);sources=tmp.sources;
                            end
                            [antenna,sources,scan,parameter]=vie_mod(antenna, sources,scan,parameter);
                            fil=[pthDALE 'DATA/LEVEL1/' runp.mod_path '/' session];
                            if ~exist([pthDALE 'DATA/LEVEL1/' runp.mod_path], 'dir')
                                mkdir([pthDALE 'DATA/LEVEL1/' runp.mod_path])
                            end
                            savestruct(fil,parameter,antenna,scan,sources);
                        else
                            fprintf('You need to run VIE_INIT for session %s before you can run VIE_MOD\n',session);
                        end
                    end

                    %% VIE_SIM

                    if runp.sim
                        fil=[pthDALE 'DATA/LEVEL1/' runp.mod_path '/' session];
                        if (exist([fil '_scan.mat'], 'file')&&...
                                exist([fil '_antenna.mat'], 'file'))

                            if ~runp.mod
                                tmp=load([fil '_scan.mat']);scan=tmp.scan;
                                tmp=load([fil '_antenna.mat']);antenna=tmp.antenna;
                                tmp=load([fil '_sources.mat']);sources=tmp.sources;
                            end
                            [antenna,scan,sources,session,parameter] = vie_sim(antenna,scan,sources,session,runp.mod_path,parameter);   % Jing SUN, Jan 10, 2012
                        else
                            fprintf('You need to run VIE_MOD for session %s before you can run VIE_SIM\n',session);
                        end
                    end
                    %         if you want to use LEVEL1 data provide the path to the data
                    %         vie_lsm(['../DATA/LEVEL1/',ngsfile(6:19)])

                    %% VIE_LSM

                    %         if you want to take the current structure arrays from the WORK directory use vie_lsm without input parameter
                    if runp.lsm
                        fil=[pthDALE 'DATA/LEVEL1/' runp.mod_path '/' session];
                        if (exist([fil '_parameter.mat'], 'file')&&...
                                exist([fil '_scan.mat'], 'file')&&...
                                exist([fil '_antenna.mat'], 'file')&&...
                                exist([fil '_sources.mat'], 'file'))

                            if ~runp.mod
                                tmp=load([fil '_parameter.mat']);
                                parameter.vie_init=tmp.parameter.vie_init;
                                parameter.vie_mod=tmp.parameter.vie_mod;
                                parameter.eop=tmp.parameter.eop;
                                tmp=load([fil '_scan.mat']);scan=tmp.scan;
                                tmp=load([fil '_antenna.mat']);antenna=tmp.antenna;
                                tmp=load([fil '_sources.mat']);sources=tmp.sources;
                            end
                            vie_lsm(antenna,sources,scan,parameter,runp.lsm_path,runp.glob_path)
                        else
                            fprintf('You need to run VIE_MOD for session %s before you can run VIE_LSM\n', session);
                        end
                        %end % Claudia 22/10/2012

                    elseif runp.lsm_scanwise %Claudia 22/10/2012

                        fil=[pthDALE 'DATA/LEVEL1/' runp.mod_path '/' session];
                        if (exist([fil '_parameter.mat'], 'file')&&...
                                exist([fil '_scan.mat'], 'file')&&...
                                exist([fil '_antenna.mat'], 'file')&&...
                                exist([fil '_sources.mat'], 'file'))

                            if ~runp.mod
                                tmp=load([fil '_parameter.mat']);
                                parameter.vie_init=tmp.parameter.vie_init;
                                parameter.vie_mod=tmp.parameter.vie_mod;
                                parameter.eop=tmp.parameter.eop;
                                tmp=load([fil '_scan.mat']);scan=tmp.scan;
                                tmp=load([fil '_antenna.mat']);antenna=tmp.antenna;
                                tmp=load([fil '_sources.mat']);sources=tmp.sources;
                            end
                            vie_lsm_scanwise(antenna,sources,scan,parameter,runp.lsm_path,runp.glob_path)
                        else
                            fprintf('You need to run VIE_MOD for session %s before you can run VIE_LSM\n', session);
                        end
                    end


                catch exception
                    if runp.error
                        error(getReport(exception, 'extended'))
                    else
                        fprintf(2, getReport(exception, 'extended', 'hyperlinks', 'off')); % display error message
                    end
                    sess_err{isess} = session_name; % save failed session name
                    fail_temp = fopen('failed_sess_temp.txt', 'a');
                    fprintf(fail_temp, '%s\n', sess_err{isess});
                    fclose(fail_temp);
                end
            end % for isess=1:num
        else %--> parallel
            %% parallel processing
            

            parfor isess = 1:number_of_sessions
                %% parameters for each session 
%                 tmp = sessionCounter;
                % Get one instance of the parameter structure for each loop instance:
                parameter = parameter_tmp;

%                 dispSessionNumber(tmp);
                
                fprintf('%s %4.0f %s %4.0f\n','session ',isess,' of ',number_of_sessions)
                session_name = process_list(isess,:);

                while session_name(length(session_name))==' '
                    session_name=session_name(1:length(session_name)-1);
                end

                % Get file format of input data file:
                parameter.data_type = 'ngs';
                if strfind(session_name, ' [vgosDB]')
                    parameter.data_type = 'vgosdb'; 
                elseif strfind(session_name, ' [VSO]')
                    parameter.data_type = 'vso'; 
                end


                %% Get the session name, the fielpath and the year (string)
                switch(parameter.data_type)
                    case 'ngs'
                        flag_absolut_path = false;
                        ind_tmp = [strfind(session_name, '\'), strfind(session_name, '/')];
                        if isempty(ind_tmp)
                            % Format: <session_name>
                            parameter.session_name  = session_name;
                            parameter.filepath      =  '../DATA/NGS/';
                        else
                            parameter.session_name  = session_name(max(ind_tmp) + 1 : end);
                            if (length(ind_tmp) == 1) && (length(session_name(1:min(ind_tmp-1))) == 4)
                                % Format: yyy/<session_name>
                                if ~isnan(str2double(session_name(1:4)))
                                    if (str2double(session_name(1:4)) < 2025) && (str2double(session_name(1:4)) > 1979)
                                        parameter.filepath =  ['../DATA/NGS/', session_name(1:4), '/'];
                                    else
                                        flag_absolut_path = true;
                                    end
                                else
                                    flag_absolut_path = true;
                                end
                            else
                                flag_absolut_path = true;
                            end
                        end

                        % Format: <absolut path>/<session_name>
                        if flag_absolut_path
                            parameter.filepath =  session_name(1 : max(ind_tmp));
                        end

                        year_tmp = str2double(parameter.session_name(1:2));
                        % Check if converversion was sucessfull:
                        if isnan(year_tmp)
                            error('Invalid NGS file name: The first two letters have to represent the year of the session! Please keept the standard naming convention, e.g. "16AUG26XU_N004"!');
                        else
                            % Convert two digit year to four digits!
                            if year_tmp < 79
                                year_tmp = year_tmp + 2000;
                            elseif year_tmp >= 79
                                year_tmp = year_tmp + 1900;
                            end
                        end
                        parameter.year = num2str(year_tmp); % year has to be saved as string!
                        fprintf(' Input file format: NGS\n');

                    case 'vso'
                        flag_absolut_path = false;
                        ind_tmp = [strfind(session_name, '\'), strfind(session_name, '/')];
                        if isempty(ind_tmp)
                            % Format: <session_name>
                            parameter.session_name  = session_name(1 : (strfind(session_name, ' [VSO]')-1));
                            parameter.filepath      =  '../DATA/VSO/';
                            parameter.year          = [];       % year is not included in the session name
                            % NOTE: The name of the folder containing the vso file has to be the session year!
                            % => parameter.year is required to look get the paths to OPT and OUITLIER files in vie_init!
                            error('Invalid filepath of input vso file(%s)! NOTE: The name of the folder containing the vso file has to be the session year!', parameter.filepath)
                        else
                            parameter.session_name  = session_name(max(ind_tmp) + 1 : (strfind(session_name, ' [VSO]')-1));
                            if (length(ind_tmp) == 1) && (length(session_name(1:min(ind_tmp-1))) == 4)
                                % Format: yyy/<session_name>
                                if ~isnan(str2double(session_name(1:4)))
                                    if (str2double(session_name(1:4)) < 2025) && (str2double(session_name(1:4)) > 1979)
                                        parameter.filepath  =  ['../DATA/VSO/', session_name(1:4), '/'];
                                        parameter.year      = session_name(1:4);
                                    else
                                        flag_absolut_path = true;
                                    end
                                else
                                    flag_absolut_path = true;
                                end
                            else
                                flag_absolut_path = true;
                            end
                        end

                        % Format: <absolut path>/<year>/<session_name>
                        % NOTE: The name of the folder containing the vso file has to be the session year!
                        % => parameter.year is required to look get the paths to OPT and OUITLIER files in vie_init!
                        if flag_absolut_path
                            parameter.filepath  = session_name(1 : max(ind_tmp));
                            % Get year:
                            if (length(ind_tmp) >= 2) && (length(session_name(ind_tmp(end-1)+1 : ind_tmp(end)-1)) == 4)
                                if ~isnan(str2double(session_name(ind_tmp(end-1)+1 : ind_tmp(end)-1)))
                                    if (str2double(session_name(ind_tmp(end-1)+1 : ind_tmp(end)-1)) < 2025) && (str2double(session_name(ind_tmp(end-1)+1 : ind_tmp(end)-1)) > 1979)
                                        parameter.year = session_name(ind_tmp(end-1)+1 : ind_tmp(end)-1);
                                    else
                                        error('Invalid filepath of input vso file(%s)! NOTE: The name of the folder containing the vso file has to be the session year!', parameter.filepath)
                                    end
                                else
                                    error('Invalid filepath of input vso file(%s)! NOTE: The name of the folder containing the vso file has to be the session year!', parameter.filepath)
                                end
                            else
                                error('Invalid filepath of input vso file(%s)! NOTE: The name of the folder containing the vso file has to be the session year!', parameter.filepath)
                            end
                        end
                        fprintf(' Input file format: VSO\n');

                    case 'vgosdb'
                        parameter.session_name  = session_name(6 : (strfind(session_name, ' [vgosDB]')-1));
                        parameter.year          = session_name(1:4);
                        parameter.filepath      = ['../DATA/vgosDB/', parameter.year, '/'];
                        fprintf(' Input file format: vgosDB\n');
                end % switch(parameter.data_type)
                session = parameter.session_name;

                try
                    %% VIE_INIT

                    %parameter.vie_init.ngsfile = sess;
                    %save('currentf_pa.mat','parameter');
                    if runp.init
                        fil=[pthDALE 'DATA/LEVEL0/' runp.init_path '/' session];
                        fprintf('Current file: %s\n', fil);
                        if isfield(runp,'ngsdir')
                            [antenna,sources,scan,parameter]=vie_init(session_name,parameter,runp.init_path,runp.ngsdir,trf,crf);
                        else
                            [antenna,sources,scan,parameter]=vie_init(session_name,parameter,runp.init_path,[],trf,crf);
                        end
                        if ~exist([pthDALE 'DATA/LEVEL0/' runp.init_path], 'dir')
                            mkdir([pthDALE 'DATA/LEVEL0/' runp.init_path])
                        end
                        savestruct(fil,parameter,antenna,scan,sources);
                    end

                    %% VIE_MOD
                    if runp.mod

                        fil=[pthDALE 'DATA/LEVEL0/' runp.init_path '/' session];
                        if (exist([fil '_parameter.mat'], 'file')&&...
                                exist([fil '_scan.mat'], 'file')&&...
                                exist([fil '_antenna.mat'], 'file')&&...
                                exist([fil '_sources.mat'], 'file'))

                            if ~runp.init
                                tmp=load([fil '_parameter.mat']);parameter.vie_init=tmp.parameter.vie_init;
                                tmp=load([fil '_scan.mat']);scan=tmp.scan;
                                tmp=load([fil '_antenna.mat']);antenna=tmp.antenna;
                                tmp=load([fil '_sources.mat']);sources=tmp.sources;
                            end
                            [antenna,sources,scan,parameter]=vie_mod(antenna, ...
                                sources,scan,parameter);
                            fil=[pthDALE 'DATA/LEVEL1/' runp.mod_path '/' session];
                            if ~exist([pthDALE 'DATA/LEVEL1/' runp.mod_path], 'dir')
                                mkdir([pthDALE 'DATA/LEVEL1/' runp.mod_path])
                            end
                            savestruct(fil,parameter,antenna,scan,sources);
                        else
                            fprintf('You need to run VIE_INIT for session %s before you can run VIE_MOD\n',session);
                        end
                    end

                    %% VIE_SIM

                    if runp.sim
                        fil=[pthDALE 'DATA/LEVEL1/' runp.mod_path '/' session];
                        if (exist([fil '_scan.mat'], 'file')&&...
                                exist([fil '_antenna.mat'], 'file'))

                            if ~runp.mod
                                tmp=load([fil '_scan.mat']);scan=tmp.scan;
                                tmp=load([fil '_antenna.mat']);antenna=tmp.antenna;
                                tmp=load([fil '_sources.mat']);sources=tmp.sources;
                            end
                            [antenna,scan,sources,session,parameter] = vie_sim(antenna,scan,sources,session,runp.mod_path,parameter);   % Jing SUN, Jan 10, 2012
                        else
                            fprintf('You need to run VIE_MOD for session %s before you can run VIE_SIM\n',session);
                        end
                    end
                    %         if you want to use LEVEL1 data provide the path to the data
                    %         vie_lsm(['../DATA/LEVEL1/',ngsfile(6:19)])

                    %% VIE_LSM
                    % if you want to take the current structure arrays from the WORK directory use vie_lsm without input parameter
                    if runp.lsm
                        fil=[pthDALE 'DATA/LEVEL1/' runp.mod_path '/' session];
                        if (exist([fil '_parameter.mat'], 'file')&&...
                                exist([fil '_scan.mat'], 'file')&&...
                                exist([fil '_antenna.mat'], 'file')&&...
                                exist([fil '_sources.mat'], 'file'))

                            if ~runp.mod
                                tmp=load([fil '_parameter.mat']);
                                parameter.vie_init=tmp.parameter.vie_init;
                                parameter.vie_mod=tmp.parameter.vie_mod;
                                parameter.eop=tmp.parameter.eop;
                                tmp=load([fil '_scan.mat']);scan=tmp.scan;
                                tmp=load([fil '_antenna.mat']);antenna=tmp.antenna;
                                tmp=load([fil '_sources.mat']);sources=tmp.sources;
                            end
                            vie_lsm(antenna,sources,scan,parameter,runp.lsm_path,runp.glob_path)
                        else
                            fprintf('You need to run VIE_MOD for session %s before you can run VIE_LSM\n', session);
                        end

                    elseif runp.lsm_scanwise %Claudia 22/10/2012

                        fil=[pthDALE 'DATA/LEVEL1/' runp.mod_path '/' session];
                        if (exist([fil '_parameter.mat'], 'file')&&...
                                exist([fil '_scan.mat'], 'file')&&...
                                exist([fil '_antenna.mat'], 'file')&&...
                                exist([fil '_sources.mat'], 'file'))

                            if ~runp.mod
                                tmp=load([fil '_parameter.mat']);
                                parameter.vie_init=tmp.parameter.vie_init;
                                parameter.vie_mod=tmp.parameter.vie_mod;
                                parameter.eop=tmp.parameter.eop;
                                tmp=load([fil '_scan.mat']);scan=tmp.scan;
                                tmp=load([fil '_antenna.mat']);antenna=tmp.antenna;
                                tmp=load([fil '_sources.mat']);sources=tmp.sources;
                            end
                            vie_lsm_scanwise(antenna,sources,scan,parameter,runp.lsm_path,runp.glob_path)
                        else
                            fprintf('You need to run VIE_MOD for session %s before you can run VIE_LSM\n', session);
                        end
                    end

                catch exception
                    if runp.error
                        error(getReport(exception, 'extended')) % display error message and abort
                    else
                        fprintf(2, getReport(exception, 'extended', 'hyperlinks', 'off')); % display error message
                    end
                    sess_err{isess} = session_name; % save failed session name

                    fail_temp = -1;
                    while fail_temp == -1 % check if file is open in other worker
                        fail_temp = fopen('failed_sess_temp.txt', 'a');
                    end
                    fprintf(fail_temp, '%s\n', sess_err{isess});
                    fclose(fail_temp);
                end
                
                if ~verLessThan('Matlab','9.2')
                    send(D, isess);
                end
            end % parfor isess = 1:num
            


            % ##### Close parallel computing #####
            if flag_release_r2013b_or_later == 0
                matlabpool close
            elseif flag_release_r2013b_or_later == 1
%                 delete(poolobj)
            end
            if ~verLessThan('Matlab','9.2')
                if isvalid(h)
                    delete(h);
                end
            end
        end % if parallel
        
    end
    
    %% VIE_GLOB

    if runp.glob
        vie_glob
    end

    if exist('process_list', 'var')
        process_list_orig = process_list;
		clear process_list;
        count = 1;
        flag = 0;
        for i_err = 1: length(sess_err) % display failed sessions
            if ~isempty(sess_err{i_err})
                fprintf(2, 'sessions %s produced an error\n', sess_err{i_err});
                process_list(count,:) = sess_err{i_err};
                count = count + 1;
                flag = 1;
            end
        end
        if flag
            save('PROCESSLIST/failed_sessions.mat','process_list');
            fprintf(2, 'a process list with all failed session was saved in PROCESSLIST/failed_sessions.mat\n');
        end
    end
    if exist('failed_sess_temp.txt', 'file') == 2
        delete 'failed_sess_temp.txt'; % delete temporary file
    end
end

if runp.lsm && runp.init && runp.sim && runp.mod
    [t_mean_sig, t_rep] = analyse_simulations(runp.lsm_path);
    
    load ('../DATA/LEVEL4/simparam.mat','simparam');
    if(~isempty(simparam.pathToStatisticsFile))
        try
            updateVieSchedppStatisticsFile(simparam.pathToStatisticsFile, t_mean_sig, t_rep);
        catch exception
            error(getReport(exception, 'extended')) % display error message and abort
        end
    end
end

if isfield(runp, 'sched') && runp.sched && exist('process_list_orig', 'var') && size(process_list_orig,1)>1 && runp.lsm && runp.init && runp.sim && runp.mod
    try
        writeLEVEL3toXLSX( runp );
    catch exception
        error(getReport(exception, 'extended')) % display error message and abort
    end
end

t = clock;
tstr1 = sprintf('%4d/%02d/%02d', t(1), t(2), t(3));
tstr2 = sprintf('%02d:%02d:%02d', t(4), t(5), round(t(6)));
fprintf('VieVS processing ends at %s, %s\n', tstr2, tstr1);
cleanVgosDB();

runtime = toc(startTime);
fprintf('VieVS runtime: %d seconds (%.2f hours)\n', uint64(runtime), runtime/3600)

end % function vie_batch


%% ##### local functions #####

function savestruct(fil,parameter,antenna,scan,sources)
    save([fil '_antenna.mat'],'antenna');
    save([fil '_sources.mat'],'sources');
    save([fil '_scan.mat'],'scan');
    save([fil '_parameter.mat'],'parameter');
end

function dispSessionNumber(number, total)
    
    t = number;
    v = [];
    while t > 0
        v = [mod(t,10) v];
        t = floor(t/10);
    end

    t = total;
    tot = [];
    while t > 0
        tot = [mod(t,10) tot];
        t = floor(t/10);
    end


    string0 = '';
    string1 = '';
    string2 = '';
    string3 = '';
    string4 = '';
    for i = 1:length(v)
        d = v(i);
        
        switch d
            case 0
                string0 = [string0 '  ___  '];
                string1 = [string1 ' / _ \\ '];
                string2 = [string2 '| | | |'];
                string3 = [string3 '| |_| |'];
                string4 = [string4 ' \\___/ '];

            case 1
                string0 = [string0 ' _ '];
                string1 = [string1 '/ |'];
                string2 = [string2 '| |'];
                string3 = [string3 '| |'];
                string4 = [string4 '|_|'];
            case 2
                string0 = [string0 ' ____  '];
                string1 = [string1 '|___ \\ '];
                string2 = [string2 '  __) |'];
                string3 = [string3 ' / __/ '];
                string4 = [string4 '|_____|'];
            case 3
                string0 = [string0 ' _____ '];
                string1 = [string1 '|___ / '];
                string2 = [string2 '  |_ \\ '];
                string3 = [string3 ' ___) |'];
                string4 = [string4 '|____/ '];
            case 4
                string0 = [string0 ' _  _   '];
                string1 = [string1 '| || |  '];
                string2 = [string2 '| || |_ '];
                string3 = [string3 '|__   _|'];
                string4 = [string4 '   |_|  '];
            case 5
                string0 = [string0 ' ____  '];
                string1 = [string1 '| ___| '];
                string2 = [string2 '|___ \\ '];
                string3 = [string3 ' ___) |'];
                string4 = [string4 '|____/ '];
            case 6
                string0 = [string0 '  __   '];
                string1 = [string1 ' / /_  '];
                string2 = [string2 '| ''_ \\ '];
                string3 = [string3 '| (_) |'];
                string4 = [string4 ' \\___/ '];
            case 7
                string0 = [string0 ' _____ '];
                string1 = [string1 '|___  |'];
                string2 = [string2 '   / / '];
                string3 = [string3 '  / /  '];
                string4 = [string4 ' /_/   '];
            case 8
                string0 = [string0 '  ___  '];
                string1 = [string1 ' ( _ ) '];
                string2 = [string2 ' / _ \\ '];
                string3 = [string3 '| (_) |'];
                string4 = [string4 ' \\___/ '];
            case 9
                string0 = [string0 '  ___  '];
                string1 = [string1 ' / _ \\ '];
                string2 = [string2 '| (_) |'];
                string3 = [string3 ' \\__, |'];
                string4 = [string4 '   /_/ '];
        end
    end
    
    string0 = [string0 '     __ '];
    string1 = [string1 '    / / '];
    string2 = [string2 '   / /  '];
    string3 = [string3 '  / /   '];
    string4 = [string4 ' /_/    '];

    for i = 1:length(tot)
        d = tot(i);
        
        switch d
            case 0
                string0 = [string0 '  ___  '];
                string1 = [string1 ' / _ \\ '];
                string2 = [string2 '| | | |'];
                string3 = [string3 '| |_| |'];
                string4 = [string4 ' \\___/ '];

            case 1
                string0 = [string0 ' _ '];
                string1 = [string1 '/ |'];
                string2 = [string2 '| |'];
                string3 = [string3 '| |'];
                string4 = [string4 '|_|'];
            case 2
                string0 = [string0 ' ____  '];
                string1 = [string1 '|___ \\ '];
                string2 = [string2 '  __) |'];
                string3 = [string3 ' / __/ '];
                string4 = [string4 '|_____|'];
            case 3
                string0 = [string0 ' _____ '];
                string1 = [string1 '|___ / '];
                string2 = [string2 '  |_ \\ '];
                string3 = [string3 ' ___) |'];
                string4 = [string4 '|____/ '];
            case 4
                string0 = [string0 ' _  _   '];
                string1 = [string1 '| || |  '];
                string2 = [string2 '| || |_ '];
                string3 = [string3 '|__   _|'];
                string4 = [string4 '   |_|  '];
            case 5
                string0 = [string0 ' ____  '];
                string1 = [string1 '| ___| '];
                string2 = [string2 '|___ \\ '];
                string3 = [string3 ' ___) |'];
                string4 = [string4 '|____/ '];
            case 6
                string0 = [string0 '  __   '];
                string1 = [string1 ' / /_  '];
                string2 = [string2 '| ''_ \\ '];
                string3 = [string3 '| (_) |'];
                string4 = [string4 ' \\___/ '];
            case 7
                string0 = [string0 ' _____ '];
                string1 = [string1 '|___  |'];
                string2 = [string2 '   / / '];
                string3 = [string3 '  / /  '];
                string4 = [string4 ' /_/   '];
            case 8
                string0 = [string0 '  ___  '];
                string1 = [string1 ' ( _ ) '];
                string2 = [string2 ' / _ \\ '];
                string3 = [string3 '| (_) |'];
                string4 = [string4 ' \\___/ '];
            case 9
                string0 = [string0 '  ___  '];
                string1 = [string1 ' / _ \\ '];
                string2 = [string2 '| (_) |'];
                string3 = [string3 ' \\__, |'];
                string4 = [string4 '   /_/ '];
        end
    end
    
    string = [string0 '\n' string1 '\n' string2 '\n' string3 '\n' string4 '\n'];
    
    fprintf(string);
end


