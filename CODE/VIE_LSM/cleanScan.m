% ************************************************************************
%   Description:
%   This function cleans the scan structure because of possible exclusions
%   of stations/sources/baselines.
%
%   Input:	
%
%   Output:
%      Scan, sources, antenna structures.
% 
%   External calls: 	
%       
%   Coded for VieVS: 
%   Jul 2012 by Matthias Madzak
%
%   Revision: 
%   dd mmm yyyy by FIRSTNAME SECONDNAME:
%   18 Jul 2016 by D. Mayer: include fistObsMjd and lastObsMjd into source
%   struct
%   13 Oct 2016 by A. Girdiuk: bug-fix: sources can be excluded when time span is written down to OPT-file
%   13 Nov 2017 by J. Gruber: Update to use sessions where several stations are excluded
%   06 Jan 2018 by J. Gruber: bug fix in the outlier section (it is now possible to remove outliers)
%   16 Jan 2018 by J. Gruber: bug fix in "No cable calibration" section
%   16 Jan 2018 by J. Gruber: Ion code flag is also considered, only ion code zero values are taken
%   17 Jan 2018 by J. Gruber: bug fix in "Excluded baselines" section
%   28 Aug 2018 by D. Landskron: shape of output slightly changed
%   28 Nov 2018 by D. Landskron: observations restrictions clarified
%   05 Dec 2018 by D. Landskron: clarification quality code / quality flag

% ************************************************************************
function [scan, sources, antenna]=cleanScan(scan, sources, antenna, parameter)

%% Options
flag_pring_debug_info = false;

%% Precalculations
nScans = size(scan,2);

% Stations:
nStats = length(antenna);
allStationNames = {antenna.name}';
allStationNames = strtrim(allStationNames);

% Sources:
% - Quasars
nSources_q = size(sources.q, 2);
if ~isempty(sources.q)
    allSourceNames_q  = {sources.q.name}';
    allSourceNames_q  = strtrim(allSourceNames_q);
else
    allSourceNames_q = {};
end

% - Spacecrafts/Satellites:
nSources_s = size(sources.s, 2);
if ~isempty(sources.s)
    allSourceNames_s  = {sources.s.name}';
    allSourceNames_s  = strtrim(allSourceNames_s);
else
    allSourceNames_s = {};
end

oneSecInDays = 1/60/60/24;

% get end-epoch of exclude stations
if ~isempty(parameter.opt.options.sta_excl)
    if ~isfield(parameter.opt.options, 'sta_excl_end')
        if isempty(parameter.opt.options.sta_excl_start)
            parameter.opt.options.sta_excl_end = '';
        else
            % if there is a start epoch
            parameter.opt.options.sta_excl_end = ones(size(parameter.opt.options.sta_excl_start))*9999999;
        end
    else 
        % there is an end-epoch-field
        % make 0 -> 99999
        parameter.opt.options.sta_excl_end(parameter.opt.options.sta_excl_end == 0) = 9999999;
    end
    exclStaCellStr = cellstr(parameter.opt.options.sta_excl);
end



%% Clean scan struct

% (1) No cable calibration
if ~isempty(parameter.opt.options.no_cab)
    
    % get antenna numbers
    antIndNoCab = find(ismember(strtrim({antenna.name}), cellstr(parameter.opt.options.no_cab)));
    
    if ~isempty(antIndNoCab)
        % for all scans
        for iS = 1 : nScans
            for kAntNoCab = 1 : length(antIndNoCab)
                curI = antIndNoCab(kAntNoCab); % current index of antenna with no cable-cal
                % if that station acutally "exists" (if the index can be accessed
                if length(scan(iS).stat) >= curI
                    curCab = scan(iS).stat(curI).cab * (1e-9);

                    obsI1NoCab = find([scan(iS).obs.i1]==curI);
                    obsI2NoCab = find([scan(iS).obs.i2]==curI);
                    
                    % Remove cable calibration which was already applied in VIe_INIT when loading the observation file:
                    if ~isempty(curCab)
                        for ionc = 1 : length(obsI1NoCab)
                            scan(iS).obs(obsI1NoCab(ionc)).obs = scan(iS).obs(obsI1NoCab(ionc)).obs + curCab;
                        end
                        for ionc = 1 : length(obsI2NoCab)
                            scan(iS).obs(obsI2NoCab(ionc)).obs = scan(iS).obs(obsI2NoCab(ionc)).obs - curCab;
                        end
                    end

                    % set cable cal to zero!
                    scan(iS).stat(curI).cab = 0;
                end
            end
        end
    end
end
    

% (2) Excluded baselines
sum_del_baseline = 0;
if ~isempty(parameter.opt.options.bas_excl)
   
    % Loop over all baselines
    for iBasel = 1 : length(parameter.opt.options.bas_excl)
        
        % get current baseline        
        curStatLog = strcmpi(allStationNames, strtrim(parameter.opt.options.bas_excl(iBasel).sta1)) | strcmpi(allStationNames, strtrim(parameter.opt.options.bas_excl(iBasel).sta2));

        
        % get indices of both current stations
        curStatInd = find(curStatLog);
        
        % go through all scans and remove stuff
        for iScan = 1 : nScans
            
            % only if both station-structs exist
            if max(curStatInd) <= size(scan(iScan).stat, 2)
                
                % if both stations take part in the session
                if length([scan(iScan).stat(curStatInd).temp]) == 2
                    
                    % get observations to delete (might not be the case)
                    A = [[scan(iScan).obs.i1]', [scan(iScan).obs.i2]'] == curStatInd(1) | [[scan(iScan).obs.i1]', [scan(iScan).obs.i2]'] == curStatInd(2);
                    ia1 = A(:,1) == 1;
                    ia2 = A(:,2) == 1;
                    obs2Delete = ia1 & ia2 ;
                    
                    if sum(obs2Delete) > 0
                        
                        sum_del_baseline = sum_del_baseline + sum(obs2Delete);
                        
                        % then we have to delete something
                        scan(iScan).obs(obs2Delete) = [];
                        scan(iScan).nobs = scan(iScan).nobs - 1;
                        
                        if sum(obs2Delete) > 1 % ACHTUNG! Test...
                            disp('Kann nicht sein!')
                            keyboard
                        end
                        
                        statStructs2Delete_1 = [];
                        statStructs2Delete_2 = [];
                        
                        % delete stat struct if needed (if that station has no observation anymore (for both stats separate)
                        if sum(sum([[scan(iScan).obs.i1]', [scan(iScan).obs.i2]'] == curStatInd(1))) == 0
                            statStructs2Delete_1 = curStatInd(1);
                        end
                        if sum(sum([[scan(iScan).obs.i1]', [scan(iScan).obs.i2]'] == curStatInd(2))) == 0
                            statStructs2Delete_2 = curStatInd(2);
                        end
                        empty_struct = get_empty_fields_struct(scan(iScan).stat(statStructs2Delete_1));
                        scan(iScan).stat(statStructs2Delete_1) = empty_struct;
                        empty_struct = get_empty_fields_struct(scan(iScan).stat(statStructs2Delete_2));
                        scan(iScan).stat(statStructs2Delete_2) = empty_struct;
                                                
                    end
                end
            end
        end
    end    
end
if flag_pring_debug_info
    fprintf(' - Observations removed due excluded baselines: %d\n', sum_del_baseline);
end


% (3) Excluded stations
nStat2excl = size(parameter.opt.options.sta_excl, 1);
if ~isempty(parameter.opt.options.sta_excl)
    
    % get station logicals to be deleted
    curStatLog = zeros(nStats, 1);
    for iStatExcl = 1 : nStat2excl
        foundLogsCurAnt = strcmpi({antenna.name}, parameter.opt.options.sta_excl(iStatExcl, :));
        if sum(foundLogsCurAnt) ~= 1
            fprintf('%s (to be excluded) not found in antenna struct\n-> not excluded!\n',parameter.opt.options.sta_excl(iStatExcl,:));
        else
            curStatLog = curStatLog(:) | foundLogsCurAnt(:);
        end
    end
    curStatInd = find(curStatLog);
    
    for iScan = 1 : nScans
        % get station indices which exist in current scan
        curStatIndInCurScan = curStatInd(curStatInd <= size(scan(iScan).stat, 2));
        
        if ~isempty(curStatIndInCurScan)
            
            % go through all antennas separately (due to time of exclusion!)
            for iStatInCurScan = 1 : length(curStatIndInCurScan)
                curSingleStatIndInCurScan = curStatIndInCurScan(iStatInCurScan);
                curAnt = strtrim(antenna(curSingleStatIndInCurScan).name); % antenna name
                curSTART = parameter.opt.options.sta_excl_start(strcmpi(exclStaCellStr,curAnt));
                curEND = parameter.opt.options.sta_excl_end(strcmpi(exclStaCellStr,curAnt)); % mjd
                
                
                for iMultStat = 1:length(curSTART) % if same station is excluded multiple times
                    curStart    = curSTART(iMultStat);
                    curEnd      = curEND(iMultStat);
                    % check if scan is in "exclude-time"
                    if (scan(iScan).mjd > curStart) && (scan(iScan).mjd < curEnd)

                        % make scan-substruct .stat empty
                        empty_struct = get_empty_fields_struct(scan(iScan).stat(curSingleStatIndInCurScan));
                        scan(iScan).stat(curSingleStatIndInCurScan) = empty_struct;

                        % delete observations
                        obs2Delete = logical(sum(ismember([scan(iScan).obs.i1; scan(iScan).obs.i2], curSingleStatIndInCurScan)));
                        scan(iScan).obs(obs2Delete) = [];

                        % decrease number of observations
                        scan(iScan).nobs = scan(iScan).nobs - sum(obs2Delete);
                        % still missing: if stats in current scan have no obs left -> make also empty => This is done at the end of this script!!
                    end
                end
            end
            
        end
    end
end


% (4) Excluded sources (only get scan indices => scans are deleted later

% preallocate
exclSourcesInd_q = zeros(size(parameter.opt.options.sour_excl,1), 1); % quasars
exclSourcesInd_s = zeros(size(parameter.opt.options.sour_excl,1), 1); % sc

exclSourcesInd_byTime_q = false(1, nScans); % Number of scans, quasars
exclSourcesInd_byTime_s = false(1, nScans); % Number of scans, sc

if ~isempty(parameter.opt.options.sour_excl)
    
    % get source logicals to be deleted
    for iSource2BeExcl = 1 : size(parameter.opt.options.sour_excl, 1)
        
        % Get source IDs:
        if ~isempty(allSourceNames_q) && (sum(strcmp(allSourceNames_q, strtrim(parameter.opt.options.sour_excl(iSource2BeExcl,:)))) > 0)
            exclSourcesInd_q(iSource2BeExcl) = find(strcmp(allSourceNames_q, strtrim(parameter.opt.options.sour_excl(iSource2BeExcl,:))));
        end
        if ~isempty(allSourceNames_s) && (sum(strcmp(allSourceNames_s, strtrim(parameter.opt.options.sour_excl(iSource2BeExcl,:)))) > 0)
            exclSourcesInd_s(iSource2BeExcl) = find(strcmp(allSourceNames_s, strtrim(parameter.opt.options.sour_excl(iSource2BeExcl,:))));
        end
        
        % Get IDs of scans in which the current source was observed:
        scansToExcludedSources_q = ismember([scan.iso], exclSourcesInd_q) & strcmp({scan.obs_type}, 'q'); % quasar
        scansToExcludedSources_s = ismember([scan.iso], exclSourcesInd_s) & strcmp({scan.obs_type}, 's'); % sc
       
        
        % Get scan IDs at which a source is partly excluded (by time):
        if logical(parameter.opt.options.sour_excl_start(iSource2BeExcl)) % If parameter.opt.options.sour_excl_start ~= 0
            
            if ~isempty(allSourceNames_q) && (sum(strcmp(allSourceNames_q, strtrim(parameter.opt.options.sour_excl(iSource2BeExcl,:)))) > 0)
                exclSourcesInd_byTime_q = exclSourcesInd_byTime_q | (scansToExcludedSources_q & (([scan.mjd] >= parameter.opt.options.sour_excl_start(iSource2BeExcl)) & ([scan.mjd]<=parameter.opt.options.sour_excl_end(iSource2BeExcl))));
                exclSourcesInd_q(iSource2BeExcl) = [];
            end
            if ~isempty(allSourceNames_s) && (sum(strcmp(allSourceNames_s, strtrim(parameter.opt.options.sour_excl(iSource2BeExcl,:)))) > 0)
                exclSourcesInd_byTime_s = exclSourcesInd_byTime_s | (scansToExcludedSources_s & (([scan.mjd] >= parameter.opt.options.sour_excl_start(iSource2BeExcl)) & ([scan.mjd]<=parameter.opt.options.sour_excl_end(iSource2BeExcl))));
                exclSourcesInd_s(iSource2BeExcl) = [];
            end
        end
    end
end


% (5) Outliers
sum_del_outliers = 0;
if parameter.outlier.flag_remove_outlier
    if ~isempty(parameter.outlier.obs2remove)
        n_outlier = length(parameter.outlier.obs2remove); % Number of outliers 
        if ~isfield(parameter.outlier.obs2remove,'sou')
            parameter.outlier.obs2remove(1).sou = '';
        end
        
        for iOutlier = 1 : n_outlier
            
            % Get indices of both stations in baseline
            curStatLog = strcmpi(allStationNames, strtrim({parameter.outlier.obs2remove(iOutlier).sta1})) | strcmpi(allStationNames,strtrim({parameter.outlier.obs2remove(iOutlier).sta2})); 
            curStatInd = find(curStatLog);

            % get scan of cur outlier (the one which is close by 1/10 second!)
            curScanLog = abs([scan.mjd] - parameter.outlier.obs2remove(iOutlier).mjd) < (oneSecInDays/10);
            curSouLog = parameter.outlier.obs2remove(iOutlier).sou;
            
            % Check, if only one scan was found!
            % - If more than one scan was found by matching the scan reference times, the stations have to be considered in addition
            flag_found_scan = true; 
            if sum(curScanLog) > 1 % More than one scan found?
                curScanLog_ids = find(curScanLog);
                flag_found_scan = false;
                if isempty(curSouLog) % for old outlier files (without source names)
                    % Check stations in scan:
                    for i_scan = 1 : sum(curScanLog)
                        scan_id = curScanLog_ids(i_scan);
                        stat_ids_in_scan = unique([scan(scan_id).obs.i1, scan(scan_id).obs.i2]);
                        if ( ismember(curStatInd(1), stat_ids_in_scan) && ismember(curStatInd(2), stat_ids_in_scan) )
                           flag_found_scan = true; 
                           break
                        end
                    end
                    if ~flag_found_scan % No observation with the current baseline found in the considered scans
                        fprintf('WARNING (outlier removal in cleanScan.m): No observation on the baseline %s - %s found at epoch %s!\n', parameter.outlier.obs2remove(iOutlier).sta1, parameter.outlier.obs2remove(iOutlier).sta2, mjd2datestr(parameter.outlier.obs2remove(iOutlier).mjd) )
                    else
                        cur_scan_id = scan_id;
                    end
                else % new outlier files with source names
                    for i_scan = 1:size(curScanLog_ids,2)
                       scan_id = curScanLog_ids(i_scan);
                       if strcmp (allSourceNames_q(scan(scan_id).iso),curSouLog)
                            cur_scan_id = scan_id;
                            flag_found_scan = true; 
                       end
                    end
                end
            else
                cur_scan_id = find(curScanLog);
            end

            if flag_found_scan
                % check if there is (still) an obs of the two stations in cur scan
                obs2Delete = sum([scan(cur_scan_id).obs.i1; scan(cur_scan_id).obs.i2] == curStatInd(1) | [scan(cur_scan_id).obs.i1; scan(cur_scan_id).obs.i2] == curStatInd(2))==2;
            else
                obs2Delete = [];
            end

            if sum(obs2Delete) == 1
                sum_del_outliers = sum_del_outliers + 1;
                % then we have to delete something
                scan(cur_scan_id).obs(obs2Delete) = [];
                scan(cur_scan_id).nobs = scan(cur_scan_id).nobs - 1;

                statStructs2Delete_1 = [];
                statStructs2Delete_2 = [];
                        
                % delete stat struct if needed (if that station has no observation anymore (for both stats separate)
                if sum(sum([scan(cur_scan_id).obs.i1; scan(cur_scan_id).obs.i2]==curStatInd(1)))==0 % If no observations for station with ID "curStatInd(1)" left => Delete station
                    statStructs2Delete_1 = curStatInd(1);
                end
                if sum(sum([scan(cur_scan_id).obs.i1; scan(cur_scan_id).obs.i2]==curStatInd(2)))==0 % If no observations for station with ID "curStatInd(1)" left => Delete station
                    statStructs2Delete_2 = curStatInd(2);
                end
                
                empty_struct = get_empty_fields_struct(scan(cur_scan_id).stat(statStructs2Delete_1));
                scan(cur_scan_id).stat(statStructs2Delete_1) = empty_struct;
                empty_struct = get_empty_fields_struct(scan(cur_scan_id).stat(statStructs2Delete_2));
                scan(cur_scan_id).stat(statStructs2Delete_2) = empty_struct;

            elseif obs2Delete > 1
                error('%d observations matche the following outlier: baseline %s - %s at MJD %f', obs2Delete, parameter.outlier.obs2remove(iOutlier).sta1, parameter.outlier.obs2remove(iOutlier).sta2, parameter.outlier.obs2remove(iOutlier).mjd)
            end
        end
    end
end
if flag_pring_debug_info
	fprintf(' - Observations removed due to outliers: %d\n', sum_del_outliers)
end


% (6) Quality of observation
sum_del_q_limit = 0;
for iScan = 1 : nScans
    
    qualOfObsOfCurScan_deluflag = [scan(iScan).obs.q_flag];    
    qualOfObsOfCurScan_IonFlag = [scan(iScan).obs.q_flag_ion];
    obs2Delete_deluflag = qualOfObsOfCurScan_deluflag > parameter.obs_restrictions.Qlim; 
    obs2Delete_IonFlag =  qualOfObsOfCurScan_IonFlag ~= 0; % fixed to zero
    
    % due to a bug by GSFC, for all observations before December 2000 also the quality code (not only the quality flag) must be checked. As long as they do not fix this bug and upload ALL vgosDB files before this date anew, this following step must be done
    % as soon as GSFC has fixed the bug and updated ALL vgosDB sessions from 1979-2000, please remove the section inside the "if" again and leave only the part after "else". Also remove the fprintf part with the variable "code_print_info" below
    if scan(1).tim(1) < 2001   &&   isfield(scan(iScan).obs,'q_code_X')   &&   isfield(scan(iScan).obs,'q_code_S')
        qualOfObsOfCurScan_q_code_X = str2double({scan(iScan).obs.q_code_X});
        qualOfObsOfCurScan_q_code_S = str2double({scan(iScan).obs.q_code_S});
        obs2Delete_q_code_X = ~(qualOfObsOfCurScan_q_code_X > 0);
        obs2Delete_q_code_S = ~(qualOfObsOfCurScan_q_code_S > 0);
        
        obs2Delete = obs2Delete_deluflag | obs2Delete_IonFlag | obs2Delete_q_code_X | obs2Delete_q_code_S;
        code_print_info = 1;
    else
        obs2Delete = obs2Delete_deluflag | obs2Delete_IonFlag;
    end
    
    if sum(obs2Delete) > 0
        sum_del_q_limit = sum_del_q_limit + sum(obs2Delete);
        scan(iScan).obs(obs2Delete) = [];
        scan(iScan).nobs=scan(iScan).nobs - sum(obs2Delete);
        
        % if there are no observations left - scan will be deleted later
        
        % delete stat (sub-)structs which are not needed anymore
        indStatStruct = 1 : size(scan(iScan).stat,2);

        % delete stat structs also if needed
        indStatStruct([[scan(iScan).obs.i1]'; [scan(iScan).obs.i2]']) = []; % those which (still) have index in i1|i2 should not be deleted
        empty_struct = get_empty_fields_struct(scan(iScan).stat(indStatStruct));
        scan(iScan).stat(indStatStruct) = empty_struct;
    end   
end
if flag_pring_debug_info
    fprintf(' - Observations removed due to quality flag limit (%d): %d\n', parameter.obs_restrictions.Qlim, sum_del_q_limit)
end

if exist('code_print_info','var')
    fprintf('%s\n%s\n\n','For all vgosDB sessions before 2001, a special handling is necessary due to a bug by GSFC. Please use NGS files instead.','Please remove this hardcoded part as soon as GSFC has updated all vgosDB files from 1979-2000 and the chis-squared are in a reasonable range..');
end        

% (7) Minimum elevation angle
sum_del_cut_off_elev = 0;
if parameter.obs_restrictions.cut_off_elev > 0
    
    for iScan = 1 : nScans
        % get stations in current scan
        allStationsInCurScan = unique([[scan(iScan).obs.i1], [scan(iScan).obs.i2]]);
        nStatsInCurScan = size(allStationsInCurScan, 2);
        
        % preallocate logicals if station(s) should be deleted in current scan
        deleteStatInCurScanLogicals = zeros(nStatsInCurScan, 1);
        
        % go through all stations and calc elevations
        for iStat = 1 : length(allStationsInCurScan)
            curStat = allStationsInCurScan(iStat);
            
            if ~isempty(scan(iScan).stat(curStat).zd)
               
                % Get elevation calculated in VIE_MOD:
                curElevation = pi/2 - scan(iScan).stat(curStat).zd; % in [rad]
                if curElevation < parameter.obs_restrictions.cut_off_elev
                    % delete this observation
                    deleteStatInCurScanLogicals(iStat) = 1;
                end
            end
        end % for all stations in cur scan
        
        % now delete bad observations
        if sum(deleteStatInCurScanLogicals) > 0
            % get observations to be deleted
            obsWithBadStation = ismember([[scan(iScan).obs.i1];[scan(iScan).obs.i2]], allStationsInCurScan(logical(deleteStatInCurScanLogicals)));
            obs2Delete = sum(obsWithBadStation,1) > 0;
            
            sum_del_cut_off_elev = sum_del_cut_off_elev + sum(obs2Delete);
            
            % delete observations which include one of stations with bad elevation 
            scan(iScan).obs(obs2Delete) = [];
            scan(iScan).nobs=scan(iScan).nobs - sum(obs2Delete);

            % if there are no observations left - scan will be deleted later anyway
            if scan(iScan).nobs > 0
                indStatStruct = 1:size(scan(iScan).stat,2);
                % delete stat structs also if needed
                indStatStruct([[scan(iScan).obs.i1]';[scan(iScan).obs.i2]']) = []; % those which have index in i1|i2 should not be deleted
                empty_struct = get_empty_fields_struct( scan(iScan).stat(indStatStruct));
                scan(iScan).stat(indStatStruct) = empty_struct;
            end    
        end
    end % all scans
end
if flag_pring_debug_info
    fprintf(' - Observations removed due to cut-off elevation (%f deg): %d\n', parameter.obs_restrictions.cut_off_elev * 180/pi, sum_del_cut_off_elev)
end



%% Delete empty scans and those to an excluded source
scansToExcludedSources_q = ismember([scan.iso], exclSourcesInd_q) & strcmp({scan.obs_type}, 'q'); % quasar
scansToExcludedSources_s = ismember([scan.iso], exclSourcesInd_s) & strcmp({scan.obs_type}, 's'); % sc
emptyScans = [scan.nobs] <= 0; % < ... just to be sure!

if flag_pring_debug_info
    fprintf(' - Observations removed due to empty scans: %d\n', sum([scan(emptyScans).nobs])); % should be = 0!
    fprintf(' - Observations removed due to sources: %d\n', sum([scan(scansToExcludedSources_q).nobs]));
end

% Remove scans:
scan( (scansToExcludedSources_q | exclSourcesInd_byTime_q | scansToExcludedSources_s | exclSourcesInd_byTime_s) | emptyScans) = [];




%% UPDATE INDICES OF ANTENNAS AND SOURCES (in scan struct)
% only when something was deleted
if sum([sum(scansToExcludedSources_s | exclSourcesInd_byTime_s), sum(scansToExcludedSources_q | exclSourcesInd_byTime_q), sum(emptyScans), nStat2excl]) > 0 
    
    % number of observations per station (one line=one station)
    nObsPerStat     = zeros(size(antenna, 2), 1);
    nObsPerSource_q   = zeros(size(sources.q, 2), 1);
    nObsPerSource_s   = zeros(size(sources.s, 2), 1);
    
    % Get number of observations per source and per station
    % => Delete those with no observations left!
    for iScan = 1 : size(scan, 2)
        
        % Stations:
        [hista,histb] = hist([scan(iScan).obs.i1, scan(iScan).obs.i2], 1:nStats);
        nObsPerStat = nObsPerStat(:) + hista(:);
        % nObsPerStat should be equal [antenna.numobs]
        
        % Sources:
        if strcmp(scan(iScan).obs_type, 'q') % quasar scan
            nObsPerSource_q(scan(iScan).iso) = nObsPerSource_q(scan(iScan).iso) + scan(iScan).nobs;

        elseif strcmp(scan(iScan).obs_type, 's') % sc scan
            nObsPerSource_s(scan(iScan).iso) = nObsPerSource_s(scan(iScan).iso) + scan(iScan).nobs;
        else
            error('Unknown soruce type: %s', scan(iScan).obs_type);
        end
    end
    
    % get new indices for antennas: e.g. antennas [1,2,3,4,5,6,7,8] -> 2,4,5,7 excluded --> [1,x,2,x,x,3,x,4]
    newIndForStat = 1 : nStats;
    % for all stations
    for iStat=1:nStats
        if nObsPerStat(iStat) == 0
            newIndForStat(iStat) = 0;
            % only if we are not at last station
            if iStat < nStats
                newIndForStat(iStat + 1:end) = newIndForStat(iStat + 1:end) - 1;
            end
        end
    end

    % get new indices for sources (quasars)
    if nSources_q > 0
        newIndForSource_q = 1 : nSources_q;
        for iSource = 1 : nSources_q
            if nObsPerSource_q(iSource)  == 0
                newIndForSource_q(iSource) = 0;
                % only if we are not at last station
                if iSource < nSources_q
                    newIndForSource_q(iSource + 1:end) = newIndForSource_q(iSource + 1:end) - 1;
                end
            end
        end
    end
    
    % get new indices for sources (sc)
    if nSources_s > 0
        newIndForSource_s = 1 : nSources_s;
        for iSource = 1 : nSources_s
            if nObsPerSource_s(iSource)  == 0
                newIndForSource_s(iSource) = 0;
                % only if we are not at last station
                if iSource < nSources_s
                    newIndForSource_s(iSource + 1:end) = newIndForSource_s(iSource + 1:end) - 1;
                end
            end
        end
    end
        
    % if there is an antenna with no obs
    % => Update station IDs
    if sum(nObsPerStat == 0) > 0       
        statWithNoObs = (nObsPerStat == 0);
 
        % go through all scans (and obs) and update to correct station index
        for iScan = 1 : size(scan, 2)
            
            % Update IDs
            % for all obs entries (with deal I don't know how)
            for iObs = 1 : size(scan(iScan).obs, 2)
                scan(iScan).obs(iObs).i1 = newIndForStat(scan(iScan).obs(iObs).i1);
                scan(iScan).obs(iObs).i2 = newIndForStat(scan(iScan).obs(iObs).i2);
            end
            
            % delete also scan(iScan).stat(...) for stations to be deleted
            scan(iScan).stat(statWithNoObs(1:size(scan(iScan).stat,2))) = [];
           
        end % for all scans
        
        % delete antenna entries (antennas with no obs)
        antenna(statWithNoObs) = [];
    end          
    
    % update number of observations per station
    newNumObsPerStat = nObsPerStat(nObsPerStat ~= 0);
    for iStat = 1 : length(antenna)   
        antenna(iStat).numobs = newNumObsPerStat(iStat);
    end
            
    % Delete sources and update IDs
    % Quasars:
    if nSources_q > 0
        if sum(nObsPerSource_q == 0) > 0
            sourWithNoObs = nObsPerSource_q == 0;

            % Update
            for iScan = 1 : size(scan, 2)
                scan(iScan).iso = newIndForSource_q(scan(iScan).iso);
            end

            % delete sources to which no observations are taken
            sources.q(sourWithNoObs) = [];   
        end

        % update number of observations per source
        newNumObsPerSou_q = nObsPerSource_q(nObsPerSource_q ~= 0);
        for iSou = 1 : length(sources.q)
            sources.q(iSou).numobs = newNumObsPerSou_q(iSou);
        end
    end
    
    % Satellites:
    if nSources_s > 0
        if sum(nObsPerSource_s == 0) > 0
            sourWithNoObs = nObsPerSource_s == 0;

            % Update
            for iScan = 1 : size(scan, 2)
                scan(iScan).iso = newIndForSource_s(scan(iScan).iso);
            end

            % delete sources to which no observations are taken
            sources.s(sourWithNoObs) = [];   
        end

        % update number of observations per source
        newNumObsPerSou_s = nObsPerSource_s(nObsPerSource_s ~= 0);
        for iSou = 1 : length(sources.s)
            sources.s(iSou).numobs = newNumObsPerSou_s(iSou);
        end
    end
    
end

%% update 'firstObsMjd' 'lastObsMjd' of antenna structs and add both fields to sources struct
stat_firstObsMjd_checked = zeros(length(antenna), 1);
src_q_firstObsMjd_checked = zeros(length(sources.q), 1);
src_s_firstObsMjd_checked = zeros(length(sources.s), 1);

for iScan = 1 : length(scan)
    statIndOfCurScan = unique([[scan(iScan).obs.i1], [scan(iScan).obs.i2]]);
    
    % Stations:
    [antenna(statIndOfCurScan(stat_firstObsMjd_checked(statIndOfCurScan) == 0)).firstObsMjd] =  deal(scan(iScan).mjd);
    stat_firstObsMjd_checked(statIndOfCurScan) = 1;
    [antenna(statIndOfCurScan).lastObsMjd] = deal(scan(iScan).mjd);
    
    % Sources:
    srcIndOfCurScan = scan(iScan).iso;
    srcTypeOfCurScan = scan(iScan).obs_type;
    
    if strcmp(srcTypeOfCurScan, 'q') % quasar scan
        if src_q_firstObsMjd_checked(srcIndOfCurScan) == 0
            sources.q(srcIndOfCurScan).firstObsMjd = scan(iScan).mjd;
            src_q_firstObsMjd_checked(srcIndOfCurScan) = 1;
        end
        sources.q(srcIndOfCurScan).lastObsMjd = scan(iScan).mjd;
        
    elseif strcmp(srcTypeOfCurScan, 's') % sc scan
        if src_s_firstObsMjd_checked(srcIndOfCurScan) == 0
            sources.s(srcIndOfCurScan).firstObsMjd = scan(iScan).mjd;
            src_s_firstObsMjd_checked(srcIndOfCurScan) = 1;
        end
        sources.s(srcIndOfCurScan).lastObsMjd = scan(iScan).mjd;
    else
        error('Unknown soruce type: %s', srcTypeOfCurScan);
    end
end


%% ############### subroutines #####################



function empty_struct = get_empty_fields_struct(in_struct)
    field_names = fields(in_struct);
    empty_struct = cell2struct( cell(1, length(field_names)), field_names, 2);
return

