% ************************************************************************
%   Description:
%   This function checks the minimum number of observations required to
%   estimate this parameter. If the given thresholds are not met, the
%   observations are removed from the vievs data structures
% 
%   Currently, only sources (quasars) are considered here.
%
%   Input:	
%      Scan, sources, antenna, parameter
%   Output:
%      Scan, sources, antenna
% 
%   External calls: 	
%   - update_source_antenna_indices.m
%       
%   Coded for VieVS: 
%   19 Dec 2019 by A. Hellerschmied

function [scan, sources, antenna] = check_num_of_obs_per_parameter(scan, sources, antenna, parameter)

%% Options
flag_pring_debug_info = true;


%% Precallocations
% Stations:
nStats = length(antenna);

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





%% Clean scan struct

% (1) ### Excluded sources (only get scan indices => scans are deleted later ###

% Get list of sources for which the threchold is not met:
q_numobs = [sources.q.numobs]'; % Number of obs per source (quasar)
q_flag = q_numobs < parameter.lsmopt.min_num_obs_per_est_source;
exclSourcesInd_q = find(q_flag);



%% Delete empty scans and those to an excluded source
scansToExcludedSources_q = ismember([scan.iso], exclSourcesInd_q) & strcmp({scan.obs_type}, 'q'); % quasar
emptyScans = [scan.nobs] <= 0; % < ... just to be sure!

if flag_pring_debug_info
    fprintf(' - Observations removed due to empty scans: %d\n', sum([scan(emptyScans).nobs])); % should be = 0!
    fprintf(' - Observations removed due to sources: %d\n', sum([scan(scansToExcludedSources_q).nobs]));
end

fprintf(1, ' - %d sources have less than %d observations (%d obs will be removed)\n', sum(q_flag), parameter.lsmopt.min_num_obs_per_est_source, sum(q_numobs(q_flag)));
if flag_pring_debug_info 
    % Print list of removed sources to CW: 
    names_of_removed_sources = allSourceNames_q(q_flag);
    q_numobs_removed_sources = q_numobs(q_flag);
    for i_src = 1 : length(names_of_removed_sources)
       fprintf(1, '   #%3d - %s (obs: %d)\n', i_src, names_of_removed_sources{i_src}, q_numobs_removed_sources(i_src));
    end
end

% Remove scans:
scan(scansToExcludedSources_q | emptyScans) = [];



%% UPDATE INDICES OF ANTENNAS AND SOURCES (in scan struct)
% only when something was deleted
if sum([sum(scansToExcludedSources_q), sum(emptyScans)]) > 0 
    
    % number of observations per station (one line=one station)
    nObsPerStat       = zeros(size(antenna, 2), 1);
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
            error('Unknown source type: %s', scan(iScan).obs_type);
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
