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

% ************************************************************************
function [scan, sources, antenna]=cleanScan(scan, sources, antenna, out_structFieldnames, allSourceNames, ini_opt, bas_excl, qualityLimit, minElevation)

%% Precalculations
vector1ToN=1:50;
nScans=size(scan,2);
nStats=length(antenna);
nSources=size(sources,2);
% define "empty" scan-substruct .stat
% substructs
subStruct_stat=struct('x', [], 'temp', [], 'pres', [], 'e', [], 'az', ...
    [], 'zd', [], 'zdry', [], 'cab', [], 'axkt', [], 'therm', [], ...
    'pantd', [], 'trop', []);

% get end-epoch of exclude stations
if ~isempty(ini_opt.sta_excl)
    if ~isfield(ini_opt, 'sta_excl_end')
        if isempty(ini_opt.sta_excl_start)
            ini_opt.sta_excl_end='';
        else
            % if there is a start epoch
            ini_opt.sta_excl_end=ones(size(ini_opt.sta_excl_start))*9999999;
        end
    else 
        % there is an end-epoch-field
        % make 0 -> 99999
        ini_opt.sta_excl_end(ini_opt.sta_excl_end==0)=9999999;
    end
    exclStaCellStr=cellstr(ini_opt.sta_excl);
end



%% Clean scan struct

% (1) No cable calibration
%fprintf('No cable-cal for stations: %1.0f\n', size(ini_opt.no_cab,1))

if isempty(ini_opt.no_cab)
else
    for k=1:size(ini_opt.no_cab,1)
        %fprintf('%s\n', ini_opt.no_cab(k,:));
    end
    
    % get antenna numbers
    antIndNoCab=find(ismember(strtrim({antenna.name}),cellstr(ini_opt.no_cab)));
    
    if ~isempty(antIndNoCab)
        % for all scans
        for iS=1:nScans
            for kAntNoCab=1:length(antIndNoCab)
                curI=antIndNoCab(kAntNoCab); % current index of antenna with no cable-cal
                % if that station acutally "exists" (if the index can be accessed
                if length(scan(iS).stat)>=curI
                    curCab=scan(iS).stat(curI).cab*(1e-9);

                    obsI1NoCab=find([scan(iS).obs.i1]==curI);
                    obsI2NoCab=find([scan(iS).obs.i2]==curI);
                    
                    for ionc = 1:length(obsI1NoCab)
                        scan(iS).obs(obsI1NoCab(ionc)).obs = scan(iS).obs(obsI1NoCab(ionc)).obs + curCab;
                    end
                    for ionc = 1:length(obsI2NoCab)
                        scan(iS).obs(obsI2NoCab(ionc)).obs = scan(iS).obs(obsI2NoCab(ionc)).obs - curCab;
                    end

                    % set cable cal to zero!
                    scan(iS).stat(curI).cab=0;
                end
            end
            
        end
    end
end
    

% (2) Excluded baselines
%fprintf('Baselines to be excluded: %1.0f\n', size(bas_excl,1))
if isempty(bas_excl)
else
    % write user info to command window
    %for k=1:size(bas_excl,1)
        %fprintf('%s\n', bas_excl(k,:))
    %end
    % for all baselines
    for iBasel=1:size(bas_excl,1)
        
        % get current baseline
        curBasel(1:2)={strtrim(bas_excl(iBasel,1:8)), ...
            strtrim(bas_excl(iBasel,10:17))};
        curStatLog=strcmpi(out_structFieldnames, curBasel(1))|...
            strcmpi(out_structFieldnames, curBasel(2));
        
        % get indices of both current stations
        curStatInd=vector1ToN(curStatLog);
        
        % go through all scans and remove stuff
        for iScan=1:nScans
            
            % only if both station-structs exist
            if max(curStatInd)<=size(scan(iScan).stat,2)
                
                % if both stations take part in the session
                if length([scan(iScan).stat(curStatInd).temp])==2
                    
                    % get observations to delete (might not be the case)
%                     obs2Delete=sum([scan(iScan).obs.i1; scan(iScan).obs.i2]==curStatInd(1) | ...
%                         [scan(iScan).obs.i1; scan(iScan).obs.i2]==curStatInd(2))==[1,1];
                    
                    A = [[scan(iScan).obs.i1]', [scan(iScan).obs.i2]']==curStatInd(1) | ... % Jakob new 17.01.2018
                        [[scan(iScan).obs.i1]', [scan(iScan).obs.i2]']==curStatInd(2);
                    ia1 = A(:,1)==1;
                    ia2 = A(:,2)==1;
                    obs2Delete = ia1 & ia2 ;
                    
                    if sum(obs2Delete)>0
                        
                        % then we have to delete something
                        scan(iScan).obs(obs2Delete)=[];
                        scan(iScan).nobs=scan(iScan).nobs-1;
                        
                        
                        % Jakob new 17.01.2018
                        statStructs2Delete_1=[];
                        statStructs2Delete_2=[];
                        % delete stat struct if needed (if that station has
                        % no observation anymore (for both stats separate)
                        if sum(sum([[scan(iScan).obs.i1]', [scan(iScan).obs.i2]']==curStatInd(1)))==0
                            statStructs2Delete_1=curStatInd(1);
                        end
                        if sum(sum([[scan(iScan).obs.i1]', [scan(iScan).obs.i2]']==curStatInd(2)))==0
                            statStructs2Delete_2=curStatInd(2);
                        end
                        
                        scan(iScan).stat(statStructs2Delete_1)=...
                            subStruct_stat;
                        
                        scan(iScan).stat(statStructs2Delete_2)=...
                            subStruct_stat;
                        
                    end
                end
            end
        end
    end
    
    
end

% fprintf('WARNING: FORTLEZA and wettzell EXCLUDED!!!\n');
% ini_opt.sta_excl=['FORTLEZA';'WETTZELL'];

% (3) Excluded stations
nStat2excl=size(ini_opt.sta_excl,1);
%fprintf('Stations to be excluded: %1.0f\n',nStat2excl)
obs2delPerStat=zeros(nStats,1); % for all stations: obs (nobs) have to be deleted!
if isempty(ini_opt.sta_excl)
else

    for k=1:nStat2excl
        %fprintf('%s ', ini_opt.sta_excl(k,:))
        %if ini_opt.sta_excl_start(k)~=0
            %fprintf('(%1.1f-%1.1f)', ini_opt.sta_excl_start(k), ini_opt.sta_excl_end(k))
        %end
        %fprintf('\n');
    end
    % get station logicals to be deleted
    curStatLog=zeros(nStats,1);
    for iStatExcl=1:size(ini_opt.sta_excl,1)
        foundLogsCurAnt=strcmpi({antenna.name},ini_opt.sta_excl(iStatExcl,:));
        if sum(foundLogsCurAnt)~=1
            fprintf('%s (to be excluded) not found in antenna struct\n-> not excluded!\n',ini_opt.sta_excl(iStatExcl,:));
        else
            curStatLog=curStatLog(:) | foundLogsCurAnt(:);
        end
    end
    
    curStatInd=vector1ToN(curStatLog);
    
    
    for iScan=1:nScans
        % get station indices which exist in current scan
        curStatIndInCurScan=curStatInd(curStatInd<=size(scan(iScan).stat,2));
        
        if isempty(curStatIndInCurScan)
        else
            % go through all antennas separately (due to time of
            % exclusion!)
            for iStatInCurScan=1:length(curStatIndInCurScan)
                curSingleStatIndInCurScan=curStatIndInCurScan(iStatInCurScan);
                curAnt=strtrim(antenna(curSingleStatIndInCurScan).name); % antenna name
                curSTART=ini_opt.sta_excl_start(strcmpi(exclStaCellStr,curAnt));
                curEND=ini_opt.sta_excl_end(strcmpi(exclStaCellStr,curAnt)); % mjd
                
                
                for iMultStat = 1:length(curSTART) % if same station is excluded multiple times
                    curStart    = curSTART(iMultStat);
                    curEnd      = curEND(iMultStat);
                    % check if scan is in "exclude-time"
                    if scan(iScan).mjd>curStart && scan(iScan).mjd<curEnd

                        % make scan-substruct .stat empty
                        scan(iScan).stat(curSingleStatIndInCurScan)=subStruct_stat;


                        % delete observations
                        obs2Delete=logical(sum(ismember(...
                            [scan(iScan).obs.i1;scan(iScan).obs.i2],curSingleStatIndInCurScan)));



                        % decrease number of observations of antennas
                        for iStat=1:nStats
                            % if we have at least the index of stat2del in current scan   
                            obs2delPerStat(iStat,1)=obs2delPerStat(iStat,1)+...
                                sum( ([[scan(iScan).obs(obs2Delete).i1],...
                                [scan(iScan).obs(obs2Delete).i2]])==iStat ); 
            %                 obs2delPerStat(curStatInd(iStatExcl))+...
            %                     sum(sum( ismember(,curStatIndInCurScan(iStatExcl)) ));
                        end

                        scan(iScan).obs(obs2Delete)=[];


                        % decrease number of observations
                        scan(iScan).nobs=scan(iScan).nobs-sum(obs2Delete);
                        % still missing: if stats in current scan have no obs left ->
                        % make also empty
                    end
                end
            end
        end
    end
end

for iStat=1:nStats
    antenna(iStat).numobs=antenna(iStat).numobs-...
        obs2delPerStat(iStat);
end


% (4) Excluded sources (only get indices - scans are deleted later)
%fprintf('Sources to be excluded: %1.0f\n', size(ini_opt.sour_excl,1))

% preallocate
exclSourcesInd=zeros(size(ini_opt.sour_excl,1),1);



exclSourcesInd_byTime = false(1,length([scan.iso]));
if ~isempty(ini_opt.sour_excl)
    % transpose source names for finding
    allSourceNamesTransposed=allSourceNames';
    
    
    % get source logicals to be deleted
    for iSource2BeExcl=1:size(ini_opt.sour_excl,1)
        foundDigitOfCurSource=strfind(allSourceNamesTransposed(:)', ini_opt.sour_excl(iSource2BeExcl,:));
        exclSourcesInd(iSource2BeExcl)=(foundDigitOfCurSource-1)/8+1;
        scansToExcludedSources=ismember([scan.iso], exclSourcesInd(iSource2BeExcl));
        if logical(ini_opt.sour_excl_start(iSource2BeExcl))
            exclSourcesInd_byTime = exclSourcesInd_byTime | (scansToExcludedSources & (([scan.mjd]>=ini_opt.sour_excl_start(iSource2BeExcl)) & ([scan.mjd]<=ini_opt.sour_excl_end(iSource2BeExcl))));
            exclSourcesInd(iSource2BeExcl)=[];
            %fprintf('%s %f %f\n', ini_opt.sour_excl(iSource2BeExcl,:), ini_opt.sour_excl_start(iSource2BeExcl),ini_opt.sour_excl_end(iSource2BeExcl));
        else
            %fprintf('%s\n', ini_opt.sour_excl(iSource2BeExcl,:))
        end
    end
	fprintf('\n')
end


% (5) Outliers
if isempty(ini_opt.scan_excl)
else
    for iOutlier=1:size(ini_opt.scan_excl,2)
        % get station indices of current outlier
%         [row_fieldnames,col_fieldnames] = size(out_structFieldnames);
%         if length(strtrim(ini_opt.scan_excl(iOutlier).sta1))<col_fieldname
            
            
%         curStatLog=strcmpi(out_structFieldnames,strtrim(ini_opt.scan_excl(iOutlier).sta1)) | ...
%             strcmpi(out_structFieldnames,strtrim(ini_opt.scan_excl(iOutlier).sta2));
        
        curStatLog = strcmpi(out_structFieldnames,strtrim({ini_opt.scan_excl(iOutlier).sta1})) | strcmpi(out_structFieldnames,strtrim({ini_opt.scan_excl(iOutlier).sta2})); % Indices of both stations in baseline
        
        curStatInd=vector1ToN(curStatLog);
        
        % get scan of cur outlier (the one which is close by 1/10 second!)
        oneSecInDays=1/60/60/24;
        curScanLog=abs([scan.mjd]-ini_opt.scan_excl(iOutlier).mjd)<(oneSecInDays/10);
        
        % Check, if only one scan was found!
        % - If more than one scan was found by matching the scan reference times, the stations have to be considered in addition
        flag_found_scan = true; 
        if sum(curScanLog) > 1 % More than one scan found?
            curScanLog_ids = find(curScanLog);
            flag_found_scan = false;
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
                fprintf('WARNING (outlier removal in cleanScan.m): No observation on the baseline %s - %s found at epoch %s!\n', ini_opt.scan_excl(iOutlier).sta1, ini_opt.scan_excl(iOutlier).sta2, mjd2datestr(ini_opt.scan_excl(iOutlier).mjd) )
            else
                cur_scan_id = scan_id;
            end
        else
            cur_scan_id = find(curScanLog);
        end
        
        if flag_found_scan
            % check if there is (still) an obs of the two stations in cur scan
            obs2Delete=sum([scan(cur_scan_id).obs.i1; scan(cur_scan_id).obs.i2]==curStatInd(1) | [scan(cur_scan_id).obs.i1; scan(cur_scan_id).obs.i2]==curStatInd(2))==2;
        else
            obs2Delete = [];
        end


        if sum(obs2Delete)>0
            % then we have to delete something
            scan(cur_scan_id).obs(obs2Delete)=[];
            scan(cur_scan_id).nobs=scan(cur_scan_id).nobs-1;

            statStructs2Delete=[];
            % delete stat struct if needed (if that station has
            % no observation anymore (for both stats separate)
            if sum(sum([scan(cur_scan_id).obs.i1; scan(cur_scan_id).obs.i2]==curStatInd(1)))==0
                statStructs2Delete(1)=curStatInd(1);
            end
            if sum(sum([scan(cur_scan_id).obs.i1; scan(cur_scan_id).obs.i2]==curStatInd(2)))==0
                statStructs2Delete(2)=curStatInd(2);
            end

            scan(cur_scan_id).stat(statStructs2Delete)=subStruct_stat;

        end
    end
end


% (6) Quality of observation
for iScan=1:nScans
    
    qualOfObsOfCurScan_deluflag = [scan(iScan).obs.q_code];    
    qualOfObsOfCurScan_IonCode = [scan(iScan).obs.q_code_ion];
    obs2Delete_deluflag = qualOfObsOfCurScan_deluflag > qualityLimit; 
    obs2Delete_IonCode =  qualOfObsOfCurScan_IonCode ~= 0; % fixed to zero
    
    obs2Delete = obs2Delete_deluflag | obs2Delete_IonCode; % 2018/01/16 Jakob, Ion code is also considered
    
    
    if sum(obs2Delete)>0
        scan(iScan).obs(obs2Delete)=[];
        scan(iScan).nobs=scan(iScan).nobs-sum(obs2Delete);
        
        % if there are no observations left - scan will be deleted later
        
        % delete stat (sub-)structs which are not needed anymore
        indStatStruct=1:size(scan(iScan).stat,2);

        % delete stat structs also if needed
        indStatStruct([[scan(iScan).obs.i1]';[scan(iScan).obs.i2]'])=[]; % those which (still) have index in i1|i2 should not be deleted

        scan(iScan).stat(indStatStruct)=subStruct_stat;
             
    end   
end

v_obs2 = [];
for i_scan = 1:length(scan)
    a = [scan(i_scan).obs.obs];
    v_obs2(length(v_obs2)+1:length(v_obs2)+length(a)) = a;
end

% (7) Minimum elevation angle
if minElevation>0
    % first get ra/de for all sources of this session
%     nSources=size(allSourceNames,1);
%     temp_sources=zeros(nSources,2);
%     for iSource=1:nSources
%         % find source in crf
%         foundSourceLog=strcmp({crf.commonName}, allSourceNames(iSource,:)) | ...
%             strcmp({crf.IERSname}, allSourceNames(iSource,:)) | ...
%             strcmp({crf.designation}, allSourceNames(iSource,:));
%         
%         % if one station is found -> write ra/de to matrix
%         if sum(foundSourceLog)==1
%             indSourceInCrf=find(foundSourceLog);
%             % always take vievsCrf (is not so critical)
%             if ~isempty(crf(indSourceInCrf).icrf2)
%                 temp_sources(iSource,1)=crf(indSourceInCrf).icrf2.ra;
%                 temp_sources(iSource,2)=crf(indSourceInCrf).icrf2.de;
%             elseif ~isempty(crf(indSourceInCrf).vievsCrf)
%                 temp_sources(iSource,1)=crf(indSourceInCrf).vievsCrf.ra;
%                 temp_sources(iSource,2)=crf(indSourceInCrf).vievsCrf.de;
%           
%             elseif ~isempty(crf(indSourceInCrf).icrf1Ext2)
%                 temp_sources(iSource,1)=crf(indSourceInCrf).icrf1Ext2.ra;
%                 temp_sources(iSource,2)=crf(indSourceInCrf).icrf1Ext2.de;
%             elseif ~isempty(crf(indSourceInCrf).VieCRF10a)
%                 temp_sources(iSource,1)=crf(indSourceInCrf).VieCRF10a.ra;
%                 temp_sources(iSource,2)=crf(indSourceInCrf).VieCRF10a.de;
%             end
%                 
%         else
%             fprintf('ERROR: Source %s not found in CRF.\n',...
%                 allSourceNames(iSource,:))
%         end
%  
%     end
    
    for iScan=1:nScans
        % get stations in current scan
        allStationsInCurScan=unique([[scan(iScan).obs.i1], [scan(iScan).obs.i2]]);
        nStatsInCurScan=size(allStationsInCurScan,2);
        
        % preallocate logicals if station(s) should be deleted in current scan
        deleteStatInCurScanLogicals=zeros(nStatsInCurScan,1);
        
        % go through all stations and calc elevations
        for iStat=1:length(allStationsInCurScan)
            curStat=allStationsInCurScan(iStat);
            
            if ~isempty(scan(iScan).stat(curStat).temp)
                % calc elevation
                curElevation=elev(scan(iScan).mjd, ...
                    [antenna(curStat).x, antenna(curStat).y, antenna(curStat).z], ...
                    sources(scan(iScan).iso).de2000, sources(scan(iScan).iso).ra2000);
                
                if curElevation<minElevation
                    % delete this observation
                    deleteStatInCurScanLogicals(iStat)=1;
                end
            end
        end % for all stations in cur scan
        
        % now delete bad observations
        if sum(deleteStatInCurScanLogicals)>0
            % get observations to be deleted
            obsWithBadStation=ismember([[scan(iScan).obs.i1];[scan(iScan).obs.i2]], ...
                allStationsInCurScan(logical(deleteStatInCurScanLogicals)));
            obs2Delete=sum(obsWithBadStation,1)>0;
            
            % delete observations which include one of stations with bad elevation 
            scan(iScan).obs(obs2Delete)=[];
            scan(iScan).nobs=scan(iScan).nobs-sum(obs2Delete);

            % if there are no observations left - scan will be deleted later anyway
            if scan(iScan).nobs>0
                indStatStruct=1:size(scan(iScan).stat,2);

                % delete stat structs also if needed
                indStatStruct([[scan(iScan).obs.i1]';[scan(iScan).obs.i2]'])=[]; % those which have index in i1|i2 should not be deleted

                scan(iScan).stat(indStatStruct)=subStruct_stat;

            end    
        end
    end % all scans
end


%% Delete empty scans and those to an excluded source
scansToExcludedSources=ismember([scan.iso], exclSourcesInd);
emptyScans=[scan.nobs]<=0; % < ... just to be sure!

scan((scansToExcludedSources | exclSourcesInd_byTime)|emptyScans)=[];

%% UPDATE INDICES OF ANTENNAS AND SOURCES (in scan struct)
% only when something was deleted
if sum([sum(scansToExcludedSources | exclSourcesInd_byTime), sum(emptyScans), nStat2excl])>0 
    
    % number of observations per station (one line=one station)
    nObsPerStat=zeros(size(antenna,2),1);
    nObsPerSource=zeros(size(sources,2),1);
    
    % get antennas with no observations, delete those from antenna and update indices
    % same with sources
    for iScan=1:size(scan,2)
        % add for each observation the entry by one
        [hista,histb]=hist([scan(iScan).obs.i1, scan(iScan).obs.i2],...
            1:nStats);
        nObsPerStat=nObsPerStat(:)+hista(:);
        % nObsPerStat should be equal [antenna.numobs]
        % same for source
        nObsPerSource(scan(iScan).iso)=nObsPerSource(scan(iScan).iso)+scan(iScan).nobs;
    end
    
    %% get new indices for antennas: e.g. antennas [1,2,3,4,5,6,7,8] -> 2,4,5,7 excluded --> [1,x,2,x,x,3,x,4]
    newIndForStat=1:nStats;
    % for all stations
    for iStat=1:nStats
        if nObsPerStat(iStat)==0
            newIndForStat(iStat)=0;
            % only if we are not at last station
            if iStat<nStats
                newIndForStat(iStat+1:end)=newIndForStat(iStat+1:end)-1;
            end
        end
    end

    %% get new indices for sources (same as for antenna)
    newIndForSource=1:nSources;
    for iSource=1:nSources
        if nObsPerSource(iSource)==0
            newIndForSource(iSource)=0;
            % only if we are not at last station
            if iSource<nSources
                newIndForSource(iSource+1:end)=newIndForSource(iSource+1:end)-1;
            end
        end
    end
        
    % if there is an antenna with no obs
    if sum(nObsPerStat==0)>0       
        statWithNoObs=nObsPerStat==0;
 
        %% go through all scans (and obs) and update to correct station index
        for iScan=1:size(scan,2)
            % for all obs entries (with deal I don't know how)
            for iObs=1:size(scan(iScan).obs,2)
                scan(iScan).obs(iObs).i1=newIndForStat(scan(iScan).obs(iObs).i1);
                scan(iScan).obs(iObs).i2=newIndForStat(scan(iScan).obs(iObs).i2);
            end
            
           
            
            % delete also scan(iScan).stat(...) for stations to be deleted
            scan(iScan).stat(statWithNoObs(1:size(scan(iScan).stat,2)))=[];
           
        end % for all scans
        
        % delete antenna entries (antennas with no obs)
        antenna(statWithNoObs)=[];
    end
    
     % source
            
            
    % if there is a source to delete
    if sum(nObsPerSource==0)>0
        sourWithNoObs=nObsPerSource==0;
 
%         matthias: this should be done in next step:...
%         % update source index in scan struct
%         % 1/2: get new indices
%         newSouInd=1:nSources; % this is (the wrong) source index for all 
%         for iSou=1:nSources
%             if sourWithNoObs(iSou)==1
%                 newSouInd(iSou:end)=newSouInd(iSou:end)-1;
%             end
%         end
%         % 2/2: update indices
%         for iScan=1:length(scan)
%             scan(iScan).iso=newSouInd(scan(iScan).iso);
%         end

        for iScan=1:size(scan,2)
            scan(iScan).iso=newIndForSource(scan(iScan).iso);
        end
        
        % delete sources to which no observations are taken
        sources(sourWithNoObs)=[];   
        
    end
    newNumObsPerSou=nObsPerSource(nObsPerSource~=0);
    % update number of observations per source
    for iSou=1:length(sources)
        sources(iSou).numobs=newNumObsPerSou(iSou);
    end
    newNumObsPerStat=nObsPerStat(nObsPerStat~=0);
    for iStat=1:length(antenna)   
        antenna(iStat).numobs=newNumObsPerStat(iStat);
    end
end

%% update 'firstObsMjd' 'lastObsMjd' of antenna structs and add both fields to sources struct
stat_firstObsMjd_checked=zeros(length(antenna),1);
src_firstObsMjd_checked=zeros(length(sources),1);

for iScan=1:length(scan)
    statIndOfCurScan=unique([[scan(iScan).obs.i1],[scan(iScan).obs.i2]]);
    % put mjd of cur scan to those stations
    [antenna(...
        statIndOfCurScan(stat_firstObsMjd_checked(statIndOfCurScan)==0)).firstObsMjd]=...
        deal(scan(iScan).mjd);
    stat_firstObsMjd_checked(statIndOfCurScan)=1;
    
    [antenna(statIndOfCurScan).lastObsMjd] = deal(scan(iScan).mjd);
    
    srcIndOfCurScan = scan(iScan).iso;
    
    if src_firstObsMjd_checked(srcIndOfCurScan)==0
        sources(srcIndOfCurScan).firstObsMjd = scan(iScan).mjd;
        src_firstObsMjd_checked(srcIndOfCurScan)=1;
    end
        sources(srcIndOfCurScan).lastObsMjd = scan(iScan).mjd;
    
end
    
% %% update 'firstObsMjd' of antenna structs
% firstObsMjd_checked=zeros(length(antenna),1);
% for iScan=1:length(scan)
%     % if all antennas are checked -> get out of for loop
%     if sum(firstObsMjd_checked)==length(firstObsMjd_checked)
%         break        
%     end
%     statIndOfCurScan=unique([[scan(iScan).obs.i1],[scan(iScan).obs.i2]]);
%     % put mjd of cur scan to those stations
%     [antenna(...
%         statIndOfCurScan(firstObsMjd_checked(statIndOfCurScan)==0)).firstObsMjd]=...
%         deal(scan(iScan).mjd);
%     firstObsMjd_checked(statIndOfCurScan)=1;
% end

