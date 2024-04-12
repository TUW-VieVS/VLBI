


function scan = changeAMB(scan, antenna, sources, parameter)



allStationNames = {antenna.name}';
allStationNames = strtrim(allStationNames);

if ~isempty(sources.q)
    allSourceNames_q  = {sources.q.name}';
    allSourceNames_q  = strtrim(allSourceNames_q);
else
    allSourceNames_q = {};
end
oneSecInDays = 1/60/60/24;

if parameter.amb.flag_change_amb
    if ~isempty(parameter.amb.obs2change)
        n_amb = length(parameter.amb.obs2change); % Number of amb 

        
        for iAmb = 1 : n_amb
            
            % Get indices of both stations in baseline
            curStatLog = strcmpi(allStationNames, strtrim({parameter.amb.obs2change(iAmb).sta1})) | strcmpi(allStationNames,strtrim({parameter.amb.obs2change(iAmb).sta2})); 
            curStatInd = find(curStatLog);

            % get scan of cur outlier (the one which is close by 1/10 second!)
            curScanLog = abs([scan.mjd] - parameter.amb.obs2change(iAmb).mjd) < (oneSecInDays/10);
            curSouLog = deblank(parameter.amb.obs2change(iAmb).sou);
            
            % Check, if only one scan was found!
            % - If more than one scan was found by matching the scan reference times, the stations have to be considered in addition
            flag_found_scan = true; 
            if sum(curScanLog) > 1 % More than one scan found?
                curScanLog_ids = find(curScanLog);
                flag_found_scan = false;

                    for i_scan = 1:size(curScanLog_ids,2)
                       scan_id = curScanLog_ids(i_scan);
                       if strcmp (allSourceNames_q(scan(scan_id).iso),curSouLog)
                            cur_scan_id = scan_id;
                            flag_found_scan = true; 
                       end
                    end
            else
                cur_scan_id = find(curScanLog);
            end

            iobs = find([scan(cur_scan_id).obs.i1] == curStatInd(1) & [scan(cur_scan_id).obs.i2] == curStatInd(2)) ;
            if isempty(iobs)
                iobs = find([scan(cur_scan_id).obs.i1] == curStatInd(2) & [scan(cur_scan_id).obs.i2] == curStatInd(1)) ;
            end
            if isempty(iobs)
                'obs not found'
                parameter.amb.obs2change(iAmb)
            end
            scan(cur_scan_id).obs(iobs).obs = scan(cur_scan_id).obs(iobs).obs + parameter.amb.obs2change(iAmb).amb*(1e-9); %s


        end
    end
end
