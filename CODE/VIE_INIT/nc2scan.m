% ************************************************************************
%   Description:
%   This function puts the data from the netCDF files into our standarized
%   format of scan structures.
%
%   Input:	
%      Both variables from the netCDF files - out_struct and nc_info.
%
%   Output:
%      The scan structure array.
% 
%   External calls: 	
%       
%   Coded for VieVS: 
%   Jul 2012 by Matthias Madzak
%
%   Revision: 
%   yyyy-mm-dd,FIRSTNAME SECONDNAME:
%   2016-06-16, M. Madzak: Errors when loading intensive sesssions: - "obs2Baseline" only provided for one for one baseline
%                                                                   - "delayFlagLikeNGS" is a skalar value
%                                                                   => Provided values are duplicated as required using repmat and a warning msg is printed to the CW
%   2016-07-06, A. Hellerschmied: Added checks for availability and validity of met. data. => Error code (error_code_invalid_met_data) in case of problems
%   2016-08-09, A. Hellerschmied: Field scan(i_scan).obs_type added.
%   2017-11-13, J.Gruber: General update. +It is now possible to choose a certain frequency band
%   2017-12-13, A. Hellerschmied: function modjuldat.m instead of date2mjd.m used.
%   2018-16-01, J.Gruber: Second redundant exception for Ion code Flag removed
%   2018-08-01, D. Landskron: bug corrected with storing the "Second" in vgosDB files
%   2018-12-05, D. Landskron: reading also the quality codes; clarification quality code / quality flag
%   2019-07-25, D. Landskron: zwet parameter added to scan structure

% ************************************************************************
function scan=nc2scan(out_struct, nc_info, fband, ioncorr, ambcorr, wrapper_data, parameter)
% fprintf('nc2scan started\n')

% ##### Options #####
error_code_invalid_met_data = -999; % Error corde for missing met. data in NGS file (numerical)

%% PREALLOCATING
nScans=out_struct.head.NumScan.val; % get number of scans, stored in Head.nc
nObs=out_struct.head.NumObs.val; % actually not required, however: nice for test
                           
% substructs for internal vievs struct
subStruct_stat=struct('x', [], 'temp', [], 'pres', [], 'e', [], 'az', ...
    [], 'zd', [], 'zdry', [], 'zwet', [], 'cab', [], 'axkt', [], 'therm', [], ...
    'pantd', [], 'trop', []);
subStruct_obs=struct('i1', [], 'i2', [], 'obs', [], 'sig', [], 'com', ...
    [], 'delion', [], 'sgdion', [], 'q_flag', [], 'q_flag_ion', [], 'q_code_X', [], 'q_code_S', []);
scan(nScans+1)=struct('mjd', [], 'stat', [], 'tim', [], ...
    'nobs', [], 'space', [], 'obs', [], 'iso', []); % +1: not working otherwise - is deleted after loop
space0.source = zeros(3,3);
space0.xp=0; space0.yp=0; space0.era=0; space0.xnut=0; space0.ynut=0;
space0.t2c=zeros(3,3);


% ##### Get cross referencing indices from folder CrossReference #####
% obs2Baseline:
nc_filename = get_nc_filename('ObsCrossRef', wrapper_data.Observation.CrossReference.files, 1);
obs2Baseline=double(out_struct.CrossReference.(nc_filename).Obs2Baseline.val)'; % it is also saved in netCDF file (->take it from there) nRows=nObs, nCols=2
% for (at least some) intensives: only one baseline given (as they are all equal) --> repmat
if size(obs2Baseline,1)==1
    obs2Baseline=repmat(obs2Baseline,nObs,1);
    fprintf('WARNING: obs2baseline only scalar! The provided values for one baseline are duplicated using repmat.\n');
end
obs2BaselineCell=num2cell(obs2Baseline);
% obs2Scan:
obs2Scan=double(out_struct.CrossReference.(nc_filename).Obs2Scan.val); % vector (lenght = nScans), giving scan number (as integer) nObs x 1
% scan2Station:
nc_filename = get_nc_filename('StationCrossRef', wrapper_data.Session.CrossReference.files, 1);
scan2Station=double(out_struct.CrossReference.(nc_filename).Scan2Station.val)'; % matrix (nRows=nScans), giving stations (also the number of observations per station is counted!), nCols = nStations
% scan2Source
nc_filename = get_nc_filename('SourceCrossRef', wrapper_data.Session.CrossReference.files, 1);
scan2Source=double(out_struct.CrossReference.(nc_filename).Scan2Source.val); % vector (nRows=nScans), giving integer source index, nCols=1

% split fband (e.g. fband = 'GroupDelayFull_bX' to observa='GroupDelayFull' and frqband='bX')
newStr = split(fband,'_');
observation = newStr{1};
freqband = newStr{2};

switch observation
    
    case 'GroupDelayFull'
        
        % Observation (tau) and sigma tau
        tau_folder = 'ObsEdit';
        tau_file = get_nc_filename({ observation , '_', freqband }, wrapper_data.Observation.ObsEdit.files, 1);
        tau_field = 'GroupDelayFull';
        sigma_tau_folder = 'Observables';
        sigma_tau_file = get_nc_filename({'GroupDelay', '_', freqband}, wrapper_data.Observation.Observables.files, 1);
        sigma_tau_field = 'GroupDelaySig';
        
        if strcmp(ambcorr,'on')
            amb_k = 0;
        elseif strcmp(ambcorr,'off')
            amb_k = -1;
        end

        
    case 'GroupDelay'

        % Observation (tau) and sigma tau
        tau_folder = 'Observables';
        tau_file = get_nc_filename({observation, '_', freqband}, wrapper_data.Observation.Observables.files, 1);
        tau_field = 'GroupDelay';
        sigma_tau_folder = 'Observables';
        sigma_tau_file = get_nc_filename({'GroupDelay', '_', freqband}, wrapper_data.Observation.Observables.files, 1);
        sigma_tau_field = 'GroupDelaySig';        

        
        if strcmp(ambcorr,'on')
            amb_k = 1;
        elseif strcmp(ambcorr,'off')
            amb_k = 0;
        end
end

fprintf('vgosdb data source:\n')
fprintf('\t observation:\t %s/%s, nc field: %s\n', tau_folder, tau_file,tau_field )
fprintf('\t sigma:\t\t %s/%s, nc field: %s\n', sigma_tau_folder, sigma_tau_file,sigma_tau_field )

% ionosphere delay correction

flag_ion_corr_available = false;
if isfield(wrapper_data.Observation, 'ObsDerived')
    nc_filename = get_nc_filename({'Cal-SlantPathIonoGroup', '_', freqband}, wrapper_data.Observation.ObsDerived.files);
    if ~isempty(nc_filename)
        flag_ion_corr_available = true;
    else
        fprintf('Cannot find ionosphere correction data: ObsDerived/%s', nc_filename)
    end
end
if flag_ion_corr_available
    tau_ion_folder = 'ObsDerived';
    tau_ion_file = strrep(nc_filename,'-','_');
    tau_ion_field = 'Cal_SlantPathIonoGroup';
    sigma_tau_ion_folder = 'ObsDerived';
    sigma_tau_ion_file = strrep(nc_filename,'-','_');
    sigma_tau_ion_field = 'Cal_SlantPathIonoGroupSigma';
    
    fprintf('\t ionosphere delay:\t %s/%s, nc field: %s\n', tau_ion_folder,tau_ion_file,tau_ion_field)
    fprintf('\t ionosphere delay sigma: %s/%s, nc field: %s\n', sigma_tau_ion_folder,sigma_tau_ion_file,sigma_tau_ion_field)
    
else % Ionosphere corrections not available!
    tau_ion_folder = {}; % ionospheric correction won't be used
    tau_ion_file = {}; % ionospheric correction won't be used
    tau_ion_field = {};
    sigma_tau_ion_folder = {};
    sigma_tau_ion_file = {};
    sigma_tau_ion_field = {};
    fprintf(' - No nc file with ionosphere corrections defined in the selected wrapper file!\n')
end

%  ionosphere correction
if strcmp(ioncorr,'on')
    fprintf('\t Iono corr: on\n')
else
    fprintf('\t Iono corr: off\n')
end
% ambiguity
if amb_k ~= 0
       
    ambN_folder = 'ObsEdit';
    
    % check if ObsEdit folder exists (mandatory for ambig values)
    if isfield(wrapper_data.Observation, 'ObsEdit')
        nc_filename = get_nc_filename({'NumGroupAmbig', '_bX'}, wrapper_data.Observation.ObsEdit.files, 1);
        ambN_file = nc_filename;
    else
         fprintf(' - No nc file with ambiguity correction defined in the selected wrapper file!\n')
    end
    
    ambN_field = 'NumGroupAmbig';
    
    % ambiguity size
    ambS_folder = 'Observables';
    ambS_file = ['AmbigSize_',freqband];
    ambS_field = 'AmbigSize';
    
    fprintf('\t ambiguity (integer number): \t %s/%s, nc field: %s\n',ambN_folder, ambN_file, ambN_field)
    fprintf('\t ambiguity (length): \t\t %s/%s, nc field: %s\n', ambS_folder, ambS_file, ambS_field)
    
end

%% DELAY:
groupDelayWAmbigCell = num2cell(out_struct.(tau_folder).(tau_file).(tau_field).val);

%% SIGMA DELAY:
groupDelaySigCell = num2cell(out_struct.(sigma_tau_folder).(sigma_tau_file).(sigma_tau_field).val);

%% IONOSPHERIC DELAY, SIGMA IONOSPHERIC DELAY and DELAY FLAG IONOSPHERIC DELAY::
ionoDelayInternalFlag = 1;
if strcmp(parameter.vie_init.iono, 'observation_database')
    if isempty(tau_ion_folder)
        ionoDelayInternalFlag = 0;
        fprintf('With this setting ionospheric delay will not be used\n')
    else
        if isfield(out_struct.(tau_ion_folder),tau_ion_file)            
            ionoDelCell = num2cell(1e9*out_struct.(tau_ion_folder).(tau_ion_file).(tau_ion_field).val(1,:)); % cell: 1 x nObs            
            ionoDelSigCell = num2cell(1e9*out_struct.(sigma_tau_ion_folder).(sigma_tau_ion_file).(sigma_tau_ion_field).val(1,:)); % cell: 1 x nObs
            if length(ionoDelCell) == 1
                if ionoDelCell{:}==0
                    ionoDelayInternalFlag = 0;
                    fprintf('Cal-SlantPathIonoGroup_bX.nc exists but no valid values, ionospheric delay will not be applied\n')
                end
            end
            if length(ionoDelSigCell) == 1
                if ionoDelSigCell{:}==0
                    ionoDelayInternalFlag = 0;
                    fprintf('Cal-SlantPathIonoGroup_bX.nc exists but no valid values for sigma, ionospheric delay will not be applied\n')
                end 
            end
            if isfield(out_struct.(tau_ion_folder).(tau_ion_file), 'Cal_SlantPathIonoGroupDataFlag') % if iono flag is given
                ionoDelFlagcell = num2cell(double(out_struct.(tau_ion_folder).(tau_ion_file).Cal_SlantPathIonoGroupDataFlag.val));
                if length(ionoDelFlagcell) == 1
                    fprintf(' - Same ionospheric delay flag (= %1.0f) used for all scans!\n', ionoDelFlagcell{1})
                end
            else
                ionoDelFlagcell = num2cell(zeros(1,length(groupDelayWAmbigCell)));
            end
            fprintf('Ionospheric delay will be used\n')
        else
            fprintf('Can find Ionospheric Delay File\n')
            warning('Ionospheric delay can not be used because was not found\n')
        end
    end
else % Take ionosphere corrections from external (ion) file:
    ionoDelayInternalFlag = 0;
    fprintf('Ionospheric delay corrections will be taken from external source.\n')
end

% in case of a zero ionoshperic delay (all ionospheric parameters will be
% set to zero, but with a correcto vector size to allow for vector addtion)
if ionoDelayInternalFlag == 0
    ionoDelCell = num2cell(zeros(1, length(groupDelayWAmbigCell)));
    ionoDelSigCell = num2cell(zeros(1, length(groupDelayWAmbigCell)));
    ionoDelFlagcell = num2cell(zeros(1, length(groupDelayWAmbigCell)));
end

%% AMBIGUITY CORRECTION
% amb_k = 0: no ambiguitues will be applied
% amb_k = 1: ambiguigies will be added
% amb_k = -1: ambiguities will be subtracted
if amb_k ~= 0
    % check for integer number ambiguities
    if isfield(out_struct.(ambN_folder),ambN_file)
        ambN= double(out_struct.(ambN_folder).(ambN_file).(ambN_field).val); % cell: nObs x 1
    else
        fprintf('Ambiguity data not available: %s is missing\n',[ambN_folder,'/',ambN_file])
    end

    % check for ambiguity size
    if isfield(out_struct.(ambS_folder),ambS_file)
        ambS = double(out_struct.(ambS_folder).(ambS_file).(ambS_field).val); % cell: nObs x 1 (sec)
    else
        fprintf('Ambiguity data not available: %s is missing\n',[ambS_folder,'/',ambS_file])        
    end
    
    % calculate ambiguity spacing and apply observation type factor
    tau_ambCell = num2cell(amb_k*ambN.*ambS);
else
    tau_ambCell = num2cell(zeros(1, length(groupDelayWAmbigCell)));
end

%% DELAY FLAG DELAY:
if isfield(wrapper_data.Observation,'ObsEdit')
    nc_filename = get_nc_filename({'Edit'}, wrapper_data.Observation.ObsEdit.files, 0);
    if ~isempty(nc_filename) % not match found in wrapper data
        delayQualityFlag = num2cell(out_struct.ObsEdit.(nc_filename).DelayFlag.val);
    else
        fprintf(' - No delay flags defined in wrapper file: delay flag is set to "0" for all observations!\n')
        delayQualityFlag = {0};
    end
else
    fprintf('No delay flags from ObsEdit/Edit*.nc file loaded: delay flag is set to "0" for all observations\n')
    delayQualityFlag = {0};
end

%% QUALITY CODES FOR X-BAND and S-BAND: 
% only used for sessions prior to 2001 in cleanScan.m 
nc_filename = get_nc_filename({'QualityCode_bX'}, wrapper_data.Observation.Observables.files, 0);
if ~isempty(nc_filename) % not mathc found in wrapper data
    qualityCode_X = num2cell(out_struct.Observables.(nc_filename).QualityCode.val);
else
    fprintf(' - No quality codes for X-band defined in wrapper file: Quality code is set to "0" for all X-band observations!\n')
    qualityCode_X = {0};
end

nc_filename = get_nc_filename({'QualityCode_bS'}, wrapper_data.Observation.Observables.files, 0);
if ~isempty(nc_filename) % not mathc found in wrapper data
    qualityCode_S = num2cell(out_struct.Observables.(nc_filename).QualityCode.val);
else
    fprintf(' - No quality codes for S-band defined in wrapper file: Quality code is set to "0" for all S-band observations!\n')
    qualityCode_S = {0};
end

%% SIGMA FINAL DELAY:  
delaySigmaTimesIonoSigma=num2cell(sqrt([groupDelaySigCell{:}].^2 + ([ionoDelSigCell{:}]*1e-9).^2)); % [sec]    cell: 1 x nObs

%% FILL SCAN STRUCT
% get number of antennas per scan
nAntPerScan=sum(scan2Station>0,2); % vector (nRows=nScans), giving number of participatin (in Scan) stations
% simply create vector from one to 40
oneToN=1:40;
% scan.mjd /.iso
% yr; mo; day; hr; minute; sec. (num cols = num scans)

nc_filename = get_nc_filename({'TimeUTC'}, wrapper_data.Scan.Scan.files, 1);
if length(out_struct.Scan.TimeUTC.Second.val) == 1
    out_struct.Scan.(nc_filename).Second.val = repmat(out_struct.Scan.(nc_filename).Second.val,length(out_struct.Scan.TimeUTC.YMDHM.val),1);
end
tim=[double(out_struct.Scan.(nc_filename).YMDHM.val); out_struct.Scan.(nc_filename).Second.val'];

scanMjd =  modjuldat(double(out_struct.Scan.(nc_filename).YMDHM.val(1,:)'), double(out_struct.Scan.(nc_filename).YMDHM.val(2,:)'), double(out_struct.Scan.(nc_filename).YMDHM.val(3,:)')) + double(out_struct.Scan.(nc_filename).YMDHM.val(4,:))'./24 + double(out_struct.Scan.(nc_filename).YMDHM.val(5,:))'./60./24 + out_struct.Scan.(nc_filename).Second.val/60/60/24;

scanMjdCell=num2cell(scanMjd);
scanSouCell=num2cell(scan2Source);
[scan(1:end-1).mjd]=deal(scanMjdCell{:});
[scan(1:end-1).iso]=deal(scanSouCell{:});
[itim,doy]=dday(tim(1,:),tim(2,:),tim(3,:),tim(4,:),tim(5,:));

% substructs need a for-loop
obsI1Index=1;
for iScan=1:nScans
    % scan.nobs
    scan(iScan).nobs=sum(obs2Scan==iScan);
    scan(iScan).tim=[tim(:,iScan);doy(iScan)];
    % +++ scan.stat +++
    
    % "preallocate"
    scan(iScan).stat(nAntPerScan(iScan))=subStruct_stat;

    % for each station in current scan get the observation number
    stationsInCurScan=logical(scan2Station(iScan,:));
    stationIndices=oneToN(stationsInCurScan);
    
    % for all stations in current scan
    for iStat=1:length(stationIndices)
        
        % Get nc file-list for the current station from the wrapper:
        stat_id = stationIndices(iStat);
        stat_name = deblank(out_struct.head.StationList.val(:,stat_id)');
        
        % Exchange '-' with '_', because '-' are not valid for field names of MATLAB structs (...and therefor not used in "wrapper_data"):
        stat_name(strfind(stat_name, '-')) = '_';
        
        wrapper_stat_file_list = wrapper_data.Station.(stat_name).files;
        
        tmp=scan2Station(scan2Station(:,stat_id) ~= 0, stat_id);
        num_of_scans_per_stat = tmp(end);
        
        
        %% #### Met. data ####
        
        % ### scan.stat.temp ###
        % Check data availability:
        nc_filename = get_nc_filename({'Met'}, wrapper_stat_file_list);
        
        if ~isempty(nc_filename)% Check, if met. data is available for this station and scan
            if isfield(out_struct.stat(stationIndices(iStat)), nc_filename)
                if ~isempty(out_struct.stat(stationIndices(iStat)).(nc_filename))
                    if isfield(out_struct.stat(stationIndices(iStat)).(nc_filename), 'TempC')
                        if size(out_struct.stat(stationIndices(iStat)).(nc_filename).TempC.val, 1) == num_of_scans_per_stat
                           tdry = out_struct.stat(stationIndices(iStat)).(nc_filename).TempC.val(scan2Station(iScan,stationIndices(iStat)));
                        else
                            tdry = error_code_invalid_met_data;
                        end
                    else
                        tdry = error_code_invalid_met_data;
                    end
                else
                    tdry = error_code_invalid_met_data;
                end
                % Check data value:
                if tdry ~= error_code_invalid_met_data
                    if (tdry < -99)
                        tdry = error_code_invalid_met_data;
                    end
                end
                scan(iScan).stat(stationIndices(iStat)).temp = tdry;

                % ### scan.stat.pres ###
                % Check data availability:
                if ~isempty(out_struct.stat(stationIndices(iStat)).(nc_filename))
                    if isfield(out_struct.stat(stationIndices(iStat)).(nc_filename), 'AtmPres')
                        if size(out_struct.stat(stationIndices(iStat)).(nc_filename).AtmPres.val, 1) == num_of_scans_per_stat
                            pres = out_struct.stat(stationIndices(iStat)).(nc_filename).AtmPres.val(scan2Station(iScan,stationIndices(iStat)));
                        else
                            pres = error_code_invalid_met_data;
                        end
                    else
                        pres = error_code_invalid_met_data;
                    end
                else
                    pres = error_code_invalid_met_data;
                end
                % Check data value:
                if pres ~= error_code_invalid_met_data
                    if (pres < 0)
                        pres = error_code_invalid_met_data;
                    end
                end
                scan(iScan).stat(stationIndices(iStat)).pres = pres;

                % ### scan.stat.e ###
                % Check data availability:
                if ~isempty(out_struct.stat(stationIndices(iStat)).(nc_filename))
                    if isfield(out_struct.stat(stationIndices(iStat)).(nc_filename), 'RelHum')
                        if size(out_struct.stat(stationIndices(iStat)).(nc_filename).RelHum.val, 1) == num_of_scans_per_stat
                            relHum = out_struct.stat(stationIndices(iStat)).(nc_filename).RelHum.val(scan2Station(iScan,stationIndices(iStat)));
                        else
                            relHum = error_code_invalid_met_data;
                        end
                    else
                        relHum = error_code_invalid_met_data;
                    end
                else
                    relHum = error_code_invalid_met_data;
                end
                % Check data value:
                if (relHum ~= error_code_invalid_met_data) && (tdry ~= error_code_invalid_met_data) 
                    if (tdry > -99) && (relHum > 0)
                        e = 6.1078 * exp((17.1 * tdry) / (235 + tdry)) * relHum;   % formula by Magnus * relative humidity
                    else
                        e = error_code_invalid_met_data;
                    end
                else
                    e = error_code_invalid_met_data;
                end
                scan(iScan).stat(stationIndices(iStat)).e = e;
            else
                scan(iScan).stat(stationIndices(iStat)).temp = error_code_invalid_met_data;
                scan(iScan).stat(stationIndices(iStat)).pres = error_code_invalid_met_data;
                scan(iScan).stat(stationIndices(iStat)).e = error_code_invalid_met_data;
            end
        else % No met. data available
            scan(iScan).stat(stationIndices(iStat)).temp = error_code_invalid_met_data;
            scan(iScan).stat(stationIndices(iStat)).pres = error_code_invalid_met_data;
            scan(iScan).stat(stationIndices(iStat)).e = error_code_invalid_met_data;
        end

            
        
        % #### Cable Cal. ####
        nc_filename = get_nc_filename({'Cal-Cable'}, wrapper_stat_file_list);
        nc_filename = strrep(nc_filename,'-','_');
        % scan.stat.cab
        if isfield(out_struct.stat(stationIndices(iStat)), nc_filename)
            if length(out_struct.stat(stationIndices(iStat)).(nc_filename).Cal_Cable.val)>1
                scan(iScan).stat(stationIndices(iStat)).cab = 1e9*out_struct.stat(stationIndices(iStat)).(nc_filename).Cal_Cable.val(scan2Station(iScan,stationIndices(iStat))); % [nano-sec]
            else
                scan(iScan).stat(stationIndices(iStat)).cab = 0; % [nano-sec]
            end
        else
            scan(iScan).stat(stationIndices(iStat)).cab = 0; % [nano-sec]            
        end
    end
    %% --- scan.stat ---
    
    % +++ scan.obs +++

    
    % "preallocate"
    scan(iScan).obs(scan(iScan).nobs)=subStruct_obs;
    [scan(iScan).obs.i1] = deal(obs2BaselineCell{obsI1Index:obsI1Index+scan(iScan).nobs-1,1});
    [scan(iScan).obs.i2] = deal(obs2BaselineCell{obsI1Index:obsI1Index+scan(iScan).nobs-1,2});

    [scan(iScan).obs.obs]=   deal(groupDelayWAmbigCell{obsI1Index:obsI1Index+scan(iScan).nobs-1}); % [sec]
    [scan(iScan).obs.sig]=deal(delaySigmaTimesIonoSigma{obsI1Index:obsI1Index+scan(iScan).nobs-1}); % [sec]
    [scan(iScan).obs.delion]=   deal(ionoDelCell{obsI1Index:obsI1Index+scan(iScan).nobs-1}); % [nano-sec]
    [scan(iScan).obs.sgdion]=   deal(ionoDelSigCell{obsI1Index:obsI1Index+scan(iScan).nobs-1}); % [nano-sec]
    
    [scan(iScan).obs.amb]=   deal(tau_ambCell{obsI1Index:obsI1Index+scan(iScan).nobs-1}); % [sec]
    
    if length(delayQualityFlag)==1 % check length of delay flag vector, if it is only 1 value for the whole session, this value will be assigned to all observations
        [scan(iScan).obs.q_flag] = deal(double(delayQualityFlag{1}));
    else
        [scan(iScan).obs.q_flag] = deal(delayQualityFlag{obsI1Index:obsI1Index+scan(iScan).nobs-1});   
    end
    
    if length(ionoDelFlagcell)==1 % check length of delay flag vector, if it is only 1 value for the whole session, this value will be assigned to all observations
        [scan(iScan).obs.q_flag_ion] = deal(ionoDelFlagcell{1}); 
    else
        [scan(iScan).obs.q_flag_ion] = deal(ionoDelFlagcell{obsI1Index:obsI1Index+scan(iScan).nobs-1});  
    end
    
    if length(qualityCode_X)==1 % check length of quality code X vector, if it is only 1 value for the whole session, this value will be assigned to all observations
        [scan(iScan).obs.q_code_X] = deal(double(qualityCode_X{1}));          
    else
        [scan(iScan).obs.q_code_X] = deal(qualityCode_X{obsI1Index:obsI1Index+scan(iScan).nobs-1});   
    end
    if length(qualityCode_S)==1 % check length of quality code S vector, if it is only 1 value for the whole session, this value will be assigned to all observations
        [scan(iScan).obs.q_code_S] = deal(double(qualityCode_S{1}));          
    else
        [scan(iScan).obs.q_code_S] = deal(qualityCode_S{obsI1Index:obsI1Index+scan(iScan).nobs-1});   
    end
    

    obsI1Index=obsI1Index+scan(iScan).nobs;

    % add cable delay
    cableCalibration        = 1;
    if cableCalibration == 1
        for iObs = 1 : length(scan(iScan).obs)
            corcab = scan(iScan).stat(scan(iScan).obs(iObs).i2).cab - scan(iScan).stat(scan(iScan).obs(iObs).i1).cab; % [ns]
            scan(iScan).obs(iObs).obs = scan(iScan).obs(iObs).obs + corcab*(1e-9);
        end
    end
    
    %  ionosphere correction
    if strcmp(ioncorr,'on')
        for iObs = 1 : length(scan(iScan).obs)
            scan(iScan).obs(iObs).obs = scan(iScan).obs(iObs).obs - scan(iScan).obs(iObs).delion*(1e-9); % [sec]
        end
    end
    
    % ambiguity correction
    if amb_k ~= 0
        for iObs = 1 : length(scan(iScan).obs)
            scan(iScan).obs(iObs).obs = scan(iScan).obs(iObs).obs + scan(iScan).obs(iObs).amb;
        end
    end
    
    % --- scan.obs ---
    
    % +++ scan.space +++
    scan(iScan).space=space0;
    % --- scan.space ---
    
    % +++ scan.obs_type +++
    scan(iScan).obs_type = 'q';
end

% delete last scan entry (which was never needed)
scan(end)=[];
% fprintf('nc2scan finished\n')





