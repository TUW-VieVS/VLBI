% Matthias Schartner, Hana Krasna

function [scans] = getScanList(PATH_in,typ)

in =[    {'IVS'}    {'TUW'}    {'GSFC'}    {'GSF'}    {'USNO'}    {'BKG'}    {'GSI'}    {'OSO'}    {'noInst'}];
wrapper_k =  'all'
wrapper_v =  'highest_version'
observation = 'GroupDelayFull'
freqband = 'bX'

wrapper_data = read_vgosdb_wrapper([PATH_in '/'], 'ses', in, wrapper_k, wrapper_v);
[out_struct, nc_info]=read_nc([PATH_in '/']);


%% get raw data

nc_filename = get_nc_filename('ObsCrossRef', wrapper_data.Observation.CrossReference.files, 1);
obs2baseline=out_struct.CrossReference.(nc_filename).Obs2Baseline.val; 
obs2scan=out_struct.CrossReference.(nc_filename).Obs2Scan.val; 

nc_filename = get_nc_filename('SourceCrossRef', wrapper_data.Session.CrossReference.files, 1);
scan2source=out_struct.CrossReference.(nc_filename).Scan2Source.val; 
srcid2source_=out_struct.CrossReference.(nc_filename).CrossRefSourceList.val'; 


srcid2source = cell(0);
for i = 1:size(srcid2source_,1)
    srcid2source{end+1} = srcid2source_(i,:);
end

nc_filename = get_nc_filename('StationCrossRef', wrapper_data.Session.CrossReference.files, 1);
staid2station_=out_struct.CrossReference.(nc_filename).CrossRefStationList.val'; 


staid2station = cell(0);
for i = 1:size(staid2station_,1)
    staid2station{end+1} = staid2station_(i,:);
end

nc_filename = get_nc_filename('DelayTheoretical', wrapper_data.Observation.ObsTheoretical.files, 1);
delayTheoretical=out_struct.ObsTheoretical.(nc_filename).DelayTheoretical.val; 

if typ==1 % baseline delay
    tau_file = get_nc_filename({ observation , '_', freqband }, wrapper_data.Observation.ObsEdit.files, 1);
    groupDelayFull = out_struct.ObsEdit.(tau_file).GroupDelayFull.val;
else % typ ==2, geocentric delay
    tau_file = get_nc_filename({ 'CorrInfo' , '_', freqband }, wrapper_data.Observation.Observables.files, 1);
    tau_file(1,9)='_';
    groupDelayFull = out_struct.Observables.(tau_file).GeocMBD.val; % Tot geocenter group delay (sec). (does not exist for VLBA)
end 

sigma_tau_file = get_nc_filename({'GroupDelay', '_', freqband}, wrapper_data.Observation.Observables.files, 1);    
sigmaGroupDelayFull = out_struct.Observables.(sigma_tau_file).GroupDelaySig.val;

rate_file = get_nc_filename({'GroupRate', '_', freqband}, wrapper_data.Observation.Observables.files, 1);    
groupRate = out_struct.Observables.(rate_file).GroupRate.val;
sigmaGroupRate = out_struct.Observables.(rate_file).GroupRateSig.val; % sec/sec

nc_filename = get_nc_filename({'Edit'}, wrapper_data.Observation.ObsEdit.files, 0);
if ~isempty(nc_filename) % not match found in wrapper data
    delayFlag = out_struct.ObsEdit.(nc_filename).DelayFlag.val;
end


PATH_TimeUTC = fullfile(PATH_in, 'Observables', 'TimeUTC.nc');
ymdhm = ncread(PATH_TimeUTC,'YMDHM');
s = ncread(PATH_TimeUTC,'Second');
ep = [double(ymdhm') s];
mjd=date2mjd(ep); % errorbar is not possible with datetime
if ymdhm(1,1)<70
    ymdhm(1,:) = ymdhm(1,:) + 2000;
elseif ymdhm(1,1)> 69 && ymdhm(1,1)< 100
    ymdhm(1,:) = ymdhm(1,:) + 1900;
end

ref_time = datetime([ymdhm' s]);

%% create struct
scans = struct();
scans(length(scan2source)).obs = struct('sta1',{},'sta2',{},'src',{},'time',{},'mjd',{},'delay',{},'delay_theoretical',{},'rate',{},'sigmaDelay',{},'sigmaRate',{},'delayFlag',{});
for i = 1:length(obs2scan)
    t_obs2scan = obs2scan(i);
    stas = cell(0);
    obs = struct();
    obs.sta1 = staid2station(obs2baseline(1,i));
    obs.sta2 = staid2station(obs2baseline(2,i));
    obs.src  = srcid2source(scan2source(t_obs2scan));
    obs.time = ref_time(i);
    obs.mjd = mjd(i);
    obs.delay = groupDelayFull(i);
    obs.delay_theoretical = delayTheoretical(i);
    obs.sigmaDelay = sigmaGroupDelayFull(i);
    obs.rate = groupRate(i);    
    obs.sigmaRate = sigmaGroupRate(i);
    obs.delayFlag = delayFlag(i); 
    scans(t_obs2scan).obs(end+1) = obs;
end
fprintf('%9d observations found\n',length(obs2scan))




%%

% PATH_obsCrossRef = fullfile(PATH_in, 'CrossReference', 'ObsCrossRef.nc');
% obs2baseline = ncread(PATH_obsCrossRef,'Obs2Baseline');
% obs2scan = ncread(PATH_obsCrossRef,'Obs2Scan');

% PATH_SourceCrossRef = fullfile(PATH_in, 'CrossReference', 'SourceCrossRef.nc');
% scan2source = ncread(PATH_SourceCrossRef,'Scan2Source');
% srcid2source_ = ncread(PATH_SourceCrossRef,'CrossRefSourceList')';
% srcid2source = cell(0);
% for i = 1:size(srcid2source_,1)
%     srcid2source{end+1} = srcid2source_(i,:);
% end

% PATH_StationCrossRef = fullfile(PATH_in, 'CrossReference', 'StationCrossRef.nc');
% staid2station_ = ncread(PATH_StationCrossRef,'CrossRefStationList')';
% staid2station = cell(0);
% for i = 1:size(staid2station_,1)
%     staid2station{end+1} = staid2station_(i,:);
% end

% PATH_DelayTheoretical = fullfile(PATH_in, 'ObsTheoretical', 'DelayTheoretical.nc');
% delayTheoretical = ncread(PATH_DelayTheoretical,'DelayTheoretical');

% if typ==1 % baseline delay
%      PATH_GroupDelay = fullfile(PATH_in, 'ObsEdit', 'GroupDelayFull_bX.nc');
% %    PATH_GroupDelay = fullfile(PATH_in, 'ObsEdit', 'GroupDelayFull_iIVS_bX.nc');
%     groupDelayFull = ncread(PATH_GroupDelay,'GroupDelayFull');   
% else % typ ==2, geocentric delay
%     PATH_GroupDelayG = fullfile(PATH_in, 'Observables', 'CorrInfo-difx_bX.nc');
%     groupDelayFull = ncread(PATH_GroupDelayG,'GeocMBD'); % Tot geocenter group delay (sec).       
% end 

% PATH_GroupDelay1 = fullfile(PATH_in, 'Observables', 'GroupDelay_bX.nc');
% sigmaGroupDelayFull = ncread(PATH_GroupDelay1,'GroupDelaySig'); % sigma for the baseline and geocentric delay is the same

% PATH_GroupRate = fullfile(PATH_in, 'Observables', 'GroupRate_bX.nc');
% groupRate = ncread(PATH_GroupRate,'GroupRate');
% sigmaGroupRate = ncread(PATH_GroupRate,'GroupRateSig'); % sec/sec

% % ncdisp(PATH_DelayFlag)
% if exist(fullfile(PATH_in, 'ObsEdit', 'Edit_iIVS.nc'))
%     PATH_DelayFlag = fullfile(PATH_in, 'ObsEdit', 'Edit_iIVS.nc');
% elseif exist(fullfile(PATH_in, 'ObsEdit', 'Edit_iGSFC.nc'))
%     PATH_DelayFlag = fullfile(PATH_in, 'ObsEdit', 'Edit_iGSFC.nc');
% else
%     PATH_DelayFlag = fullfile(PATH_in, 'ObsEdit', 'Edit.nc');
% end
% delayFlag = ncread(PATH_DelayFlag,'DelayFlag'); %DelayFlag

