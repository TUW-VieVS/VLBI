function [scans] = getScanList(PATH_in)
%% get raw data
PATH_obsCrossRef = fullfile(PATH_in, 'CrossReference', 'ObsCrossRef.nc');
obs2baseline = ncread(PATH_obsCrossRef,'Obs2Baseline');
obs2scan = ncread(PATH_obsCrossRef,'Obs2Scan');

PATH_SourceCrossRef = fullfile(PATH_in, 'CrossReference', 'SourceCrossRef.nc');
scan2source = ncread(PATH_SourceCrossRef,'Scan2Source');
srcid2source_ = ncread(PATH_SourceCrossRef,'CrossRefSourceList')';
srcid2source = cell(0);
for i = 1:size(srcid2source_,1)
    srcid2source{end+1} = srcid2source_(i,:);
end

PATH_StationCrossRef = fullfile(PATH_in, 'CrossReference', 'StationCrossRef.nc');
staid2station_ = ncread(PATH_StationCrossRef,'CrossRefStationList')';
staid2station = cell(0);
for i = 1:size(staid2station_,1)
    staid2station{end+1} = staid2station_(i,:);
end

PATH_DelayTheoretical = fullfile(PATH_in, 'ObsTheoretical', 'DelayTheoretical.nc');
delayTheoretical = ncread(PATH_DelayTheoretical,'DelayTheoretical');

PATH_GroupDelay = fullfile(PATH_in, 'ObsEdit', 'GroupDelayFull_iIVS_bX.nc');
%PATH_GroupDelay = fullfile(PATH_in, 'ObsEdit', 'GroupDelayFull_bX.nc');
groupDelayFull = ncread(PATH_GroupDelay,'GroupDelayFull');

PATH_GroupRate = fullfile(PATH_in, 'Observables', 'GroupRate_bX.nc');
groupRate = ncread(PATH_GroupRate,'GroupRate');

PATH_TimeUTC = fullfile(PATH_in, 'Observables', 'TimeUTC.nc');
ymdhm = ncread(PATH_TimeUTC,'YMDHM');
s = ncread(PATH_TimeUTC,'Second');
ref_time = datetime([ymdhm' s]);

%% create struct
scans = struct();
scans(length(scan2source)).obs = struct('sta1',{},'sta2',{},'src',{},'time',{},'delay',{},'delay_theoretical',{},'rate',{});
for i = 1:length(obs2scan)
    t_obs2scan = obs2scan(i);
    stas = cell(0);
    obs = struct();
    obs.sta1 = staid2station(obs2baseline(1,i));
    obs.sta2 = staid2station(obs2baseline(2,i));
    obs.src  = srcid2source(scan2source(t_obs2scan));
    obs.time = ref_time(i);
    obs.delay = groupDelayFull(i);
    obs.delay_theoretical = delayTheoretical(i);
    obs.rate = groupRate(i);
    scans(t_obs2scan).obs(end+1) = obs;
end
fprintf('%9d observations found\n',length(obs2scan))



