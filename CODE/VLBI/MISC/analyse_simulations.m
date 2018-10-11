%% summarize results from VieVS simulations
% written by: Matthias Schartner
clc, clear, close all;

% change sub-direcotry name here:
name = 'dummy';


path = pwd;
if strcmp(path(end-3:end),'WORK')
    path = '../DATA/LEVEL3/';
elseif strcmp(path(end-3:end),'MISC')
    path = '../../../VLBI/DATA/LEVEL3/';
end


if ~isfolder([path name])
    error('check folder name and path to LEVEL3 directory!')
end

%% read x_ variable
% get list of files from level 3 output
l3 = dir([path name '/x_*.mat']);
names = {l3.name};
sessionIds = unique(cellfun(@(x){x(1:end-9)}, names))';
nSessions = length(sessionIds);

%% preallocate memory
nSta = zeros(nSessions,1);
nSrc = zeros(nSessions,1);
nSim = zeros(nSessions,1);

stations = cell(nSessions,1);

dut1_mx = cell(nSessions,1);
dut1_val = cell(nSessions,1);

xpol_mx = cell(nSessions,1);
xpol_val = cell(nSessions,1);
 
ypol_mx = cell(nSessions,1);
ypol_val = cell(nSessions,1);

nutdx_mx = cell(nSessions,1);
nutdx_val = cell(nSessions,1);

nutdy_mx = cell(nSessions,1);
nutdy_val = cell(nSessions,1);

source_mx = cell(nSessions,1);
source_val = cell(nSessions,1);

station_mx = cell(nSessions,1);
station_val = cell(nSessions,1);

%% read x_ variables
for iSession=1:nSessions
    
    thisId = sessionIds{iSession};
    fprintf('Processing file %s (%d of %d)\n',thisId(3:end), iSession, nSessions);
    idx = cellfun(@(x) [any(x)], strfind(names, thisId));
    thisNSim = sum(idx);
    fprintf('    %4d simulations found\n', thisNSim);
    
    files = l3(idx);
    
    dut1_mx{iSession}  = zeros(thisNSim,1);
    dut1_val{iSession} = zeros(thisNSim,1);

    xpol_mx{iSession}  = zeros(thisNSim,1);
    xpol_val{iSession} = zeros(thisNSim,1);

    ypol_mx{iSession}  = zeros(thisNSim,1);
    ypol_val{iSession} = zeros(thisNSim,1);

    nutdx_mx{iSession}  = zeros(thisNSim,1);
    nutdx_val{iSession} = zeros(thisNSim,1);

    nutdy_mx{iSession}  = zeros(thisNSim,1);
    nutdy_val{iSession} = zeros(thisNSim,1);

    load([path name '/' files(1).name])
    nSta(iSession) = length(x_.coorx);
    nSrc(iSession) = length(x_.soura);
    nSim(iSession) = thisNSim;
    stations{iSession} = {x_.antenna.name};

    station_mx{iSession}  = zeros(thisNSim,nSta(iSession));
    station_val{iSession} = zeros(thisNSim,nSta(iSession));

    
    for iSim=1:thisNSim
        load([path name '/' files(iSim).name])
        
        if(~isempty(x_.dut1.col))
            dut1_mx{iSession}(iSim) = mean(x_.dut1.mx) * 1000;
            dut1_val{iSession}(iSim) = mean(x_.dut1.val) * 1000;
        end
        
        if(~isempty(x_.xpol.col))
            xpol_mx{iSession}(iSim) = mean(x_.xpol.mx) * 1000;
            xpol_val{iSession}(iSim) = mean(x_.xpol.val) * 1000;
        end
        
        if(~isempty(x_.ypol.col))
            ypol_mx{iSession}(iSim) = mean(x_.ypol.mx) * 1000;
            ypol_val{iSession}(iSim) = mean(x_.ypol.val) * 1000;
        end
        
        if(~isempty(x_.nutdx.col))
            nutdx_mx{iSession}(iSim) = mean(x_.nutdx.mx) * 1000;
            nutdx_val{iSession}(iSim) = mean(x_.nutdx.val) * 1000;
        end
        
        if(~isempty(x_.nutdy.col))
            nutdy_mx{iSession}(iSim) = mean(x_.nutdy.mx) * 1000;
            nutdy_val{iSession}(iSim) = mean(x_.nutdy.val) * 1000;
        end
        
        
        if(~isempty([x_.coorx.col]))
            coorx_mx = [x_.coorx.mx];
            coory_mx = [x_.coory.mx];
            coorz_mx = [x_.coorz.mx];
            
            coorx_val = [x_.coorx.val];
            coory_val = [x_.coory.val];
            coorz_val = [x_.coorz.val];
            
            % sum coordinates up to 3d vector
            station_mx{iSession}(iSim,:) = ...
                sqrt(coorx_mx.^2+coory_mx.^2+coorz_mx.^2) * 10;
            station_val{iSession}(iSim,:) = ...
                sqrt(coorx_val.^2+coory_val.^2+coorz_val.^2) * 10;
        end
        
    end
end
fprintf('\n\n');

%% calculate repeatabilities and mean sigmas
dut1_mean_sig = cellfun(@(x) mean(x), dut1_mx);
dut1_rep      = cellfun(@(x) std(x),  dut1_val);

xpol_mean_sig = cellfun(@(x) mean(x), dut1_mx);
xpol_rep      = cellfun(@(x) std(x),  dut1_val);

ypol_mean_sig = cellfun(@(x) mean(x), dut1_mx);
ypol_rep      = cellfun(@(x) std(x),  dut1_val);

nutdx_mean_sig = cellfun(@(x) mean(x), dut1_mx);
nutdx_rep      = cellfun(@(x) std(x),  dut1_val);

nutdy_mean_sig = cellfun(@(x) mean(x), dut1_mx);
nutdy_rep      = cellfun(@(x) std(x),  dut1_val);


% do the same for stations (sessions can have varying station networks)
allStations = unique([stations{:}]);
nStations = length(allStations);

station_mean_sig = zeros(nSessions, nStations);
station_rep      = zeros(nSessions, nStations);

for iSession=1:nSessions
    thisStas = stations{iSession};
    this_station_mean_sig = mean(station_mx{iSession},1);
    this_station_rep = std(station_val{iSession},0,1);
    
    for iSta=1:length(thisStas)
        thisSta = thisStas(iSta);
        idx = strcmp(thisSta,allStations);
        station_mean_sig(iSession, idx) = this_station_mean_sig(iSta);
        station_rep(iSession, idx) = this_station_rep(iSta);
    end
end


%% summarize results for mean sigma
t_mean_sig = table(nSim, dut1_mean_sig, xpol_mean_sig, ypol_mean_sig, nutdx_mean_sig, nutdy_mean_sig,'RowNames',sessionIds);
t_mean_sig.Properties.VariableNames = {'nr_sim','dUT1','x_pol','y_pol','x_nut','y_nut'};
t_mean_sig.Properties.VariableDescriptions = {'number of simulations', 'dUT1', 'x polar motion', 'y polar motion', 'x nutation', 'y nutation'};
t_mean_sig.Properties.VariableUnits = {'','mus','muas','muas','muas','muas'};
for iSta = 1:nStations
    tmp = table(station_mean_sig(:,iSta),'VariableNames',replace(strtrim(allStations(iSta)),'-','_'));
    tmp.Properties.VariableUnits = {'mm'};
    tmp.Properties.VariableDescriptions = {'3d coordinates'};
    t_mean_sig = [t_mean_sig tmp];
end


%% summarize results for repeatability
t_rep = table(nSim, dut1_rep, xpol_rep, ypol_rep, nutdx_rep, nutdy_rep,'RowNames',sessionIds);
t_rep.Properties.VariableNames = {'nr_sim','dUT1','x_pol','y_pol','x_nut','y_nut'};
t_rep.Properties.VariableDescriptions = {'number of simulations', 'dUT1', 'x polar motion', 'y polar motion', 'x nutation', 'y nutation'};
t_rep.Properties.VariableUnits = {'','mus','muas','muas','muas','muas'};
for iSta = 1:nStations
    tmp = table(station_rep(:,iSta),'VariableNames',replace(strtrim(allStations(iSta)),'-','_'));
    tmp.Properties.VariableUnits = {'mm'};
    tmp.Properties.VariableDescriptions = {'3d coordinates'};
    t_rep = [t_rep tmp];
end

%% displayresults
summary(t_mean_sig)
fprintf('\nmean formal error:\n')
disp(t_mean_sig)
fprintf('\nrepeatability:\n')
disp(t_rep)