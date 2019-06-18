%% summarize results from VieVS simulations
% written by: Matthias Schartner
function [t_mean_sig, t_rep] = analyse_simulations(name)


% change sub-direcotry name here:
if isempty(name)
    name = 'DUMMY';
end

path = pwd;
if strcmp(path(end-3:end),'WORK')
    path = '../DATA/LEVEL3/';
elseif strcmp(path(end-3:end),'MISC')
    path = '../../DATA/LEVEL3/';
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
    
    load([path name '/' files(1).name])
    nSta(iSession) = length(x_.coorx);
    nSrc(iSession) = length(x_.soura);
    nSim(iSession) = 0;
    stations{iSession} = {x_.antenna.name};

        
    if(~isempty(x_.dut1.col))
        dut1_mx{iSession} = mean(x_.dut1.mx) * 1000;
        dut1_val{iSession} = mean(x_.dut1.val) * 1000;
        nSim(iSession) = length(dut1_mx{iSession});
    end

    if(~isempty(x_.xpol.col))
        xpol_mx{iSession} = mean(x_.xpol.mx) * 1000;
        xpol_val{iSession} = mean(x_.xpol.val) * 1000;
        nSim(iSession) = length(xpol_mx{iSession});
    end

    if(~isempty(x_.ypol.col))
        ypol_mx{iSession} = mean(x_.ypol.mx) * 1000;
        ypol_val{iSession} = mean(x_.ypol.val) * 1000;
        nSim(iSession) = length(ypol_mx{iSession});
    end

    if(~isempty(x_.nutdx.col))
        nutdx_mx{iSession} = mean(x_.nutdx.mx) * 1000;
        nutdx_val{iSession} = mean(x_.nutdx.val) * 1000;
        nSim(iSession) = length(nutdx_mx{iSession});
    end

    if(~isempty(x_.nutdy.col))
        nutdy_mx{iSession} = mean(x_.nutdy.mx) * 1000;
        nutdy_val{iSession} = mean(x_.nutdy.val) * 1000;
        nSim(iSession) = length(nutdy_mx{iSession});
    end


    if(~isempty([x_.coorx.col]))
        nant = length(x_.coorx);
        nsim = length(x_.coorx(1).mx);
        
        coorx_mx = zeros(nant,nsim);
        coory_mx = zeros(nant,nsim);
        coorz_mx = zeros(nant,nsim);
        coorx_val = zeros(nant,nsim);
        coory_val = zeros(nant,nsim);
        coorz_val = zeros(nant,nsim);
        
        for s = 1:nant
        	coorx_mx(s,:) = x_.coorx(s).mx;
            coory_mx(s,:) = x_.coory(s).mx;
            coorz_mx(s,:) = x_.coorz(s).mx;

            coorx_val(s,:) = x_.coorx(s).val;
            coory_val(s,:) = x_.coory(s).val;
            coorz_val(s,:) = x_.coorz(s).val;
        end

        % sum coordinates up to 3d vector
        station_mx{iSession} = sqrt(coorx_mx.^2+coory_mx.^2+coorz_mx.^2) * 10;
        station_val{iSession} = sqrt(coorx_val.^2+coory_val.^2+coorz_val.^2) * 10;
        nSim(iSession) = length(station_mx{iSession});
    end
        
end
fprintf('\n\n');

%% calculate repeatabilities and mean sigmas
dut1_mean_sig = cellfun(@(x) mean(x), dut1_mx);
dut1_rep      = cellfun(@(x) std(x),  dut1_val);

xpol_mean_sig = cellfun(@(x) mean(x), xpol_mx);
xpol_rep      = cellfun(@(x) std(x),  xpol_val);

ypol_mean_sig = cellfun(@(x) mean(x), ypol_mx);
ypol_rep      = cellfun(@(x) std(x),  ypol_val);

nutdx_mean_sig = cellfun(@(x) mean(x), nutdx_mx);
nutdx_rep      = cellfun(@(x) std(x),  nutdx_val);

nutdy_mean_sig = cellfun(@(x) mean(x), nutdy_mx);
nutdy_rep      = cellfun(@(x) std(x),  nutdy_val);


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
t_mean_sig.Properties.VariableNames = {'n_sim','dUT1','x_pol','y_pol','x_nut','y_nut'};
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
t_rep.Properties.VariableNames = {'n_sim','dUT1','x_pol','y_pol','x_nut','y_nut'};
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

save(sprintf('../OUT/%s.mat',name),'t_mean_sig','t_rep', 'dut1_mx', 'dut1_val', 'xpol_mx', 'xpol_val', 'ypol_mx', 'ypol_val', 'nutdx_mx', 'nutdx_val', 'nutdy_mx', 'nutdy_val', 'stations', 'station_mx', 'station_val')

end