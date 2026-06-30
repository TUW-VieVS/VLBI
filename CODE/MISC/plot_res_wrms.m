% called by L3plot.m

function [fig1, fig2, fig3, fig4, fig5, fig6] = plot_res_wrms(res, opt_,tableout)

x_all  = double(res.mainVal(:));
mx_all = double(res.sigma_residuals_aposteriori(:));

sources_all   = double(res.source(:));
baselines_all = double(res.baselineOfObs);

sourceNames = cellstr(strtrim(res.allSourceNames));
statNames   = cellstr(strtrim(res.allStatNames));

%%
colors1 =[55,  126, 184 %  #377eb8 'blue':  
     255, 127, 0 %    #ff7f00 'orange':
     77,  175, 74 %,   #4daf4a 'green': 
      247, 129, 191 %,  #f781bf 'pink': 
      166, 86,  40 %,   #a65628 'brown':
     152, 78,  163 %,  #984ea3 'purple':
     153, 153, 153 %,  #999999 'gray':  
      228, 26,  28 %,   #e41a1c 'red':  
    222, 222, 0] ./255; %    #dede00 'yellow': 


colors2 = ([55,126,184;
          255,127,0;
          77,175,74;
          247,129,191;
          166,86,40;
          152,78,163;
          153,153,153;
          228,26,28;
          222,222,0] +...
          [50 50 50
          0 50 50
          50 50 50
          0 50 50
          50 50 50
          50 50 50
          50 50 50
          10 50 50
          10 10 50])./ 255;

clrs=[colors1;colors2]; % 18 ant max

%% Diagnostics of all observations

nObsTotal = numel(x_all);

isProblem = ~isfinite(x_all) | ~isfinite(mx_all) | mx_all <= 0;

Tdiagnostics = table( ...
    nObsTotal, ...
    sum(isnan(x_all)), ...
    sum(isinf(x_all)), ...
    sum(isnan(mx_all)), ...
    sum(isinf(mx_all)), ...
    sum(mx_all == 0), ...
    sum(mx_all < 0), ...
    sum(isProblem), ...
    'VariableNames', { ...
    'Nobs_total', ...
    'NaN_residuals', ...
    'Inf_residuals', ...
    'NaN_sigma', ...
    'Inf_sigma', ...
    'Zero_sigma', ...
    'Negative_sigma', ...
    'Problematic_observations'});

disp('Observation diagnostics')
disp(Tdiagnostics)

if tableout
    writetable(Tdiagnostics,'observation_diagnostics.csv');
end

%% WRMS-valid observations only for mathematical WRMS computation

validWRMS = isfinite(x_all) & isfinite(mx_all) & mx_all > 0;

x         = x_all(validWRMS);
mx        = mx_all(validWRMS);
sources   = sources_all(validWRMS);
baselines = baselines_all(validWRMS,:);

%% Global statistics

globalWRMS = calcWRMS(x,mx);
z = x ./ mx;

Tglobal = table( ...
    numel(x), ...
    nObsTotal, ...
    sum(isProblem), ...
    round(mean(x),2), ...
    round(median(x),2), ...
    round(std(x),2), ...
    round(min(x),2), ...
    round(max(x),2), ...
    round(sqrt(mean(x.^2)),2), ...
    round(globalWRMS,2), ...
    round(mean(z),2), ...
    round(median(z),2), ...
    round(std(z),2), ...
    round(sqrt(mean(z.^2)),2), ...
    round(100 * sum(abs(z) > 3) / numel(z),2), ...
    round(100 * sum(abs(x) > 1.0) / numel(x),2), ...
    'VariableNames', { ...
    'Nobs_used_for_WRMS', ...
    'Nobs_total', ...
    'Problematic_observations', ...
    'MeanResidual_cm', ...
    'MedianResidual_cm', ...
    'StdResidual_cm', ...
    'MinResidual_cm', ...
    'MaxResidual_cm', ...
    'RMS_cm', ...
    'WRMS_cm', ...
    'MeanNormResidual', ...
    'MedianNormResidual', ...
    'StdNormResidual', ...
    'RMSNormResidual', ...
    'PercentAbsNormResidualLarger3', ...
    'PercentAbsResidualLarger1cm'});

disp('Global residual statistics')
disp(Tglobal)

if tableout
    writetable(Tglobal,'residual_global_statistics.csv');
end
fprintf('\nGlobal WRMS: %.2f cm\n',globalWRMS);

%% WRMS per source

nSources = numel(sourceNames);

srcName = {};
srcWRMS = [];
srcNobs = [];
srcNproblem = [];

srcEstNNR= [];
srcEstPWLO= [];
srcNNR = [];
srcPWLO = [];

for i = 1:nSources

    idxAll = sources_all == i;
    idx    = sources == i;

    if any(idxAll)

        srcName{end+1,1} = sourceNames{i};
        srcNobs(end+1,1) = sum(idxAll);
        srcNproblem(end+1,1) = sum(idxAll & isProblem);
        
        srcEstNNR(end+1,1) = opt_.est_sourceNNR;
        srcEstPWLO(end+1,1) = opt_.pw_sou;
        srcNNR(end+1,1) = opt_.source(i).nnr_inc;
        srcPWLO(end+1,1) = opt_.source(i).rade_inc;

        if any(idx)
            srcWRMS(end+1,1) = calcWRMS(x(idx),mx(idx));
        else
            srcWRMS(end+1,1) = NaN;
        end

    end

end

Tsrc = table(srcName,srcNobs,srcNproblem,srcWRMS, srcEstNNR, srcNNR, srcEstPWLO, srcPWLO,...
    'VariableNames',{'Source','Nobs','ProblematicObs','WRMS_cm','estNNR','incNNR','estPwlo','incPwlo'});

Tsrc.WRMS_cm = round(Tsrc.WRMS_cm,2);
Tsrc = sortrows(Tsrc,'WRMS_cm','descend','MissingPlacement','last');






souFIXED = false;
if Tsrc.estNNR(1)==0 & Tsrc.estPwlo(1)==0
    souFIXED = true;
end
souNNR = false;
if Tsrc.estNNR(1)==1
    souNNR = true;
    idso0 = find(Tsrc.incNNR ==0);
end
souPWLO = false;
if Tsrc.estPwlo(1)==1
    souPWLO = true;
    idso0 = find(Tsrc.incPwlo ==1);
end


%% WRMS per station

nSta = numel(statNames);

staName = {};
staWRMS = [];
staNobs = [];
staNproblem = [];

for i = 1:nSta

    idxAll = baselines_all(:,1) == i | baselines_all(:,2) == i;
    idx    = baselines(:,1) == i | baselines(:,2) == i;

    if any(idxAll)

        staName{end+1,1} = statNames{i};
        staNobs(end+1,1) = sum(idxAll);
        staNproblem(end+1,1) = sum(idxAll & isProblem);

        if any(idx)
            staWRMS(end+1,1) = calcWRMS(x(idx),mx(idx));
        else
            staWRMS(end+1,1) = NaN;
        end

    end

end

Tsta = table(staName,staNobs,staNproblem,staWRMS, ...
    'VariableNames',{'Station','Nobs','ProblematicObs','WRMS_cm'});

Tsta.WRMS_cm = round(Tsta.WRMS_cm,2);
Tsta = sortrows(Tsta,'WRMS_cm','descend','MissingPlacement','last');

%% WRMS per baseline

baselineID_all = strings(size(baselines_all,1),1);

for i = 1:size(baselines_all,1)

    s1 = baselines_all(i,1);
    s2 = baselines_all(i,2);

    if s1 > 0 && s2 > 0 && ...
       s1 <= numel(statNames) && s2 <= numel(statNames)

        name1 = statNames{s1};
        name2 = statNames{s2};

        pair = sort({name1,name2});
        baselineID_all(i) = string(pair{1}) + "-" + string(pair{2});

    else
        baselineID_all(i) = "invalid-baseline";
    end

end

baselineID = baselineID_all(validWRMS);
uniqueBase = unique(baselineID_all);

baseName = {};
baseWRMS = [];
baseNobs = [];
baseNproblem = [];

for i = 1:numel(uniqueBase)

    idxAll = baselineID_all == uniqueBase(i);
    idx    = baselineID == uniqueBase(i);

    baseName{end+1,1} = char(uniqueBase(i));
    baseNobs(end+1,1) = sum(idxAll);
    baseNproblem(end+1,1) = sum(idxAll & isProblem);

    if any(idx)
        baseWRMS(end+1,1) = calcWRMS(x(idx),mx(idx));
    else
        baseWRMS(end+1,1) = NaN;
    end

end

Tbase = table(baseName,baseNobs,baseNproblem,baseWRMS, ...
    'VariableNames',{'Baseline','Nobs','ProblematicObs','WRMS_cm'});

Tbase.WRMS_cm = round(Tbase.WRMS_cm,2);
Tbase = sortrows(Tbase,'WRMS_cm','descend','MissingPlacement','last');

%% Save CSV tables
if tableout
    writetable(Tsrc,'wrms_sources.csv');
    writetable(Tsta,'wrms_stations.csv');
    writetable(Tbase,'wrms_baselines.csv');
end

disp('WRMS per source')
disp(Tsrc)

disp('WRMS per station')
disp(Tsta)

disp('WRMS per baseline')
disp(Tbase)

%% PDF Page 1: Overview

fig1 = figure('Color','w','Position',[50 50 1600 1000]);

tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile
histogram(x,50)
xlabel('Residual [cm]')
ylabel('Number of observations')
title('Residual distribution')
grid on
box on

nexttile
histogram(z,50)
xlabel('Normalized residual x / mx [-]')
ylabel('Number of observations')
title('Normalized residual distribution')
xline(-3,'--','-3\sigma')
xline( 3,'--','+3\sigma')
grid on
box on




nexttile
scatter(Tsrc.Nobs,Tsrc.WRMS_cm,35,'filled')
hold on
if ~souFIXED
    scatter(Tsrc.Nobs(idso0),Tsrc.WRMS_cm(idso0),35,'filled')
    if souPWLO
        legend('fixed','pwlo')
    else
        legend('NNR','wo NNR')
    end
end
hold off
xlabel('Number of observations')
ylabel('WRMS [cm]')
title('Sources: WRMS vs Nobs')
grid on
box on

nexttile
scatter(Tbase.Nobs,Tbase.WRMS_cm,35,'filled')
xlabel('Number of observations')
ylabel('WRMS [cm]')
title('Baselines: WRMS vs Nobs')
grid on
box on

sgtitle('VLBI post-fit delay residual statistics - overview','FontWeight','bold')


%% PDF Page 2: all sources

fig2 = figure('Color','w','Position',[50 50 2200 1000]);

tiledlayout(1,1,'TileSpacing','compact','Padding','compact');

nexttile
nShow = height(Tsrc);

souorder = [1:nShow];

yyaxis left
% b1 = bar(Tsrc.WRMS_cm);
b1 = plot(souorder,Tsrc.WRMS_cm,'-o','LineWidth',1.1,'MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',clrs(9,:));
if ~souFIXED
    hold on
    b2 = plot(souorder(idso0),Tsrc.WRMS_cm(idso0),'o','LineWidth',1.1,'MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',clrs(8,:));
end
hold off
ylabel('WRMS [cm]')

yyaxis right
p1 = plot(souorder,Tsrc.Nobs,':o','LineWidth',0.5,'MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',clrs(9,:));
if ~souFIXED
    hold on
    b3 = plot(souorder(idso0),Tsrc.Nobs(idso0),'o','LineWidth',1.1,'MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',clrs(8,:));
end

ylabel('Number of observations')

set(gca, ...
    'XLim',[0 nShow+1], ...
    'XTick',1:nShow, ...
    'XTickLabel',Tsrc.Source, ...
    'XTickLabelRotation',90, ...
    'FontSize',10)

title('All sources sorted by WRMS')
if ~souFIXED
    if souPWLO
       legend([b1 b2 p1],{'WRMS fixed','WRMS pwlo','Nobs'},'Location','best')
    else
        legend([b1 b2 p1],{'WRMS NNR','WRMS wo NNR','Nobs'},'Location','best')
    end
else
    legend([b1 p1],{'WRMS','Nobs'},'Location','best')
end
grid on
box on

sgtitle('WRMS and number of observations - all sources','FontWeight','bold')


%% PDF Page 3: all stations and diagnostic bars

fig3 = figure('Color','w','Position',[50 50 1600 1000]);

tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile([1 2])
nShow = height(Tsta);
staorder = [1:nShow];

yyaxis left
% b2 = bar(Tsta.WRMS_cm);
b2 = plot(staorder,Tsta.WRMS_cm,'-o','LineWidth',1.1,'MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',clrs(1,:));
ylabel('WRMS [cm]')

yyaxis right
p2 = plot(staorder,Tsta.Nobs,':o','LineWidth',0.5,'MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',clrs(2,:));
ylabel('Number of observations')

set(gca, ...
    'XLim',[0 nShow+1], ...
    'XTick',1:nShow, ...
    'XTickLabel',Tsta.Station, ...
    'XTickLabelRotation',45, ...
    'FontSize',12)

title('All stations sorted by WRMS')
legend([b2 p2],{'WRMS','Nobs'},'Location','best')
grid on
box on

nexttile
bar([sum(~isProblem), sum(isProblem)])
set(gca,'XTickLabel',{'Usable for WRMS','Problematic'})
ylabel('Number of observations')
title('Usable and problematic observations')
grid on
box on

nexttile
bar([sum(isnan(x_all)), sum(isinf(x_all)), sum(isnan(mx_all)), ...
     sum(isinf(mx_all)), sum(mx_all == 0), sum(mx_all < 0)])
set(gca, ...
    'XTickLabel',{'NaN x','Inf x','NaN mx','Inf mx','mx=0','mx<0'}, ...
    'XTickLabelRotation',45)
ylabel('Number of observations')
title('Types of problematic observations')
grid on
box on

sgtitle('Stations and observation diagnostics','FontWeight','bold')


%% PDF Page 4: all baselines

fig4 = figure('Color','w','Position',[50 50 2400 1000]);

tiledlayout(1,1,'TileSpacing','compact','Padding','compact');

nexttile
nShow = height(Tbase);
basorder = [1:nShow];

yyaxis left
% b3 = bar(Tbase.WRMS_cm);
b3 = plot(basorder,Tbase.WRMS_cm,'-o','LineWidth',1.1,'MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',clrs(1,:));

ylabel('WRMS [cm]')

yyaxis right
p3 = plot(1:nShow,Tbase.Nobs,':o','LineWidth',0.5,'MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor',clrs(2,:));
ylabel('Number of observations')

set(gca, ...
    'XLim',[0 nShow+1], ...
    'XTick',1:nShow, ...
    'XTickLabel',Tbase.Baseline, ...
    'XTickLabelRotation',90, ...
    'FontSize',14)

title('All baselines sorted by WRMS')
legend([b3 p3],{'WRMS','Nobs'},'Location','best')
grid on
box on

sgtitle('WRMS and number of observations - all baselines','FontWeight','bold')


%% PDF Page 5: problematic observations by source and station

fig5 = figure('Color','w','Position',[50 50 1800 1000]);

tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile
bar(Tsrc.ProblematicObs)
set(gca, ...
    'XLim',[0 height(Tsrc)+1], ...
    'XTick',1:height(Tsrc), ...
    'XTickLabel',Tsrc.Source, ...
    'XTickLabelRotation',90, ...
    'FontSize',6)
ylabel('Problematic observations')
title('Problematic observations per source')
grid on
box on

nexttile
bar(Tsta.ProblematicObs)
set(gca, ...
    'XTick',1:height(Tsta), ...
    'XTickLabel',Tsta.Station, ...
    'XTickLabelRotation',45, ...
    'FontSize',8)
ylabel('Problematic observations')
title('Problematic observations per station')
grid on
box on

sgtitle('Problematic observations by source and station','FontWeight','bold')


%% PDF Page 6: boxplot residuals per station

stationGroup = {};
stationResiduals = [];

for i = 1:nSta

    idx = baselines(:,1) == i | baselines(:,2) == i;

    if any(idx)

        stationResiduals = [stationResiduals; x(idx)];
        stationGroup = [stationGroup; repmat(statNames(i),sum(idx),1)];

    end

end

fig6 = figure('Color','w','Position',[50 50 1800 900]);

boxplot(stationResiduals,stationGroup)

%xlabel('Station')
ylabel('Residual [cm]')
title('Residual distribution per station')
xtickangle(45)
grid on
box on

if tableout
    fprintf('\nSaved CSV files:\n')
    fprintf('  observation_diagnostics.csv\n')
    fprintf('  residual_global_statistics.csv\n')
    fprintf('  wrms_sources.csv\n')
    fprintf('  wrms_stations.csv\n')
    fprintf('  wrms_baselines.csv\n')
end

%% Local functions

function wrms = calcWRMS(x,mx)

    x  = double(x(:));
    mx = double(mx(:));

    valid = isfinite(x) & isfinite(mx) & mx > 0;

    x  = x(valid);
    mx = mx(valid);

    if isempty(x)
        wrms = NaN;
        return
    end

    P = 1 ./ (mx.^2);

    wrms = sqrt(sum((x.^2).*P) / sum(P));

end

end