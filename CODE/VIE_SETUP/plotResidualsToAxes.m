% #########################################################################
% #     plotResidualsToAxes
% #########################################################################
%
% DESCRITPION
% This function plots the residuals as chosen in the interface
% (stationwise, baseline wise, all) to the axes in the panel of the
% residuals.
%
% AUTHOR
%   Matthias Madzak
%
% INPUT
%   handles     structure from the GUI (also containing data, eg residuals)
%
% OUTPUT
%   handles     structure from the GUI (also containing data, eg residuals)
%
% CHANGES
%   2014-06-04, A. Hellerschmied: Use "hours since session
%       start" as time scale in the plotting residuals window.
%   2016-09-10, A. Girdiuk: function was suited for 'Select Data' panel usage
%   2016-09-21, A. Girdiuk: typos-fix
%   2016-10-10, A. Hellerschmied: Enabled satellite sources
%

function handles = plotResidualsToAxes(handles)

% ##### Options #####
% Set plotstyle for residual values; 1 = lines, 2 = lines + markers, 3 = scatterplot
plotstyle = 2;


% ##### Init #####
% "clear" outlier stuff
set(handles.axes_plot_residuals,                'ButtonDownFcn', '')
set(handles.figure_vievs2,              'WindowButtonUpFcn', '')
set(handles.togglebutton_plot_residuals_selectOutliers, 'Value', 0)
set(handles.pushbutton_plot_residuals_removeOutliers, 'Enable', 'Off')
set(handles.pushbutton_plot_residuals_removeOutliers, 'String', 'Remove Outliers')

set(handles.togglebutton_plot_residuals_selectData,'Value',0)
% set(handles.togglebutton_plot_residuals_SelectedData_output, 'Enable', 'Off')
set(handles.pushbutton_selectedData_writeOPT,'Enable','Off');
% set(handles.radiobutton_unit_plot, 'Enable', 'Off')
% set(handles.radiobutton_unit_UTC, 'Enable', 'Off')
% set(handles.radiobutton_unit_MJD, 'Enable', 'Off')
% set(handles.edit_plot_residuals_interval_input_1_unit, 'Enable', 'Off')
% set(handles.edit_plot_residuals_interval_input_2_unit, 'Enable', 'Off')
% set(handles.edit_plot_residuals_interval_show,'String','xx.xx - xx.xx')

% set(handles.edit_plot_residuals_interval_input_1_unit, 'String', 'hh')
% set(handles.edit_plot_residuals_interval_input_2_unit, 'String', 'hh')

% remove all black object (box and crosses)
allLineHandles=findobj(handles.axes_plot_residuals,'Type','line', 'color', 'k');
delete(allLineHandles);
allLineHandles=findobj(handles.axes_plot_residuals,'Type','line', ...
    'color', [0 0 0.04]);
delete(allLineHandles);

curSession=get(handles.popupmenu_plot_residuals_session, 'Value');
curSessionName = get(handles.popupmenu_plot_residuals_session, 'String');
curSessionName = curSessionName{curSession};

SessionStartTimeMJD =  handles.data.plot.res(curSession).mjd(1);
SessionEndTimeMJD   =  (handles.data.plot.res(curSession).mjd(end)-SessionStartTimeMJD)*24;

mfc = [31,120,180
    51,160,44
    227,26,28
    255,127,0
    202,178,214
    106,61,154
    177,89,40
    166,206,227
    178,223,138
    251,154,153
    253,191,111
    202,178,214
    255,255,153
    ]/255;

mfc = repmat(mfc,8,1);

mec = [repmat([0 0 0],13,1)
    repmat([1 1 0],13,1)
    repmat([1 0 1],13,1)
    repmat([0 1 1],13,1)
    repmat([1 1 1],13,1)
    repmat([1 0 0],13,1)
    repmat([0 1 0],13,1)
    repmat([0 0 1],13,1)
    ];

% ##### Choose between First / Main Solution #####

% if first solution is chosen
if get(handles.radiobutton_plot_residuals_firstSolution, 'Value')
    val=handles.data.plot.res(curSession).firstVal;
    if ~isempty(handles.data.plot.res(curSession).sigma_residuals_aposteriori)
        val_sigma_cm = handles.data.plot.res(curSession).sigma_from_fringe_fitting'*physconst('LightSpeed')*1e2;
    else
        fprintf('sigma of observation does not exist for this session, re-run the current session to store sigmas for plotting\n')
        val_sigma_cm = zeros(length(val),1);
    end
    % if main solution is chosen
else
    val=handles.data.plot.res(curSession).mainVal;
    if ~isempty(handles.data.plot.res(curSession).sigma_residuals_aposteriori)
        val_sigma_cm = handles.data.plot.res(curSession).sigma_residuals_aposteriori;
    else
        fprintf('sigma of observation does not exist for this session, re-run the current session to store sigmas for plotting\n')
        val_sigma_cm = zeros(length(val),1);        
    end
end

% Residuals limits for plot
ResLim(1) = floor(min(val)); ResLim(2) = ceil(max(val));

% Get data from handles-structure:
outlier     = handles.data.plot.res(curSession).outlier;
antennas    = handles.data.plot.res(curSession).allStatNames;
baselines   = handles.data.plot.res(curSession).allStatNames(handles.data.plot.res(curSession).baselineOfObs);
% - Name of observed source for each observation in session:
obs_type    = handles.data.plot.res(curSession).obs_type;
sat_ind     = strcmp(cellstr(obs_type), 's');
quasar_ind  = strcmp(cellstr(obs_type), 'q');
source_id_list = handles.data.plot.res(curSession).source;
if sum(sat_ind) > 0
    sources(sat_ind) = handles.data.plot.res(curSession).allSatelliteNames(source_id_list(sat_ind));
end
if  sum(quasar_ind) > 0
    sources(quasar_ind) = handles.data.plot.res(curSession).allSourceNames(source_id_list(quasar_ind));
end
sources = sources';

mjd             = handles.data.plot.res(curSession).mjd;
DurationHours   = (mjd - SessionStartTimeMJD) * 24;


% ##### make axes empty #####
cla(handles.axes_plot_residuals);
% set axes as current axes and hold on
set(gcf,'CurrentAxes',handles.axes_plot_residuals)
%hold(handles.axes_plot_residuals, 'on');

% ##### Listener for Axes Changes ####
if isfield(handles, 'axesListener')
    delete(handles.axesListener);
end
hax1 = handles.axes_plot_residuals;
prop1 = findPropertyHandle( hax1 , 'XLim' ) ;
handles.axesListener = addlistener( hax1, prop1 , 'PostSet', @() [] ) ;
pl.props2link = 'XLim' ;   % name of properties to link
pl.listeners  = handles.axesListener ;         % listener handle
pl.metaprops2link = prop1; % metaproperty handles
pl.handles = handles;
handles.axesListener.Callback = @(h,e) AxesLimitsChanged(h,e,pl) ;


plotOutliers=0;

background = [val; -val];
backgroundTime = [DurationHours; DurationHours];

% ##### Selection of Residual Values for Plotting #####

plotName = '';
% #### 1.) ALL RESIDUALS ####
% if all residuals should be plotted
if get(handles.radiobutton_plot_residuals_perAll, 'Value')
    plotName = 'all';
    % all residuals should be plotted when "per All" is chosen!
    valsOfCurSelection=[val; -val];
    sigsOfCurSelection=[val_sigma_cm; -val_sigma_cm];
    horAxis = [DurationHours; DurationHours];
    
    % see if there are outliers
    allOutliersLog=zeros(length(val),1);
    allOutliersLog(outlier)=1;
    
    % get all (no real selection as in station/source/baseline wise)
    if sum(allOutliersLog) > 0
        plotOutliers=1; % Flag
        outlierIndIWantToFind=1:length(outlier);
        indOutliersOfCurSelection=outlier;
        
        % preallocate
        baselinesForOutlier=cell(length(outlier), 1);
        
        for k=1:length(outlier)
            baselinesForOutlier(k,1)=baselines(outlier(k),1);
            baselinesForOutlier(k,2)=baselines(outlier(k),2);
        end
    end
    set(handles.togglebutton_plot_residuals_selectData, 'Enable', 'Off');
    
    
    % #### 2.) STATION-WISE ####
    % if station-wise residuals should be plotted
elseif get(handles.radiobutton_plot_residuals_perStat, 'Value')
    % get chosen station
    allStationsInMenu=get(handles.popupmenu_plot_residuals_station, 'String');
    curStation=allStationsInMenu{get(handles.popupmenu_plot_residuals_station, 'Value')};
    plotName = curStation;
    % get values where chosen station takes part
    obsWithCurStation=sum(~cellfun(@isempty, strfind(baselines, curStation)),2);
    valsOfCurSelection=val(logical(obsWithCurStation));
    sigsOfCurSelection=val_sigma_cm(logical(obsWithCurStation));
    horAxis = DurationHours(logical(obsWithCurStation));
    
    % for values where station is second: *(-1)
    stationIsSecond=~cellfun(@isempty, strfind(baselines(logical(obsWithCurStation),2), curStation));
    valsOfCurSelection(stationIsSecond)=valsOfCurSelection(stationIsSecond)*-1;
    sigsOfCurSelection(stationIsSecond)=sigsOfCurSelection(stationIsSecond)*-1;
    
    
    % see if this station has outliers
    allOutliersLog=zeros(length(val),1);
    allOutliersLog(outlier)=1;
    
    % if there are outliers at chosen station
    if sum(allOutliersLog & obsWithCurStation) > 0
        plotOutliers=1;
        
        % preallocate
        indOutliersOfCurSelection=zeros(sum(allOutliersLog & obsWithCurStation),1); % those indices show values in chosenStation-values
        baselinesForOutlier=cell(sum(allOutliersLog & obsWithCurStation), 1);
        
        % get index (in outlier) of outlier I want to find (ie all which include curStat)
        outlierIndIWantToFind=find(ismember(outlier, find(allOutliersLog & obsWithCurStation)));
        
        % find them
        for k=1:sum(allOutliersLog & obsWithCurStation)
            indOutliersOfCurSelection(k)=sum(obsWithCurStation(1:outlier(outlierIndIWantToFind(k))));
            baselinesForOutlier(k,1)=baselines(outlier(outlierIndIWantToFind(k)),1);
            baselinesForOutlier(k,2)=baselines(outlier(outlierIndIWantToFind(k)),2);
        end
    end
    
    set(handles.togglebutton_plot_residuals_selectData, 'Enable', 'On');
    
    % #### 3.) BASELINE-WISE ####
    % baseline wise plot
elseif get(handles.radiobutton_plot_residuals_perBasel, 'Value')
    % get chosen baseline
    allBaselinesInMenu=get(handles.popupmenu_plot_residuals_baseline, 'String');
    curBaseline=allBaselinesInMenu{get(handles.popupmenu_plot_residuals_baseline, 'Value')};
    plotName = curBaseline;
    
    % get values of chosen baseline
    obsWithCurSelection=sum(~cellfun(@isempty, strfind(baselines, curBaseline(1:8))),2) & ...
        sum(~cellfun(@isempty, strfind(baselines, curBaseline(10:17))),2);
    valsOfCurSelection=val(logical(obsWithCurSelection));
    sigsOfCurSelection=val_sigma_cm(logical(obsWithCurSelection));
    horAxis = DurationHours(logical(obsWithCurSelection));
    
    % see if this station has outliers
    allOutliersLog=zeros(length(val),1);
    allOutliersLog(outlier)=1;
    
    if sum(allOutliersLog & obsWithCurSelection) > 0
        plotOutliers=1;
        
        % get station names of baseline (which is equal for all
        % observationsin baseline-wise plots!)
        baselinesForOutlier=...
            repmat([{curBaseline(1:8)}, {curBaseline(10:17)}], length(outlier), 1);
        
        % preallocate
        indOutliersOfCurSelection=zeros(sum(allOutliersLog & obsWithCurSelection),1);
        
        % get index (in outlier) of outlier I want to find (ie all which include curStat)
        outlierIndIWantToFind=find(ismember(outlier, find(allOutliersLog & obsWithCurSelection)));
        
        % find them
        for k=1:sum(allOutliersLog & obsWithCurSelection)
            indOutliersOfCurSelection(k)=sum(obsWithCurSelection(1:outlier(outlierIndIWantToFind(k))));
        end
    end
    
    set(handles.togglebutton_plot_residuals_selectData, 'Enable', 'Off');
    
    % #### 4.) SOURCE-WISE ####
    % source wise residuals
else
    % get chosen source
    allSourcesInMenu=get(handles.popupmenu_plot_residuals_source, 'String');
    curSource=allSourcesInMenu{get(handles.popupmenu_plot_residuals_source, 'Value')};
    plotName = curSource;
    
    % get values of chosen source
    obsWithCurSelection=strcmp(sources, curSource);
    valsOfCurSelection=val(logical(obsWithCurSelection));
    sigsOfCurSelection=val_sigma_cm(logical(obsWithCurSelection));
    horAxis = DurationHours(logical(obsWithCurSelection));
    
    % see if this station has outliers
    allOutliersLog=zeros(length(val),1);
    allOutliersLog(outlier)=1;
    
    if sum(allOutliersLog & obsWithCurSelection) > 0
        plotOutliers=1;
        
        % preallocate
        indOutliersOfCurSelection=zeros(sum(allOutliersLog & obsWithCurSelection),1);
        baselinesForOutlier=cell(sum(allOutliersLog & obsWithCurSelection), 1);
        
        % get index (in outlier) of outlier I want to find (ie all which include curSource)
        outlierIndIWantToFind=find(ismember(outlier, find(allOutliersLog & obsWithCurSelection)));
        
        % find them
        for k=1:sum(allOutliersLog & obsWithCurSelection)
            indOutliersOfCurSelection(k)=sum(obsWithCurSelection(1:outlier(outlierIndIWantToFind(k))));
            baselinesForOutlier(k,1)=baselines(outlier(outlierIndIWantToFind(k)),1);
            baselinesForOutlier(k,2)=baselines(outlier(outlierIndIWantToFind(k)),2);
        end
    end
    
    set(handles.togglebutton_plot_residuals_selectData, 'Enable', 'On');
    
end % end - all/sourcewise/stationwise/baselinewise


% ##### PLOTTING #####

% #### plot residuals ####

hold(handles.axes_plot_residuals, 'on');
handles.data.plot.plottedResiduals=[];

% Set X-Lim Mode to 'auto'
% set(handles.axes_plot_residuals, 'XLimMode', 'auto')
set(handles.axes_plot_residuals, 'XLim', [-SessionEndTimeMJD*0.1 SessionEndTimeMJD*1.1]);
% set(handles.axes_plot_residuals, 'XTick',0:6:24);
% xTicksLabel = datetime(get(handles.axes_plot_residuals,'XTick')/24+SessionStartTimeMJD,'ConvertFrom','ModifiedJulianDate');
% hhmm = datestr(xTicksLabel,'hh:MM');
% date = unique(datestr(xTicksLabel,'dd.mm.yyyy'),'rows');
% set(handles.axes_plot_residuals, 'XTickLabel',hhmm);
% dateLab = '';
% for i=1:size(date,1)
%     if i == 1
%         dateLab = [dateLab date(i,:)];
%     else
%         dateLab = [dateLab ' - ' date(i,:)];
%     end
% end
% xlabel(dateLab);

switch plotstyle
    case 1 % lines only
        handles.data.plot.plottedResiduals = plot(handles.axes_plot_residuals, horAxis, valsOfCurSelection, 'b-','HitTest','off');
    case 2 % lines and markers
        hold on
        plot(handles.axes_plot_residuals, backgroundTime, background, 'Marker','.', 'LineStyle','none','MarkerSize',5,'Color',[.7 .7 .7],'HitTest','off');
        %         handles.data.plot.plottedResiduals = plot(handles.axes_plot_residuals, horAxis, valsOfCurSelection, 'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor',[55,126,184]/255,'MarkerSize',8,'HitTest','off');
        handles.data.plot.plottedResiduals = errorbar(handles.axes_plot_residuals, horAxis, valsOfCurSelection,sigsOfCurSelection, 'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor',[55,126,184]/255,'Color',[55,126,184]/255,'MarkerSize',8,'HitTest','off','CapSize',0);
        
        % ### plot all ###
        if get(handles.radiobutton_plot_residuals_perAll, 'Value')
            h = [];
            allStationsInMenu=get(handles.popupmenu_plot_residuals_station, 'String');
            for i=1:length(allStationsInMenu)
                curStation=allStationsInMenu{i};
                tmp = contains(baselines, curStation);
                pos = tmp(:,1);
                neg = tmp(:,2);
                obsWithCurStation = [val(pos); val(neg)*-1];
                sigWithCurStation = [val_sigma_cm(pos); val_sigma_cm(neg)*-1];
                timWithCurStation = [DurationHours(pos); DurationHours(neg)];
                %                 h(i) = plot(handles.axes_plot_residuals, timWithCurStation, obsWithCurStation,'LineStyle','none','Marker','o','MarkerEdgeColor',mec(i,:),'MarkerFaceColor',mfc(i,:),'MarkerSize',8,'HitTest','off');
                h(i) = errorbar(handles.axes_plot_residuals, timWithCurStation, obsWithCurStation,sigWithCurStation,'LineStyle','none','Marker','o','MarkerEdgeColor',mec(i,:),'MarkerFaceColor',mfc(i,:),'Color',mfc(i,:),'MarkerSize',8,'HitTest','off','CapSize',0);
                
            end
            legend(h,allStationsInMenu);
            % get values where chosen station takes part
            
            % ### plot station ###
        elseif get(handles.radiobutton_plot_residuals_perStat, 'Value')
            
            idx=get(handles.popupmenu_plot_residuals_station, 'Value');
            name = allStationsInMenu{idx};
            handles.data.plot.plottedResiduals.MarkerFaceColor = mfc(idx,:);
            handles.data.plot.plottedResiduals.MarkerEdgeColor = mec(idx,:);
            
            if isempty(valsOfCurSelection)
                legend({'all'});
            else
                legend({'all',name});
            end
            
            % ### plot baseline ###
        elseif get(handles.radiobutton_plot_residuals_perBasel, 'Value')
            name=allBaselinesInMenu{get(handles.popupmenu_plot_residuals_baseline, 'Value')};
            
            if isempty(valsOfCurSelection)
                legend({'all'});
            else
                legend({'all',name});
            end
            
            % ### plot source ###
        else
            obsWithCurSelection=strcmp(sources, curSource);
            names = {};
            h = [];
            allStationsInMenu=get(handles.popupmenu_plot_residuals_station, 'String');
            for i=1:length(allStationsInMenu)
                curStation=allStationsInMenu{i};
                tmp = contains(baselines, curStation);
                pos = tmp(:,1);
                neg = tmp(:,2);
                obsWithCurStation = [val(pos & obsWithCurSelection); val(neg & obsWithCurSelection)*-1];
                sigWithCurStation = [val_sigma_cm(pos & obsWithCurSelection); val_sigma_cm(neg & obsWithCurSelection)*-1];
                timWithCurStation = [DurationHours(pos & obsWithCurSelection); DurationHours(neg & obsWithCurSelection)];
                if ~isempty(obsWithCurStation)
                    %                     h(end+1) = plot(handles.axes_plot_residuals, timWithCurStation, obsWithCurStation,'LineStyle','none','Marker','o','MarkerEdgeColor',mec(i,:),'MarkerFaceColor',mfc(i,:),'MarkerSize',8,'HitTest','off');
                    h(end+1) = errorbar(handles.axes_plot_residuals, timWithCurStation, obsWithCurStation,sigWithCurStation,'LineStyle','none','Marker','o','MarkerEdgeColor',mec(i,:),'MarkerFaceColor',mfc(i,:),'Color',mfc(i,:),'MarkerSize',8,'HitTest','off','CapSize',0);
                    names{end+1} = curStation;
                end
            end
            legend(h,names);
            
        end
        
    case 3 % scatterplot
        handles.data.plot.plottedResiduals = plot(handles.axes_plot_residuals, horAxis, valsOfCurSelection, 'b*','HitTest','off');
end

% Label X- and Y-Axis:
ylabel('Residuals [cm]');
% title(strrep(sprintf('%s (%s)',curSessionName,plotName),'_','\_'));
% title(sprintf('residual limits [%5.0f: %5.0f] [cm]',floor(min(valsOfCurSelection)),ceil(max(valsOfCurSelection))));


% #### plot outliers ####
if plotOutliers==1 && ...
        get(handles.radiobutton_plot_residuals_mainSolution, 'Value')
    
    if get(handles.radiobutton_plot_residuals_perStat, 'Value') || get(handles.radiobutton_plot_residuals_perBasel, 'Value')
        handles.data.plot.plottedOutlierBoxes=...
            plot(horAxis(indOutliersOfCurSelection), valsOfCurSelection(indOutliersOfCurSelection), 'x', 'color', 'k', 'markersize', 10, 'DisplayName', 'outliers','LineWidth',3);
    else
        handles.data.plot.plottedOutlierBoxes=...
            plot([horAxis(indOutliersOfCurSelection); horAxis(indOutliersOfCurSelection)], ...
            [valsOfCurSelection(indOutliersOfCurSelection); -valsOfCurSelection(indOutliersOfCurSelection)], ...
            'x', 'color', 'k', 'markersize', 10, 'DisplayName', 'outliers','LineWidth',3);
    end
    
    % plot station indices / source names to outliers
    for k=1:length(indOutliersOfCurSelection)
        curStr='';
        % if station numbers should be printed
        %         if get(handles.checkbox_plot_residuals_showStatNumbers, 'Value')
        %             ant1=find(~cellfun(@isempty, strfind(antennas, baselinesForOutlier{k,1})));
        %             ant2=find(~cellfun(@isempty, strfind(antennas, baselinesForOutlier{k,2})));
        %             curStr=[curStr, num2str(ant1),'-',num2str(ant2)];
        %         end
        
        % if source names should be printed
        if get(handles.checkbox_plot_residuals_showSourceNames, 'Value')
            if ~isempty(curStr)
                curStr=[curStr, ' '];
            end
            curStr=[curStr, char(sources(outlier(outlierIndIWantToFind(k))))];
        end
        text(horAxis(indOutliersOfCurSelection(k)), valsOfCurSelection(indOutliersOfCurSelection(k)), ...
            ['\fontsize{14}', '\color[rgb]{0 0 0}',...
            curStr], 'HitTest', 'off')
    end
end

ylim=get(gca,'ylim');
xlim=get(gca,'xlim');


txt = '';
if get(handles.radiobutton_plot_residuals_firstSolution, 'Value')
    if get(handles.radiobutton_plot_residuals_perAll, 'Value')
        nObs = length(valsOfCurSelection)/2;
    else
        nObs = length(valsOfCurSelection);
    end
    txt = [txt sprintf('#obs: %d of %d\n',nObs,length(val))];
    
    if ~isempty(handles.data.plot.res(curSession).mo_first)
        txt = [txt sprintf('chi^2: %.3f',handles.data.plot.res(curSession).mo_first^2)];
    end
else
    if get(handles.radiobutton_plot_residuals_perAll, 'Value')
        nObs = length(valsOfCurSelection)/2;
    else
        nObs = length(valsOfCurSelection);
    end
    txt = [txt sprintf('#obs: %d of %d\n',nObs,length(val))];
    if ~isempty(handles.data.plot.res(curSession).wrms)
        txt = [txt sprintf('wrms: %.3f cm (%.3f ps)\n', handles.data.plot.res(curSession).wrms, handles.data.plot.res(curSession).wrms*100/2.99792458)];
    end
    
    if ~isempty(handles.data.plot.res(curSession).mo)
        txt = [txt sprintf('chi^2: %.3f',handles.data.plot.res(curSession).mo^2)];
    end
end

text(xlim(1)+(xlim(2)-xlim(1))*0.025,ylim(2)-(ylim(2)-ylim(1))*0.025,txt,'VerticalAlignment','top');

% txt2 = datestr(datetime(SessionStartTimeMJD,'ConvertFrom','modifiedJulianDate'));
% text(xlim(1)+(xlim(2)-xlim(1))*0.025,ylim(1)+(ylim(2)-ylim(1))*0.05,['first observation: ' txt2],'VerticalAlignment','top');

hold off

hold(handles.axes_plot_residuals, 'off')