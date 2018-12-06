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
set(handles.togglebutton_plot_residuals_SelectedData_output, 'Enable', 'Off')
set(handles.pushbutton_selectedData_writeOPT,'Enable','Off');
set(handles.radiobutton_unit_plot, 'Enable', 'Off')
set(handles.radiobutton_unit_UTC, 'Enable', 'Off')
set(handles.radiobutton_unit_MJD, 'Enable', 'Off')
set(handles.edit_plot_residuals_interval_input_1_unit, 'Enable', 'Off')
set(handles.edit_plot_residuals_interval_input_2_unit, 'Enable', 'Off')
set(handles.edit_plot_residuals_interval_show,'String','xx.xx - xx.xx')

set(handles.edit_plot_residuals_interval_input_1_unit, 'String', 'hh')
set(handles.edit_plot_residuals_interval_input_2_unit, 'String', 'hh')
    
% remove all black object (box and crosses)
allLineHandles=findobj(handles.axes_plot_residuals,'Type','line', 'color', 'k');
delete(allLineHandles);
allLineHandles=findobj(handles.axes_plot_residuals,'Type','line', ...
    'color', [0 0 0.04]);
delete(allLineHandles);

curSession=get(handles.popupmenu_plot_residuals_session, 'Value');

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

mec = [repmat([0 0 0],12,1)
       repmat([1 1 1],12,1)
       repmat([1 0 0],12,1)
       repmat([0 1 0],12,1)
       repmat([0 0 1],12,1)
       repmat([1 1 0],12,1)
       repmat([1 0 1],12,1)
       repmat([0 1 1],12,1)
       ];

% ##### Choose between First / Main Solution #####

% if first solution is chosen
if get(handles.radiobutton_plot_residuals_firstSolution, 'Value')
    val=handles.data.plot.res(curSession).firstVal;
% if main solution is chosen
else
    val=handles.data.plot.res(curSession).mainVal;
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

plotOutliers=0;

background = [val; -val];
backgroundTime = [DurationHours; DurationHours];

% ##### Selection of Residual Values for Plotting #####
    
% #### 1.) ALL RESIDUALS ####
% if all residuals should be plotted
if get(handles.radiobutton_plot_residuals_perAll, 'Value')
    % all residuals should be plotted when "per All" is chosen!
    valsOfCurSelection=[val; -val];
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
            
    
% #### 2.) STATION-WISE ####
% if station-wise residuals should be plotted    
elseif get(handles.radiobutton_plot_residuals_perStat, 'Value')
    % get chosen station
    allStationsInMenu=get(handles.popupmenu_plot_residuals_station, 'String');
    curStation=allStationsInMenu{get(handles.popupmenu_plot_residuals_station, 'Value')};
    % get values where chosen station takes part
    obsWithCurStation=sum(~cellfun(@isempty, strfind(baselines, curStation)),2);
    valsOfCurSelection=val(logical(obsWithCurStation));
    horAxis = DurationHours(logical(obsWithCurStation));

    % for values where station is second: *(-1)
    stationIsSecond=~cellfun(@isempty, strfind(baselines(logical(obsWithCurStation),2), curStation));
    valsOfCurSelection(stationIsSecond)=valsOfCurSelection(stationIsSecond)*-1;

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
        
    
% #### 3.) BASELINE-WISE ####
% baseline wise plot
elseif get(handles.radiobutton_plot_residuals_perBasel, 'Value') 
    % get chosen baseline
    allBaselinesInMenu=get(handles.popupmenu_plot_residuals_baseline, 'String');
    curBaseline=allBaselinesInMenu{get(handles.popupmenu_plot_residuals_baseline, 'Value')};

    % get values of chosen baseline
    obsWithCurSelection=sum(~cellfun(@isempty, strfind(baselines, curBaseline(1:8))),2) & ...
        sum(~cellfun(@isempty, strfind(baselines, curBaseline(10:17))),2);
    valsOfCurSelection=val(logical(obsWithCurSelection));
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
    
    
% #### 4.) SOURCE-WISE ####    
% source wise residuals
else 
    % get chosen source
    allSourcesInMenu=get(handles.popupmenu_plot_residuals_source, 'String');
    curSource=allSourcesInMenu{get(handles.popupmenu_plot_residuals_source, 'Value')};

    % get values of chosen source
    obsWithCurSelection=strcmp(sources, curSource);
    valsOfCurSelection=val(logical(obsWithCurSelection));
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

end % end - all/sourcewise/stationwise/baselinewise
    

% ##### PLOTTING #####

% #### plot residuals ####

hold(handles.axes_plot_residuals, 'on');
handles.data.plot.plottedResiduals=[];

% Set X-Lim Mode to 'auto'
set(handles.axes_plot_residuals, 'XLimMode', 'auto')
set(handles.axes_plot_residuals, 'XLim', [-1 SessionEndTimeMJD+1]);
set(handles.axes_plot_residuals, 'XTick',0:6:24);

switch plotstyle
    case 1 % lines only
        handles.data.plot.plottedResiduals = plot(handles.axes_plot_residuals, horAxis, valsOfCurSelection, 'b-','HitTest','off');
    case 2 % lines and markers
        hold on
        plot(handles.axes_plot_residuals, backgroundTime, background, 'Marker','.', 'LineStyle','none','MarkerSize',5,'Color',[.7 .7 .7],'HitTest','off');
        handles.data.plot.plottedResiduals = plot(handles.axes_plot_residuals, horAxis, valsOfCurSelection, 'LineStyle','none','Marker','o','MarkerEdgeColor','k','MarkerFaceColor',[55,126,184]/255,'MarkerSize',8,'HitTest','off');

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
                timWithCurStation = [DurationHours(pos); DurationHours(neg)];
                h(i) = plot(handles.axes_plot_residuals, timWithCurStation, obsWithCurStation,'LineStyle','none','Marker','o','MarkerEdgeColor',mec(i,:),'MarkerFaceColor',mfc(i,:),'MarkerSize',8,'HitTest','off');
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
            name=allSourcesInMenu{get(handles.popupmenu_plot_residuals_source, 'Value')};
            
            if isempty(valsOfCurSelection)
                legend({'all'});
            else
                legend({'all',name});
            end
        end
        
        hold off
    case 3 % scatterplot
        handles.data.plot.plottedResiduals = plot(handles.axes_plot_residuals, horAxis, valsOfCurSelection, 'b*','HitTest','off');
end

% Label X- and Y-Axis:
xlabel('Hours from session start');
ylabel('Residuals [cm]');
title(sprintf('* Residuals in limits [%5.0f: %5.0f] [cm]',floor(min(valsOfCurSelection)),ceil(max(valsOfCurSelection))));

% #### plot outliers ####
if plotOutliers==1 && ...
        get(handles.radiobutton_plot_residuals_mainSolution, 'Value')

    handles.data.plot.plottedOutlierBoxes=...
        plot(horAxis(indOutliersOfCurSelection), valsOfCurSelection(indOutliersOfCurSelection), 'rs', 'HitTest', 'off');

    % plot station indices / source names to outliers
    for k=1:length(indOutliersOfCurSelection)
        curStr='';
        % if station numbers should be printed
        if get(handles.checkbox_plot_residuals_showStatNumbers, 'Value')
            ant1=find(~cellfun(@isempty, strfind(antennas, baselinesForOutlier{k,1})));
            ant2=find(~cellfun(@isempty, strfind(antennas, baselinesForOutlier{k,2})));
            curStr=[curStr, num2str(ant1),'-',num2str(ant2)];
        end

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


hold(handles.axes_plot_residuals, 'off')   