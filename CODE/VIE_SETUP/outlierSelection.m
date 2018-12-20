% #########################################################################
% #     outlierSelection
% #########################################################################
%
% DESCRITPION
% This function creates a box for the selected area and marks the selected
% values in the box.
%
% AUTHOR 
%   Matthias Madzak
%
% INPUT
%   handles     structure from the GUI (also containing data, eg residuals)
%   hObject
%   eventdata
%
% OUTPUT
%
% CHANGES
%   2014-06-04, A. Hellerschmied: Adaption to use "hours since session
%       start" as time scale in the plotting residuals window.
%   2016-09-10, A. Girdiuk: function was suited for 'Select Data' panel usage
%   2017-04-04, A. Girdiuk: errors were caught when default 'hh' was not
%                           replaced with desired values

function outlierSelection(hObject, eventdata, handles)

% set "Remove Outleir" button to enable
%set(handles.pushbutton_plot_residuals_removeOutliers, 'Enable', 'On')

handles.data.plot.intervalSelectedatPlot = 1;

% if ~get(handles.togglebutton_plot_residuals_SelectedData_output,'Value')
% get current position
curPoint=get(get(hObject, 'CurrentAxes'), 'CurrentPoint');
curPoint=curPoint(1,1:2);
% draw lines ...
x=[min(curPoint(1), handles.data.plot.outlierBoxStart(1)), ...
    max(curPoint(1), handles.data.plot.outlierBoxStart(1))];
y=[min(curPoint(2), handles.data.plot.outlierBoxStart(2)), ...
    max(curPoint(2), handles.data.plot.outlierBoxStart(2))];
% end

% get current position
%curPoint=get(get(hObject, 'CurrentAxes'), 'CurrentPoint');
%curPoint=curPoint(1,1:2);


% if old lines exist - delete them
allLineHandles=findobj(handles.axes_plot_residuals,'Type','line', 'color', 'k');
delete(allLineHandles);
allLineHandles=findobj(handles.axes_plot_residuals,'Type','line', 'color', [0 0 0.04]);
delete(allLineHandles);
% if ~isempty(allLineHandles)
%     
% end

% draw lines ...
%x=[min(curPoint(1), handles.data.plot.outlierBoxStart(1)), ...
 %   max(curPoint(1), handles.data.plot.outlierBoxStart(1))];
%y=[min(curPoint(2), handles.data.plot.outlierBoxStart(2)), ...
 %   max(curPoint(2), handles.data.plot.outlierBoxStart(2))];
%handles.data.plot.lineHandle=line(...
 %   [x(1), x(2), x(2), x(1), x(1)], [y(2), y(2), y(1), y(1), y(2)],...
  %  'color', 'k');

% ... and mark values within box

% get values which are currently plotted
% get index of currently selected session
curSession=get(handles.popupmenu_plot_residuals_session, 'Value');

% #### Per Station ####
if get(handles.radiobutton_plot_residuals_perStat, 'Value')
    curStation=handles.data.plot.res(curSession).allStatNames{...
        get(handles.popupmenu_plot_residuals_station, 'Value')};
    baselineNames=handles.data.plot.res(curSession).allStatNames(handles.data.plot.res(curSession).baselineOfObs);
    obsWithCurSelection=sum(~cellfun(@isempty, strfind(...
        baselineNames, ...
        curStation)),2);
    
    % if station-wise is chosen - again: multiply second-antenna-observations
    % with -1 (do the multiplication later)
    secondAntObservations=~cellfun(@isempty, strfind(baselineNames(logical(obsWithCurSelection),2), curStation));
    handles.data.plot.currentStation = curStation;
% #### Per Baseline ####    
elseif get(handles.radiobutton_plot_residuals_perBasel, 'Value')
    allBaselinesInMenu=get(handles.popupmenu_plot_residuals_baseline, 'String');
    curBaseline=allBaselinesInMenu{get(handles.popupmenu_plot_residuals_baseline, 'Value')}; % chosen (in popupmenu) baseline
    baselineNames=handles.data.plot.res(curSession).allStatNames(handles.data.plot.res(curSession).baselineOfObs); % baseline names for all observations
    obsWithCurSelection=sum(~cellfun(@isempty, strfind(...
        baselineNames, curBaseline(1:8))),2) & ...
        sum(~cellfun(@isempty, strfind(...
        baselineNames, curBaseline(10:17))),2);
    
% #### Per Source ####
elseif get(handles.radiobutton_plot_residuals_perSource, 'Value') % sourcewise residuals
    allSourcesInPopupmenu=get(handles.popupmenu_plot_residuals_source, 'String');
    curSource=allSourcesInPopupmenu{get(handles.popupmenu_plot_residuals_source, 'Value')};
    sourceNames=handles.data.plot.res(curSession).allSourceNames(handles.data.plot.res(curSession).source); % source names for each observation
    obsWithCurSelection=strcmp(sourceNames, curSource);
    handles.data.plot.currentSources = curSource;
% #### All Values ####    
else
    obsWithCurSelection=ones(size(handles.data.plot.res(curSession).mjd,1),1);
end

% get first or main residuals
if get(handles.radiobutton_plot_residuals_firstSolution, 'Value')
    curVals=handles.data.plot.res(curSession).firstVal(logical(obsWithCurSelection));
else
    curVals=handles.data.plot.res(curSession).mainVal(logical(obsWithCurSelection));
end

% now do the *(-1)
if get(handles.radiobutton_plot_residuals_perStat, 'Value')
    curVals(secondAntObservations)=curVals(secondAntObservations)*-1;
end

% ##### get values which lie in the box #####

% Get X-Axis Values for current selection:
SessionStartTimeMJD =  handles.data.plot.res(curSession).mjd(1);
mjd = handles.data.plot.res(curSession).mjd;
DurationHours = (mjd(logical(obsWithCurSelection)) - SessionStartTimeMJD) * 24; 

% if get(handles.togglebutton_plot_residuals_SelectedData_output,'Value')
%     if ~strcmp(get(handles.edit_plot_residuals_interval_input_1_unit,'String'),'hh') && ...
%         ~strcmp(get(handles.edit_plot_residuals_interval_input_2_unit,'String'),'hh') 
%             
%         x(1) = sscanf(get(handles.edit_plot_residuals_interval_input_1_unit,'String'),'%f');
%         x(2) = sscanf(get(handles.edit_plot_residuals_interval_input_2_unit,'String'),'%f');
%         y(2)=max(curVals);y(1)=min(curVals);
%     else

% get current position
curPoint=get(get(hObject, 'CurrentAxes'), 'CurrentPoint');
curPoint=curPoint(1,1:2);
% draw lines ...
x=[min(curPoint(1), handles.data.plot.outlierBoxStart(1)), ...
    max(curPoint(1), handles.data.plot.outlierBoxStart(1))];
y=[min(curPoint(2), handles.data.plot.outlierBoxStart(2)), ...
    max(curPoint(2), handles.data.plot.outlierBoxStart(2))];

%     end
% end

handles.data.plot.lineHandle=line(...
    [x(1), x(2), x(2), x(1), x(1)], [y(2), y(2), y(1), y(1), y(2)],...
        'color', 'k', 'DisplayName', 'selection');
set(get(get(handles.data.plot.lineHandle,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

handles.data.plot.lineHandle=line(...
    [x(1), x(2), x(2), x(1), x(1)], [-y(2), -y(2), -y(1), -y(1), -y(2)],...
        'color', 'k', 'DisplayName', 'selection');
 set(get(get(handles.data.plot.lineHandle,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


% Check Y-Values:
valsWithinValues = (curVals<=y(2) &  curVals>=y(1)) | (curVals<=-y(1) &  curVals>=-y(2));
% Check X-Values:
temp=DurationHours(valsWithinValues)<=x(2) & DurationHours(valsWithinValues)>=x(1);

valsWithinValues(valsWithinValues==1)=temp;

indToPlot=find(valsWithinValues);
hold(handles.axes_plot_residuals, 'on');

if get(handles.togglebutton_plot_residuals_selectOutliers,'Value')
    % write number of outliers to button where they can be removed
    set(handles.pushbutton_plot_residuals_removeOutliers, 'String','Remove')
    % plotting
	handles.data.plot.outlierMarksHandle=plot(handles.axes_plot_residuals, [DurationHours(indToPlot) DurationHours(indToPlot)], ...
        [curVals(indToPlot) -curVals(indToPlot)], 'x', 'color', 'k', 'markersize', 10, 'DisplayName', 'outlier','LineWidth',3);
	for i = 2:length(handles.data.plot.outlierMarksHandle)
        set(get(get(handles.data.plot.outlierMarksHandle(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
	end
    hold(handles.axes_plot_residuals, 'off');
else
    % if old lines exist - delete them
    allLineHandles=findobj(handles.axes_plot_residuals,'Type','line', 'color', 'k');
    delete(allLineHandles);
    allLineHandles=findobj(handles.axes_plot_residuals,'Type','line', 'color', [0 0 0.04]);
    delete(allLineHandles);
    
	[~, ~, ~, hour, minu, sec] = mjd2date(SessionStartTimeMJD);
    SessionStartTimeUTC = hour + minu/60 + sec/3600;    
    handles.data.plot.foget_unit =0;

    % for non-intensive
%     if get(handles.radiobutton_unit_UTC,'Value') && ...
%             get(handles.togglebutton_plot_residuals_SelectedData_output,'Value')    
%         if DurationHours(end)>5
%             ind=x<SessionStartTimeUTC;
%             if sum(x<SessionStartTimeUTC)>0
%                 x(ind) = x(ind)+24-SessionStartTimeUTC;
%                 x(~ind) = x(~ind) - SessionStartTimeUTC;
%             else
%                 x = x - SessionStartTimeUTC;
%             end
%         else
%                 x = x - SessionStartTimeUTC;
%         end
%     end
    yLimits = get(handles.axes_plot_residuals,'YLim');

    % plotting
    p = plot(handles.axes_plot_residuals, [x(1) x(1)],yLimits, ...
         'color', 'k', 'markersize', 10, 'DisplayName', 'selection');
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

    p = plot(handles.axes_plot_residuals, [x(2) x(2)],yLimits, ...
         'color', 'k', 'markersize', 10, 'DisplayName', 'selection');
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

    hold(handles.axes_plot_residuals, 'off');

end

% Save the change you made to the structure
guidata(hObject,handles)
