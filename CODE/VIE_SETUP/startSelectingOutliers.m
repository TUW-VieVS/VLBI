% #########################################################################
% #     startSelectingOutliers
% #########################################################################
%
% DESCRITPION
% This function is called when the mouse is moved (while clicking) over the
% axes. togglebutton_plot_residuals_selectOutliers_Callback must be
% selected!
%
% AUTHOR 
%   Matthias Madzak
%
% INPUT
%   handles     structure from the GUI (also containing data, eg residuals)
%   hObject     handle to checkbox_plot_residuals_showStatNumbers (see GCBO)
%   eventdata   reserved - to be defined in a future version of MATLAB
%
% OUTPUT
%   handles     structure from the GUI (also containing data, eg residuals)
%
% CHANGES
%   2016-09-10, A. Girdiuk: function was suited for 'Select Data' panel usage
%

function startSelectingOutliers(hObject, eventdata, handles)

% if isfield(handles.data.plot, 'outlierMarksHandle')
%     get(handles.data.plot.outlierMarksHandle)
% end

% save "start" point of box to handles structure
startPointBox=get(get(handles.figure_vievs2, 'CurrentAxes'), 'CurrentPoint');
%handles.data.plot.startPointBox=startPointBox;
startPointBox=startPointBox(1,1:2);
handles.data.plot.outlierBoxStart=startPointBox;

% if ~get(handles.togglebutton_plot_residuals_selectOutliers,'Value')
%     set(handles.radiobutton_unit_plot, 'Enable', 'On')
%     set(handles.radiobutton_unit_UTC, 'Enable', 'On')
%     set(handles.radiobutton_unit_MJD, 'Enable', 'On')
% end

% ##### IF MOUSE IS MOVED WHILE CLICKING - BOX! #####
handles.data.plot.foget_unit =0;
handles.data.plot.intervalSelectedatPlot=0;

set(handles.figure_vievs2, 'WindowButtonMotionFcn', {@outlierSelection,handles})


% ##### OTHERWiSE - SINGLE OUTLIER SELECTION: #####

% if only one point is chosen - select (black cross) this one
% find closest point
curXLim=get(handles.axes_plot_residuals, 'XLim'); % needed for "scaling" for distances
curYLim=get(handles.axes_plot_residuals, 'YLim');
xExtent=curXLim(2)-curXLim(1);
yExtent=curYLim(2)-curYLim(1);
axesPos=get(handles.axes_plot_residuals, 'Position');
xScaling=(xExtent/axesPos(4))/(yExtent/axesPos(3));
% get plotted values
plottedX=get(handles.data.plot.plottedResiduals, 'XData');
plottedY=get(handles.data.plot.plottedResiduals, 'YData');

% calc distances (scaled for not axis equal)
distances=sqrt((plottedX-startPointBox(1)).^2.+...
    (plottedY-startPointBox(2)).^2.*xScaling);

% plot
% clear black (or nearly black) values
allLineHandles=findobj(handles.axes_plot_residuals,'Type','line', 'color', 'k');
delete(allLineHandles);
allLineHandles=findobj(handles.axes_plot_residuals,'Type','line', 'color', [0 0 0.04]);
delete(allLineHandles);

hold(handles.axes_plot_residuals, 'on');
handles.data.plot.outlierMarksHandle=plot(handles.axes_plot_residuals,...
    [plottedX(abs(distances)==min(abs(distances))) plottedX(abs(distances)==min(abs(distances)))],...
    [plottedY(abs(distances)==min(abs(distances))) -plottedY(abs(distances)==min(abs(distances)))], 'x',...
    'color', [0 0 0.04], 'markersize', 10, 'LineWidth', 3, 'DisplayName', 'outlier');
for i = 2:length(handles.data.plot.outlierMarksHandle)
    set(get(get(handles.data.plot.outlierMarksHandle(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

hold(handles.axes_plot_residuals, 'off');

if get(handles.togglebutton_plot_residuals_selectOutliers,'Value')
	
    % write number of outliers to button where they can be removed
%     set(handles.pushbutton_plot_residuals_removeOutliers, 'String', ...
%         sprintf('Remove %1.0f Outliers', 1))

else
         %
        curSession=get(handles.popupmenu_plot_residuals_session, 'Value');
        SessionStartTimeMJD =  handles.data.plot.res(curSession).mjd(1);
        [year, month, day, hour, minu, sec] = mjd2date(SessionStartTimeMJD);
        %
        SessionTimeUTC = hour + minu/60 + sec/3600 + startPointBox(1);
        if sum(SessionTimeUTC>24)>0
            SessionTimeUTC(SessionTimeUTC>24)=SessionTimeUTC(SessionTimeUTC>24)-24;
        end
        %
        handles.data.plot.unit_MJD  = startPointBox(1)/24+SessionStartTimeMJD;
        handles.data.plot.unit_UTC  = SessionTimeUTC;
        handles.data.plot.unit_plot = startPointBox(1);
        
    % get values which are currently plotted
    % get index of currently selected session

% 	if get(handles.radiobutton_unit_MJD,'Value')
%         set(handles.edit_plot_residuals_interval_show, 'String', ...
%             sprintf(num2str(startPointBox(1)/24+SessionStartTimeMJD)))
%     elseif get(handles.radiobutton_unit_UTC,'Value')
%          set(handles.edit_plot_residuals_interval_show, 'String', ...
%             sprintf(num2str(SessionTimeUTC)))
% 	elseif get(handles.radiobutton_unit_plot,'Value')
%          set(handles.edit_plot_residuals_interval_show, 'String', ...
%             sprintf(num2str(startPointBox(1))))
%     end
end

% Save the change you made to the structure
guidata(hObject,handles)