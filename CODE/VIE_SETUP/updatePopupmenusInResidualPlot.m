% #########################################################################
% #     updatePopupmenusInResidualPlot
% #########################################################################
%
% DESCRITPION
% This function updates the popupmenus in the residual plot panel after a
% new session is chosen.
%
% AUTHOR 
%   moved in separate function by A. Girdiuk
%
% INPUT
%   handles      structure from the GUI (also containing data, eg residuals)
%   hObject      ...unused so far.
%
% OUTPUT
%   handles      structure from the GUI (also containing data, eg residuals)
%
% CHANGES
% 2016-10-10, A. Hellerschmied: Satellite sources enabled (names are loaded to popup menu)
% 2017-04-03, A. Hellerschmied: Minor bug fix (satellite source name vector transposed) 
% 2017-05-08, A. Hellerschmied: Minor bug fix (satellite source names displayed correctly) 
%
function handles=updatePopupmenusInResidualPlot(hObject, handles)


allsubfoldersInMenu=get(handles.popupmenu_plot_residuals_folder, 'String');

allSessionsInMenu=get(handles.popupmenu_plot_residuals_session, 'String');
curSelectedSessionNr=get(handles.popupmenu_plot_residuals_session, 'Value');
if ~isempty(allSessionsInMenu)
    
    % get res-file data from selected session
    res=handles.data.plot.res(curSelectedSessionNr);

    % update menus
    % (a) station
    set(handles.popupmenu_plot_residuals_station, 'Value', 1)
    set(handles.popupmenu_plot_residuals_station, 'String', res.allStatNames)
    % (b) baselines
    baselines=getAllStringCombinations(res.allStatNames);
    set(handles.popupmenu_plot_residuals_baseline, 'Value', 1)
    set(handles.popupmenu_plot_residuals_baseline, 'String', ...
        strcat(baselines(:,1), repmat('-', length(baselines),1), baselines(:,2)))
    % (c) sources
    % Create list of all sources
    % => First, all quasars, followed by satellites (if available)!
    if ~isfield(res, 'allSatelliteNames')
        res.allSatelliteNames = {};
    end
    set(handles.popupmenu_plot_residuals_source, 'Value', 1)
    if size(res.allSatelliteNames, 2) > 1
        res.allSatelliteNames = res.allSatelliteNames';
    end
    set(handles.popupmenu_plot_residuals_source, 'String', [res.allSourceNames; res.allSatelliteNames])
    % (d) statinos for clock reference
    set(handles.popupmenu_plot_residuals_refClock, 'Value', 1)
    set(handles.popupmenu_plot_residuals_refClock, 'String', res.allStatNames)
    
    % update first/main radiobuttons (if 
    if isempty(res.firstVal)
        set(handles.radiobutton_plot_residuals_firstSolution, 'Enable', 'off')
        set(handles.radiobutton_plot_residuals_mainSolution, 'Value', 1)        % one of the two must have been chosen...
    else
        set(handles.radiobutton_plot_residuals_firstSolution, 'Enable', 'on')
    end
    if isempty(res.mainVal)
        set(handles.radiobutton_plot_residuals_mainSolution, 'Enable', 'Off')
        set(handles.radiobutton_plot_residuals_firstSolution, 'Value', 1)        % one of the two must have been chosen...
    else
        set(handles.radiobutton_plot_residuals_mainSolution, 'Enable', 'On')
    end
end

