% #########################################################################
% #     loadResFiles
% #########################################################################
%
% This function loads all res_ files from DATA/LEVEL3 and updates the
% popupmenu in the interface
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
%
%   2016-09-10, A. Girdiuk: function was suited for 'Select Data' panel usage
%   2016-10-10, A. Hellerschmied: Enabled satellite sources ("res.obs_type" and "res.allSatelliteNames" loaded)
%

function handles=loadResFiles(hObject, handles)


% make everything empty/first values in menu items +clear those
% make axes empty
cla(handles.axes_plot_residuals);
set(handles.popupmenu_plot_residuals_session, 'Value', 1)
set(handles.popupmenu_plot_residuals_station, 'Value', 1)
set(handles.popupmenu_plot_residuals_baseline, 'Value', 1)
set(handles.popupmenu_plot_residuals_source, 'Value', 1)
set(handles.popupmenu_plot_residuals_refClock, 'Value', 1)
set(handles.popupmenu_plot_residuals_session, 'String', ' ')
set(handles.popupmenu_plot_residuals_station, 'String', ' ')
set(handles.popupmenu_plot_residuals_baseline, 'String',' ')
set(handles.popupmenu_plot_residuals_source, 'String', ' ')
set(handles.popupmenu_plot_residuals_refClock, 'String', ' ')



% make everything disable at first
set(handles.radiobutton_plot_residuals_firstSolution, 'Enable', 'off')
set(handles.radiobutton_plot_residuals_mainSolution, 'Enable', 'off')
set(handles.radiobutton_plot_residuals_perStat, 'Enable', 'off')
set(handles.popupmenu_plot_residuals_station, 'Enable', 'off')
set(handles.radiobutton_plot_residuals_perBasel, 'Enable', 'off')
set(handles.popupmenu_plot_residuals_baseline, 'Enable', 'off')
set(handles.radiobutton_plot_residuals_perSource, 'Enable', 'off')
set(handles.popupmenu_plot_residuals_source, 'Enable', 'off')
set(handles.radiobutton_plot_residuals_perAll, 'Enable', 'off')
% set(handles.checkbox_plot_residuals_showStatNumbers, 'Enable', 'off')
set(handles.checkbox_plot_residuals_showSourceNames, 'Enable', 'off')
set(handles.popupmenu_plot_residuals_refClock, 'Enable', 'off')
set(handles.pushbutton_plot_residuals_clockRef_get, 'Enable', 'off')
set(handles.pushbutton_plot_residuals_clockRef_set, 'Enable', 'off')
set(handles.pushbutton_plot_residuals_clockBreak, 'Enable', 'off')
set(handles.togglebutton_plot_residuals_selectOutliers, 'Enable', 'off')
set(handles.pushbutton_plot_residuals_removeOutliers, 'Enable', 'off')
set(handles.togglebutton_plot_residuals_selectData, 'Enable', 'off')
% set(handles.togglebutton_plot_residuals_SelectedData_output, 'Enable', 'off')

% set(handles.radiobutton_unit_plot, 'Enable', 'Off')
% set(handles.radiobutton_unit_UTC, 'Enable', 'Off')
% set(handles.radiobutton_unit_MJD, 'Enable', 'Off')

% get all res_ files in folder
allSubfolders=get(handles.popupmenu_plot_residuals_folder, 'string');
curSubfolder=allSubfolders{get(handles.popupmenu_plot_residuals_folder, 'Value')};

% get all res_ files in folder
allResFiles=dir(['../DATA/LEVEL3/', curSubfolder, '/res_*.mat']);

if ~isempty(allResFiles)
    
    % enable parts of the interface
    set(handles.popupmenu_plot_residuals_session, 'Enable', 'On')
    set(handles.radiobutton_plot_residuals_firstSolution, 'Enable', 'On')
    set(handles.radiobutton_plot_residuals_mainSolution, 'Enable', 'On')
    set(handles.radiobutton_plot_residuals_perStat, 'Enable', 'On')
    set(handles.radiobutton_plot_residuals_perBasel, 'Enable', 'On')
    set(handles.radiobutton_plot_residuals_perSource, 'Enable', 'On')
    set(handles.radiobutton_plot_residuals_perAll, 'Enable', 'On')
%     set(handles.checkbox_plot_residuals_showStatNumbers, 'Enable', 'On')
    set(handles.checkbox_plot_residuals_showSourceNames, 'Enable', 'On')
    set(handles.popupmenu_plot_residuals_refClock, 'Enable', 'On')
    set(handles.pushbutton_plot_residuals_clockRef_get, 'Enable', 'On')
    set(handles.pushbutton_plot_residuals_clockRef_set, 'Enable', 'On')
    set(handles.pushbutton_plot_residuals_clockBreak, 'Enable', 'On')
    set(handles.togglebutton_plot_residuals_selectOutliers, 'Enable', 'On')
    set(handles.togglebutton_plot_residuals_selectData, 'Enable', 'On')
    
    % create short names for popupmenu
    resNamesShort=strrep({allResFiles.name}, 'res_', '');
    resNamesShort=strrep(resNamesShort, '.mat', '');
    
    % create waitbar
    h_waitbar = waitbar(0,'Please wait...');
    
    % for all res files
    nResFiles=size(allResFiles,1);
    handles.data.plot.res=struct('allStatNames',[], 'allSourceNames', [],...
        'baselineOfObs', [],...
        'source', [],...
        'mjd', [],...
        'mainVal', [],...
        'outlier', [],...
        'allSatelliteNames', [],...
        'obs_type', [],...
        'firstVal', [],...
        'mo',[],...
        'mo_first',[],...
        'wrms',[],...
        'sigma_from_fringe_fitting',[],...
        'sigma_residuals_aposteriori',[]);
    
    for iFile=1:nResFiles


        % load res file
        load(['../DATA/LEVEL3/', curSubfolder, '/', allResFiles(iFile).name]);

        % add fields which do ot exist
        if ~isfield(res, 'firstVal')
            res.firstVal=[];
        end
        if ~isfield(res, 'mainVal')
            res.mainVal=[];
        end
        if ~isfield(res, 'outlier')
            res.outlier=[];
        end
        if ~isfield(res, 'allSatelliteNames')
            res.allSatelliteNames={};
        end
        if ~isfield(res, 'obs_type')
            res.obs_type = repmat('q', length(res.mjd), 1);
        end
        if ~isfield(res, 'sigma_from_fringe_fitting')
            res.sigma_from_fringe_fitting=[];
        end
        if ~isfield(res, 'sigma_residuals_aposteriori')
            res.sigma_residuals_aposteriori=[];
        end
        fields=fieldnames(handles.data.plot.res);
        for iF=1:length(fields)
            if isfield (res, fields{iF})
                handles.data.plot.res(iFile).(fields{iF})=res.(fields{iF}); % this manual adding of fields is required when the order of the fields is differnt (Matthias)
            else
                handles.data.plot.res(iFile).(fields{iF})=[]; % this manual adding of fields is required when the order of the fields is differnt (Matthias)
            end
        end
        % update waitbar
        waitbar(iFile/nResFiles,h_waitbar,sprintf('res_ file %1.0f / %1.0f is loaded.', iFile, nResFiles))
    end % for allres files
    
    % close waitbar
    close(h_waitbar);

    % update popupmenus
    % (a) sessions
    set(handles.popupmenu_plot_residuals_session, 'String', resNamesShort);

    % update popupmenus in panel
    handles=updatePopupmenusInResidualPlot(hObject, handles);
    
    % plot right away
    handles=plotResidualsToAxes(handles);
end

% Update handles structure
guidata(hObject, handles);
