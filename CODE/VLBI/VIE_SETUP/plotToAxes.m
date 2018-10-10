% #########################################################################
% #     plotToAxes
% #########################################################################
%
% DESCRITPION
% This function plots the chosen values (parameters) to the main axes in
% the parameter plotting window.
%
% AUTHOR 
%   Matthias Madzak
%
% INPUT
%   handles      structure from the GUI (also containing data, eg residuals)
%   hObject      ...unused so far.
%
% COUPLING
%   mjd2datestr.m
%
% OUTPUT
%   handles      structure from the GUI (also containing data, eg residuals)
%
% CHANGES
%   2014-08-18, A. Hellerschmied: Adaptions to plot source coordinates in parameter plotting window.
%   2016-09-10, A.Girdiuk: plot review
%   2017-05-09, A. Hellerschmied: Satellite sources are considered
%   2017-06-20, A. Hellerschmied: Solved problem with invisible figure used for saving plots in vie_setup
%

function handles=plotToAxes(hObject, handles)

% define parameters whcih are estimated per session (and not per
% stations or source), e.g. dut1, nutdx, nutxy
paramsEstimPerSess={'xpol', 'ypol', 'dut1', 'nutdx', 'nutdy'};

% define parameters whcih are estimated per source (and not per
% stations or session), e.g. soura, soude 
paramsEstimPerSource={'soura', 'soude', 'sat_pos1', 'sat_pos2', 'sat_pos3'};


% create an invisbile figure (for saving) if not done already
if isfield(handles, 'figure_save') && isvalid(handles.figure_save)
else
    handles.figure_save=figure('name', 'Figure for saving', 'visible', 'off');
end

% make axes empty
cla(handles.axes_plot);
cla(get(handles.figure_save, 'currentaxes'));

% get panels to be plotted
plotPanelLogical=zeros(3,1);
for iPanel=1:length(plotPanelLogical)
    if ~isempty(get(eval(['handles.popupmenu_plot_folder', num2str(iPanel), '_param']), 'String'))
        plotPanelLogical(iPanel)=1;
    end
end

% if nothing should be plotted -> return
if sum(plotPanelLogical)==0
    return;
end

% index of panels that should be plotted
panelsToBePlotted=find(plotPanelLogical);

% define colors
myColors={'k', 'r', 'b'};

% Loop over all panels that should be plotted
% => Load data
for iPanelToPlot=1:sum(plotPanelLogical)
   
    % current panel to be plotted
    indPanel=panelsToBePlotted(iPanelToPlot);
    
    % get antenna/parameter to be plotted
    allAntennas=get(eval(['handles.popupmenu_plot_folder', num2str(indPanel), '_stat']), 'String');
    allParameters=get(eval(['handles.popupmenu_plot_folder', num2str(indPanel), '_param']), 'String');
    curAnt=allAntennas{get(eval(['handles.popupmenu_plot_folder', num2str(indPanel), '_stat']), 'Value')};
    curParam=allParameters{get(eval(['handles.popupmenu_plot_folder', num2str(indPanel), '_param']), 'Value')};
    
    % if paremeter = "soude" or "soura" => Get current source name!
    if strcmp(curParam, 'soude') || strcmp(curParam, 'soura') || strcmp(curParam, 'sat_pos1') || strcmp(curParam, 'sat_pos2') || strcmp(curParam, 'sat_pos3')
        allSources=get(eval(['handles.popupmenu_plot_folder', num2str(indPanel), '_sources']), 'String');
        curSrc=allSources{get(eval(['handles.popupmenu_plot_folder', num2str(indPanel), '_sources']), 'Value')};
    end
    
    % get number of sessions
    if get(eval(['handles.radiobutton_plot_folder', num2str(indPanel), '_oneSess']), 'Value')
        nSessions=1;
    else
        nSessions=size(handles.plot.x_files(indPanel).x_,2);
    end

    % get data to be plotted
    % preallocate
    dummyVal=99999;
    tmpVal=ones(10000, nSessions)*dummyVal;
    tmpMjd=tmpVal;
    tmpMx=tmpVal;
    
    % for all sessions -> get data
    for iSession=1:nSessions
        % make empty
        valsOfSession=[];
        mjdsOfSession=[];
        mxsOfSession=[];
        
        % provide index of session to be plotted (needed when just one
        % session should be plotted
        if get(eval(['handles.radiobutton_plot_folder', num2str(indPanel), '_oneSess']), 'Value')
           %allSessionsInPopupmenu=get(eval(['handles.popupmenu_plot_folder', num2str(indPanel), '_sessions']), 'String');
            curSessionInd=get(eval(['handles.popupmenu_plot_folder', num2str(indPanel), '_sessions']), 'Value');
        else
            curSessionInd=iSession;
        end
          
        % if parameter was estimated at all ( - it could be empty when just
        % "preallocated" or the field val might not exist(no lsm performed)
        if ~isempty(eval(['handles.plot.x_files(indPanel).x_(curSessionInd).', curParam])) && isfield(eval(['handles.plot.x_files(indPanel).x_(curSessionInd).', curParam]), 'val')
                    
            % if a parameter was selected which is estimated per session (not
            % per station or source (eg dut1, nutdx, nutdy, xpol)
            if ismember(curParam, paramsEstimPerSess)
                
                % Disable popupmenus and ralated text which are not used
                set(eval(['handles.popupmenu_plot_folder', num2str(indPanel), '_stat']), 'Enable', 'off');
                set(eval(['handles.text_plot_folder', num2str(indPanel), '_stat']), 'Enable', 'off');
                set(eval(['handles.popupmenu_plot_folder', num2str(indPanel), '_sources']), 'Enable', 'off');
                set(eval(['handles.text_plot_folder', num2str(indPanel), '_sources']), 'Enable', 'off');
                
                % only take values which are not NaN
                if ~isnan(eval(['handles.plot.x_files(indPanel).x_(curSessionInd).', curParam,'.val(1)']))
                    % get values (mjd, val, mx) of current session
                    valsOfSession=eval(['handles.plot.x_files(indPanel).x_(curSessionInd).', curParam,'.val']);
                    mjdsOfSession=eval(['handles.plot.x_files(indPanel).x_(curSessionInd).', curParam,'.mjd;']);
                    mxsOfSession=eval(['handles.plot.x_files(indPanel).x_(curSessionInd).', curParam,'.mx;']);
                end
                
            % If the selected parameter was estimated per source    
            elseif ismember(curParam, paramsEstimPerSource)
                
                % Disable/enable popupmenus and ralated text which are not used/used
                set(eval(['handles.popupmenu_plot_folder', num2str(indPanel), '_stat']), 'Enable', 'off');
                set(eval(['handles.text_plot_folder', num2str(indPanel), '_stat']), 'Enable', 'off');
                set(eval(['handles.popupmenu_plot_folder', num2str(indPanel), '_sources']), 'Enable', 'on');
                set(eval(['handles.text_plot_folder', num2str(indPanel), '_sources']), 'Enable', 'on');
                
                % get index of chosen source in current session
                indChosenSrcInCurSess = [];
                % - Check quasar sources:
                if strcmp(curParam, 'soude') || strcmp(curParam, 'soura')
                    indChosenSrcInCurSess = find(strcmp({handles.plot.x_files(indPanel).sources(curSessionInd).source.name}, curSrc));
                end
                % - Check satellite sources:
                if strcmp(curParam, 'sat_pos1') || strcmp(curParam, 'sat_pos2') || strcmp(curParam, 'sat_pos3')
                    indChosenSrcInCurSess = find(strcmp({handles.plot.x_files(indPanel).sources(curSessionInd).satellite.name}, curSrc));
                end

                % if we have found the chosen station in current session
                if ~isempty(indChosenSrcInCurSess)
                    % get values (mjd, val, mx) of current session
                    valsOfSession=eval(['handles.plot.x_files(indPanel).x_(curSessionInd).', curParam,'(indChosenSrcInCurSess).val;']);
                    mjdsOfSession=eval(['handles.plot.x_files(indPanel).x_(curSessionInd).', curParam,'(indChosenSrcInCurSess).mjd;']);
                    mxsOfSession=eval(['handles.plot.x_files(indPanel).x_(curSessionInd).', curParam,'(indChosenSrcInCurSess).mx;']);
                end
            % If a parameter was chosen which is estimated per station
            else
                
                % Disable/enable popupmenus and ralated text which are not used/used
                set(eval(['handles.popupmenu_plot_folder', num2str(indPanel), '_stat']), 'Enable', 'on');
                set(eval(['handles.text_plot_folder', num2str(indPanel), '_stat']), 'Enable', 'on');
                set(eval(['handles.popupmenu_plot_folder', num2str(indPanel), '_sources']), 'Enable', 'off');
                set(eval(['handles.text_plot_folder', num2str(indPanel), '_sources']), 'Enable', 'off');                
                
                % get index of chosen station (EG TIGO) in current session
                indChosenStatInCurSess=find(strcmp({handles.plot.x_files(indPanel).x_(curSessionInd).antenna.name}, curAnt));

                % if we have found the chosen station in current session
                if ~isempty(indChosenStatInCurSess)


                    % get values (mjd, val, mx) of current session
                    valsOfSession=eval(['handles.plot.x_files(indPanel).x_(curSessionInd).', curParam,'(indChosenStatInCurSess).val;']);
                    mjdsOfSession=eval(['handles.plot.x_files(indPanel).x_(curSessionInd).', curParam,'(indChosenStatInCurSess).mjd;']);
                    mxsOfSession=eval(['handles.plot.x_files(indPanel).x_(curSessionInd).', curParam,'(indChosenStatInCurSess).mx;']);

                end
            end

            % if the parameter was estimated in the current session
            if ~isempty(valsOfSession)
                % put the current values to the big value-matrix
                tmpVal(1:length(valsOfSession), iSession)=valsOfSession;
                tmpMjd(1:length(mjdsOfSession), iSession)=mjdsOfSession;
                tmpMx (1:length(mxsOfSession ), iSession)=mxsOfSession;
            end

        end
        
    end % for all sessions -> get data
    
    % make one vector out of matrix
    data.vals{indPanel}=tmpVal(:);
    data.mjds{indPanel}=tmpMjd(:);
    data.mxs{indPanel}=tmpMx(:);
    
    % sort chronological
    [data.mjds{indPanel}, sortInd]=sort(data.mjds{indPanel});
    data.vals{indPanel}=data.vals{indPanel}(sortInd);
    data.mxs{indPanel} =data.mxs{indPanel} (sortInd);
    
    % delete dummy entries
    data.vals{indPanel}(data.vals{indPanel}==dummyVal)=[];
    data.mjds{indPanel}(data.mjds{indPanel}==dummyVal)=[];
    data.mxs{indPanel}(data.mxs{indPanel}==dummyVal)=[];
    
end

mjd_min = min(vertcat(data.mjds{:,:})); % time ref.
    
% Loop over all panels that should be plotted
% => Plot data
for iPanelToPlot=1:sum(plotPanelLogical)
    
    % current panel to be plotted
    indPanel=panelsToBePlotted(iPanelToPlot);
    
    % Get time reference
    tmp_cellStr = get(handles.popupmenu_plot_select_time_ref_format, 'String');
%     mjd_min = 
    switch(tmp_cellStr{get(handles.popupmenu_plot_select_time_ref_format, 'Value')})
        case 'MJD'
            t_ref = data.mjds{indPanel};
            xlabel_str = 'MJD';
        case 'Hours since session start'
            if ~isempty(data.mjds{indPanel})
                t_ref = (data.mjds{indPanel} - mjd_min) * 24;
                xlabel_str = sprintf('Hours since session start (%s)', mjd2datestr(mjd_min));
            else
                t_ref = data.mjds{indPanel};
                xlabel_str = sprintf('Hours since session start');
            end
        case 'Minutes since session start'
            if ~isempty(data.mjds{indPanel})
                t_ref = (data.mjds{indPanel} - mjd_min) * 24 * 60;
                xlabel_str = sprintf('Minutes since session start (%s)', mjd2datestr(mjd_min));
            else
                t_ref = data.mjds{indPanel};
                xlabel_str = sprintf('Minutes since session start');
            end
            
        case 'Seconds since session start'
            if ~isempty(data.mjds{indPanel})
                t_ref = (data.mjds{indPanel} - mjd_min) * 24 * 60 * 60;
                xlabel_str = sprintf('Seconds since session start (%s)', mjd2datestr(mjd_min));
            else
                t_ref = data.mjds{indPanel};
                xlabel_str = sprintf('Seconds since session start');
            end
    end


    
    % plot
    if get(handles.checkbox_plot_show_errorbars, 'Value')
        handles.plotting_parameter_data{indPanel}                       = errorbar(handles.axes_plot, t_ref, data.vals{indPanel}, data.mxs{indPanel}, ['-', myColors{indPanel}], 'marker','*');
        handles.plotting_parameter_xlabel{indPanel}                     = xlabel(xlabel_str);
        handles.plotting_parameter_current_figure_save{indPanel}        = errorbar(get(handles.figure_save, 'CurrentAxes'), t_ref, data.vals{indPanel}, data.mxs{indPanel}, ['-', myColors{indPanel}], 'marker','*');
        handles.plotting_parameter_current_figure_save_xlabel{indPanel} = xlabel(xlabel_str);
    else
        handles.plotting_parameter_data{indPanel}                       = plot(handles.axes_plot, t_ref, data.vals{indPanel}, ['-', myColors{indPanel}], 'marker','*');
        handles.plotting_parameter_xlabel{indPanel}                     = xlabel(xlabel_str);
        handles.plotting_parameter_current_figure_save{indPanel}        = plot(get(handles.figure_save, 'CurrentAxes'), t_ref, data.vals{indPanel}, ['-', myColors{indPanel}], 'marker','*');
        handles.plotting_parameter_current_figure_save_xlabel{indPanel} = xlabel(xlabel_str);
    end
    hold(handles.axes_plot, 'on');
    hold(get(handles.figure_save, 'CurrentAxes'), 'on');
    
end

% hold axes off
hold(handles.axes_plot, 'off');
hold(get(handles.figure_save, 'CurrentAxes'), 'off');

% save plotted data to handles struct
handles.plot.data=data;