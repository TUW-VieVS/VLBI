% #########################################################################
% #     update_antennaParameterPopupmenus_plotParameters
% #########################################################################
%
% DESCRITPION
% This function updates the popupmenus "parameter", "sources" and "station" according
% to the chosen options. The selected sessions only consists of several
% parameters and stations and only those should be visible.
%
% AUTHOR 
%   Matthias Madzak
%
% INPUT
%   handles      structure from the GUI (also containing data, eg residuals)
%   chosenPanel  Number of the chosen panel/folder (1, 2, or 3)
%
% OUTPUT
%   handles     structure from the GUI (also containing data, eg residuals)
%
% CHANGES
%   2014-08-18, A. Hellerschmied: Adaptions to plot source coordinates in parameter plotting window.
%   2016-09-10, A.Girdiuk: antenna name is written shorter
%   2016-10-11, A. Hellerschmied: Check, if "handles.plot.x_files(chosenPanel).sources(iFile).source" is empty. Required, if only satellites were observed in a session. 
%   2017-05-09, A. Hellerschmied: Satellite sources are considered
%   2017-06-30, A. Hellerschmied: Minor bug fixed (source names)
%

function handles=update_antennaParameterPopupmenus_plotParameters(handles,chosenPanel)

% get current selection of antenna, source and parameter popupmenus (so that the
% selection can be the same after update)
allParamBeforeUpdate=get(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_param']), 'String');
if isempty(allParamBeforeUpdate) || max(size(allParamBeforeUpdate))==1 % if the popupmenu is either empty or just a char ' '
    selectedParamBeforeUpdate='';
else
    selectedParamBeforeUpdate=allParamBeforeUpdate{get(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_param']), 'Value')};
end
allAntBeforeUpdate=get(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_stat']), 'String');
if isempty(allAntBeforeUpdate) || max(size(allAntBeforeUpdate))==1
    selectedAntBeforeUpdate='';
else
    selectedAntBeforeUpdate=allAntBeforeUpdate{get(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_stat']), 'Value')};
end
allSrcBeforeUpdate=get(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_sources']), 'String');
if isempty(allSrcBeforeUpdate) || max(size(allSrcBeforeUpdate))==1
    selectedSrcBeforeUpdate='';
else
    selectedSrcBeforeUpdate=allSrcBeforeUpdate{get(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_sources']), 'Value')};
end


% reset popupmenus
set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_param']), 'Value', 1)
set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_stat']), 'Value', 1)
set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel),'_sources']), 'Value', 1)
set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_param']), 'String', {''})
set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_stat']), 'String', {''})
set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel),'_sources']), 'String', ' ')
set(eval(['handles.text_plot_folder', num2str(chosenPanel),'_unit']), 'String', '')

nSess=size(handles.plot.x_files(chosenPanel).x_,2);

% get "selected" session(s) - one session, or all sessions, dependent on
% the radiobutton "handles.radiobutton_plot_folder_oneSess"
if eval(['get(handles.radiobutton_plot_folder', num2str(chosenPanel), '_oneSess, ''Value'')'])
    indOfSelSessions=eval(['get(handles.popupmenu_plot_folder', num2str(chosenPanel), '_sessions, ''Value'')']);
else
    indOfSelSessions=1:nSess; % get the vector from one to nSessions
end

% ##### (1) update parameter popupmenu #####
% update parameters to be chosen
%fieldsToBeRemoved={'antenna', 'qclk', 'rclk', 'rqclk', 'soude', 'soura'};
fieldsToBeRemoved={'antenna', 'qclk', 'rclk', 'rqclk'};

% get estimated fields (also not estimated parames do have a field)
estimatedFields=getEstimatedFields(handles, indOfSelSessions, chosenPanel);

% get "durchschnitt".. when a field of those does not exist -> it can't be
% deleted!
fieldsThatCanBeRemoved=intersect(fieldnames(handles.plot.x_files(chosenPanel).x_(1)), fieldsToBeRemoved);
set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_param']), 'String', estimatedFields(~ismember(estimatedFields,fieldsThatCanBeRemoved)));

% ##### (2) update antenna popupmenu #####
if eval(['get(handles.radiobutton_plot_folder', num2str(chosenPanel), '_oneSess, ''Value'')'])
    allAntNames={handles.plot.x_files(chosenPanel).x_(indOfSelSessions).antenna.name};
else
    % all sessions are selected
    allAntNames={''};
    for iFile=1 : nSess
        allAntNames=unique([allAntNames, {handles.plot.x_files(chosenPanel).x_(iFile).antenna.name}]);
    end
end
set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_stat']), 'String', allAntNames(~strcmp(allAntNames, '')));


% ##### (3) update sources popupmenu #####

% Check if source coordinates have been estimated:
FlagSrcEstimated = 0;
for i_fiel = 1: length(estimatedFields)
    if (strcmp(estimatedFields(i_fiel), 'soude')) || (strcmp(estimatedFields(i_fiel), 'soura')) || (strcmp(estimatedFields(i_fiel), 'sat_pos1')) || (strcmp(estimatedFields(i_fiel), 'sat_pos2')) || (strcmp(estimatedFields(i_fiel), 'sat_pos3'))
       FlagSrcEstimated = 1;
       break;
    end
end

if FlagSrcEstimated
    % If only one session is chosen:
    if eval(['get(handles.radiobutton_plot_folder', num2str(chosenPanel), '_oneSess, ''Value'')'])
        if ~isempty(handles.plot.x_files(chosenPanel).sources(indOfSelSessions).source)
            allQuasarNames  = {handles.plot.x_files(chosenPanel).sources(indOfSelSessions).source.name};
        else
            allQuasarNames  = {''};
        end
        if ~isempty(handles.plot.x_files(chosenPanel).sources(indOfSelSessions).satellite)
            allSatNames     = {handles.plot.x_files(chosenPanel).sources(indOfSelSessions).satellite.name};
        else
            allSatNames     = {''};
        end
        allSrcNames     = [allQuasarNames, allSatNames];
    
    % If all sessions in subfolder are chosen:    
    else
        allQuasarNames  = {''};
        allSatNames     = {''};
        % Loop over all sessions in list:
        for iFile = 1 : nSess
            if ~isempty(handles.plot.x_files(chosenPanel).sources(iFile).source)
                allQuasarNames  = unique([allQuasarNames, {handles.plot.x_files(chosenPanel).sources(iFile).source.name}]);
            end
            if ~isempty(handles.plot.x_files(chosenPanel).sources(iFile).satellite)
                allSatNames     = unique([allSatNames, {handles.plot.x_files(chosenPanel).sources(iFile).satellite.name}]);
            end
        end
        allSrcNames     = [allQuasarNames, allSatNames];
    end
else
    allSrcNames={' '};
    set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_sources']), 'Enable', 'off');
end

set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_sources']), 'String', allSrcNames(~strcmp(allSrcNames, '')));


% ##### if the popupmenus include the ant/param/source which was selected before -> take that one again (convenience)#####
previousParamInNewMenu=ismember(estimatedFields(~ismember(estimatedFields,fieldsThatCanBeRemoved)), selectedParamBeforeUpdate);
previousAntInNewMenu=ismember(allAntNames(~strcmp(allAntNames, '')), selectedAntBeforeUpdate);
previousSrcInNewMenu=ismember(allSrcNames(~strcmp(allSrcNames, '')), selectedSrcBeforeUpdate);

if sum(previousParamInNewMenu)>0
    set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_param']), 'Value', find(previousParamInNewMenu))
end
if sum(previousAntInNewMenu)>0
    set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_stat']), 'Value', find(previousAntInNewMenu))
end
if sum(previousSrcInNewMenu)>0
    set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_sources']), 'Value', find(previousSrcInNewMenu))
end

% write unit to text (as info for the user)
unitPerParam={'pwclk', 'cm'; 'rclk', 'cm/day'; 'qclk', 'cm/(day^2)'; ...
    'rqclk', 'cm/(day^2) warning: rate+quadr'; 'zwd', 'cm'; 'ngr', 'cm'; 'egr', 'cm';...
    'xpol', 'mas'; 'ypol', 'mas'; 'dut1', 'ms'; 'nutdx', 'mas'; 'nutdy',...
    'mas'; 'coorx', 'cm'; 'coory', 'cm'; 'coorz', 'cm'; 'soude', 'mas'; 'soura', 'mas';...
    'sat_pos1', 'cm'; 'sat_pos2', 'cm'; 'sat_pos3', 'cm'};

% get selected parameter
allParamInPopupmenu=get(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_param']), 'String');
curSelParam=allParamInPopupmenu{get(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_param']), 'Value')};
strForText=unitPerParam{strcmp(unitPerParam(:,1), curSelParam),2};
if isempty(strForText)
    set(eval(['handles.text_plot_folder', num2str(chosenPanel), '_unit']), 'String', '')
else
    set(eval(['handles.text_plot_folder', num2str(chosenPanel), '_unit']), 'String', ['[', strForText, ']'])
end




