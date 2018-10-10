% #########################################################################
% #     setAllPanelsToInvisible
% #########################################################################
%
% DESCRITPION
% This function set panels invisible
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
function setAllPanelsToInvisible(hObject, handles)

% set bar for plotting options to visible
set(handles.uitoolbar_plotOptions, 'Visible', 'Off');

% add to y-position: +27, is changed when plotting-tab is chosen
if (~strcmp(get(hObject, 'Tag'), handles.menu_plotting_parameters) && ...
        strcmp(get(handles.uipanel_plot_parameters, 'Visible'), 'on')) ||...
        (~strcmp(get(hObject, 'Tag'), handles.menu_plotting_residuals) && ...
        strcmp(get(handles.uipanel_plot_residuals, 'Visible'), 'on'))
    
    curPosition=get(handles.figure_vievs2, 'Position');
    curPosition(2)=curPosition(2)+27;
    set(handles.figure_vievs2, 'Position', curPosition)
end


for iPanel=1:length(handles.allMainUiPanels)
    set(handles.allMainUiPanels(iPanel), 'Visible', 'Off');
end

% Update handles structure
guidata(hObject, handles);
