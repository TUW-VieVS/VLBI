% #########################################################################
% #     endSelectingOutliers
% #########################################################################
%
% DESCRITPION
% This function consists of the end of the outlier selection, ie. the
% function of mouse motion is stopped, and the values are shown.
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
%
% CHANGES

function endSelectingOutliers(hObject, eventdata, handles)
% This function consists of the end of the outlier selection, ie. the
% function of mouse motion is stopped, and the values are shown.
set(handles.figure_vievs2, 'WindowButtonMotionFcn', '')