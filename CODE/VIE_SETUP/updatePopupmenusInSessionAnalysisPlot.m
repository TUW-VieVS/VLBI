% #########################################################################
% #     updatePopupmenusInSessionAnalysisPlot
% #########################################################################
%
% DESCRITPION
% This function updates the popupemenues in the session analysis plot
% panel after a new sessino is selected.
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
%   2016-09-21, A. Girdiuk: bug-fix, estimated parameters can be chosen only  to plot correlation matrix
%   2017-06-14, A. Hellerschmied: updated for plotting sat. positions estimates in correlation matrix
%
function handles=updatePopupmenusInSessionAnalysisPlot(hObject, handles)


% if there is at least one session available
if ~isempty(get(handles.popupmenu_plot_sessionAnalysis_session, 'String'))
    % get previously selected popupmenu entries (to get this one again
    % later if exists...)
    allFirstPar=get(handles.popupmenu_plot_sessionAnalysis_corrMat_firstPar, 'String');
    allLastPar=get(handles.popupmenu_plot_sessionAnalysis_corrMat_secondPar, 'String');
    if strcmp(allFirstPar,' ')
        preSelStrings{1}=' ';
    else
        preSelStrings{1}=allFirstPar{get(handles.popupmenu_plot_sessionAnalysis_corrMat_firstPar, 'Value')};
    end
    if strcmp(allLastPar, ' ')
        preSelStrings{2}=' ';
    else
        preSelStrings{2}=allLastPar{get(handles.popupmenu_plot_sessionAnalysis_corrMat_secondPar, 'Value')};
    end
    
    curSelSessionIndex=get(handles.popupmenu_plot_sessionAnalysis_session, 'Value');
    
       
    possibleParameters={'pwclk'; 'rqclk'; 'zwd'; 'ngr'; 'egr'; 'xpol'; 'ypol'; 'dut1'; 'nutdx'; 'nutdy'; 'soura'; 'soude'; 'coorx'; 'coory'; 'coorz'; 'sat_pos1'; 'sat_pos2'; 'sat_pos3'; 'scale'; 'bdclko'}; 
       
    % get fields of current session (only those which were estimated)
    fieldsOfCurSes=cell(length(possibleParameters),1);
    for iPossiblePar=1:length(possibleParameters)
        if isfield(handles.data.plot.sessionAnalysis.x_files(1).x_(curSelSessionIndex), possibleParameters(iPossiblePar))
            if isfield(handles.data.plot.sessionAnalysis.x_files(1).x_(curSelSessionIndex).(possibleParameters{iPossiblePar}), 'val')
                partials = {handles.data.plot.sessionAnalysis.x_files(1).x_(curSelSessionIndex).(possibleParameters{iPossiblePar}).val};
                ind_partials=find(~cellfun(@isempty,partials));
                if ~isempty(ind_partials)
                    fieldsOfCurSes{iPossiblePar}=possibleParameters{iPossiblePar};
                end
            end
        end
    end
      
    % delete cell entries of parameters which were not estimated
    fieldsOfCurSes(cellfun(@isempty, fieldsOfCurSes))=[];
       %     fieldnames(handles.data.plot.sessionAnalysis.x_(curSelSessionIndex));
    
    possibleFieldsOfCurSes=fieldsOfCurSes(...
        ismember(fieldsOfCurSes, possibleParameters));
       
    % update popupmenus
    set(handles.popupmenu_plot_sessionAnalysis_corrMat_firstPar, 'Value', 1)
    set(handles.popupmenu_plot_sessionAnalysis_corrMat_secondPar, 'Value', 1)
    set(handles.popupmenu_plot_sessionAnalysis_corrMat_firstPar, 'String', possibleFieldsOfCurSes)
    set(handles.popupmenu_plot_sessionAnalysis_corrMat_secondPar, 'String', possibleFieldsOfCurSes)
    
    % if previously selected fields exist - use them, otherwise first/last
    if sum(strcmp(possibleFieldsOfCurSes, preSelStrings{1}))>0
        set(handles.popupmenu_plot_sessionAnalysis_corrMat_firstPar, 'Value', find(strcmp(possibleFieldsOfCurSes, preSelStrings{1})));
    end
    if sum(strcmp(possibleFieldsOfCurSes, preSelStrings{2}))>0
        set(handles.popupmenu_plot_sessionAnalysis_corrMat_secondPar, 'Value', find(strcmp(possibleFieldsOfCurSes, preSelStrings{2})));
    else
        set(handles.popupmenu_plot_sessionAnalysis_corrMat_secondPar, 'Value', size(possibleFieldsOfCurSes,1))
    end
end

