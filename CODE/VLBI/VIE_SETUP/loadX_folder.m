% #########################################################################
% #     loadX_folder
% #########################################################################
%
% DESCRITPION
% This function is used to load x_ files from the subfolder given in the
% plotting panel
%
% AUTHOR 
%   Matthias Madzak
%
% INPUT
%   handles     structure from the GUI (also containing data, eg residuals)
%   hObject
%
% OUTPUT
%   handles     structure from the GUI (also containing data, eg residuals)
%
% CHANGES
%   2014-08-18, A. Hellerschmied: Adaptions to plot source coordinates in
%       parameter plotting window.
%	2016-08-31, A. Girdiuk: code optimized -> working time improved
%   2016-09-10, A. Girdiuk: antenna names are taken from x_ file -> no need to load them
%                           and there is no need to check data in x_
%                           because all of fields are written and exist even 
%                           they are not estimated you just do not have data
%   2016-09-21, A. Girdiuk: bug-fix
%   2016-10-11, A. Hellerschmied: Fields (soude, soura) added, if missing in loaded x_ file.
%   2017-05-09, A. Hellerschmied: satellite sources are loaded to "handles.plot.x_files(chosenPanel).sources(iFile).satellite"
%   2017-06-20, A. Hellerschmied: Solved problem with invisible figure used for saving plots in vie_setup
%

function handles=loadX_folder(hObject, handles)

% get panel (1|2|3) where the load button was pushed
if strcmp(get(hObject, 'Tag'), get(handles.pushbutton_plot_folder1_load, 'Tag'))
    chosenPanel=1;
elseif strcmp(get(hObject, 'Tag'), get(handles.pushbutton_plot_folder2_load, 'Tag'))
    chosenPanel=2;
elseif strcmp(get(hObject, 'Tag'), get(handles.pushbutton_plot_folder3_load, 'Tag'))
    chosenPanel=3;
end

% clear figure
cla(handles.axes_plot)
if isfield(handles, 'figure_save') && isvalid(handles.figure_save)
	cla(get(handles.figure_save, 'CurrentAxes'))
end

% delete all popupmenus (+text showing the unit)
set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel),'_sessions']), 'Value', 1)
set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel),'_param']), 'Value', 1)
set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel),'_stat']), 'Value', 1)
set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel),'_sources']), 'Value', 1)
set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel),'_sessions']), 'String', ' ')
set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel),'_param']), 'String', ' ')
set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel),'_stat']), 'String', ' ')
set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel),'_sources']), 'String', ' ')
set(eval(['handles.text_plot_folder', num2str(chosenPanel),'_unit']), 'String', '')

% get chosen subfolder
allSubfolders   = get(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel),'_subfolder']), 'string');
curSubfolder    = allSubfolders{get(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel),'_subfolder']), 'Value')};

% get all x_ files in folder
allX_Files      = dir(['../DATA/LEVEL3/', curSubfolder, '/x_*.mat']);

% create short names for popupmenu
x_NamesShort    = strrep({allX_Files.name}, 'x_', '');
x_NamesShort    = strrep(x_NamesShort, '.mat', '');

% some files (eg simulated) don't have antenna files
antennaFileNotExist = zeros(size(allX_Files,1),1);

% get number of x_ files
nX_Files        = size(allX_Files,1);

% create waitbar
h_waitbar = waitbar(0,'Please wait while x_ files are loaded');

% remove previous data if it exists 
if isfield(handles,'plot') && length(handles.plot.x_files) >= chosenPanel
    handles.plot.x_files(chosenPanel).x_(1:end)         = [];
    handles.plot.x_files(chosenPanel).sources(1:end)    = [];
end
% for all x_ files/sessions in choosen subfolder
for iFile = 1 : nX_Files
    
    % Load parameter file (opt_)
    load(['../DATA/LEVEL3/', curSubfolder, '/', strrep(allX_Files(iFile).name, 'x_', 'opt_')]);
    
    % get LEVEL1 subfolder (for antenne file) from parameter file in LEVEL3
    if ~isfield(opt_, 'level1OutDir')
        curLevel1SubDir = curSubfolder;
    else
        curLevel1SubDir = opt_.level1OutDir;
    end    
    
    
    
    % define antenna file
    curAntFile = ['../DATA/LEVEL1/', curLevel1SubDir, '/', strrep(strrep(allX_Files(iFile).name, 'x_', ''), '.mat', '_antenna.mat')];
    
    if ~exist(curAntFile, 'file')
        antennaFileNotExist(iFile) = 1;
    else
        % load x_ file
        load(['../DATA/LEVEL3/', curSubfolder, '/', allX_Files(iFile).name]);
        
        % Check, if all fields in x_ are available.
        % => Init. empty fields, id required:
        if ~isfield(x_, 'soura')
            x_.soura = [];
        end
        if ~isfield(x_, 'soude')
            x_.soude = [];
        end
        
        % Save current x_ file to "handles.plot.x_files." structure:
        if iFile > 1
        	x_ = orderfields(x_, handles.plot.x_files(chosenPanel).x_(iFile-1));
            handles.plot.x_files(chosenPanel).x_(iFile) = x_;
        else
        	handles.plot.x_files(chosenPanel).x_ = x_;
        end
        
        % Get all source parameters from opt_ structure and save it to to "handles.plot.x_files." structure:
        %  - e.g. in case of satellite observations only, opt_.source is missing => Init. emopty field in handes struct.
        if isfield(opt_, 'source')
            handles.plot.x_files(chosenPanel).sources(iFile).source = opt_.source;
        else
            handles.plot.x_files(chosenPanel).sources(iFile).source = [];
        end
        if isfield(opt_, 'satellite')
            handles.plot.x_files(chosenPanel).sources(iFile).satellite = opt_.satellite;
        else
            handles.plot.x_files(chosenPanel).sources(iFile).satellite = [];
        end
        
    end

    waitbar(iFile/nX_Files, h_waitbar, sprintf('File %1.0f / %1.0f is loaded.', iFile, nX_Files))
end % for all x_ files

% delete entries which have no antenna file
handles.plot.x_files(chosenPanel).x_(logical(antennaFileNotExist))=[];

x_NamesShort(logical(antennaFileNotExist))=[];
allX_Files(logical(antennaFileNotExist))=[];

% update popupmenu
% if there is no session to be plotted
if isempty(x_NamesShort)
    set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_sessions']), 'String', ' ');
    msgbox('There is no session where data may be visualized!','No session has both x_ and antenna file!','warn')
else
    set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_sessions']), 'String', x_NamesShort);

    % update "parameter" and "station" popupmenus and put the fields allAnt and
    % allPar to the handles struct (for gettting the index for plotting)
    handles=update_antennaParameterPopupmenus_plotParameters(handles,chosenPanel);

    % enable objects in panel
    set(eval(['handles.radiobutton_plot_folder', num2str(chosenPanel), '_oneSess']), 'Enable', 'On')
    set(eval(['handles.radiobutton_plot_folder', num2str(chosenPanel), '_allSess']), 'Enable', 'On')
    set(eval(['handles.text_plot_folder', num2str(chosenPanel), '_param']), 'Enable', 'On')
    set(eval(['handles.text_plot_folder', num2str(chosenPanel), '_stat']), 'Enable', 'On')
    set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_param']), 'Enable', 'On')
    set(eval(['handles.popupmenu_plot_folder', num2str(chosenPanel), '_stat']), 'Enable', 'On')

end
% close waitbar
close(h_waitbar)

% plot right away
handles=plotToAxes(hObject, handles);

% Update handles structure
guidata(hObject, handles);
