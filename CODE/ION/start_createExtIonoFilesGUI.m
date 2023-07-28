function varargout = start_createExtIonoFilesGUI(varargin)
% START_CREATEEXTIONOFILESGUI MATLAB code for start_createExtIonoFilesGUI.fig
%      START_CREATEEXTIONOFILESGUI, by itself, creates a new START_CREATEEXTIONOFILESGUI or raises the existing
%      singleton*.
%
%      H = START_CREATEEXTIONOFILESGUI returns the handle to a new START_CREATEEXTIONOFILESGUI or the handle to
%      the existing singleton*.
%
%      START_CREATEEXTIONOFILESGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in START_CREATEEXTIONOFILESGUI.M with the given input arguments.
%
%      START_CREATEEXTIONOFILESGUI('Property','Value',...) creates a new START_CREATEEXTIONOFILESGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before start_createExtIonoFilesGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to start_createExtIonoFilesGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help start_createExtIonoFilesGUI

% Last Modified by GUIDE v2.5 09-Aug-2011 09:42:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @start_createExtIonoFilesGUI_OpeningFcn, ...
    'gui_OutputFcn',  @start_createExtIonoFilesGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before start_createExtIonoFilesGUI is made visible.
function start_createExtIonoFilesGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to start_createExtIonoFilesGUI (see VARARGIN)

% Choose default command line output for start_createExtIonoFilesGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes start_createExtIonoFilesGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% check if we are in the WORK dir, just as a reminder-write an error if not
curPath=pwd;
if ~strcmp(curPath(end-3:end), 'WORK')
    fprintf('Probably need to change to ''VieVS/WORK'' directory...\n')
end


% UIWAIT makes start_createExtTropoFilesGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% add auxiliary folder to path
addpath('../ION/PROGRAM/auxiliary/'); % still needed for gpt (rest is in VIE_MOD_V1d)
addpath('../OUT/'); % for azel_out_auto.m



% add auxiliary folder to path
% addpath('..\auxiliary\');


% % go to folder of this .m file and make it to cd
% slashIdx=strfind(mfilename('fullpath'), '\');
% curMFilePath=mfilename('fullpath');
% cd(curMFilePath(1:slashIdx(end)));

% fill listboxes
% define folder where ngs year folders are
ngsdir='../DATA/NGS/';

% get files and folders of this path
ngsdirContent=dir(ngsdir);

% go through all entries of ngsdir to find folders
ngsYrs=cell(size(ngsdirContent,1),1);
for k=1:size(ngsdirContent,1)
    % if is a directory, add it to cell array
    if exist([ngsdir, ngsdirContent(k).name], 'dir')
        ngsYrs{k}=ngsdirContent(k).name;
    end
end

% delete empty rows
ngsYrs(cellfun(@isempty, ngsYrs))=[];


% update listbox 1 (first from left)
set(handles.listbox_1, 'String', ngsYrs);

% fill listbox 4 (first from right, containing process lists)
% get content of process list folder
%pldir='../../../WORK/PROCESSLIST/';
pldir='../WORK/PROCESSLIST/';
pldirContent=dir(pldir);

% preallocating
process_lists=cell(size(pldirContent,1),1);
for k=1:size(pldirContent,1)
    % find out if is .mat file
    if size(pldirContent(k).name,2)>=5 % we need min "1.mat"
        if strcmp(pldirContent(k).name(end-3:end), '.mat')
            process_lists{k}= pldirContent(k).name(1:end-4);
        end
    end
end

% delete empty cells
process_lists(cellfun(@isempty, process_lists))=[];

% update listbox
set(handles.listbox_4, 'String', process_lists);

% get old content for listbox3
if exist('GUIsavings.mat', 'file')
    load('GUIsavings.mat');
    set(handles.listbox_3, 'String', listbox3_content);
    
    %listbox3_content=get(handles.listbox_3, 'String');
end


% Save the change you made to the structure
guidata(hObject,handles)



% --- Outputs from this function are returned to the command line.
function varargout = start_createExtIonoFilesGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_createTropoFiles.
function button_createTropoFiles_Callback(hObject, eventdata, handles)
% hObject    handle to button_createTropoFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set button to disable
% set(handles.button_createTropoFiles, 'Enable', 'Off');

tic

% ##### get chosen sessions in listbox 3 #####
chosenSessions = get(handles.listbox_3, 'String');

'##### GUI does not work with vgosDB. Define the process list in start_createExtIonoFilesGUI.m line 174:'
% chosenSessions = ['2022/22FEB11KL [vgosDB]'
%                   '2022/22FEB11KR [vgosDB]']
% OR
% load('PROCESSLIST/pl.mat')
% chosenSessions = process_list;


% Check, if session was selected:
if isempty(chosenSessions)
    text = 'Chose session from listbox';
    msgbox(text,'Help','warn');
    
    % set button to disable
    set(handles.button_createTropoFiles, 'Enable', 'On');
else
    
    % ##### create iono files #####
    
    % check if there are azel files for alle sessions
    % preallocate
    noAzelAvailable     = ones(size(chosenSessions, 1), 1);
    wantedAzelPath      = cell(size(chosenSessions, 1), 1); % path for all files where azel is located
    
    % get all files in azel directory
    azelPath = '../OUT/AZEL/';
    
    % if AZEL folder does not exist, create it
    if ~exist(azelPath, 'dir')
        mkdir(azelPath);
    end
    
    % get files in AZEL dir
    filesInAzelDir = dirr([azelPath, '*azel']);
    
    % for all chosen sessions
    for k = 1 : size(chosenSessions, 1)
        
        curSession = chosenSessions(k, 6:end);
        
        % Check, is the session name (in process list) contains the tags " [vgosDB]" or " [VSO]":
        % => This is required to get the actual filename in case of vgosDB or vso files.
        % => Just remove everything after the first blank!
        blank_ind = min(strfind(curSession, ' '));
        if ~isempty(blank_ind)
            curSession = curSession(1 : blank_ind-1);
        end
        
        % if there is an azel file with the same name as the current
        % session in the main directory AZEL
        if sum(~cellfun(@isempty, strfind({filesInAzelDir.name}, curSession)))~=0
            noAzelAvailable(k)  = 0;
            wantedAzelPath{k}   = azelPath;
        else % look for azel-file in subfolder of year (AZEL/2008/)
            % find index of current year in filesInAzelDir (where
            % azel-file could/should be
            wantedAzelPath{k} = [azelPath, chosenSessions(k,1:4), '/'];
            
            % if folder of current year does not exist, create it
            if ~exist(wantedAzelPath{k}, 'dir')
                mkdir(wantedAzelPath{k});
            end
            
            % get file in year-directory
            %filesInAzelYearDir=dirr(wantedAzelPath{k}, ['azel_', curSession]);
            
            % check if found or not
            %if ~isempty(filesInAzelYearDir)
            %noAzelAvailable(k)=0;
            %wantedAzelPath{k}=azelPath;
            %end
            
        end
    end
    
    % DO THIS IN ANY CASE | if there is at least one with available azel
    % get iono model (one of CODE, IGS,...)
    if get(handles.button_model_code, 'Value') == 1
        ionoModel = 'CODE';
    elseif get(handles.button_model_IGS, 'Value') == 1
        ionoModel = 'IGS';
    elseif get(handles.button_model_gnss, 'Value') == 1
        ionoModel = 'GNSS';
    elseif get(handles.button_model_GNSSAltimetry, 'Value') == 1
        ionoModel = 'GNSSAltimetry';
    elseif get(handles.button_model_GNSSAltimetryFC, 'Value') == 1
        ionoModel = 'GNSSAltimetryFC';
    end

'You can specify the model in start_createExtIonoFilesGUI.m line 265'
% ionoModel = 'IGS'
% ionoModel = 'UQR'
% ionoModel = 'EMR'
% ionoModel = 'CODE'
% ionoModel = 'JPL'
% ionoModel = 'ESA'

    
    % delete empty entries for all sessions where no azel file is available
    %noAzelAvailablePath=noAzelAvailablePath(~noAzelAvailable);
    
    % get user defined subdirectory
    subdirectory = get(handles.edit_subdir, 'String');
    
    % get files where we have azel files
    %azelAvailInd=find(noAzelAvailable==0);
    
    % preallocate array for error messages
    %errorMsgs=zeros(sum(noAzelAvailable==0),1);
    errorMsgs = zeros(size(chosenSessions, 1), 1);
    
    % ##### create external errorMsgstropospheric files for all sessions #####
    for k = 1 : size(chosenSessions, 1)
        
        curSession  = chosenSessions(k,6:end);
        % Check, is the session name (in process list) contains the tags " [vgosDB]" or " [VSO]":
        % => This is required to get the actual filename in case of vgosDB or vso files.
        % => Just remove everything after the first blank!
        blank_ind = min(strfind(curSession, ' '));
        if ~isempty(blank_ind)
            curSession = curSession(1 : blank_ind-1);
        end
        
        % Get name and path of AZEL file:
        curAzelFile = [wantedAzelPath{k}, 'azel_', curSession, '.txt'];
        
        % call the function createExtIonoFiles.m
        errorMsgs(k, 1) = createExtIonoFiles(curSession, subdirectory, curAzelFile, ionoModel);
    end
    
    if ~exist('errorMsgs', 'var')
        errorMsgs = 999;
    end
    
    % write processing log
    fprintf('\n\n ======== LOG =========\n\n');
    % ERROR=0: sucessful
    fprintf('%1.0f file(s): successfully created\n', sum(errorMsgs == 0))
    if sum(errorMsgs == 0) > 0
        idxNoError=find(errorMsgs == 0);
        
        for k = 1 : sum(errorMsgs == 0)
            fprintf('     session: %s\n', chosenSessions(idxNoError(k), :));
        end
    end
    fprintf('\n');
    
    % without error message (error is already found in this code)
    %     if sum(noAzelAvailable==1)>0
    %         fprintf('%1.0f file(s): No azel-file found\n', sum(noAzelAvailable));
    %
    %         % get indices of all stations with no azel available
    %         indNoAzelAvailable=find(noAzelAvailable);
    %
    %         % for all stations with this error
    %         for k=1:sum(noAzelAvailable)
    %             fprintf('     %s\n', chosenSessions(indNoAzelAvailable(k),:));
    %         end
    %     end
    %     fprintf('\n');
    
    % for all other error messages (6 currently)
    errorTexts={'Name of subfolder not allowed', ...
        'azelFile not found (probably need to process session in VieVS)', ...
        'iono map does not exist (was downloaded - see above)', ...
        'GNSS/GNSS+Altimetry/GNSS+Altimetry+F/C not yet implemented', ...
        'antenna struct not found (probably need to process session in VieVS)', ...
        'station (in azel) was not found in antenna struct'};
    
    % for all possible error messages
    for k = 1 : 6
        
        % if this error message was found
        if sum(errorMsgs == k) > 0
            fprintf('%1.0f file(s): %s\n', sum(errorMsgs == k), errorTexts{k})
            % write each session with this error message into one row
            
            %chosenSessionsAzelAvailable=chosenSessions(noAzelAvailable==0,:);
            
            % get indices with current error message
            idxCurError = find(errorMsgs == k);
            for m = 1 : sum(errorMsgs == k)
                fprintf('     %s\n', chosenSessions(idxCurError(m), :));
            end
            
            % write spacing line
            fprintf('\n');
        end
        
    end
    
    fprintf(' ====== END LOG =======\n\n');
    
    
end
toc

% set(handles.button_createTropoFiles, 'Enable', 'On');



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_subdir_Callback(hObject, eventdata, handles)
% hObject    handle to edit_subdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_subdir as text
%        str2double(get(hObject,'String')) returns contents of edit_subdir as a double


% --- Executes during object creation, after setting all properties.
function edit_subdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_subdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_1.
function listbox_1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_1

% get folder that was selected
contents = cellstr(get(hObject,'String'));
curFolder = contents{get(hObject,'Value')};

folderContent = dir(['../DATA/NGS/', curFolder]);

% update listbox 2
set(handles.listbox_2,'Value',1);
set(handles.listbox_2, 'String', {folderContent(3:end).name});

% save last selected year as var
handles.data.lastYearSelected=curFolder;

% Save the change you made to the structure
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function listbox_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_2.
function listbox_2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_2

% get selected file
contents = cellstr(get(hObject,'String'));
curFile=contents{get(hObject,'Value')};

% get content of listbox 3
curListbox3content=get(handles.listbox_3, 'String');

% if no file selected -> don't do anything
if ~isempty(curFile)
    
    if ~isempty(curListbox3content)
        
        % make cell out of this
        curListbox3contentCell=cell(size(curListbox3content,1),1);
        for k=1:size(curListbox3content,1)
            curListbox3contentCell{k,1}=curListbox3content(k,:);
        end
        
        % update listbox 3 if not already found in listbox 3
        if cellfun(@isempty,strfind(curListbox3contentCell, curFile))
            newEntries={curListbox3content; [handles.data.lastYearSelected, '/', curFile]};
            set(handles.listbox_3, 'String', char(newEntries));
        end
        
        % if listbox three is empty
    else
        set(handles.listbox_3,'Value',1);
        set(handles.listbox_3, 'String', [curListbox3content; handles.data.lastYearSelected, '/', curFile]);
    end
    
    % save new content to file on disk for loading at next startup
    listbox3_content=get(handles.listbox_3, 'String');
    save('GUIsavings.mat', 'listbox3_content');
    
    % Save the change you made to the structure
    guidata(hObject,handles)
end



% --- Executes during object creation, after setting all properties.
function listbox_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_3.
function listbox_3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_3

% get current content of listbox 3
curListbox3=get(handles.listbox_3, 'String');

if ~isempty(curListbox3)
    
    % get listbox entries into one cell array
    curListbox3_cell=cell(size(curListbox3, 1), 1);
    for k=1:size(curListbox3, 1)
        curListbox3_cell{k}=curListbox3(k,:);
    end
    
    % get selected item
    contents = cellstr(get(hObject,'String'));
    curSelectedItem=contents{get(hObject,'Value')};
    
    
    
    % find ind of selected file in listbox content
    idx2delete=find(~cellfun(@isempty, strfind(curListbox3_cell, curSelectedItem)));
    
    % delete entry
    curListbox3(idx2delete,:)=[];
    
    % update listbox 3
    set(handles.listbox_3,'Value',1);
    set(handles.listbox_3, 'String', curListbox3);
    
    % save new content to file on disk for loading at next startup
    listbox3_content=get(handles.listbox_3, 'String');
    save('GUIsavings.mat', 'listbox3_content');
    
    % Save the change you made to the structure
    guidata(hObject,handles)
    
end



% --- Executes during object creation, after setting all properties.
function listbox_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_4.
function listbox_4_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_4


% get selected process list
contents = cellstr(get(hObject,'String'));
curProcList=contents{get(hObject,'Value')};

% load(['../../../WORK/PROCESSLIST/', curProcList, '.mat']);
load(['./PROCESSLIST/', curProcList, '.mat']);


% for all entries in process list
for k = 1 : size(process_list, 1)% - 1
    curName=process_list(k,:);
    % get content of listbox 3
    curListbox3content=get(handles.listbox_3, 'String');
    
    if ~isempty(curListbox3content)
        
        % make cell out of this
        curListbox3contentCell=cell(size(curListbox3content,1),1);
        for m=1:size(curListbox3content,1)
            curListbox3contentCell{m,1}=curListbox3content(m,:);
        end
        % update listbox 3 if not already found in listbox 3
        if cellfun(@isempty,strfind(curListbox3contentCell, curName))
            set(handles.listbox_3, 'String', [curListbox3content; curName]);
        end
        % if listbox 3 is empty
    else
        set(handles.listbox_3,'Value',1);
        set(handles.listbox_3, 'String', [curListbox3content; curName]);
    end
end

% save content of listbox 3 to file on disk
listbox3_content = get(handles.listbox_3, 'String');
save('GUIsavings.mat', 'listbox3_content');

% Save the change you made to the structure
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function listbox_4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_model_gnss.
function button_model_gnss_Callback(hObject, eventdata, handles)
% hObject    handle to button_model_gnss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button_model_gnss


% --- Executes on button press in button_model_GNSSAltimetry.
function button_model_GNSSAltimetry_Callback(hObject, eventdata, handles)
% hObject    handle to button_model_GNSSAltimetry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button_model_GNSSAltimetry


% --- Executes on button press in button_model_GNSSAltimetryFC.
function button_model_GNSSAltimetryFC_Callback(hObject, eventdata, handles)
% hObject    handle to button_model_GNSSAltimetryFC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button_model_GNSSAltimetryFC


% --- Executes on button press in button_model_code.
function button_model_code_Callback(hObject, eventdata, handles)
% hObject    handle to button_model_code (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button_model_code


% --- Executes on button press in button_model_IGS.
function button_model_IGS_Callback(hObject, eventdata, handles)
% hObject    handle to button_model_IGS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of button_model_IGS


% --- Executes on key press with focus on button_model_code and none of its controls.
function button_model_code_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to button_model_code (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
