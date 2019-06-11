function varargout = createSupersourcesFile(varargin)
% CREATESUPERSOURCESFILE MATLAB code for createSupersourcesFile.fig
%      CREATESUPERSOURCESFILE, by itself, creates a new CREATESUPERSOURCESFILE or raises the existing
%      singleton*.
%
%      H = CREATESUPERSOURCESFILE returns the handle to a new CREATESUPERSOURCESFILE or the handle to
%      the existing singleton*.
%
%      CREATESUPERSOURCESFILE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CREATESUPERSOURCESFILE.M with the given input arguments.
%
%      CREATESUPERSOURCESFILE('Property','Value',...) creates a new CREATESUPERSOURCESFILE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before createSupersourcesFile_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to createSupersourcesFile_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help createSupersourcesFile

% Last Modified by GUIDE v2.5 05-Dec-2018 09:57:17
% Revision:
%   26 Sep 2016 by Hana Krasna: VieCRF13 instead of VieCRF10a, IERS or IVS
%   names can be choosen in UserCRF

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @createSupersourcesFile_OpeningFcn, ...
                   'gui_OutputFcn',  @createSupersourcesFile_OutputFcn, ...
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


% --- Executes just before createSupersourcesFile is made visible.
function createSupersourcesFile_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to createSupersourcesFile (see VARARGIN)

superstatFolder='../CRF/';

% define names of all files and save it in handles struct
handles.data.allFileDescriptions={...
    'ICRF2', 'icrf2VcsOnly', 'ICRF3sx', 'VieCRF10a', 'vievsCrf', 'gsf2015b',...
    'userCrf', 'translation'};

handles.data.allFileNames={'icrf2-non-vcs.dat', 'icrf2-vcs-only.dat',...
    'icrf3sx.txt', 'VieCRF13.txt', 'vievsCrf.txt','gsf2015b_astro.sou.txt', '''userdefined''', ...
    'IVS_SrcNamesTable.txt'};

% load gui state if file exist
guiStateFile=[superstatFolder, '/SavedGuiData_supersources.txt'];
if exist(guiStateFile, 'file')
    fid=fopen(guiStateFile);
    nFiles=size(handles.data.allFileDescriptions, 2);
    
    guiOptLoaded=textscan(fid, '%s %s', nFiles*3+1, 'delimiter', '|', ...
        'commentstyle', '#');
    
    fclose(fid);
    
    for iFile=1:nFiles
        eval(['set(handles.radiobutton_', handles.data.allFileDescriptions{iFile}, '_download, ''Value'', str2double(guiOptLoaded{2}{3*iFile-2}))'])
        if str2double(guiOptLoaded{2}{3*iFile-2})
            eval(['set(handles.edit_', handles.data.allFileDescriptions{iFile},', ''Enable'', ''Off'');']);
            eval(['set(handles.pushbutton_', handles.data.allFileDescriptions{iFile},'_browse, ''Enable'', ''Off'');']);
        else
            eval(['set(handles.edit_', handles.data.allFileDescriptions{iFile},', ''Enable'', ''On'');']);
            eval(['set(handles.pushbutton_', handles.data.allFileDescriptions{iFile},'_browse, ''Enable'', ''On'');']);
        end
        eval(['set(handles.radiobutton_', handles.data.allFileDescriptions{iFile}, '_specify, ''Value'', str2double(guiOptLoaded{2}{3*iFile-1}))'])
        eval(['set(handles.edit_', handles.data.allFileDescriptions{iFile},', ''String'', guiOptLoaded{2}{3*iFile});']);

    end
    set(handles.edit_createSupersource, 'String', guiOptLoaded{2}{3*nFiles+1})
end
    

% Choose default command line output for createSupersourcesFile
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes createSupersourcesFile wait for user response (see UIRESUME)
% uiwait(handles.figure_main);

function saveGuiDataToDisk(handles)
% This function saves the current state of the gui to a textfile.

outFile='../CRF/SavedGuiData_supersources.txt';

fid=fopen(outFile, 'w');

% for all files
for iFile=1:size(handles.data.allFileDescriptions,2)
    % "header"
    fprintf(fid, '# %s --------------\n', handles.data.allFileNames{iFile});
    % download
    fprintf(fid, '%21s, download | %1.0f\n', handles.data.allFileDescriptions{iFile},...
        get(...
        eval(['handles.radiobutton_', handles.data.allFileDescriptions{iFile}, '_download']), 'Value'));
    % specify
    fprintf(fid, '%22s, specify | %1.0f\n', handles.data.allFileDescriptions{iFile},...
        get(...
        eval(['handles.radiobutton_', handles.data.allFileDescriptions{iFile}, '_specify']), 'Value'));
    % edit
    fprintf(fid, '%18s, editTextbox | %s\n', handles.data.allFileDescriptions{iFile},...
        get(...
        eval(['handles.edit_', handles.data.allFileDescriptions{iFile}]), 'String'));
    % end
    fprintf(fid, '#\n');
    
end

% write outfile
fprintf(fid, '                        outFile | %s\n',...
    get(handles.edit_createSupersource, 'String'));

% close txt file
fclose(fid);

% --- Outputs from this function are returned to the command line.
function varargout = createSupersourcesFile_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_searchForFiles.
function pushbutton_searchForFiles_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_searchForFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

searchPath='../CRF/data/';

nFiles=size(handles.data.allFileDescriptions,2);

foundLogicals=zeros(nFiles,1);

for k=1:nFiles
    
    % search
    foundFiles=dir([searchPath, handles.data.allFileNames{k}]);
    
    if ~isempty(foundFiles)
        eval(['set(handles.edit_', handles.data.allFileDescriptions{k},', ''String'', [searchPath, foundFiles(1).name]);']);
        
        % enable textbox and browse button
        eval(['set(handles.edit_', handles.data.allFileDescriptions{k},', ''Enable'', ''On'');']);
        eval(['set(handles.pushbutton_', handles.data.allFileDescriptions{k},'_browse, ''Enable'', ''On'');']);
        
        % ... set 1 to specify and 0 to download ...
        eval(['set(handles.radiobutton_', handles.data.allFileDescriptions{k},'_download, ''Value'', 0);']);
        eval(['set(handles.radiobutton_', handles.data.allFileDescriptions{k},'_specify, ''Value'', 1);']);
        
        foundLogicals(k,1)=1;
    end
end

h=msgbox(sprintf('%1.0f of %1.0f files found and written to textboxes',...
    sum(foundLogicals), nFiles), 'Automatic Search Done', 'warn');

saveGuiDataToDisk(handles);


function edit_vievsCrf_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vievsCrf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vievsCrf as text
%        str2double(get(hObject,'String')) returns contents of edit_vievsCrf as a double


% --- Executes during object creation, after setting all properties.
function edit_vievsCrf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vievsCrf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_vievsCrf_browse.
function pushbutton_vievsCrf_browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_vievsCrf_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseForFile(hObject,handles);


function browseForFile(hObject, handles)
% This function lets the user select a file
% get current entry of strings
tagOfHObject=get(hObject,'Tag');
curEditString=get(eval(['handles.edit_', tagOfHObject(12:end-7)]), 'String');

if isempty(curEditString)
    curEditString='../CRF/data/';
end

% get file from explorer
[FileName, PathName] = uigetfile('*.*','Select an input file', curEditString);

if ~isempty(FileName)
    if ~strcmp('0', num2str(FileName))
        
        % set edit textbox
        set(eval(['handles.edit_', tagOfHObject(12:end-7)]), 'String', [PathName, FileName])
        %set(handles.edit_nsCodesFile, 'String', [PathName, FileName]);

        % save GUI data
        saveGuiDataToDisk(handles);        
    end
end




function edit_ICRF3sx_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ICRF3sx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ICRF3sx as text
%        str2double(get(hObject,'String')) returns contents of edit_ICRF3sx as a double


% --- Executes during object creation, after setting all properties.
function edit_ICRF3sx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ICRF3sx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_ICRF3sx_browse.
function pushbutton_ICRF3sx_browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ICRF3sx_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseForFile(hObject,handles);



function edit_ICRF2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ICRF2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ICRF2 as text
%        str2double(get(hObject,'String')) returns contents of edit_ICRF2 as a double


% --- Executes during object creation, after setting all properties.
function edit_ICRF2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ICRF2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_ICRF2_browse.
function pushbutton_ICRF2_browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ICRF2_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseForFile(hObject,handles);



% --- Executes on button press in pushbutton_browse.
function pushbutton_browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get current entry of strings
curEditString=get(handles.edit_createSupersource, 'String');

if isempty(curEditString)
    curEditString='../CRF/';
end

% get file from explorer
[FileName, PathName] = uiputfile('*.mat','Select an output .mat-file', curEditString);

if ~isempty(FileName)
    if ~strcmp('0', num2str(FileName))
        
        % set edit textbox
        set(handles.edit_createSupersource, 'String', [PathName, FileName])
        %set(handles.edit_nsCodesFile, 'String', [PathName, FileName]);

        % save GUI data
        saveGuiDataToDisk(handles);        
    end
end


% --- Executes on button press in pushbutton_create.
function pushbutton_create_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_create (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nFiles=size(handles.data.allFileNames,2);
inFiles(nFiles)=struct('name', []);

% for all files
for k=1:nFiles
    inFiles(k).name=eval(['get(handles.edit_', handles.data.allFileDescriptions{k},', ''String'');']);
end

% save the information about the source names in the user crf
UserCrfNames.IERS=get(handles.radiobutton_UserIERS,'Value');
UserCrfNames.IVS=get(handles.radiobutton_UserIVS,'Value');

message=mk_supersourceFile(inFiles, get(handles.edit_createSupersource, 'String'),UserCrfNames);

% write user info message
msgbox(message, 'Information about creation of supersource file', 'warn');


function edit_createSupersource_Callback(hObject, eventdata, handles)
% hObject    handle to edit_createSupersource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_createSupersource as text
%        str2double(get(hObject,'String')) returns contents of edit_createSupersource as a double


% --- Executes during object creation, after setting all properties.
function edit_createSupersource_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_createSupersource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_VieCRF10a_Callback(hObject, eventdata, handles)
% hObject    handle to edit_VieCRF10a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_VieCRF10a as text
%        str2double(get(hObject,'String')) returns contents of edit_VieCRF10a as a double


% --- Executes during object creation, after setting all properties.
function edit_VieCRF10a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_VieCRF10a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_VieCRF10a_browse.
function pushbutton_VieCRF10a_browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_VieCRF10a_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseForFile(hObject,handles);


% --- Executes when selected object is changed in uipanel2.
function uipanel2_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel2 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


if get(handles.radiobutton_ICRF2_download, 'Value')
    state='Off';
else
    state='On';
end
set(handles.edit_ICRF2, 'Enable', state);
set(handles.pushbutton_ICRF2_browse, 'Enable', state)

saveGuiDataToDisk(handles)


% --- Executes when selected object is changed in uipanel3.
function uipanel3_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel3 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

if get(handles.radiobutton_ICRF3sx_download, 'Value')
    state='Off';
else
    state='On';
end
set(handles.edit_ICRF3sx, 'Enable', state)
set(handles.pushbutton_ICRF3sx_browse, 'Enable', state)

saveGuiDataToDisk(handles)


% --- Executes when selected object is changed in uipanel6.
function uipanel6_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel6 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

if get(handles.radiobutton_VieCRF10a_download, 'Value')
    state='Off';
else
    state='On';
end
set(handles.edit_VieCRF10a, 'Enable', state)
set(handles.pushbutton_VieCRF10a_browse, 'Enable', state)

saveGuiDataToDisk(handles)


% --- Executes when selected object is changed in uipanel4.
function uipanel4_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel4 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

if get(handles.radiobutton_vievsCrf_download, 'Value')
    state='Off';
else
    state='On';
end
set(handles.edit_vievsCrf, 'Enable', state)
set(handles.pushbutton_vievsCrf_browse, 'Enable', state)

saveGuiDataToDisk(handles)



function edit_userCrf_Callback(hObject, eventdata, handles)
% hObject    handle to edit_userCrf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_userCrf as text
%        str2double(get(hObject,'String')) returns contents of edit_userCrf as a double


% --- Executes during object creation, after setting all properties.
function edit_userCrf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_userCrf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_userCrf_browse.
function pushbutton_userCrf_browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_userCrf_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseForFile(hObject,handles);



function edit_blokq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_blokq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_blokq as text
%        str2double(get(hObject,'String')) returns contents of edit_blokq as a double


% --- Executes during object creation, after setting all properties.
function edit_blokq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_blokq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_blokq_browse.
function pushbutton_blokq_browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_blokq_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_translation_Callback(hObject, eventdata, handles)
% hObject    handle to edit_translation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_translation as text
%        str2double(get(hObject,'String')) returns contents of edit_translation as a double


% --- Executes during object creation, after setting all properties.
function edit_translation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_translation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_translation_browse.
function pushbutton_translation_browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_translation_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


browseForFile(hObject,handles);

function edit_icrf2VcsOnly_Callback(hObject, eventdata, handles)
% hObject    handle to edit_icrf2VcsOnly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_icrf2VcsOnly as text
%        str2double(get(hObject,'String')) returns contents of edit_icrf2VcsOnly as a double


% --- Executes during object creation, after setting all properties.
function edit_icrf2VcsOnly_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_icrf2VcsOnly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_icrf2VcsOnly_browse.
function pushbutton_icrf2VcsOnly_browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_icrf2VcsOnly_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


browseForFile(hObject,handles);


% --- Executes on button press in pushbutton_gsf2015b_browse.
function pushbutton_gsf2015b_browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_gsf2015b_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
browseForFile(hObject,handles);



function edit_gsf2015b_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gsf2015b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gsf2015b as text
%        str2double(get(hObject,'String')) returns contents of edit_gsf2015b as a double


% --- Executes during object creation, after setting all properties.
function edit_gsf2015b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gsf2015b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
