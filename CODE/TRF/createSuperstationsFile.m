function varargout = createSuperstationsFile(varargin)
% This interface creates superstations files for VieVS: These .mat files
% include all the information available for all stations in ns_codes.txt
% (IVS stations). These include TRF (ITRF, VTRF) for different epochs
% (2005, 2008), deformation, loading coefficients, "our" new coordinates
% after an earthquake...
% created: Matthias Madzak, 6.10.2011

% updates:
% 05 June 2014 by Hana Krasna, VieTRF10a changed to VieTRF13
% 18 June 2014 by Hana Krasna, Vienna tidal APL added
% 20 Oct  2017 by Hana Krasna, DTRF2014 added


% CREATESUPERSTATIONSFILE MATLAB code for createSuperstationsFile.fig
%      CREATESUPERSTATIONSFILE, by itself, creates a new CREATESUPERSTATIONSFILE or raises the existing
%      singleton*.
%
%      H = CREATESUPERSTATIONSFILE returns the handle to a new CREATESUPERSTATIONSFILE or the handle to
%      the existing singleton*.
%
%      CREATESUPERSTATIONSFILE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CREATESUPERSTATIONSFILE.M with the given input arguments.
%
%      CREATESUPERSTATIONSFILE('Property','Value',...) creates a new CREATESUPERSTATIONSFILE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before createSuperstationsFile_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to createSuperstationsFile_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help createSuperstationsFile

% Last Modified by GUIDE v2.5 20-Dec-2019 13:16:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @createSuperstationsFile_OpeningFcn, ...
                   'gui_OutputFcn',  @createSuperstationsFile_OutputFcn, ...
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


% --- Executes just before createSuperstationsFile is made visible.
function createSuperstationsFile_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to createSuperstationsFile (see VARARGIN)


% define all panel tags (names) to be made visible when other panel is selected
handles.allPanels={'uipanel_mainInfoFiles', 'uipanel_trf', ...
    'uipanel_atmoLoading', 'uipanel_oceanLoading'};


handles.data.allFileDescriptions={'nsCodes', 'antennaInfo', 'eccdat',...
        'gravdef',...
        'itrf2014','dtrf2014', 'vtrf2014dat', 'ivstrf2014b'...
        'vieTrf13', 's12vienna','s12vandam', 's12gsfc', 's12userOwn',...
        'oceanLoadingTPXO72', 'oceanLoadingGOT00', 'oceanLoadingEOT08a', 'oceanLoadingFES2004', 'oceanLoadingAG06', 'oceanLoadingUserOwn' ...
        'opoleloadcoefcmcor', 'opoleloadUserOwn', 'vievsTrf', 'userOwnTrf'};
    

handles.data.allFileNames={'ns-codes.txt', 'antenna-info.txt', 'ECCDAT.ecc',...
    'gravitationalDeformation.txt',...
    'ITRF2014-IVS-TRF.snx', 'DTRF2014_VLBI.snx', 'VTRF2014_final.snx', ...
    'IVS_TRF2014b.SSC.txt', 'VieTRF13.txt',...
    's1_s2_s3_cm_noib_vlbi.dat','s12_cm_noib_vandam.dat', 's12_cm_noib_gsfc.dat',  's12userOwn.dat',...
    'ocean_loading_TPXO72.TXT', 'ocean_loading_GOT00.TXT', 'ocean_loading_EOT08a.TXT', 'ocean_loading_FES2004.TXT', 'ocean_loading_AG06.TXT', 'ocean_loading_userOwn.TXT',...
    'opoleloadcoefcmcor.txt', 'opoleloadUserOwn.txt', 'vievsTrf.txt', 'userOwnTrf.txt'};
    
% load gui options if "SavedGuiData_superstations.txt" exist
guiSateFile='../TRF/SavedGuiData_superstations.txt';
nFiles=size(handles.data.allFileDescriptions,2); % thats the desired length of the gui-state file
if exist(guiSateFile, 'file')
    fid=fopen(guiSateFile);
    % load
    guiOptLoaded=textscan(fid, '%s %s', nFiles*3+1, 'delimiter', '|', 'commentstyle', '#'); % +1: also outFile should be read
    
    % close
    fclose(fid);
    
    % if enough files were read in
    if length(guiOptLoaded{1})==(nFiles+1) % then I have enough files from default file
        for k=1:nFiles
            eval(['set(handles.edit_', handles.data.allFileDescriptions{k},'File, ''String'', guiOptLoaded{2}{', num2str(k),'});']);
        end
        
        % set output file and path to textbox
        set(handles.edit_createSuperstations, 'String', guiOptLoaded{2}{k+1});
    else
        % seems like: old gui-state txt file (not enough input options)
        delete(guiSateFile);
    end
end    

% Choose default command line output for createSuperstationsFile
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = createSuperstationsFile_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function makeAllPanelsInvisible(handles)
% This function makes all panels (main files, trf, ocean loading, atmo
% loading) invisible, so that one (the one selected in the menu) can be
% made visible

a=5;
for k=1:size(handles.allPanels,2)
    set(handles.(handles.allPanels{k}), 'Visible', 'Off')
end


function edit_nsCodesFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nsCodesFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nsCodesFile as text
%        str2double(get(hObject,'String')) returns contents of edit_nsCodesFile as a double


% --- Executes during object creation, after setting all properties.
function edit_nsCodesFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nsCodesFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function saveGuiDataToDisk(handles)
% This function saves the current state of the interface to a file
% (SavedGuiData.txt) to be loaded at next time when the gui opens.

guiSateFile='../TRF/SavedGuiData_superstations.txt';
% open file
fid=fopen(guiSateFile, 'w');

text={'ns-codes, editTextbox', ...
    'antenna-info, editTextBox',...
    'eccdat, editTextBox',...
    'gravdef, editTextBox',...
    'itrf2014, editTextBox',...
    'dtrf2014, editTextBox',...
    'vtrf2014dat, editTextBox',...
    'ivstrf2014b, editTextBox',...
    'vieTrf13, editTextBox',...
    's12vienna, editTextBox',...    
    's12vandam, editTextBox',...
    's12gsfc, editTextBox',...
    's12userOwn, editTextbox',...
    'oceanLoadTPXO72, editTextBox',...
    'oceanLoadGOT00, editTextBox',...
    'oceanLoadEOT08a, editTextBox',...
    'oceanLoadFES2004, editTextBox',...
    'oceanLoadAG06, editTextBox',...
    'oceanLoadUserOwn, editTextBox',...
    'opoleloadcoefcmcor, editTextBox',...
    'opoleloadUserOwn, editTextBox',...
    'vievsTrf, editTextBox',...
    'userOwnTrf, editTextBox'};

% get max length of text
l=num2str(max(cellfun('length', text)) + 2);

nFiles=size(handles.data.allFileNames,2);

% write gui options to file
for k=1:nFiles
    fprintf(fid, ['# ', handles.data.allFileNames{k},' --------------\n']);
    fprintf(fid, ['%',l,'s | %s\n'], text{k},  eval(['get(handles.edit_', handles.data.allFileDescriptions{k},'File, ''String'')']));
    fprintf(fid, '#\n');
end

% write output path+file to file
fprintf(fid, ['\n#\n%', l,'s | %s\n'], 'outFile', get(handles.edit_createSuperstations, 'String'));

fclose(fid);

% --- Executes on button press in pushbutton_nsCodesBrowse.
function pushbutton_nsCodesBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_nsCodesBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)



function browseFunction(hObject, handles)
% This button performs the browsing (Browse) of the wanted files, ie it
% creates an open dialog and writes (if chosen) the chosen file into the
% wanted textbox.

% get handle of edit textbox I want to write the chosen file into
        
% 1. get all children of panel where the current Browse button is inside
allChildren=get(get(hObject, 'Parent'), 'children');

% 2. get the handle of the edit textbox (there are sometimes to edit -> so
% also check for the Tag "File" at end of Tag!
for k=1:size(allChildren,1)
    curTag=get(allChildren(k),'Tag');
    if strcmp(get(allChildren(k), 'Style'), 'edit') && ...
            strcmp(curTag(end-3:end), 'File')
        break;
    end
end    

% get file from explorer
[FileName, PathName] = uigetfile('*.*',['Select an ', get(get(hObject, 'Parent'), 'title'),' file'], get(allChildren(k), 'String'));

if ~isempty(FileName)
    if ~strcmp('0', num2str(FileName))
        
        % set edit textbox
        set(allChildren(k), 'String', [PathName, FileName])
        %set(handles.edit_nsCodesFile, 'String', [PathName, FileName]);


        % save GUI data
        saveGuiDataToDisk(handles);        
    end
end
    
% --- Executes on button press in pushbutton_antennaInfoBrowse.
function pushbutton_antennaInfoBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_antennaInfoBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)


function edit_antennaInfoFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_antennaInfoFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_antennaInfoFile as text
%        str2double(get(hObject,'String')) returns contents of edit_antennaInfoFile as a double


% --- Executes during object creation, after setting all properties.
function edit_antennaInfoFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_antennaInfoFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_searchForFiles.
function pushbutton_searchForFiles_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_searchForFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This function searches the 'data' folder and puts the names of the
% found files (if more than one found, simply the first is taken) into the
% corresponding textbox

% get path to be searched
searchPath='../TRF/data/';

% get number of files
nFiles=size(handles.data.allFileNames,2);
% preallocate logical for finding files (0=not found, 1=found)
foundLogicals=zeros(nFiles,1);

% for all files
for k=1:nFiles
    % search
    foundFiles=dir([searchPath, handles.data.allFileNames{k}]);%   [curFile(1:allSlashes(size(allSlashes,2))), '*x_*']);
    
    % if found,  
    if ~isempty(foundFiles)
        % put it to textbox 
        eval(['set(handles.edit_', handles.data.allFileDescriptions{k},'File, ''String'', [searchPath, foundFiles(1).name]);']);
        
        % ... set Enables for textbox and browse button ...
%         eval(['set(handles.edit_', handles.data.allFileDescriptions{k},'File, ''Enable'', ''On'');']);
%         eval(['set(handles.pushbutton_', handles.data.allFileDescriptions{k},'Browse, ''Enable'', ''On'');']);
              
        % ... and write to foundLogicals
        foundLogicals(k,1)=1;
    else
        % put it to textbox 
        eval(['set(handles.edit_', handles.data.allFileDescriptions{k},'File, ''String'', '''');']);
    end   
end
    

% write message textbox
h=msgbox(sprintf('%1.0f of %1.0f files found and written to textboxes', sum(foundLogicals), nFiles), 'Automatic Search Done', 'warn');  

 
saveGuiDataToDisk(handles);

% --- Executes on button press in pushbutton_opoleloadcoefcmcorBrowse.
function pushbutton_opoleloadcoefcmcorBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_opoleloadcoefcmcorBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)



function edit_opoleloadcoefcmcorFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_opoleloadcoefcmcorFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_opoleloadcoefcmcorFile as text
%        str2double(get(hObject,'String')) returns contents of edit_opoleloadcoefcmcorFile as a double


% --- Executes during object creation, after setting all properties.
function edit_opoleloadcoefcmcorFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_opoleloadcoefcmcorFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_eccdatBrowse.
function pushbutton_eccdatBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_eccdatBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)



function edit_eccdatFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_eccdatFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_eccdatFile as text
%        str2double(get(hObject,'String')) returns contents of edit_eccdatFile as a double


% --- Executes during object creation, after setting all properties.
function edit_eccdatFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_eccdatFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_s12vandamBrowse.
function pushbutton_s12vandamBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_s12vandamBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)



function edit_s12vandamFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_s12vandamFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_s12vandamFile as text
%        str2double(get(hObject,'String')) returns contents of edit_s12vandamFile as a double


% --- Executes during object creation, after setting all properties.
function edit_s12vandamFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_s12vandamFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in uipanel_antenna.
function uipanel_antenna_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_antenna 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% if "download" button was chosen
if strcmp(get(eventdata.NewValue, 'Tag'), 'radiobutton_antennaDownload')
    set(handles.edit_antennaFile, 'Enable', 'Off');
    set(handles.pushbutton_antennaBrowse, 'Enable', 'Off');
else
    set(handles.edit_antennaFile, 'Enable', 'On');
    set(handles.pushbutton_antennaBrowse, 'Enable', 'On');
end

saveGuiDataToDisk(handles)


% --- Executes on button press in pushbutton_vievsTrfBrowse.
function pushbutton_vievsTrfBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_vievsTrfBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)



function edit_vievsTrfFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vievsTrfFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vievsTrfFile as text
%        str2double(get(hObject,'String')) returns contents of edit_vievsTrfFile as a double


% --- Executes during object creation, after setting all properties.
function edit_vievsTrfFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vievsTrfFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel_opoleloadcoefcmcor.
function uipanel_opoleloadcoefcmcor_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_opoleloadcoefcmcor 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% if "download" button was chosen
if strcmp(get(eventdata.NewValue, 'Tag'), 'radiobutton_opoleloadcoefcmcorDownload')
    set(handles.edit_opoleloadcoefcmcorFile, 'Enable', 'Off');
    set(handles.pushbutton_opoleloadcoefcmcorBrowse, 'Enable', 'Off');
else
    set(handles.edit_opoleloadcoefcmcorFile, 'Enable', 'On');
    set(handles.pushbutton_opoleloadcoefcmcorBrowse, 'Enable', 'On');
end

saveGuiDataToDisk(handles)


% --- Executes on button press in pushbutton_s12gsfcBrowse.
function pushbutton_s12gsfcBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_s12gsfcBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)



function edit_s12gsfcFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_s12gsfcFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_s12gsfcFile as text
%        str2double(get(hObject,'String')) returns contents of edit_s12gsfcFile as a double


% --- Executes during object creation, after setting all properties.
function edit_s12gsfcFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_s12gsfcFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_oceanLoadingFES2004Browse.
function pushbutton_oceanLoadingFES2004Browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_oceanLoadingFES2004Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)


function edit_oceanLoadingFES2004File_Callback(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingFES2004File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_oceanLoadingFES2004File as text
%        str2double(get(hObject,'String')) returns contents of edit_oceanLoadingFES2004File as a double


% --- Executes during object creation, after setting all properties.
function edit_oceanLoadingFES2004File_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingFES2004File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel_oceanLoadingFES2004.
function uipanel_oceanLoadingFES2004_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_oceanLoadingFES2004 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% if "download" button was chosen
if strcmp(get(eventdata.NewValue, 'Tag'), 'radiobutton_oceanLoadingFES2004Download')
    set(handles.edit_oceanLoadingFES2004File, 'Enable', 'Off');
    set(handles.pushbutton_oceanLoadingFES2004Browse, 'Enable', 'Off');
else
    set(handles.edit_oceanLoadingFES2004File, 'Enable', 'On');
    set(handles.pushbutton_oceanLoadingFES2004Browse, 'Enable', 'On');
end

saveGuiDataToDisk(handles)


% --- Executes on button press in pushbutton_createSuperstationsBrowse.
function pushbutton_createSuperstationsBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_createSuperstationsBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get current entry of strings
curEditString=get(handles.edit_createSuperstations, 'String');

if isempty(curEditString)
    curEditString='../TRF/';
end

% get file from explorer
[FileName, PathName] = uiputfile('*.mat','Select an output .mat-file', curEditString);

if ~isempty(FileName)
    if ~strcmp('0', num2str(FileName))
        
        % set edit textbox
        set(handles.edit_createSuperstations, 'String', [PathName, FileName])
        %set(handles.edit_nsCodesFile, 'String', [PathName, FileName]);


        % save GUI data
        saveGuiDataToDisk(handles);        
    end
end


% --- Executes on button press in pushbutton_createSuperstationsCreate.
function pushbutton_createSuperstationsCreate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_createSuperstationsCreate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This function calls the program to create the wanted superstations file.

% save gui state
saveGuiDataToDisk(handles)

% get number of input files
nInFiles=size(handles.data.allFileNames,2);

% preallocating
inFiles(nInFiles)=struct('name', []);

% for all files -> write to one struct
% Get locations of input files and fieldnames (optional)
for k=1:nInFiles
    inFiles(k).name=eval(['get(handles.edit_', handles.data.allFileDescriptions{k},'File, ''String'');']);
    if isfield(handles, eval(['''edit_', handles.data.allFileDescriptions{k},'_structFieldName''']))
        inFiles(k).fieldName=eval(['get(handles.edit_', handles.data.allFileDescriptions{k},'_structFieldName, ''String'');']);
    end
end

% put (additional) fieldnames for ocean loading to struct


% call function
message=mk_superstatFile(inFiles, get(handles.edit_createSuperstations, 'String'));

% write message box
msgbox(message, 'Information about creation of superstation file', 'warn');

function edit_createSuperstations_Callback(hObject, eventdata, handles)
% hObject    handle to edit_createSuperstations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_createSuperstations as text
%        str2double(get(hObject,'String')) returns contents of edit_createSuperstations as a double


saveGuiDataToDisk(handles);


% --- Executes during object creation, after setting all properties.
function edit_createSuperstations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_createSuperstations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menu_mainInfoFiles_Callback(hObject, eventdata, handles)
% hObject    handle to menu_mainInfoFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function menu_trf_Callback(hObject, eventdata, handles)
% hObject    handle to menu_trf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function menu_atmo_Callback(hObject, eventdata, handles)
% hObject    handle to menu_atmo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function menu_oceanLoading_Callback(hObject, eventdata, handles)
% hObject    handle to menu_oceanLoading (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





function edit_oceanLoadingAG06File_Callback(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingAG06File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_oceanLoadingAG06File as text
%        str2double(get(hObject,'String')) returns contents of edit_oceanLoadingAG06File as a double


% --- Executes during object creation, after setting all properties.
function edit_oceanLoadingAG06File_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingAG06File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_oceanLoadingAG06Browse.
function pushbutton_oceanLoadingAG06Browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_oceanLoadingAG06Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)

function edit_oceanLoadingAG06_structFieldName_Callback(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingAG06_structFieldName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_oceanLoadingAG06_structFieldName as text
%        str2double(get(hObject,'String')) returns contents of edit_oceanLoadingAG06_structFieldName as a double


% --- Executes during object creation, after setting all properties.
function edit_oceanLoadingAG06_structFieldName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingAG06_structFieldName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_oceanLoadingTPXO72File_Callback(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingTPXO72File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_oceanLoadingTPXO72File as text
%        str2double(get(hObject,'String')) returns contents of edit_oceanLoadingTPXO72File as a double


% --- Executes during object creation, after setting all properties.
function edit_oceanLoadingTPXO72File_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingTPXO72File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_oceanLoadingTPXO72Browse.
function pushbutton_oceanLoadingTPXO72Browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_oceanLoadingTPXO72Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)

function edit_oceanLoadingTPXO72_structFieldName_Callback(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingTPXO72_structFieldName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_oceanLoadingTPXO72_structFieldName as text
%        str2double(get(hObject,'String')) returns contents of edit_oceanLoadingTPXO72_structFieldName as a double


% --- Executes during object creation, after setting all properties.
function edit_oceanLoadingTPXO72_structFieldName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingTPXO72_structFieldName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_oceanLoadingEOT08aFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingEOT08aFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_oceanLoadingEOT08aFile as text
%        str2double(get(hObject,'String')) returns contents of edit_oceanLoadingEOT08aFile as a double


% --- Executes during object creation, after setting all properties.
function edit_oceanLoadingEOT08aFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingEOT08aFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_oceanLoadingEOT08aBrowse.
function pushbutton_oceanLoadingEOT08aBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_oceanLoadingEOT08aBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)

function edit_oceanLoadingEOT08a_structFieldName_Callback(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingEOT08a_structFieldName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_oceanLoadingEOT08a_structFieldName as text
%        str2double(get(hObject,'String')) returns contents of edit_oceanLoadingEOT08a_structFieldName as a double


% --- Executes during object creation, after setting all properties.
function edit_oceanLoadingEOT08a_structFieldName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingEOT08a_structFieldName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_oceanLoadingGOT00File_Callback(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingGOT00File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_oceanLoadingGOT00File as text
%        str2double(get(hObject,'String')) returns contents of edit_oceanLoadingGOT00File as a double


% --- Executes during object creation, after setting all properties.
function edit_oceanLoadingGOT00File_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingGOT00File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_oceanLoadingGOT00Browse.
function pushbutton_oceanLoadingGOT00Browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_oceanLoadingGOT00Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)

function edit_oceanLoadingGOT00_structFieldName_Callback(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingGOT00_structFieldName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_oceanLoadingGOT00_structFieldName as text
%        str2double(get(hObject,'String')) returns contents of edit_oceanLoadingGOT00_structFieldName as a double


% --- Executes during object creation, after setting all properties.
function edit_oceanLoadingGOT00_structFieldName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingGOT00_structFieldName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_oceanLoadingFES2004_structFieldName_Callback(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingFES2004_structFieldName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_oceanLoadingFES2004_structFieldName as text
%        str2double(get(hObject,'String')) returns contents of edit_oceanLoadingFES2004_structFieldName as a double


% --- Executes during object creation, after setting all properties.
function edit_oceanLoadingFES2004_structFieldName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingFES2004_structFieldName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_vieTrf13File_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vieTrf13File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vieTrf13File as text
%        str2double(get(hObject,'String')) returns contents of edit_vieTrf13File as a double


% --- Executes during object creation, after setting all properties.
function edit_vieTrf13File_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vieTrf13File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_vieTrf13Browse.
function pushbutton_vieTrf13Browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_vieTrf13Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
browseFunction(hObject, handles)


function edit_userOwnTrfFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_userOwnTrfFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_userOwnTrfFile as text
%        str2double(get(hObject,'String')) returns contents of edit_userOwnTrfFile as a double


% --- Executes during object creation, after setting all properties.
function edit_userOwnTrfFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_userOwnTrfFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_userOwnTrfBrowse.
function pushbutton_userOwnTrfBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_userOwnTrfBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)



function edit_s12userOwnFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_s12userOwnFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_s12userOwnFile as text
%        str2double(get(hObject,'String')) returns contents of edit_s12userOwnFile as a double


% --- Executes during object creation, after setting all properties.
function edit_s12userOwnFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_s12userOwnFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_s12userOwnBrowse.
function pushbutton_s12userOwnBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_s12userOwnBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)


function edit_oceanLoadingUserOwnFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingUserOwnFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_oceanLoadingUserOwnFile as text
%        str2double(get(hObject,'String')) returns contents of edit_oceanLoadingUserOwnFile as a double


% --- Executes during object creation, after setting all properties.
function edit_oceanLoadingUserOwnFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingUserOwnFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_oceanLoadingUserOwnBrowse.
function pushbutton_oceanLoadingUserOwnBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_oceanLoadingUserOwnBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)


function edit_oceanLoadingUserOwn_structFieldName_Callback(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingUserOwn_structFieldName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_oceanLoadingUserOwn_structFieldName as text
%        str2double(get(hObject,'String')) returns contents of edit_oceanLoadingUserOwn_structFieldName as a double


% --- Executes during object creation, after setting all properties.
function edit_oceanLoadingUserOwn_structFieldName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingUserOwn_structFieldName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_opoleloadUserOwnFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_opoleloadUserOwnFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_opoleloadUserOwnFile as text
%        str2double(get(hObject,'String')) returns contents of edit_opoleloadUserOwnFile as a double


% --- Executes during object creation, after setting all properties.
function edit_opoleloadUserOwnFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_opoleloadUserOwnFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_opoleloadUserOwnBrowse.
function pushbutton_opoleloadUserOwnBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_opoleloadUserOwnBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)



function edit_s12userOwn_structFieldName_Callback(hObject, eventdata, handles)
% hObject    handle to edit_s12userOwn_structFieldName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_s12userOwn_structFieldName as text
%        str2double(get(hObject,'String')) returns contents of edit_s12userOwn_structFieldName as a double


% --- Executes during object creation, after setting all properties.
function edit_s12userOwn_structFieldName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_s12userOwn_structFieldName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_opoleloadUserOwn_structFieldName_Callback(hObject, eventdata, handles)
% hObject    handle to edit_opoleloadUserOwn_structFieldName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_opoleloadUserOwn_structFieldName as text
%        str2double(get(hObject,'String')) returns contents of edit_opoleloadUserOwn_structFieldName as a double


% --- Executes during object creation, after setting all properties.
function edit_opoleloadUserOwn_structFieldName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_opoleloadUserOwn_structFieldName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% 
% 
% function edit_s12gsfcFile_Callback(hObject, eventdata, handles)
% % hObject    handle to edit_s12gsfcFile (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of edit_s12gsfcFile as text
% %        str2double(get(hObject,'String')) returns contents of edit_s12gsfcFile as a double
% 
% 
% % --- Executes during object creation, after setting all properties.
% function edit_s12gsfcFile_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edit_s12gsfcFile (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% 
% 
% % --- Executes on button press in pushbutton_s12gsfcBrowse.
% function pushbutton_s12gsfcBrowse_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton_s12gsfcBrowse (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)



function edit_s12viennaFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_s12viennaFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_s12viennaFile as text
%        str2double(get(hObject,'String')) returns contents of edit_s12viennaFile as a double


% --- Executes during object creation, after setting all properties.
function edit_s12viennaFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_s12viennaFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_s12viennaBrowse.
function pushbutton_s12viennaBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_s12viennaBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
browseFunction(hObject, handles)



function edit_vtrf2014datFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vtrf2014datFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vtrf2014datFile as text
%        str2double(get(hObject,'String')) returns contents of edit_vtrf2014datFile as a double


% --- Executes during object creation, after setting all properties.
function edit_vtrf2014datFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vtrf2014datFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_vtrf2014datBrowse.
function pushbutton_vtrf2014datBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_vtrf2014datBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)

% --------------------------------------------------------------------
function menu_atmo_set_Callback(hObject, eventdata, handles)
% hObject    handle to menu_atmo_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set all panels to invisible
makeAllPanelsInvisible(handles)

set(handles.uipanel_atmoLoading, 'Visible', 'on');

% --------------------------------------------------------------------
function menu_trf_set_Callback(hObject, eventdata, handles)
% hObject    handle to menu_trf_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set all panels to invisible
makeAllPanelsInvisible(handles)

set(handles.uipanel_trf, 'Visible', 'on');

% --------------------------------------------------------------------
function menu_mainInfoFiles_set_Callback(hObject, eventdata, handles)
% hObject    handle to menu_mainInfoFiles_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set all panels to invisible
makeAllPanelsInvisible(handles)

set(handles.uipanel_mainInfoFiles, 'Visible', 'on');

% --------------------------------------------------------------------
function menu_oceanLoading_set_Callback(hObject, eventdata, handles)
% hObject    handle to menu_oceanLoading_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set all panels to invisible
makeAllPanelsInvisible(handles)

set(handles.uipanel_oceanLoading, 'Visible', 'on');




% This function is executed on a right-click on an edit textbox
function rightClickOnEdit(hObject,eventdata,handles)

% if there is an entrie in the listbox at all
if ~isempty(get(hObject, 'String'))
    
    % create the menu at proper position
    curMousePosition=get(handles.figure_main, 'CurrentPoint');
    cmenu = uicontextmenu('Parent',handles.figure_main,'Position',curMousePosition);
%     cmenu = uicontextmenu;
    
    % create function handles
    %fHandles_openOPTfile=@{openOPTfile, handles};

    % create entries with the pointer to the callbacks
    item1 = uimenu(cmenu, 'Label', 'Open file', 'Callback', {@openfile,hObject,handles});
        
    % set visible
    set(cmenu, 'Visible', 'on');
end


function openfile(src,eventdata,hObject,handles)
% This function opens the OPT file (empty=create or open existing)
% of the currently selected session

% get current filename
curf=get(handles.(get(hObject, 'Tag')),'String');


% if exist: open
if exist(curf, 'file')
    if ispc
        openFailed=dos(['start notepad++ ',curf]);
        
        if openFailed==1
            dos(['start wordpad ',curf]);
        end
        %eval(['!wordpad ', wantedOPTfile])
    else
        system(['xterm -e ''vi ',curf '''']);
    end
% else: throw msgbox
else
    msgbox(sprintf('File %s\nnot found!',curf),'File not found','warn');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_vtrf2014datFile.
function edit_vtrf2014datFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_vtrf2014datFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rightClickOnEdit(hObject,eventdata,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_vieTrf13File.
function edit_vieTrf13File_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_vieTrf13File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rightClickOnEdit(hObject,eventdata,handles);

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_vievsTrfFile.
function edit_vievsTrfFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_vievsTrfFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rightClickOnEdit(hObject,eventdata,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_userOwnTrfFile.
function edit_userOwnTrfFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_userOwnTrfFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rightClickOnEdit(hObject,eventdata,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_oceanLoadingFES2004File.
function edit_oceanLoadingFES2004File_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingFES2004File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rightClickOnEdit(hObject,eventdata,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_oceanLoadingGOT00File.
function edit_oceanLoadingGOT00File_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingGOT00File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rightClickOnEdit(hObject,eventdata,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_oceanLoadingEOT08aFile.
function edit_oceanLoadingEOT08aFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingEOT08aFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rightClickOnEdit(hObject,eventdata,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_oceanLoadingTPXO72File.
function edit_oceanLoadingTPXO72File_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingTPXO72File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rightClickOnEdit(hObject,eventdata,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_oceanLoadingAG06File.
function edit_oceanLoadingAG06File_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingAG06File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rightClickOnEdit(hObject,eventdata,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_oceanLoadingUserOwnFile.
function edit_oceanLoadingUserOwnFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_oceanLoadingUserOwnFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rightClickOnEdit(hObject,eventdata,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_opoleloadcoefcmcorFile.
function edit_opoleloadcoefcmcorFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_opoleloadcoefcmcorFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rightClickOnEdit(hObject,eventdata,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_opoleloadUserOwnFile.
function edit_opoleloadUserOwnFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_opoleloadUserOwnFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rightClickOnEdit(hObject,eventdata,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_s12vandamFile.
function edit_s12vandamFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_s12vandamFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rightClickOnEdit(hObject,eventdata,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_s12gsfcFile.
function edit_s12gsfcFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_s12gsfcFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rightClickOnEdit(hObject,eventdata,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_s12viennaFile.
function edit_s12viennaFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_s12viennaFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rightClickOnEdit(hObject,eventdata,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_s12userOwnFile.
function edit_s12userOwnFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_s12userOwnFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rightClickOnEdit(hObject,eventdata,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_nsCodesFile.
function edit_nsCodesFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_nsCodesFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rightClickOnEdit(hObject,eventdata,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_antennaInfoFile.
function edit_antennaInfoFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_antennaInfoFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rightClickOnEdit(hObject,eventdata,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_eccdatFile.
function edit_eccdatFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_eccdatFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rightClickOnEdit(hObject,eventdata,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_maskFile.
function edit_maskFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_maskFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rightClickOnEdit(hObject,eventdata,handles);

function edit_ivstrf2014bFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ivstrf2014bFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ivstrf2014bFile as text
%        str2double(get(hObject,'String')) returns contents of edit_ivstrf2014bFile as a double


% --- Executes during object creation, after setting all properties.
function edit_ivstrf2014bFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ivstrf2014bFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_ivstrf2014bBrowse.
function pushbutton_ivstrf2014bBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ivstrf2014bBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
browseFunction(hObject, handles)


function edit_itrf2014File_Callback(hObject, eventdata, handles)
% hObject    handle to edit_itrf2014File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_itrf2014File as text
%        str2double(get(hObject,'String')) returns contents of edit_itrf2014File as a double


% --- Executes during object creation, after setting all properties.
function edit_itrf2014File_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_itrf2014File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_itrf2014Browse.
function pushbutton_itrf2014Browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_itrf2014Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)


% --- Executes on button press in pushbutton_dtrf2014Browse.
function pushbutton_dtrf2014Browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_dtrf2014Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
browseFunction(hObject, handles)

function edit_dtrf2014File_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dtrf2014File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dtrf2014File as text
%        str2double(get(hObject,'String')) returns contents of edit_dtrf2014File as a double


% --- Executes during object creation, after setting all properties.
function edit_dtrf2014File_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dtrf2014File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%Markus
% --- Executes on button press in pushbutton_gravdefBrowse.
function pushbutton_gravdefBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_gravdefBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
browseFunction(hObject, handles)



function edit_gravdefFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gravdefFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gravdefFile as text
%        str2double(get(hObject,'String')) returns contents of edit_gravdefFile as a double


% --- Executes during object creation, after setting all properties.
function edit_gravdefFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gravdefFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_gravdefFile.
function edit_gravdefFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_gravdefFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rightClickOnEdit(hObject,eventdata,handles);
