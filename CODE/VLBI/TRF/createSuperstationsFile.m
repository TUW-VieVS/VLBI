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

% Last Modified by GUIDE v2.5 20-Oct-2017 13:39:10

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

% define names of all files and save it in handles struct
% handles.data.allFileDescriptions={'nsCodes', 'antennaInfo', 'antenna', 'eccdat', 'blokq', 'itrf2005', ...
%         'vlbidiscont', 'itrf2008', 'itrf2014','dtrf2014', 'vtrf2008dat', 'vtrf2014dat', 'ivstrf2014b'...
%         'vieTrf13', 'mask', 'equip', 's12vienna','s12vandam', 's12gsfc', 's12userOwn',...
%         'oceanLoadingFES2004', 'oceanLoadingGOT00', 'oceanLoadingEOT08a', 'oceanLoadingTPXO72', 'oceanLoadingAG06', 'oceanLoadingUserOwn' ...
%         'opoleloadcoefcmcor', 'opoleloadUserOwn', 'vievsTrf', 'userOwnTrf'};
handles.data.allFileDescriptions={'nsCodes', 'antennaInfo', 'eccdat', 'itrf2005', ...
        'itrf2008', 'itrf2014','dtrf2014', 'vtrf2008dat', 'vtrf2014dat', 'ivstrf2014b'...
        'vieTrf13', 's12vienna','s12vandam', 's12gsfc', 's12userOwn',...
        'oceanLoadingFES2004', 'oceanLoadingGOT00', 'oceanLoadingEOT08a', 'oceanLoadingTPXO72', 'oceanLoadingAG06', 'oceanLoadingUserOwn' ...
        'opoleloadcoefcmcor', 'opoleloadUserOwn', 'vievsTrf', 'userOwnTrf'};
    
% handles.data.allFileNames={'ns-codes.txt', 'antenna-info.txt', 'antenna.dat', 'ECCDAT.ecc', 'blokq.dat', 'ITRF2005_VLBI.SSC.txt',...
%     'VLBI-DISCONT.txt', 'ITRF2008_VLBI.SSC.txt', 'ITRF2014-IVS-TRF.SNX', 'DTRF2014_VLBI.snx', 'VTRF2008.dat', 'VTRF2014_final.snx', ...
%     'IVS_TRF2014b.SSC.txt', 'VieTRF13.txt', 'mask.cat', 'equip.cat',...
%     's1_s2_s3_cm_noib_vlbi.dat','s12_cm_noib_vandam.dat', 's12_cm_noib_gsfc.dat',  's12userOwn.dat',...
%     'ocean_loading_FES2004.TXT', 'ocean_loading_GOT00.TXT', 'ocean_loading_EOT08a.TXT', 'ocean_loading_TPXO72.TXT', 'ocean_loading_AG06.TXT', 'ocean_loading_userOwn.TXT',...
%     'opoleloadcoefcmcor.txt', 'opoleloadUserOwn.txt', 'vievsTrf.txt', 'userOwnTrf.txt'};
handles.data.allFileNames={'ns-codes.txt', 'antenna-info.txt', 'ECCDAT.ecc', 'ITRF2005_VLBI.SSC.txt',...
    'ITRF2008_VLBI.SSC.txt', 'ITRF2014-IVS-TRF.SNX', 'DTRF2014_VLBI.snx', 'VTRF2008.dat', 'VTRF2014_final.snx', ...
    'IVS_TRF2014b.SSC.txt', 'VieTRF13.txt',...
    's1_s2_s3_cm_noib_vlbi.dat','s12_cm_noib_vandam.dat', 's12_cm_noib_gsfc.dat',  's12userOwn.dat',...
    'ocean_loading_FES2004.TXT', 'ocean_loading_GOT00.TXT', 'ocean_loading_EOT08a.TXT', 'ocean_loading_TPXO72.TXT', 'ocean_loading_AG06.TXT', 'ocean_loading_userOwn.TXT',...
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
    if length(guiOptLoaded{1})==(nFiles*3+1) % then I have enough files from default file
        for k=1:nFiles
            eval(['set(handles.radiobutton_', handles.data.allFileDescriptions{k},'Download, ''Value'', str2double(guiOptLoaded{2}{', num2str(3*k-2),'}));']);
            if str2double(guiOptLoaded{2}{3*k-2})==1
                eval(['set(handles.edit_', handles.data.allFileDescriptions{k},'File, ''Enable'', ''Off'');']);
                eval(['set(handles.pushbutton_', handles.data.allFileDescriptions{k},'Browse, ''Enable'', ''Off'');']);
            else
                eval(['set(handles.edit_', handles.data.allFileDescriptions{k},'File, ''Enable'', ''On'');']);
                eval(['set(handles.pushbutton_', handles.data.allFileDescriptions{k},'Browse, ''Enable'', ''On'');']);
            end
            eval(['set(handles.radiobutton_', handles.data.allFileDescriptions{k},'Specify, ''Value'', str2double(guiOptLoaded{2}{', num2str(3*k-1),'}));']);
            eval(['set(handles.edit_', handles.data.allFileDescriptions{k},'File, ''String'', guiOptLoaded{2}{', num2str(3*k),'});']);
        end

        % set output file and path to textbox
        set(handles.edit_createSuperstations, 'String', guiOptLoaded{2}{3*k+1});
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

% define descriptions
text={'ns-codes, download', 'ns-codes, specify', 'ns-codes, editTextbox', ...
    'antenna-info, download', 'antenna-info, specify', 'antenna-info, editTextBox',...
    'eccdat, download', 'eccdat, specify', 'eccdat, editTextBox',...
    'itrf2005, download', 'itrf2005, specify', 'itrf2005, editTextBox',...
    'itrf2008, download', 'itrf2008, specify', 'itrf2008, editTextBox',...
    'itrf2014, download', 'itrf2014, specify', 'itrf2014, editTextBox',...
    'dtrf2014, download', 'dtrf2014, specify', 'dtrf2014, editTextBox',...
    'vtrf2008dat, download', 'vtrf2008dat, specify', 'vtrf2008dat, editTextBox',...
    'vtrf2014dat, download', 'vtrf2014dat, specify', 'vtrf2014dat, editTextBox',...
    'ivstrf2014bdat, download', 'ivstrf2014b, specify', 'ivstrf2014b, editTextBox',...
    'vieTrf13, download', 'vieTrf13, specify', 'vieTrf13, editTextBox',...
    's12vienna, download', 's12vienna, specify', 's12vienna, editTextBox',...    
    's12vandam, download', 's12vandam, specify', 's12vandam, editTextBox',...
    's12gsfc, download', 's12gsfc, specify', 's12gsfc, editTextBox',...
    's12userOwn, download', 's12userOwn, specify', 's12userOwn, editTextbox',...
    'oceanLoadFES2004, download', 'oceanLoadFES2004, specify', 'oceanLoadFES2004, editTextBox',...
    'oceanLoadGOT00, download', 'oceanLoadGOT00, specify', 'oceanLoadGOT00, editTextBox',...
    'oceanLoadEOT08a, download', 'oceanLoadEOT08a, specify', 'oceanLoadEOT08a, editTextBox',...
    'oceanLoadTPXO72, download', 'oceanLoadTPXO72, specify', 'oceanLoadTPXO72, editTextBox',...
    'oceanLoadAG06, download', 'oceanLoadAG06, specify', 'oceanLoadAG06, editTextBox',...
    'oceanLoadUserOwn, download', 'oceanLoadUserOwn, specify', 'oceanLoadUserOwn, editTextBox',...
    'opoleloadcoefcmcor, download', 'opoleloadcoefcmcor, specify', 'opoleloadcoefcmcor, editTextBox',...
    'opoleloadUserOwn, download', 'opoleloadUserOwn, specify', 'opoleloadUserOwn, editTextBox',...
    'vievsTrf, download', 'vievsTrf, specify', 'vievsTrf, editTextBox',...
    'userOwnTrf, download', 'userOwnTrf, specify', 'userOwnTrf, editTextBox'};

% text={'ns-codes, download', 'ns-codes, specify', 'ns-codes, editTextbox', ...
%     'antenna-info, download', 'antenna-info, specify', 'antenna-info, editTextBox',...
%     'antenna, download', 'antenna, specify', 'antenna, editTextBox',...
%     'eccdat, download', 'eccdat, specify', 'eccdat, editTextBox',...
%     'blokq, download', 'blokq, specify', 'blokq, editTextBox',...
%     'itrf2005, download', 'itrf2005, specify', 'itrf2005, editTextBox',...
%     'vlbidiscont, download', 'vlbidiscont, specify', 'vlbidiscont, editTextBox',...
%     'itrf2008, download', 'itrf2008, specify', 'itrf2008, editTextBox',...
%     'itrf2014, download', 'itrf2014, specify', 'itrf2014, editTextBox',...
%     'dtrf2014, download', 'dtrf2014, specify', 'dtrf2014, editTextBox',...
%     'vtrf2008dat, download', 'vtrf2008dat, specify', 'vtrf2008dat, editTextBox',...
%     'vtrf2014dat, download', 'vtrf2014dat, specify', 'vtrf2014dat, editTextBox',...
%     'ivstrf2014bdat, download', 'ivstrf2014b, specify', 'ivstrf2014b, editTextBox',...
%     'vieTrf13, download', 'vieTrf13, specify', 'vieTrf13, editTextBox',...
%     'mask, download', 'mask, specify', 'mask, editTextBox',...
%     'equip, download', 'equip, specify', 'equip, editTextBox',...
%     's12vienna, download', 's12vienna, specify', 's12vienna, editTextBox',...    
%     's12vandam, download', 's12vandam, specify', 's12vandam, editTextBox',...
%     's12gsfc, download', 's12gsfc, specify', 's12gsfc, editTextBox',...
%     's12userOwn, download', 's12userOwn, specify', 's12userOwn, editTextbox',...
%     'oceanLoadFES2004, download', 'oceanLoadFES2004, specify', 'oceanLoadFES2004, editTextBox',...
%     'oceanLoadGOT00, download', 'oceanLoadGOT00, specify', 'oceanLoadGOT00, editTextBox',...
%     'oceanLoadEOT08a, download', 'oceanLoadEOT08a, specify', 'oceanLoadEOT08a, editTextBox',...
%     'oceanLoadTPXO72, download', 'oceanLoadTPXO72, specify', 'oceanLoadTPXO72, editTextBox',...
%     'oceanLoadAG06, download', 'oceanLoadAG06, specify', 'oceanLoadAG06, editTextBox',...
%     'oceanLoadUserOwn, download', 'oceanLoadUserOwn, specify', 'oceanLoadUserOwn, editTextBox',...
%     'opoleloadcoefcmcor, download', 'opoleloadcoefcmcor, specify', 'opoleloadcoefcmcor, editTextBox',...
%     'opoleloadUserOwn, download', 'opoleloadUserOwn, specify', 'opoleloadUserOwn, editTextBox',...
%     'vievsTrf, download', 'vievsTrf, specify', 'vievsTrf, editTextBox',...
%     'userOwnTrf, download', 'userOwnTrf, specify', 'userOwnTrf, editTextBox'};

% get max length of text
l=num2str(max(cellfun('length', text)));

nFiles=size(handles.data.allFileNames,2);

% write gui options to file
for k=1:nFiles
    fprintf(fid, ['# ', handles.data.allFileNames{k},' --------------\n%',l,'s | %1.0f\n%',l,'s | %1.0f\n%',l,'s | %s\n#\n'], ...
        text{3*k-2}, eval(['get(handles.radiobutton_', handles.data.allFileDescriptions{k},'Download, ''Value'')']),...
        text{3*k-1}, eval(['get(handles.radiobutton_', handles.data.allFileDescriptions{k},'Specify, ''Value'')']),...
        text{3*k},   eval(['get(handles.edit_', handles.data.allFileDescriptions{k},'File, ''String'')']) );
end

% write output path+file to file
fprintf(fid, ['\n#\n%', l,'s | %s\n'], 'outFile', get(handles.edit_createSuperstations, 'String'));

fclose(fid);


% --- Executes when selected object is changed in uipanel_ns.
function uipanel_ns_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_ns 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% if "donwload newest" is selected
if strcmp(get(eventdata.NewValue, 'Tag'), 'radiobutton_nsCodesDownload')
    set(handles.edit_nsCodesFile, 'Enable', 'off');
    set(handles.pushbutton_nsCodesBrowse, 'Enable', 'off');
else
    set(handles.edit_nsCodesFile, 'Enable', 'on');
    set(handles.pushbutton_nsCodesBrowse, 'Enable', 'on');
end

saveGuiDataToDisk(handles)


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

% --- Executes when selected object is changed in uipanel_antennaInfo.
function uipanel_antennaInfo_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_antennaInfo 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% if "download" button was chosen
if strcmp(get(eventdata.NewValue, 'Tag'), 'radiobutton_antennaInfoDownload')
    set(handles.edit_antennaInfoFile, 'Enable', 'Off');
    set(handles.pushbutton_antennaInfoBrowse, 'Enable', 'Off');
else
    set(handles.edit_antennaInfoFile, 'Enable', 'On');
    set(handles.pushbutton_antennaInfoBrowse, 'Enable', 'On');
end

saveGuiDataToDisk(handles)
    
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

% This function searches the 'neededFiles' folder and puts the names of the
% found files (if more than one found, simply the first is taken) into the
% corresponding textbox

% get path to be searched
searchPath='../TRF/neededFiles/';

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
        eval(['set(handles.edit_', handles.data.allFileDescriptions{k},'File, ''Enable'', ''On'');']);
        eval(['set(handles.pushbutton_', handles.data.allFileDescriptions{k},'Browse, ''Enable'', ''On'');']);
        
        % ... set 1 to specify and 0 to download ...
        eval(['set(handles.radiobutton_', handles.data.allFileDescriptions{k},'Download, ''Value'', 0);']);
        eval(['set(handles.radiobutton_', handles.data.allFileDescriptions{k},'Specify, ''Value'', 1);']);
        
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



% --- Executes on button press in pushbutton_vtrf2008datBrowse.
function pushbutton_vtrf2008datBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_vtrf2008datBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)



function edit_vtrf2008datFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vtrf2008datFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vtrf2008datFile as text
%        str2double(get(hObject,'String')) returns contents of edit_vtrf2008datFile as a double


% --- Executes during object creation, after setting all properties.
function edit_vtrf2008datFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vtrf2008datFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_itrf2008Browse.
function pushbutton_itrf2008Browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_itrf2008Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)


function edit_itrf2008File_Callback(hObject, eventdata, handles)
% hObject    handle to edit_itrf2008File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_itrf2008File as text
%        str2double(get(hObject,'String')) returns contents of edit_itrf2008File as a double


% --- Executes during object creation, after setting all properties.
function edit_itrf2008File_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_itrf2008File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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


% --- Executes on button press in pushbutton_itrf2005Browse.
function pushbutton_itrf2005Browse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_itrf2005Browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)



function edit_itrf2005File_Callback(hObject, eventdata, handles)
% hObject    handle to edit_itrf2005File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_itrf2005File as text
%        str2double(get(hObject,'String')) returns contents of edit_itrf2005File as a double


% --- Executes during object creation, after setting all properties.
function edit_itrf2005File_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_itrf2005File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_blokqBrowse.
function pushbutton_blokqBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_blokqBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)



function edit_blokqFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_blokqFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_blokqFile as text
%        str2double(get(hObject,'String')) returns contents of edit_blokqFile as a double


% --- Executes during object creation, after setting all properties.
function edit_blokqFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_blokqFile (see GCBO)
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


% --- Executes on button press in pushbutton_antennaBrowse.
function pushbutton_antennaBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_antennaBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)


function edit_antennaFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_antennaFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_antennaFile as text
%        str2double(get(hObject,'String')) returns contents of edit_antennaFile as a double


% --- Executes during object creation, after setting all properties.
function edit_antennaFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_antennaFile (see GCBO)
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


% --- Executes on button press in pushbutton_equipBrowse.
function pushbutton_equipBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_equipBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)



function edit_equipFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_equipFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_equipFile as text
%        str2double(get(hObject,'String')) returns contents of edit_equipFile as a double


% --- Executes during object creation, after setting all properties.
function edit_equipFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_equipFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_maskBrowse.
function pushbutton_maskBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_maskBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)


function edit_maskFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maskFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maskFile as text
%        str2double(get(hObject,'String')) returns contents of edit_maskFile as a double


% --- Executes during object creation, after setting all properties.
function edit_maskFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maskFile (see GCBO)
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


% --- Executes when selected object is changed in uipanel_eccdat.
function uipanel_eccdat_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_eccdat 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% if "download" button was chosen
if strcmp(get(eventdata.NewValue, 'Tag'), 'radiobutton_eccdatDownload')
    set(handles.edit_eccdatFile, 'Enable', 'Off');
    set(handles.pushbutton_eccdatBrowse, 'Enable', 'Off');
else
    set(handles.edit_eccdatFile, 'Enable', 'On');
    set(handles.pushbutton_eccdatBrowse, 'Enable', 'On');
end

saveGuiDataToDisk(handles)


% --- Executes when selected object is changed in uipanel_blokq.
function uipanel_blokq_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_blokq 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% if "download" button was chosen
if strcmp(get(eventdata.NewValue, 'Tag'), 'radiobutton_blokqDownload')
    set(handles.edit_blokqFile, 'Enable', 'Off');
    set(handles.pushbutton_blokqBrowse, 'Enable', 'Off');
else
    set(handles.edit_blokqFile, 'Enable', 'On');
    set(handles.pushbutton_blokqBrowse, 'Enable', 'On');
end

saveGuiDataToDisk(handles)


% --- Executes when selected object is changed in uipanel_itrf2005.
function uipanel_itrf2005_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_itrf2005 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% if "download" button was chosen
if strcmp(get(eventdata.NewValue, 'Tag'), 'radiobutton_itrf2005Download')
    set(handles.edit_itrf2005File, 'Enable', 'Off');
    set(handles.pushbutton_itrf2005Browse, 'Enable', 'Off');
else
    set(handles.edit_itrf2005File, 'Enable', 'On');
    set(handles.pushbutton_itrf2005Browse, 'Enable', 'On');
end

saveGuiDataToDisk(handles)


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


% --- Executes when selected object is changed in uipanel_itrf2008.
function uipanel_itrf2008_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_itrf2008 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% if "download" button was chosen
if strcmp(get(eventdata.NewValue, 'Tag'), 'radiobutton_itrf2008Download')
    set(handles.edit_itrf2008File, 'Enable', 'Off');
    set(handles.pushbutton_itrf2008Browse, 'Enable', 'Off');
else
    set(handles.edit_itrf2008File, 'Enable', 'On');
    set(handles.pushbutton_itrf2008Browse, 'Enable', 'On');
end

saveGuiDataToDisk(handles)


% --- Executes when selected object is changed in uipanel_vtrf2008dat.
function uipanel_vtrf2008dat_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_vtrf2008dat 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% if "download" button was chosen
if strcmp(get(eventdata.NewValue, 'Tag'), 'radiobutton_vtrf2008datDownload')
    set(handles.edit_vtrf2008datFile, 'Enable', 'Off');
    set(handles.pushbutton_vtrf2008datBrowse, 'Enable', 'Off');
else
    set(handles.edit_vtrf2008datFile, 'Enable', 'On');
    set(handles.pushbutton_vtrf2008datBrowse, 'Enable', 'On');
end

saveGuiDataToDisk(handles)


% --- Executes when selected object is changed in uipanel_mask.
function uipanel_mask_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_mask 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% if "download" button was chosen
if strcmp(get(eventdata.NewValue, 'Tag'), 'radiobutton_maskDownload')
    set(handles.edit_maskFile, 'Enable', 'Off');
    set(handles.pushbutton_maskBrowse, 'Enable', 'Off');
else
    set(handles.edit_maskFile, 'Enable', 'On');
    set(handles.pushbutton_maskBrowse, 'Enable', 'On');
end

saveGuiDataToDisk(handles)


% --- Executes when selected object is changed in uipanel_equip.
function uipanel_equip_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_equip 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% if "download" button was chosen
if strcmp(get(eventdata.NewValue, 'Tag'), 'radiobutton_equipDownload')
    set(handles.edit_equipFile, 'Enable', 'Off');
    set(handles.pushbutton_equipBrowse, 'Enable', 'Off');
else
    set(handles.edit_equipFile, 'Enable', 'On');
    set(handles.pushbutton_equipBrowse, 'Enable', 'On');
end

saveGuiDataToDisk(handles)


% --- Executes when selected object is changed in uipanel_s12vandam.
function uipanel_s12vandam_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_s12vandam 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% if "download" button was chosen
if strcmp(get(eventdata.NewValue, 'Tag'), 'radiobutton_s12vandamDownload')
    set(handles.edit_s12vandamFile, 'Enable', 'Off');
    set(handles.pushbutton_s12vandamBrowse, 'Enable', 'Off');
else
    set(handles.edit_s12vandamFile, 'Enable', 'On');
    set(handles.pushbutton_s12vandamBrowse, 'Enable', 'On');
end

saveGuiDataToDisk(handles)


% --- Executes when selected object is changed in uipanel_vlbidiscont.
function uipanel_vlbidiscont_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_vlbidiscont 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% if "download" button was chosen
if strcmp(get(eventdata.NewValue, 'Tag'), 'radiobutton_vlbidiscontDownload')
    set(handles.edit_vlbidiscontFile, 'Enable', 'Off');
    set(handles.pushbutton_vlbidiscontBrowse, 'Enable', 'Off');
else
    set(handles.edit_vlbidiscontFile, 'Enable', 'On');
    set(handles.pushbutton_vlbidiscontBrowse, 'Enable', 'On');
end

saveGuiDataToDisk(handles)


% --- Executes on button press in pushbutton_vlbidiscontBrowse.
function pushbutton_vlbidiscontBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_vlbidiscontBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

browseFunction(hObject, handles)



function edit_vlbidiscontFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vlbidiscontFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vlbidiscontFile as text
%        str2double(get(hObject,'String')) returns contents of edit_vlbidiscontFile as a double


% --- Executes during object creation, after setting all properties.
function edit_vlbidiscontFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vlbidiscontFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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


% --- Executes when selected object is changed in uipanel_s12gsfc.
function uipanel_s12gsfc_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_s12gsfc 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% if "download" button was chosen
if strcmp(get(eventdata.NewValue, 'Tag'), 'radiobutton_s12gsfcDownload')
    set(handles.edit_s12gsfcFile, 'Enable', 'Off');
    set(handles.pushbutton_s12gsfcBrowse, 'Enable', 'Off');
else
    set(handles.edit_s12gsfcFile, 'Enable', 'On');
    set(handles.pushbutton_s12gsfcBrowse, 'Enable', 'On');
end

saveGuiDataToDisk(handles)


% --- Executes when selected object is changed in uipanel_vievsTrf.
function uipanel_vievsTrf_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_vievsTrf 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% if "download" button was chosen
if strcmp(get(eventdata.NewValue, 'Tag'), 'radiobutton_vievsTrfDownload')
    set(handles.edit_vievsTrfFile, 'Enable', 'Off');
    set(handles.pushbutton_vievsTrfBrowse, 'Enable', 'Off');
else
    set(handles.edit_vievsTrfFile, 'Enable', 'On');
    set(handles.pushbutton_vievsTrfBrowse, 'Enable', 'On');
end

saveGuiDataToDisk(handles)


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
for k=1:nInFiles
    % check if should be downloaded (Check the "Values" of the according radio buttions!)
    if get(eval(['handles.radiobutton_', handles.data.allFileDescriptions{k}, 'Download']), 'Value')
        switch handles.data.allFileDescriptions{k}
            case 'nsCodes'
                url='ftp://ivscc.gsfc.nasa.gov/pub/control/ns-codes.txt';
                urlwrite(url, '../TRF/neededFiles/ns-codes.txt');
                inFiles(1).name='../TRF/neededFiles/ns-codes.txt';
                fprintf('File %s was downloaded\n', handles.data.allFileNames{k});
            case 'antennaInfo'
                url='http://vlbi.geod.uni-bonn.de/IVS-AC/Conventions/antenna-info.txt';
                urlwrite(url, '../TRF/neededFiles/antenna-info.txt');
                inFiles(2).name='../TRF/neededFiles/antenna-info.txt';
                fprintf('File %s was downloaded\n', handles.data.allFileNames{k});
            case 'antenna'
                url='http://gemini.gsfc.nasa.gov/solve_save/antenna.dat';
                urlwrite(url, '../TRF/neededFiles/antenna.dat');
                inFiles(3).name='../TRF/neededFiles/antenna.dat';
                fprintf('File %s was downloaded\n', handles.data.allFileNames{k});
            case 'eccdat'
                url='http://gemini.gsfc.nasa.gov/solve_save/ECCDAT.ecc';
                urlwrite(url, '../TRF/neededFiles/ECCDAT.ecc');
                inFiles(4).name='../TRF/neededFiles/ECCDAT.ecc';
                fprintf('File %s was downloaded\n', handles.data.allFileNames{k});
            case 'blokq'
                url='http://gemini.gsfc.nasa.gov/apriori_files/blokq.dat';
                urlwrite(url, '../TRF/neededFiles/blokq.dat');
                inFiles(5).name='../TRF/neededFiles/blokq.dat';
                fprintf('File %s was downloaded\n', handles.data.allFileNames{k});
            case 'itrf2005'

                url='http://itrf.ensg.ign.fr/ITRF_solutions/2005/doc/ITRF2005_VLBI.SSC.txt';
                urlwrite(url, '../TRF/neededFiles/ITRF2005_VLBI.SSC.txt');
                inFiles(6).name='../TRF/neededFiles/ITRF2005_VLBI.SSC.txt';
                fprintf('File %s was downloaded\n', handles.data.allFileNames{k});
            case 'vlbidiscont'
                url='http://vlbi.geod.uni-bonn.de/IVS-AC/data/VLBI-DISCONT.txt';
                urlwrite(url, '../TRF/neededFiles/VLBI-DISCONT.txt');
                inFiles(7).name='../TRF/neededFiles/VLBI-DISCONT.txt';
                fprintf('File %s was downloaded\n', handles.data.allFileNames{k});   
            case 'itrf2008'
                url='http://itrf.ensg.ign.fr/ITRF_solutions/2008/doc/ITRF2008_VLBI.SSC.txt';
                urlwrite(url, '../TRF/neededFiles/ITRF2008_VLBI.SSC.txt');
                inFiles(8).name='../TRF/neededFiles/ITRF2008_VLBI.SSC.txt';
                fprintf('File %s was downloaded\n', handles.data.allFileNames{k});
            case 'vtrf2008dat'
                fprintf('Where can the file VTRF2008.dat be downloaded?\n');
                keyboard;
            case 'vtrf2014dat'
                fprintf('Where can the file vtrf2014 be downloaded?\n');
                keyboard;
            case 'ivstrf2014b'
                url='http://ccivs.bkg.bund.de/combination/QUAT/TRF/IVS_TRF2014b.SSC.txt';
                urlwrite(url, '../TRF/neededFiles/IVS_TRF2014b.SSC.txt');
                inFiles(11).name='../TRF/neededFiles/IVS_TRF2014b.SSC.txt';
                fprintf('File %s was downloaded\n', handles.data.allFileNames{k});
            case 'vieTrf13'
                fprintf('The vieTrf file has to be created manually (Hana!) and cannout be downloaded!\n');
                keyboard;
            case 'mask'
                fprintf('Where can the file mask.cat be downloaded?\n');
                keyboard;
                %                 1 'nsCodes', 2'antennaInfo', '3antenna', 4'eccdat', 5'blokq', 6'itrf2005', ...
%         7'vlbidiscont', 8'itrf2008', 9'vtrf2008dat', 10'vtrf2014dat', 11'ivstrf2014b'...
%         12'vieTrf13', 13'mask', 14'equip', 15's12vienna','16s12vandam', '17s12gsfc', '18s12userOwn',...
%         19'oceanLoadingFES2004', 20'oceanLoadingGOT00', 21'oceanLoadingEOT08a', 22'oceanLoadingTPXO72', 23'oceanLoadingAG06', 24'oceanLoadingUserOwn' ...
%         25'opoleloadcoefcmcor', 'opoleloadUserOwn', 'vievsTrf', 'userOwnTrf'
            case 'equip'
                fprintf('Where can the file equip.cat be downloaded?\n');
                keyboard;
            case 's12vienna'
                fprintf('Can the file s1_s2_s3_cm_noib_vlbi.dat be downloaded?\nIt might be a ''VieVS internal'' file.\n');
                keyboard;
            case 's12vandam'
                fprintf('Can the file s12_cm_noib_vandam.dat be downloaded?\nIt might be a ''VieVS internal'' file.\n');
                keyboard;
            case 's12gsfc'
                fprintf('Can the file s12_cm_noib_gsfc.dat be downloaded?\nIt might be a ''VieVS internal'' file.\n');
                keyboard;
            case 'opoleloadcoefcmcor'
                url='ftp://tai.bipm.org/iers/conv2010/chapter7/opoleloadcoefcmcor.txt.gz';
                urlwrite(url, '../TRF/neededFiles/opoleloadcoefcmcor.txt.gz');
                fprintf('File ../TRF/neededFiles/opoleloadcoefcmcor.txt.gz was downloaded,\nplease unzip manually to same folder and type ''return'' and hit enter!\nThanks!\n');
                inFiles(25).name='../neededFiles/opoleloadcoefcmcor.txt';
            case 'oceanLoadingFES2004'
                fprintf('Where can the file oceanLoadingFES2004 be downloaded?\n');
                keyboard;
            case 'vievsTrf'
                fprintf('The VieVS trf file has to be created manually and cannout be downloaded!\n');
                keyboard;
        end
                
    else
        % 'Specify' radiobutton was checked -> get textbox entry for name
        inFiles(k).name=eval(['get(handles.edit_', handles.data.allFileDescriptions{k},'File, ''String'');']);
        if isfield(handles, eval(['''edit_', handles.data.allFileDescriptions{k},'_structFieldName''']))
            inFiles(k).fieldName=eval(['get(handles.edit_', handles.data.allFileDescriptions{k},'_structFieldName, ''String'');']);
        end
    end
    
    %inFiles(k).path=eval(['get(handles.edit_', handles.data.allFileDescriptions{k}, 'File, ''String'');']);
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

% --- Executes when selected object is changed in uipanel_vieTrf13.
function uipanel_vieTrf13_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_vieTrf13 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% if "donwload newest" is selected
if strcmp(get(eventdata.NewValue, 'Tag'), 'radiobutton_vieTrf13Download')
    set(handles.edit_vieTrf13File, 'Enable', 'off');
    set(handles.pushbutton_vieTrf13Browse, 'Enable', 'off');
else
    set(handles.edit_vieTrf13File, 'Enable', 'on');
    set(handles.pushbutton_vieTrf13Browse, 'Enable', 'on');
end

saveGuiDataToDisk(handles)



function edit_userOwnTrfFile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_userOwnTrfFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_userOwnTrfFile as text
%        str2double(get(hObject,'String')) returns contents of edit_userOwnTrfFile as a double


% --- Executes during object creation, after setting all properties.
function edit_userOwnTrfFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_userOwnTrfFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

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


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_itrf2005File.
function edit_itrf2005File_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_itrf2005File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rightClickOnEdit(hObject,eventdata,handles);


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
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_itrf2008File.
function edit_itrf2008File_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_itrf2008File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rightClickOnEdit(hObject,eventdata,handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_vtrf2008datFile.
function edit_vtrf2008datFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_vtrf2008datFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

rightClickOnEdit(hObject,eventdata,handles);


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
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_vlbidiscontFile.
function edit_vlbidiscontFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_vlbidiscontFile (see GCBO)
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
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_antennaFile.
function edit_antennaFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_antennaFile (see GCBO)
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
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_blokqFile.
function edit_blokqFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_blokqFile (see GCBO)
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


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over edit_equipFile.
function edit_equipFile_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to edit_equipFile (see GCBO)
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
