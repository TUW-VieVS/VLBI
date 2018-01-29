% 05.09.2016 - Matthias Schartner: creator
% 13.09.2016 - Matthias Schartner: 	station selection added
%									now uses tagalong_star
%									possibility to change starttime for manual scheduling
% 02.09.2016 - Matthias Schartner: 	bugfixes, new parameters and improvements
% 29.05.2017 - Matthias Schartner:  bugfixes


function varargout = sched_manual(varargin)

% Last Modified by GUIDE v2.5 09-Feb-2017 10:56:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sched_manual_OpeningFcn, ...
                   'gui_OutputFcn',  @sched_manual_OutputFcn, ...
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


% --- Executes just before sched_manual is made visible.
function sched_manual_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sched_manual (see VARARGIN)
handles.figure1.Name='VieVS SCHED manually';
handles.stations = varargin{1};
handles.twin = varargin{2};
handles.source = varargin{3};
% change name to commonname if it exists
commonnames = {handles.source.commoname};
boolcommon = ~strcmp(commonnames,'        ');
idx = find(boolcommon);
for i = 1:length(idx)
    handles.source(idx(i)).name = commonnames{idx(i)};
end

% save defining flag
defining = {handles.source.info};
defining = strfind(defining,'ICRF2 def');
defining = ~cellfun(@isempty,defining);
for i = 1:length(handles.source)
    handles.source(i).defining = defining(i);
end

handles.obsmode = varargin{4};
handles.sched = struct();
PARA = varargin{5};
if ~PARA.STARMODE
    handles.uipanel_STARmode.Visible = 'off';
end

handles.obs = struct();
stanum = length(handles.stations);
for ista = 1 : stanum
    handles.obs.staobs(ista).ifirst   = 1;
    handles.obs.staobs(ista).endmjd   = PARA.startmjd;
    handles.obs.staobs(ista).az       = handles.stations(ista).lim11;
    handles.obs.staobs(ista).el       = handles.stations(ista).lim12;
    handles.obs.staobs(ista).ha       = handles.stations(ista).lim11;
    handles.obs.staobs(ista).dc       = handles.stations(ista).lim12;
    handles.obs.staobs(ista).ncat     = 0;
    handles.obs.staobs(ista).cattn(1) = PARA.startmjd;
    handles.obs.staobs(ista).catsn(1) = 0;
end
% initialize 'srcobs'
srcnum = length(handles.source);
for isrc = 1 : srcnum
    handles.obs.srcobs(isrc).obsmjd = 0.0;
    handles.obs.srcobs(isrc).nscan = 0.0;
end
handles.obs.srcat = 0;

handles.startsched = PARA.startmjd;
handles.endsched = PARA.endmjd;
handles.time = PARA.startmjd;
handles.PARA = PARA;
handles.edit_SNR_X_Band.String = PARA.MIN_SNR(1);
handles.edit_SNR_S_Band.String = PARA.MIN_SNR(2);
handles.edit_MIN_SRCRP.String = PARA.MIN_SRCRP;
handles.edit_MAXSLEWTIME.String = PARA.MAXSLEWTIME;
handles.edit_MAX_WAIT.String = PARA.MAX_WAIT;
handles.edit_MAX_SCAN.String = PARA.MAX_SCAN;
handles.edit_MIN_SCAN.String = PARA.MIN_SCAN;
handles.edit_SUNDIST.String = PARA.MIN_SUNDIST*180/pi;
handles.edit_CUTOFEL.String = PARA.MIN_CUTEL*180/pi;
handles.edit_CADENCE.String = PARA.CADENCE;
handles.edit_weight_scanendtime.String = PARA.WEIGHT_SCAN_END_TIME;
handles.edit_weight_skyCoverage.String = PARA.WEIGHT_SKY_COVERAGE;
handles.edit_weight_numberOfObservations.String = PARA.WEIGHT_NUMBER_OF_OBS;
handles.edit_min_stanum.String = PARA.MIN_STANUM;
handles.edit_min_stascan.String = PARA.MIN_STASCAN;
handles.edit_min_stanum_fi.String = PARA.MIN_STANUM_FI;

if PARA.SUBNETTING == 1
    handles.radiobutton_subnetting_no.Value = 1;
    handles.radiobutton_subnetting_yes.Value = 0;
else
    handles.radiobutton_subnetting_yes.Value = 1;
    handles.radiobutton_subnetting_n0.Value = 0;
end

handles.listbox_useStation.String = {handles.stations.name};

[year, month, day, hour, minute, second] = tymdhms(handles.endsched);
[datestr] = tdatestr(year, month, day, hour, minute, second);
handles.edit_dateTime.String = datestr;
handles.edit_mjd.String = handles.endsched;
handles.edit_duration.String = PARA.duration;

[year, month, day, hour, minute, second] = tymdhms(handles.startsched);
[datestr] = tdatestr(year, month, day, hour, minute, second);
handles.text_timeDate.String = datestr;
handles.text_timeMJD.String = handles.startsched;
handles.text_timeH.String = 0;

handles.popupmenu_STRONGANT.String = {handles.stations.name};
% Choose default command line output for sched_manual
handles.output = hObject;

stanames = {handles.stations.name};
handles.popupmenu_station.String = stanames;
handles.uitable_sched.Data(1:length(stanames),1) = stanames;
handles.uitable_sched.Data(1:length(stanames),1) = stanames;

srcnames = {handles.source.name};
[~,idx] = sort(srcnames);
srcnames = srcnames(idx);
handles.whatsupSrcNameIdx = idx;
handles.listbox_sources.String = srcnames;

cnames = {handles.stations.po};
cnames(2:length(cnames)+1)=cnames;
cnames(1)={'names'};
cnames(end+1) = {'total'};
handles.uitable_whatsup.ColumnName = cnames;
cw = ones(size(cnames))*30;
cw(1)=60;
cw(end)=60;
handles.uitable_whatsup.ColumnWidth = num2cell(cw);
cf = cnames;
cf(1)={'char'};
cf(2:end-1)={'logical'};
cf(end) = {'numeric'};
handles.uitable_whatsup.ColumnFormat=cf;

handles.uitable_sched.Data(:,2)={'none'};

axes(handles.axes_skyCoverage)
polar([0:0.1:2*pi 2*pi],90*ones(1,64)/90,'k');
hold on
h = polar([0:0.1:2*pi 2*pi],30*ones(1,64)/90);
set(h,'Color',[0.8725    0.8725    0.8725])
h = polar([0:0.1:2*pi 2*pi],60*ones(1,64)/90);
set(h,'Color',[0.8725    0.8725    0.8725])
h = polar([0:0.1:2*pi 2*pi],60*ones(1,64)/90);
set(h,'Color',[0.8725    0.8725    0.8725])
hold off
h = findall(gca,'type','line');
% delete(h(end));
h = findall(gca,'type','text');
delete(h([13,14]));
view([90 -90])

jpl = load_jpl_Earth('jpl_421');
handles.jpl = jpl;
handles.listbox_scans.String = {};

handles.alreadyCalced_recommended = 0;
handles.alreadyCalced_source = 0;
handles.alreadyCalced_station = 0;
handles.alreadyCalced_whatsup = 0;

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes sched_manual wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sched_manual_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.figure1;



% --- Executes during object creation, after setting all properties.
function popupmenu_mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu_mode.
function popupmenu_mode_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
if strcmp(contents{get(hObject,'Value')},'auto schedule')
    handles.uipanel_autoSched.Visible = 'on';
    handles.uipanel_manualSched.Visible = 'off';
    handles.uipanel_continue.Visible = 'off';
else
    handles.uipanel_autoSched.Visible = 'off';
    handles.uipanel_manualSched.Visible = 'on';
    handles.uipanel_continue.Visible = 'off';
    handles.PARA.STARMODE = 0;
    updateTable(hObject, eventdata, handles);
    handles.popupmenu_mode.Enable = 'on';
    for i = 1:length(handles.listbox_scans.String)
        listbox_scans_Callback(hObject, eventdata, handles);
    end
end

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_mode


% --- Executes on button press in pushbutton_startAutoSched.
function pushbutton_startAutoSched_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_startAutoSched (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

station = handles.stations;
twin = handles.twin;
source = handles.source;
obsmode = handles.obsmode;
PARA = handles.PARA;

if strcmp(handles.uibuttongroup_strategy.SelectedObject.String,'Source-based strategy')
    PARA.OPTIMIZATION = 1;
    PARA.SRCNUM = str2double(handles.popupmenu_nSourceBased.String{handles.popupmenu_nSourceBased.Value});
else
    PARA.OPTIMIZATION = 2;
    PARA.DISTRIBUTE = handles.checkbox_distributeStationBased.Value;
end

% Fillin Modes
if handles.radiobutton_fi1.Value && handles.radiobutton_fi2.Value
    PARA.FILLINMODE = 12;
elseif handles.radiobutton_fi1.Value
    PARA.FILLINMODE = 1;
elseif handles.radiobutton_fi2.Value
    PARA.FILLINMODE = 2;
else
    PARA.FILLINMODE = 0;
end


%STAR mode
PARA.STARMODE = handles.checkbox_STARmode.Value;
PARA.STRONGANT = handles.popupmenu_STRONGANT.String{handles.popupmenu_STRONGANT.Value};
PARA.CADENCE = str2double(handles.edit_CADENCE.String);

% starttime
PARA.startmjd = handles.time;

% endtime
so = handles.uibuttongroup_schedUntil.SelectedObject.String;
if strcmp(so,'date/time')
    startT = textscan(handles.edit_dateTime.String,'%f-%f-%f %f:%f:%f');
    PARA.endmjd = date2mjd(cell2mat(startT));
    PARA.duration = (PARA.endmjd-PARA.startmjd)*24;
elseif strcmp(so,'mjd')
    PARA.endmjd = min([str2double(handles.edit_mjd.String) handles.endsched]);
    PARA.duration = (PARA.endmjd-PARA.startmjd)*24;
elseif strcmp(so,'duration [h]')
    PARA.endmjd = min([PARA.startmjd+str2double(handles.edit_duration.String)/24 handles.endsched]);
    PARA.duration = (PARA.endmjd-PARA.startmjd)*24;
else
    PARA.endmjd = handles.endsched;
    PARA.duration = (PARA.endmjd-PARA.startmjd)*24;
end

if (PARA.OPTIMIZATION == 1)
    [~, srcat, catpair] = isrcseq(source, PARA);
    
    % grid information
    if ~isequal(handles.obs.srcat,srcat)
        handles.obs.srcn(1:length(srcat)) = 0;
        handles.obs.icatpair = 0;
    end
else
    srcat = 0;
    catpair = 0;
    handles.obs.srcn   = 0.0;
    handles.obs.icatpair = 0.0;
end

for i = 1:length(handles.listbox_dontUseStation.String)
    sta = handles.listbox_dontUseStation.String(i);
    staid = find(strcmp(sta,{station.name}));
    station(staid).downum = station(staid).downum+1;
    if station(staid).downum == 1
        station(staid).downstart = PARA.startmjd;
        station(staid).downend = PARA.endmjd;
    else
        station(staid).downstart = [[station(staid).downstart] PARA.startmjd];
        station(staid).downend   = [[station(staid).downend] PARA.endmjd];
    end
end

PARA.SRCFRINGE = '********';
fprintf('\n\n----------------------------------------------\n');
fprintf('        #### START auto-sched ####\n')
[year, month, day, hour, minute, second] = tymdhms(PARA.startmjd);
[datestr] = tdatestr(year, month, day, hour, minute, second);
fprintf('start: %s (mjd: %.5f) \n',datestr,PARA.startmjd)
[year, month, day, hour, minute, second] = tymdhms(PARA.endmjd);
[datestr] = tdatestr(year, month, day, hour, minute, second);
fprintf('until: %s (mjd: %.5f) \n',datestr,PARA.endmjd)
fprintf('       %.4f hours \n',PARA.duration)
fprintf('----------------------------------------------\n');
handles.popupmenu_mode.Enable = 'off';
set( findall(handles.uipanel_autoSched, '-property', 'Enable'), 'Enable', 'off');
drawnow
[sched,obs] = vie_sched_src(station, twin, source, obsmode, srcat, catpair, PARA, handles.obs, handles.jpl);   
if PARA.STARMODE
   [sched] = vie_sched_tagalong_star(station, source, handles.obsmode, sched, PARA); 
end
handles.popupmenu_mode.Enable = 'on';
set( findall(handles.uipanel_autoSched, '-property', 'Enable'), 'Enable', 'on');
uibuttongroup_strategy_SelectionChangedFcn(handles.uibuttongroup_strategy.SelectedObject, eventdata, handles);
uibuttongroup_schedUntil_SelectionChangedFcn(handles.uibuttongroup_schedUntil.SelectedObject, eventdata, handles);
checkbox_STARmode_Callback(handles.checkbox_STARmode, eventdata, handles);
if isempty(fieldnames(handles.sched))
    handles.sched = [sched];
else
    handles.sched = [handles.sched, sched];
end
handles.obs = obs;
handles.time = max([obs.staobs.endmjd]);


[year, month, day, hour, minute, second] = tymdhms(handles.time);
[datestr] = tdatestr(year, month, day, hour, minute, second);
handles.text_timeDate.String = datestr;
handles.text_timeMJD.String = sprintf('%.3f',handles.time);
handles.text_timeH.String = (handles.time-handles.startsched)*24;
fprintf('\n\n----------------------------------------------\n')
fprintf('        #### END auto-sched ####\n')
fprintf('end: %s (mjd: %.5f) \n',datestr,handles.time)
fprintf('----------------------------------------------\n')
handles.alreadyCalced_recommended = 0;
handles.alreadyCalced_source = 0;
handles.alreadyCalced_station = 0;
handles.alreadyCalced_whatsup = 0;

% Update handles structure
guidata(hObject, handles);

if handles.time > handles.endsched
    sched = handles.sched;
    save('..\DATA\LEVEL5\sched_manually.mat','sched')
    close(handles.figure1);
end


% --- Executes during object creation, after setting all properties.
function popupmenu_nSourceBased_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_nSourceBased (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uibuttongroup_strategy.
function uibuttongroup_strategy_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup_strategy 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(hObject.String,'Source-based strategy')
    handles.popupmenu_nSourceBased.Enable = 'on';
    handles.checkbox_distributeStationBased.Enable = 'off';
    handles.text7.Enable = 'on';
    handles.checkbox_distributeStationBased.Value = 0;
    handles.checkbox_STARmode.Value = 0;
    checkbox_STARmode_Callback(handles.checkbox_STARmode, eventdata, handles)
    handles.checkbox_STARmode.Enable = 'off';
else 
    handles.popupmenu_nSourceBased.Enable = 'off';
    handles.checkbox_distributeStationBased.Enable = 'on';
    handles.text7.Enable = 'off';
    handles.checkbox_STARmode.Enable = 'on';
end


% --- Executes during object creation, after setting all properties.
function edit_SNR_S_Band_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SNR_S_Band (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_SNR_X_Band_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SNR_X_Band (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in uibuttongroup_schedUntil.
function uibuttongroup_schedUntil_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup_schedUntil 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(hObject.String,'date/time')
    handles.edit_dateTime.Enable = 'on';
    handles.edit_mjd.Enable = 'off';
    handles.edit_duration.Enable = 'off';
elseif  strcmp(hObject.String,'mjd')
    handles.edit_dateTime.Enable = 'off';
    handles.edit_mjd.Enable = 'on';
    handles.edit_duration.Enable = 'off';
elseif  strcmp(hObject.String,'duration [h]')
    handles.edit_dateTime.Enable = 'off';
    handles.edit_mjd.Enable = 'off';
    handles.edit_duration.Enable = 'on';
else
    handles.edit_dateTime.Enable = 'off';
    handles.edit_mjd.Enable = 'off';
    handles.edit_duration.Enable = 'off';
end


% --- Executes during object creation, after setting all properties.
function edit_duration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_mjd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mjd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_dateTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dateTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_MIN_SRCRP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MIN_SRCRP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_MAXSLEWTIME_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MAXSLEWTIME (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_MAX_WAIT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MAX_WAIT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_MAX_SCAN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MAX_SCAN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_MIN_SCAN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MIN_SCAN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_SUNDIST_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SUNDIST (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_CUTOFEL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_CUTOFEL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_STARmode.
function checkbox_STARmode_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_STARmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value')
    handles.text34.Enable = 'on';
    handles.popupmenu_STRONGANT.Enable = 'on';
    handles.text35.Enable = 'on';
    handles.edit_CADENCE.Enable = 'on';
    handles.radiobutton_fi1.Value = 0;
    handles.radiobutton_fi2.Value = 0;
else
    handles.text34.Enable = 'off';
    handles.popupmenu_STRONGANT.Enable = 'off';
    handles.text35.Enable = 'off';
    handles.edit_CADENCE.Enable = 'off';
end
% Hint: get(hObject,'Value') returns toggle state of checkbox_STARmode


% --- Executes during object creation, after setting all properties.
function popupmenu_STRONGANT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_STRONGANT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function edit_CADENCE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_CADENCE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes on selection change in popupmenu_table.
function popupmenu_table_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateTable(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_table contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_table


% --- Executes during object creation, after setting all properties.
function popupmenu_table_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_station.
function popupmenu_station_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_station (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
print_skyplot(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_station contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_station


% --- Executes during object creation, after setting all properties.
function popupmenu_station_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_station (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function updateTable(hObject, eventdata, handles)
if strcmp(handles.popupmenu_table.String(handles.popupmenu_table.Value),'stations')
    handles.uipanel_staobs.Visible = 'on';
    handles.uipanel_srcobs.Visible = 'off';
    handles.uipanel_whatsup.Visible = 'off';
    handles.uipanel_recommended.Visible = 'off';
    
    if ~handles.alreadyCalced_station
        handles.uitable_staobs.Data = {};
        
        n = length(handles.stations);
        storage(1:n,1)={handles.stations.name};
        az = wrapTo360(round(([handles.obs.staobs.az])*180/pi));
        storage(1:n,2)=num2cell(az);
        el = round([handles.obs.staobs.el]*180/pi);
        storage(1:n,3)=num2cell(el);
        endmjd = [handles.obs.staobs.endmjd];
    %     handles.uitable_staobs.Data(1:n,4)=num2cell(endmjd);
        for i = 1:n
            [year, month, day, hour, minute, second] = tymdhms(endmjd(i));
            [datestr] = tdatestr(year, month, day, hour, minute, second);
            storage(i,4) = {['  ' datestr]};
        end

        sched = handles.sched;
        if isfield(handles.sched,'scan')
            allScans = [handles.sched.scan];
               % save observation individually to improve performance
            allStanum = [allScans.nsta];
            allStart = [allScans.startmjd];
            allSource = [allScans.srcid];
            allnSta = [allScans.nsta];
            allnObs = (allnSta.*(allnSta-1))/2;
            allStar = [handles.source(allSource).star];
            alldefining = [handles.source(allSource).defining];
            handles.allstaobs = [];
            handles.allstaobs = [allScans.sta];
            c = 1;
            for i = 1:length(allStanum)
                for j = 1:allStanum(i)
                    handles.allstaobs(c).startmjd=allStart(i);
                    handles.allstaobs(c).srcid = allSource(i);
                    handles.allstaobs(c).star = allStar(i);
                    handles.allstaobs(c).defining = alldefining(i);
                    handles.allstaobs(c).srcname = [' ' handles.source(allSource(i)).name];
                    c = c+1;
                end
            end
            
            print_skyplot(hObject, eventdata, handles);
        else
            handles.allstaobs = [];
        end
        set(handles.uitable_staobs,'Data',storage);
        handles.alreadyCalced_station = 1;
    end
    
elseif strcmp(handles.popupmenu_table.String(handles.popupmenu_table.Value),'sources')
    handles.uipanel_staobs.Visible = 'off';
    handles.uipanel_srcobs.Visible = 'on';
    handles.uipanel_whatsup.Visible = 'off';
    handles.uipanel_recommended.Visible = 'off';
    
    if ~handles.alreadyCalced_source
    
        handles.uitable_srcobs.Data = {};
        storage = cell(length(handles.source),8);

        n = length(handles.source);
        storage(1:n,1)={handles.source.name};
        ra = round([handles.source.ra]*180/pi);
        storage(1:n,2)=num2cell(ra);
        de = round([handles.source.de]*180/pi);
        storage(1:n,3)=num2cell(de);
        storage(1:n,4)=num2cell([handles.source.defining]);
        storage(1:n,5)=num2cell(logical([handles.source.star]));
        nscan = [handles.obs.srcobs.nscan];
        storage(1:n,6)=num2cell(nscan);
        mjd = [handles.obs.srcobs.obsmjd];
        for i = 1:n
            if nscan(i)>0
                storage(i,8)={mjd(i)};
                [year, month, day, hour, minute, second] = tymdhms(mjd(i));
                [datestr] = tdatestr(year, month, day, hour, minute, second);
                storage(i,7)={['  ' datestr]};
            end
        end
        set(handles.uitable_srcobs,'Data',storage);
        handles.alreadyCalced_source = 1;
    end
    
elseif strcmp(handles.popupmenu_table.String(handles.popupmenu_table.Value),sprintf('what''s up'))
    handles.uipanel_staobs.Visible = 'off';
    handles.uipanel_srcobs.Visible = 'off';
    handles.uipanel_whatsup.Visible = 'on';
    handles.uipanel_recommended.Visible = 'off';
    
    if ~handles.alreadyCalced_whatsup
        handles.uitable_whatsup.Data = {};
        storage = cell(length(handles.source),length(handles.stations));
        srcnames = {handles.source.name};
        storage(1:length(srcnames),1)=srcnames;
        source = handles.source;
        station = handles.stations;
        endmjd = [handles.obs.staobs.endmjd];

        llh = [station.llh];
        lon = llh(1:3:end);
        lat = llh(2:3:end);
        idx = handles.whatsupSrcNameIdx;

        for isrc = 1:length(source)
            ra = source(isrc).ra;
            de = source(isrc).de;
            [az, el, ha, dc] = zazel_s(endmjd, lon, lat, ra, de);
            tot = 0;
            for ista = 1 : length(station)
                [lup] = zlup(endmjd(ista), az(ista), el(ista), ha(ista), dc, ista, station, handles.PARA.MIN_CUTEL);
                storage(idx==isrc,ista+1) = {lup};
                if lup
                    tot = tot +1;
                end
            end  
            storage(idx==isrc,ista+2) = {tot};
        end
        set(handles.uitable_whatsup,'Data',storage);
        handles.alreadyCalced_whatsup = 1;
    end
    
elseif strcmp(handles.popupmenu_table.String(handles.popupmenu_table.Value),sprintf('recommended'))
    handles.uipanel_staobs.Visible = 'off';
    handles.uipanel_srcobs.Visible = 'off';
    handles.uipanel_whatsup.Visible = 'off';
    handles.uipanel_recommended.Visible = 'on';
    
    if ~handles.alreadyCalced_recommended
        handles.uitable_recommended.Data = {};
               
        if strcmp(handles.popupmenu_recommended.String{handles.popupmenu_recommended.Value},'station based')
            handles.uitable_recommended.Data = {};
            handles.uitable_recommended.ColumnName = {'','score','start ','duration','end','obs','source 1','obs 1','source 2','obs 2'};
            subconall = bestSubcon( 'station',handles.stations,handles.twin,handles.source,...
            handles.obs.staobs,handles.obs.srcobs,handles.obsmode,handles.PARA, 0 );
        elseif strcmp(handles.popupmenu_recommended.String{handles.popupmenu_recommended.Value},'source based (1)')
            handles.uitable_recommended.Data = {};
            handles.uitable_recommended.ColumnName = {'','score','start ','duration','end','obs','source 1','obs 1'};
        elseif strcmp(handles.popupmenu_recommended.String{handles.popupmenu_recommended.Value},'source based (2)')
            handles.uitable_recommended.Data = {};
            handles.uitable_recommended.ColumnName = {'','score','start ','duration','end','obs','source 1','obs 1','source 2','obs 2'};
        elseif strcmp(handles.popupmenu_recommended.String{handles.popupmenu_recommended.Value},'source based (4)')
            handles.uitable_recommended.Data = {};
            handles.uitable_recommended.ColumnName = {'','score','start ','duration','end','obs','source 1','obs 1','source 2','obs 2','source 3','obs 3','source 4','obs 4'};
        end

        if exist('subconall','var')
            storage = cell(length(subconall),10);
            if subconall(1).nscan ~= 0 || length(subconall)>1
                [yy, xysort] = ssubsort(handles.source, handles.stations, handles.obs.staobs,...
                    handles.obs.srcobs, subconall, handles.PARA);
                score = xysort(yy);
                subconSort = subconall(yy);
                n = length(subconSort);

                for i = 0:n-1
                    storage{n-i,2} = sprintf('%.4f',score(i+1));
                    str = '';
                    for j = 1:length(subconSort(i+1).scan)
                        storage{n-i,5+2*j} = handles.source(subconSort(i+1).scan(j).srcid).name;
                        str = {handles.stations([subconSort(i+1).scan(j).sta.staid]).po};
                        str = strjoin(str);
                        storage{n-i,5+2*j+1} = str;
                    end
                    sta = [subconSort(i+1).scan.nsta];
                    storage{n-i,6} = sprintf('%d',sum((sta.*(sta-1))./2));
                    start = min([subconSort(i+1).scan.startmjd]);
                    [year, month, day, hour, minute, second] = tymdhms(start);
                    [datestr] = tdatestr(year, month, day, hour, minute, second);
                    storage{n-i,3} = datestr(12:end);

                    dur = [subconSort(i+1).scan];
                    dur = [dur.sta];
                    dur = max([dur.duration]);

                    storage{n-i,4} = sprintf('%d',dur);
                    [year, month, day, hour, minute, second] = tymdhms(start+dur/86400);
                    [datestr] = tdatestr(year, month, day, hour, minute, second);
                    storage{n-i,5} = datestr(12:end);
                end
                set(handles.uitable_recommended,'Data',storage);
            else
                msgbox('no possible subconfiguration found... Check your settings', 'Error','error');
            end
        end
        handles.alreadyCalced_recommended = 1;
    end
end
guidata(hObject, handles);

function storage = popupmenu_sourceSort_helpfunction_Callback(handles, storage)
if strcmp(handles.popupmenu_sourceSort.String(handles.popupmenu_sourceSort.Value),'name')
    d = storage(:,1);
    [~,idx] = sort(d);
    storage = storage(idx,:);
elseif strcmp(handles.popupmenu_sourceSort.String(handles.popupmenu_sourceSort.Value),'ra')
    d = storage(:,2);
    d = cell2mat(d);
    [~,idx] = sort(d);
    storage = storage(idx,:);
elseif strcmp(handles.popupmenu_sourceSort.String(handles.popupmenu_sourceSort.Value),'de')
    d = storage(:,3);
    d = cell2mat(d);
    [~,idx] = sort(d);
    storage = storage(idx,:);
elseif strcmp(handles.popupmenu_sourceSort.String(handles.popupmenu_sourceSort.Value),'def')
    d = storage(:,4);
    d = cell2mat(d);
    [~,idx] = sort(d,'descend');
    storage = storage(idx,:);
elseif strcmp(handles.popupmenu_sourceSort.String(handles.popupmenu_sourceSort.Value),'star')
    d = storage(:,5);
    d = cell2mat(d);
    [~,idx] = sort(d,'descend');
    storage = storage(idx,:);
elseif strcmp(handles.popupmenu_sourceSort.String(handles.popupmenu_sourceSort.Value),'scans')
    d = storage(:,6);
    d = cell2mat(d);
    [~,idx] = sort(d,'descend');
    storage = storage(idx,:);
elseif strcmp(handles.popupmenu_sourceSort.String(handles.popupmenu_sourceSort.Value),'time')
    d = storage(:,8);
    d = cell2mat(d);
    [~,idx] = sort(d,'descend');
    if ~isempty(idx)
        storage.Data = storage.Data(idx,:);
    end
end



% --- Executes on selection change in popupmenu_sourceSort.
function storage = popupmenu_sourceSort_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_sourceSort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_sourceSort contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_sourceSort
storage = handles.uitable_srcobs.Data;
storage = popupmenu_sourceSort_helpfunction_Callback(handles, storage);
handles.uitable_srcobs.Data = storage;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_sourceSort_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_sourceSort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in legend.
function legend_Callback(hObject, eventdata, handles)
% hObject    handle to legend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = msgbox({'red: source' 'blue: defining source' 'green: current position' '-----------' 'circle: source' 'star: STAR source'});


% --- Executes on selection change in listbox_sources.
function listbox_sources_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_sources (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idx = handles.listbox_sources.Value;
item = handles.listbox_sources.String(idx);
newString = [handles.listbox_scans.String;item];
handles.listbox_scans.String = newString;
if idx == length(handles.listbox_sources.String)
    handles.listbox_sources.Value = idx-1;
end
handles.listbox_sources.String(idx) = [];
updateSched(handles.uitable_sched, eventdata, handles)
handles.uitable_sched.Data(:,2)={'none'};
if ~isempty(handles.listbox_scans.String)
    handles.pushbutton_continue.Enable = 'on';
else
    handles.pushbutton_continue.Enable = 'off';
end
% Hints: contents = cellstr(get(hObject,'String')) returns listbox_sources contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_sources


% --- Executes during object creation, after setting all properties.
function listbox_sources_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_sources (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_scans.
function listbox_scans_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_scans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idx = handles.listbox_scans.Value;
if ~isempty(idx)
    item = handles.listbox_scans.String(idx);
    newString = [handles.listbox_sources.String;item];
    handles.listbox_sources.String = newString;
    if idx == length(handles.listbox_scans.String) && idx>1
        handles.listbox_scans.Value = idx-1;
    end
    handles.listbox_scans.String(idx) = [];
    updateSched(handles.uitable_sched, eventdata, handles)
    handles.uitable_sched.Data(:,2)={'none'};
else
    handles.listbox_scans.Value=1;
end
newList = handles.listbox_sources.String;
[~,idx] = sort(newList);
newList = newList(idx);
handles.listbox_sources.String = newList;
if ~isempty(handles.listbox_scans.String)
    handles.pushbutton_continue.Enable = 'on';
else
    handles.pushbutton_continue.Enable = 'off';
end

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_scans contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_scans


% --- Executes during object creation, after setting all properties.
function listbox_scans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_scans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function updateSched(hObject, eventdata, handles)
str = handles.listbox_scans.String';
str(2:length(str)+1)=str;
str(1)={'none'};
handles.uitable_sched.ColumnFormat{2}=str;

ncolumn = length(str);
cnames = handles.listbox_scans.String';
cnames(3:length(cnames)+2)=cnames;
cnames(1)={'station'};
cnames(2)={'scan'};
handles.uitable_sched.ColumnName = cnames;

handles.uitable_sched.ColumnFormat(3:2+ncolumn)={'logical'};

handles.uitable_sched.Data(:,3:end)=[];
sourceName = handles.listbox_scans.String;
station = handles.stations;
endmjd = [handles.obs.staobs.endmjd];

llh = [station.llh];
lon = llh(1:3:end);
lat = llh(2:3:end);

for isrc = 1:length(sourceName)
    thisSource = handles.source(strcmp(sourceName(isrc),{handles.source.name}));
    ra = thisSource.ra;
    de = thisSource.de;
    [az, el, ha, dc] = zazel_s(endmjd, lon, lat, ra, de);

    for ista = 1 : length(station)
        [lup] = zlup(endmjd(ista), az(ista), el(ista), ha(ista), dc, ista, station, handles.PARA.MIN_CUTEL);

        handles.uitable_sched.Data(ista,isrc+2) = {lup};

    end  
end


% --- Executes on selection change in popupmenu_recommended.
function popupmenu_recommended_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_recommended (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateTable(handles.popupmenu_table, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_recommended contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_recommended


% --- Executes during object creation, after setting all properties.
function popupmenu_recommended_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_recommended (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_recommended_starmode.
function checkbox_recommended_starmode_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_recommended_starmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateTable(handles.popupmenu_table, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of checkbox_recommended_starmode


% --- Executes on button press in pushbutton_submit.
function pushbutton_submit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_submit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
subcon = handles.subcon;
[flag, staid] = singleCheckSource(handles.source, handles.stations, subcon, handles.PARA, handles.jpl);
if flag ~= 1
    popupmenu_mode_Callback(handles.popupmenu_mode, eventdata, handles)
    switch flag
        case 2
            msgbox(sprintf('ERROR while checking with rigorouse model: source not visible at start of the observation (%s)',handles.stations(staid).name),'ERROR','error')
        case 3
            msgbox(sprintf('ERROR while checking with rigorouse model: source not visible at end of the observation (%s)',handles.stations(staid).name),'ERROR','error')
        case 4
            msgbox(sprintf('ERROR while checking with rigorouse model: source not visible during observation (%s)',handles.stations(staid).name),'ERROR','error')
        case 5
            msgbox(sprintf('ERROR while checking with rigorouse model: problems with cable wrap consistency (%s)',handles.stations(staid).name),'ERROR','error')
        otherwise
            msgbox('ERROR','ERROR','error')
    end
else
    if isfield(handles.sched,'scan')
        scanuminfo = length([handles.sched.scan]);
        nsched = length(handles.sched);
    else
        scanuminfo = 0;
        nsched = 0;
    end
    
    if handles.radiobutton_continueFillin1.Value
        % #### fill-in 1 scheduling mode ####   
        while true
            [subcon_f1] = sscanfi1(handles.stations, handles.twin,  handles.obs.staobs, handles.source, handles.obs.srcobs, handles.obsmode, subcon, handles.PARA, handles.jpl);
            if (subcon_f1.nscan > 0)
                fprintf('------------ fill in 1 scan -----------  (sched: %3d)\n',nsched+1);
                [handles.obs.staobs, handles.obs.srcobs, endmjdmax, scanuminfo] = supdate(handles.stations, handles.source, subcon_f1, handles.obs.staobs, handles.obs.srcobs, scanuminfo, handles.PARA);
                nsched = nsched + 1;
                if isempty(fieldnames(handles.sched))
                    handles.sched = subcon_f1;
                else
                    handles.sched(length(handles.sched)+1) = subcon_f1;
                end
            elseif (subcon_f1.nscan == 0)
                break;
            end
        end
    end

    % #### update structures ####
    fprintf('------------ %d-scan subnet ------------  (sched: %3d)\n',subcon.nscan,nsched+1);

    [handles.obs.staobs, handles.obs.srcobs, endmjdmax, scanuminfo] = supdate(handles.stations, handles.source, subcon, handles.obs.staobs, handles.obs.srcobs, scanuminfo, handles.PARA);
    if (subcon.nscan > 0)
        nsched = nsched + 1;
        if isempty(fieldnames(handles.sched))
            handles.sched = subcon;
        else
            handles.sched(length(handles.sched)+1) = subcon;
        end
    end

    if handles.radiobutton_continueFillin2.Value
        % #### fill-in 2 scheduling mode ####
        endmjdmaxfi2 = endmjdmax;
        while true
            [subcon_f2] = sscanfi2(handles.stations, handles.twin, handles.obs.staobs, handles.source, handles.obs.srcobs, handles.obsmode, endmjdmaxfi2, handles.PARA, handles.jpl);
            if (subcon_f2.nscan > 0)
                fprintf('------------ fill in 2 scan -----------  (sched: %3d)\n',nsched+1);
                [handles.obs.staobs, handles.obs.srcobs, endmjdmax, scanuminfo] = supdate(handles.stations, handles.source, subcon_f2, handles.obs.staobs, handles.obs.srcobs, scanuminfo, handles.PARA);
                nsched = nsched + 1;
                handles.sched(nsched) = subcon_f2;
            elseif (subcon_f2.nscan == 0)
                break;
            end
        end
    end



    fprintf('=======================================\n');
    handles.time = max([handles.obs.staobs.endmjd]);

    [year, month, day, hour, minute, second] = tymdhms(handles.time);
    [datestr] = tdatestr(year, month, day, hour, minute, second);
    handles.text_timeDate.String = datestr;
    handles.text_timeMJD.String = sprintf('%.3f',handles.time);
    handles.text_timeH.String = (handles.time-handles.startsched)*24;

    
    if handles.time > handles.endsched
        sched = handles.sched;
        save('..\DATA\LEVEL5\sched_manually.mat','sched')
        close(handles.figure1);
    end
    handles.alreadyCalced_recommended = 0;
    handles.alreadyCalced_source = 0;
    handles.alreadyCalced_station = 0;
    handles.alreadyCalced_whatsup = 0;
    guidata(hObject, handles);
    popupmenu_mode_Callback(handles.popupmenu_mode, eventdata, handles)
end







% --- Executes on button press in pushbutton_continue.
function pushbutton_continue_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_continue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.uipanel_manualSched.Visible = 'off';
handles.uipanel_continue.Visible = 'on';
handles.popupmenu_mode.Enable = 'off';
handles.uitable_perStation.Data = {};
handles.uitable_perBaseline.Data = {};

subcon = struct();
subcon.nscan = length(handles.listbox_scans.String);
obsSources = handles.listbox_scans.String;
counter = 1;
for i = 1:subcon.nscan
    subcon.scan(counter).srcid = find(strcmp({handles.source.name},obsSources(i)));
    subcon.scan(counter).startmjd = 0;
    subcon.scan(counter).nsta = 0;
    for j = 1:length(handles.stations)
        if strcmp(handles.uitable_sched.Data(j,2),obsSources(i))
            subcon.scan(counter).nsta = subcon.scan(counter).nsta+1;
            subcon.scan(counter).sta(subcon.scan(counter).nsta).staid = j;
            subcon.scan(counter).sta(subcon.scan(counter).nsta).slewtime = 0;
            subcon.scan(counter).sta(subcon.scan(counter).nsta).startmjd = 0;
            subcon.scan(counter).sta(subcon.scan(counter).nsta).az = 0;
            subcon.scan(counter).sta(subcon.scan(counter).nsta).el = 0;
            subcon.scan(counter).sta(subcon.scan(counter).nsta).ha = 0;
            subcon.scan(counter).sta(subcon.scan(counter).nsta).dc = 0;
            subcon.scan(counter).sta(subcon.scan(counter).nsta).duration = 0;
            subcon.scan(counter).sta(subcon.scan(counter).nsta).endmjd = 0;
        end
    end
    if subcon.scan(counter).nsta == 0
        subcon.scan(counter)=[];
        subcon.nscan = subcon.nscan -1;
    else
        counter = counter +1;
    end
end

subcon1 = sstartmjd(handles.source, handles.stations, handles.twin, handles.obs.staobs, subcon, handles.PARA);
[subcon2,actbl] = sduration(handles.source, handles.stations, handles.twin, handles.obsmode, subcon1, handles.PARA);

bool = 1;
if subcon2.nscan == subcon.nscan
    for i = 1:subcon2.nscan
        if ~isequal([subcon2.scan(i).sta.staid],[subcon.scan(i).sta.staid])
            bool = 0;
            break
        end
    end
else
    bool = 0;
end
   
if bool == 0
    waitfor(msgbox('some changes were made due to internal checks. To avoid this click "back" and try to change your settings!','WARNING','warn'))
end

subcon = subcon2;

if subcon.nscan == 0
    pushbutton_back_Callback(hObject, eventdata, handles);
    waitfor(msgbox('There was no possible observation!','NO SCAN POSSIBLE','error'))
end

storage = perStationUpdate(subcon,handles);
handles.uitable_perStation.Data = storage;

 
counter = 1;
for i = 1:length(actbl)
    bls = actbl(i).actbl;
    for j = 1:length(bls)
        storageBL(counter,1) = {handles.stations(bls(j).staid1).name};
        storageBL(counter,2) = {handles.stations(bls(j).staid2).name};
        storageBL(counter,3) = {true};
        storageBL(counter,4) = {bls(j).duration(1)};
        storageBL(counter,5) = {bls(j).duration(2)};
        counter = counter+1;
    end
    if i ~= length(actbl)
        storageBL(counter,3) = {false};
        counter = counter+1;
    end
    
end
handles.subcon = subcon;
handles.uitable_perBaseline.Data = storageBL;
guidata(hObject, handles);

function storage = perStationUpdate(subcon,handles)
counter = 1;
for i = 1:subcon.nscan
    thisScan = subcon.scan(i);
    for j = 1:thisScan.nsta
        storage(counter,1) = {handles.source(thisScan.srcid).name};
        storage(counter,2) = {handles.stations(thisScan.sta(j).staid).name};

        [year, month, day, hour, minute, second] = tymdhms(thisScan.startmjd);
        [datestr] = tdatestr(year, month, day, hour, minute, second);
        storage(counter,3) = {datestr(12:end)};

        storage(counter,4) = {thisScan.sta(j).duration};
        
        [year, month, day, hour, minute, second] = tymdhms(thisScan.startmjd+thisScan.sta(j).duration/86400);
        [datestr] = tdatestr(year, month, day, hour, minute, second);
        storage(counter,5) = {datestr(12:end)};
        counter = counter +1;
    end
    counter = counter +1;
end


% --- Executes on button press in radiobutton_continueFillin2.
function radiobutton_continueFillin2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_continueFillin2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_continueFillin2


% --- Executes on button press in radiobutton_continueFillin1.
function radiobutton_continueFillin1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_continueFillin1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_continueFillin1


% --- Executes when entered data in editable cell(s) in uitable_recommended.
function uitable_recommended_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_recommended (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
idx = eventdata.Indices(1);
handles.uitable_recommended.Data(:,1) = {logical(0)};
handles.uitable_recommended.Data(idx,1) = {logical(1)};

for i = 1:length(handles.listbox_scans.String)
    listbox_scans_Callback(hObject, eventdata, handles);
end

nsources = (length(handles.uitable_recommended.Data(1,:))-6)/2;
for i = 1:nsources
    source = handles.uitable_recommended.Data(idx,5+2*i);
    if ~isequal(source,{[]})
        srcid = find(strcmp(handles.listbox_sources.String,source));
        handles.listbox_sources.Value = srcid;
        listbox_sources_Callback(hObject, eventdata, handles);
    end
end

for i = 1:nsources
    source = handles.uitable_recommended.Data(idx,5+2*i);
    if ~isequal(source,{[]})
        sta = handles.uitable_recommended.Data(idx,5+2*i+1);
        sta = strsplit(sta{1});
        for j = 1:length(sta)
            thisSta = sta{j};
            handles.uitable_sched.Data(find(strcmp({handles.stations.po},thisSta)),2)=source;
        end
    end
end

% handles    structure with handles and user data (see GUIDATA)


function edit_PARA_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MIN_SRCRP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Minimum SNR
handles.PARA.MIN_SNR(1) = str2double(handles.edit_SNR_X_Band.String);
handles.PARA.MIN_SNR(2) = str2double(handles.edit_SNR_S_Band.String);

% Parameters
handles.PARA.MIN_SRCRP = str2double(handles.edit_MIN_SRCRP.String);
handles.PARA.MAXSLEWTIME = str2double(handles.edit_MAXSLEWTIME.String);
handles.PARA.MAX_WAIT = str2double(handles.edit_MAX_WAIT.String);
handles.PARA.MAX_SCAN = str2double(handles.edit_MAX_SCAN.String);
handles.PARA.MIN_SCAN = str2double(handles.edit_MIN_SCAN.String);
handles.PARA.MIN_SUNDIST = str2double(handles.edit_SUNDIST.String) /180*pi;
handles.PARA.MIN_CUTEL = str2double(handles.edit_CUTOFEL.String) /180*pi;

% individual weight factors
handles.PARA.WEIGHT_NUMBER_OF_OBS = str2double(handles.edit_weight_scanendtime.String);
handles.PARA.WEIGHT_SKY_COVERAGE = str2double(handles.edit_weight_skyCoverage.String);
handles.PARA.WEIGHT_SCAN_END_TIME = str2double(handles.edit_weight_numberOfObservations.String);

% minimum number of stations
handles.PARA.MIN_STANUM = str2double(handles.edit_min_stanum.String);
handles.PARA.MIN_STASCAN = str2double(handles.edit_min_stascan.String);
handles.PARA.MIN_STANUM_FI = str2double(handles.edit_min_stanum_fi.String);

%subnetting 
if handles.radiobutton_subnetting_yes.Value
    handles.PARA.SUBNETTING = 2;
else
    handles.PARA.SUBNETTING = 1;
end

guidata(hObject, handles);


% --- Executes on button press in pushbutton_changeSettings.
function pushbutton_changeSettings_Callback(hObject, eventdata, handles)
prompt = {'min source repeat','max slewtime','max wait','max scan','min scan','min sundist','cut-off elevation','min SNR X-Band','min SNR S-Band','min stations per subconfiguration','min stations per scan','min stations per fillin mode','subnetting (1 = false, 2 = true)', 'weight number of obs','weight sky coverage', 'weight scan end time'};
dlg_title = 'Settings';
num_lines = 1;
defaultans = {num2str(handles.PARA.MIN_SRCRP),num2str(handles.PARA.MAXSLEWTIME),num2str(handles.PARA.MAX_WAIT),num2str(handles.PARA.MAX_SCAN),...
    num2str(handles.PARA.MIN_SCAN),num2str(handles.PARA.MIN_SUNDIST/pi*180),num2str(handles.PARA.MIN_CUTEL/pi*180),num2str(handles.PARA.MIN_SNR(1)),num2str(handles.PARA.MIN_SNR(2)),...
    num2str(handles.PARA.MIN_STANUM),num2str(handles.PARA.MIN_STASCAN),num2str(handles.PARA.MIN_STANUM_FI),num2str(handles.PARA.SUBNETTING),...
    num2str(handles.PARA.WEIGHT_NUMBER_OF_OBS),num2str(handles.PARA.WEIGHT_SKY_COVERAGE),num2str(handles.PARA.WEIGHT_SCAN_END_TIME)};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

if ~isempty(answer)
    handles.PARA.MIN_SRCRP = str2double(answer{1});
    handles.PARA.MAXSLEWTIME = str2double(answer{2});
    handles.PARA.MAX_WAIT = str2double(answer{3});
    handles.PARA.MAX_SCAN = str2double(answer{4});
    handles.PARA.MIN_SCAN = str2double(answer{5});
    handles.PARA.MIN_SUNDIST = str2double(answer{6})/180*pi;
    handles.PARA.MIN_CUTEL = str2double(answer{7})/180*pi;
    handles.PARA.MIN_SNR(1) = str2double(answer{8});
    handles.PARA.MIN_SNR(2) = str2double(answer{9});
    handles.PARA.MIN_STANUM = str2double(answer{10});
    handles.PARA.MIN_STASCAN = str2double(answer{11});
    handles.PARA.MIN_STANUM_FI = str2double(answer{12});
    handles.PARA.SUBNETTING = str2double(answer{13});
    handles.PARA.WEIGHT_NUMBER_OF_OBS = str2double(answer{14});
    handles.PARA.WEIGHT_SKY_COVERAGE = str2double(answer{15});
    handles.PARA.WEIGHT_SCAN_END_TIME = str2double(answer{16});

    % add changes on auto sched window
    handles.edit_MIN_SRCRP.String = handles.PARA.MIN_SRCRP;
    handles.edit_MAXSLEWTIME.String = handles.PARA.MAXSLEWTIME;
    handles.edit_MAX_WAIT.String = handles.PARA.MAX_WAIT;
    handles.edit_MAX_SCAN.String = handles.PARA.MAX_SCAN;
    handles.edit_MIN_SCAN.String = handles.PARA.MIN_SCAN;
    handles.edit_SUNDIST.String = handles.PARA.MIN_SUNDIST*180/pi;
    handles.edit_CUTOFEL.String = handles.PARA.MIN_CUTEL*180/pi;
    handles.edit_SNR_X_Band.String = handles.PARA.MIN_SNR(1);
    handles.edit_SNR_S_Band.String = handles.PARA.MIN_SNR(2);
    handles.edit_min_stanum.String = handles.PARA.MIN_STANUM;
    handles.edit_min_stascan.String = handles.PARA.MIN_STASCAN;
    handles.edit_min_stanum_fi.String = handles.PARA.MIN_STANUM_FI;
    handles.edit_weight_scanendtime.String = handles.PARA.WEIGHT_SCAN_END_TIME;
    handles.edit_weight_skyCoverage.String = handles.PARA.WEIGHT_SKY_COVERAGE;
    handles.edit_weight_numberOfObservations.String = handles.PARA.WEIGHT_NUMBER_OF_OBS;
    if handles.PARA.SUBNETTING == 1
        handles.radiobutton_subnetting_no.Value = 1;
        handles.radiobutton_subnetting_yes.Value = 0;
    else
        handles.radiobutton_subnetting_yes.Value = 1;
        handles.radiobutton_subnetting_n0.Value = 0;
    end
    
    if strcmp(hObject.String,'settings')
        handles.alreadyCalced_recommended = 0;
        updateTable(hObject, eventdata, handles);
    end
    
end
guidata(hObject, handles);


% --- Executes on button press in pushbutton_back.
function pushbutton_back_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.uipanel_manualSched.Visible = 'on';
handles.uipanel_continue.Visible = 'off';
handles.popupmenu_mode.Enable = 'on';



% --- Executes on selection change in popupmenu_sortWhatsUp.
function popupmenu_sortWhatsUp_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_sortWhatsUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.popupmenu_sortWhatsUp.String(handles.popupmenu_sortWhatsUp.Value),'name')
    d = handles.uitable_whatsup.Data(:,1);
    [~,idx] = sort(d);
    handles.uitable_whatsup.Data = handles.uitable_whatsup.Data(idx,:);
elseif strcmp(handles.popupmenu_sortWhatsUp.String(handles.popupmenu_sortWhatsUp.Value),'total')
    d = handles.uitable_whatsup.Data(:,end);
    d = cell2mat(d);
    [~,idx] = sort(d,'descend');
    handles.uitable_whatsup.Data = handles.uitable_whatsup.Data(idx,:);
end

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_sortWhatsUp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_sortWhatsUp


% --- Executes during object creation, after setting all properties.
function popupmenu_sortWhatsUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_sortWhatsUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_lastMinutes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lastMinutes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function print_skyplot(hObject, eventdata, handles)

if isfield(handles.allstaobs,'staid')
    lastMinutes = str2double(handles.edit_lastMinutes.String);
    ax = handles.axes_skyCoverage;
    for i = 5:length(ax.Children)
        delete(ax.Children(1));
    end
    staid = handles.popupmenu_station.String(handles.popupmenu_station.Value);
    staid = find(strcmp(staid,{handles.stations.name}));
    obs = handles.allstaobs([handles.allstaobs.staid] == staid);
    obs = obs([obs.startmjd] > handles.time-lastMinutes/1440);
    boolStar = [obs.star];
    boolDef = [obs.defining];
    boolSource = ~(boolStar | boolDef);

    obsStar = obs(logical(boolStar));
    obsDef = obs(logical(boolDef));
    obsSource = obs(logical(boolSource));
    hold on
    pdef = polar([obsDef.az],[90-[obsDef.el]*180/pi]./90);
    set(pdef,'marker','o','MarkerEdgeColor','k','MarkerFaceColor','b','Linestyle','none');

    p1 = polar([obsSource.az],[90-[obsSource.el]*180/pi]./90);
    set(p1,'marker','o','MarkerEdgeColor','k','MarkerFaceColor','r','Linestyle','none');

    p1Star = polar([obsStar.az],[90-[obsStar.el]*180/pi]./90);
    set(p1Star,'marker','p','MarkerEdgeColor','k','markerSize',10,'MarkerFaceColor','r','Linestyle','none')
    
    pnow = polar([handles.obs.staobs(staid).az],[90-[handles.obs.staobs(staid).el]*180/pi]./90);
    set(pnow,'marker','d','MarkerEdgeColor','k','MarkerFaceColor','g','Linestyle','none');
    hold off
    
    
end

function edit_lastMinutes_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lastMinutes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
print_skyplot(hObject, eventdata, handles)


% --- Executes when entered data in editable cell(s) in uitable_sched.
function uitable_sched_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_sched (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
row = eventdata.Indices(1);
col = eventdata.Indices(2);
data = handles.uitable_sched.Data(row,col);
colSource = find(strcmp(data,handles.uitable_sched.ColumnName));
if ~strcmp(data{1},'none')
    if ~handles.uitable_sched.Data{row,colSource}
        handles.uitable_sched.Data(row,col) = {'none'};
        msgbox('source is not visible','ERROR','error');
    end
end

% --- Executes when entered data in editable cell(s) in uitable_perBaseline.
function uitable_perBaseline_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_perBaseline (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
idx = eventdata.Indices(1);
idxS = [handles.subcon.scan.nsta];
idxS = (idxS.*(idxS-1))/2;

idxBLs = [];
counter = 1;
for i = 1:length(idxS)
    idxBLs = [idxBLs counter];
    counter =  counter+idxS(i)-1;
    idxBLs = [idxBLs counter];
    counter = counter +2;
end

sta1 = handles.uitable_perBaseline.Data(idx,1);
bool = cell2mat([handles.uitable_perBaseline.Data(:,3)]);
dur1 = strcmp(handles.uitable_perBaseline.Data(:,1),sta1);
dur2 = strcmp(handles.uitable_perBaseline.Data(:,2),sta1);
dur = dur1 | dur2;
dur = dur & bool;
d1 = max([handles.uitable_perBaseline.Data{dur,4}]);
d2 = max([handles.uitable_perBaseline.Data{dur,5}]);
dTot1 = max([d1,d2]);
if isempty(dTot1) || dTot1 < handles.PARA.MIN_SCAN
    dTot1 = handles.PARA.MIN_SCAN;
end

sta2 = handles.uitable_perBaseline.Data(idx,2);
bool = cell2mat([handles.uitable_perBaseline.Data(:,3)]);
dur1 = strcmp(handles.uitable_perBaseline.Data(:,1),sta2);
dur2 = strcmp(handles.uitable_perBaseline.Data(:,2),sta2);
dur = dur1 | dur2;
dur = dur & bool;
d1 = max([handles.uitable_perBaseline.Data{dur,4}]);
d2 = max([handles.uitable_perBaseline.Data{dur,5}]);
dTot2 = max([d1,d2]);
if isempty(dTot2) || dTot2 < handles.PARA.MIN_SCAN
    dTot2 = handles.PARA.MIN_SCAN;
end

sta1idx = find(strcmp({handles.stations.name},sta1));
sta2idx = find(strcmp({handles.stations.name},sta2));

for i = 1:handles.subcon.nscan
    for j = 1:handles.subcon.scan(i).nsta
        if handles.subcon.scan(i).sta(j).staid == sta1idx
            handles.subcon.scan(i).sta(j).duration = dTot1;
            handles.subcon.scan(i).sta(j).endmjd = handles.subcon.scan(i).startmjd+dTot1/86400;
            end1 = handles.subcon.scan(i).sta(j).endmjd;
        end
        if handles.subcon.scan(i).sta(j).staid == sta2idx
            handles.subcon.scan(i).sta(j).duration = dTot2;
            handles.subcon.scan(i).sta(j).endmjd = handles.subcon.scan(i).startmjd+dTot2/86400;
            end2 = handles.subcon.scan(i).sta(j).endmjd;
        end
    end
end

for i = 1:size(handles.uitable_perStation.Data,1)
    if strcmp(sta1,handles.uitable_perStation.Data(i,2)) 
        handles.uitable_perStation.Data(i,4) = {dTot1};
        
        [year, month, day, hour, minute, second] = tymdhms(end1);
        [datestr] = tdatestr(year, month, day, hour, minute, second);
        handles.uitable_perStation.Data(i,5) = {datestr(12:end)};
    end
    if strcmp(sta2,handles.uitable_perStation.Data(i,2))
        handles.uitable_perStation.Data(i,4) = {dTot2};
        
        [year, month, day, hour, minute, second] = tymdhms(end2);
        [datestr] = tdatestr(year, month, day, hour, minute, second);
        handles.uitable_perStation.Data(i,5) = {datestr(12:end)};
    end
end
guidata(hObject, handles);

% --- Executes when entered data in editable cell(s) in uitable_perStation.
function uitable_perStation_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_perStation (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
row = eventdata.Indices(1);
col = eventdata.Indices(2);
data = handles.uitable_perStation.Data{row,col};

if col == 4
    if data < handles.PARA.MIN_SCAN
        data = handles.PARA.MIN_SCAN;
    end


    staid = find(strcmp({handles.stations.name},handles.uitable_perStation.Data{row,2}));
    subcon = handles.subcon;
    for i = 1:handles.subcon.nscan
        bool = [subcon.scan(i).sta.staid] == staid;
        if any(bool)
            subcon.scan(i).sta(bool).duration = data;
            endtime = subcon.scan(i).startmjd + data/86400;
            subcon.scan(i).sta(bool).endmjd = endtime;
        end
    end

    [year, month, day, hour, minute, second] = tymdhms(endtime);
    [datestr] = tdatestr(year, month, day, hour, minute, second);
    handles.uitable_perStation.Data(row,5) = {datestr(12:end)};

    handles.subcon = subcon;
    guidata(hObject, handles);
else
    subcon = handles.subcon;
    staid = find(strcmp({handles.stations.name},handles.uitable_perStation.Data{row,2}));
    for i = 1:subcon.nscan
        if any([subcon.scan(i).sta.staid]==staid)
            scanNr = i;
            minStartDate = handles.subcon.scan(i).startmjd;
            break;
        end
    end
    [year, month, day, hour, minute, second] = tymdhms(minStartDate);
    [datestr] = tdatestr(year, month, day, hour, minute, second);
    newTime = duration(str2double(data(1:2)),str2double(data(4:5)),str2double(data(7:8)));
    oldTime = duration(hour,minute,second);
    sec = seconds(newTime-oldTime);
    
    if sec<0 && mod(minStartDate,1)>.9
        sec = seconds(newTime);
        handles.subcon.scan(i).startmjd = floor(minStartDate)+1+sec/86400;
        for i = 1:handles.subcon.scan(scanNr).nsta
            handles.subcon.scan(scanNr).sta(i).startmjd = handles.subcon.scan(scanNr).startmjd;
            handles.subcon.scan(scanNr).sta(i).endmjd = handles.subcon.scan(scanNr).sta(i).startmjd + handles.subcon.scan(scanNr).sta(i).duration/86400;
        end
        storage = perStationUpdate(handles.subcon,handles);
        handles.uitable_perStation.Data = storage;

    elseif sec>0
        handles.subcon.scan(i).startmjd = minStartDate+sec/86400;
        for i = 1:handles.subcon.scan(scanNr).nsta
            handles.subcon.scan(scanNr).sta(i).startmjd = handles.subcon.scan(scanNr).startmjd;
            handles.subcon.scan(scanNr).sta(i).endmjd = handles.subcon.scan(scanNr).sta(i).startmjd + handles.subcon.scan(scanNr).sta(i).duration/86400;
        end
        storage = perStationUpdate(handles.subcon,handles);
        handles.uitable_perStation.Data = storage;
        
    else
        handles.uitable_perStation.Data(row,col) = {eventdata.PreviousData};
        msgbox('new start time is too far from minimal start time.','ERROR','error')
    end
end
guidata(hObject, handles);


% --- Executes on selection change in listbox_useStation.
function listbox_useStation_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_useStation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if length(handles.listbox_useStation.String)>2
    idx = handles.listbox_useStation.Value;
    item = handles.listbox_useStation.String(idx);

    newString = [handles.listbox_dontUseStation.String;item];
    handles.listbox_dontUseStation.String = newString;

    if idx == length(handles.listbox_useStation.String)
        handles.listbox_useStation.Value = idx-1;
    end
    handles.listbox_useStation.String(idx) = [];
end

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_useStation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_useStation


% --- Executes during object creation, after setting all properties.
function listbox_useStation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_useStation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_dontUseStation.
function listbox_dontUseStation_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_dontUseStation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idx = handles.listbox_dontUseStation.Value;
if ~isempty(idx)
    item = handles.listbox_dontUseStation.String(idx);
    
    newString = [handles.listbox_useStation.String;item];
    handles.listbox_useStation.String = newString;
    
    if idx == length(handles.listbox_dontUseStation.String) && idx>1
        handles.listbox_dontUseStation.Value = idx-1;
    end
    handles.listbox_dontUseStation.String(idx) = [];
else
    handles.listbox_dontUseStation.Value=1;
end
newList = handles.listbox_useStation.String;
[~,idx] = sort(newList);
newList = newList(idx);
handles.listbox_useStation.String = newList;

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_dontUseStation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_dontUseStation


% --- Executes during object creation, after setting all properties.
function listbox_dontUseStation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_dontUseStation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_weight_scanendtime_Callback(hObject, eventdata, handles)
% hObject    handle to edit_weight_scanendtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_weight_scanendtime as text
%        str2double(get(hObject,'String')) returns contents of edit_weight_scanendtime as a double


% --- Executes during object creation, after setting all properties.
function edit_weight_scanendtime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_weight_scanendtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_weight_skyCoverage_Callback(hObject, eventdata, handles)
% hObject    handle to edit_weight_skyCoverage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_weight_skyCoverage as text
%        str2double(get(hObject,'String')) returns contents of edit_weight_skyCoverage as a double


% --- Executes during object creation, after setting all properties.
function edit_weight_skyCoverage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_weight_skyCoverage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_weight_numberOfObservations_Callback(hObject, eventdata, handles)
% hObject    handle to edit_weight_numberOfObservations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_weight_numberOfObservations as text
%        str2double(get(hObject,'String')) returns contents of edit_weight_numberOfObservations as a double


% --- Executes during object creation, after setting all properties.
function edit_weight_numberOfObservations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_weight_numberOfObservations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_min_stanum_Callback(hObject, eventdata, handles)
% hObject    handle to edit_min_stanum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_min_stanum as text
%        str2double(get(hObject,'String')) returns contents of edit_min_stanum as a double


% --- Executes during object creation, after setting all properties.
function edit_min_stanum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_min_stanum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_min_stascan_Callback(hObject, eventdata, handles)
% hObject    handle to edit_min_stascan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_min_stascan as text
%        str2double(get(hObject,'String')) returns contents of edit_min_stascan as a double


% --- Executes during object creation, after setting all properties.
function edit_min_stascan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_min_stascan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_min_stanum_fi_Callback(hObject, eventdata, handles)
% hObject    handle to edit_min_stanum_fi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_min_stanum_fi as text
%        str2double(get(hObject,'String')) returns contents of edit_min_stanum_fi as a double


% --- Executes during object creation, after setting all properties.
function edit_min_stanum_fi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_min_stanum_fi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
