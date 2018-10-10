% ************************************************************************
%   Description:
%   GUIGLOB_RF M-file for guiglob_rf.fig
%   2nd Graphical User Interface for vie_glob
%   Gain informations for estimation of the TRF and CRF
%
%   Input:										
%      parGS               name of the parameter and if it will be
%                          estimated (1), fixed (0) or reduced (2)
%      paths               paths to the data
%
%   Output:                stored in VieVS/OUT/GLOB/
%      paths               paths to the data
%
%   External calls: 	
%      globind             					    											
%       
%   Coded for VieVS: 
%   01 Aug 2010 by Hana Spicakova
%
%   Revision: 
%   05 Oct 2010 by Hana Spicakova: choice between NNT and NNT+dz for source
%      coordinates
%   12 Oct 2010 by Hana Spicakova: choice between NNT/NNR and NNT/NNR/NNS
%      for station coordinates
%   20 Jan 2011 by Hana Spicakova: some sources can be fixed to apriori
%      coordinates
%   14 Feb 2011 by Hana Spicakova: between different station velocity ties
%      can be used
%   20 Feb 2011 by Hana Spicakova: some stations (with only few
%      observations) can be session-wise reduced
%   24 Jul 2011 by Hana Spicakova: some sources can be session-wise reduced
%
%%


function varargout = guiglob_rf(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guiglob_rf_OpeningFcn, ...
                   'gui_OutputFcn',  @guiglob_rf_OutputFcn, ...
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


% --- Executes just before guiglob_rf is made visible.
function guiglob_rf_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for guiglob_rf
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

load(['../../OUT/GLOB/parGS'],'parGS');
[g] = globind(parGS);

if (parGS(g.g_coord(1)).id==0 && parGS(g.g_vel(1)).id==0)
    set(handles.cb_trf_6p,'Visible','off')
    set(handles.cb_trf_7p,'Visible','off')
    set(handles.text1,'Visible','off')
    set(handles.text2,'Visible','on')
    set(handles.file_datumsta,'Visible','off')
    set(handles.file_discont,'Visible','on')
    set(handles.cb_antred,'Visible','off')
    set(handles.text14,'Visible','off')
end
set(handles.cb_trf_6p,'Value',1)
set(handles.cb_trf_7p,'Value',0)

set(handles.text5,'Visible','off')
set(handles.file_velconst,'Visible','off')
set(handles.text13,'Visible','off')
set(handles.file_velties,'Visible','off')
set(handles.file_antred,'Visible','off')
set(handles.text15,'Visible','off')


if parGS(g.g_vel(1)).id==0
    set(handles.cb_velconst,'Visible','off')
    set(handles.cb_velties,'Visible','off')
end

if parGS(g.g_srade(1)).id==0
    set(handles.text6,'Visible','off')
    set(handles.text10,'Visible','off')
    set(handles.text12,'Visible','off')
    set(handles.file_datumsou,'Visible','off')
    set(handles.file_fixedsou,'Visible','off')
    set(handles.cb_datumsou,'Visible','off','Value',0)
    set(handles.cb_nnrsou,'Visible','off','Value',0)
    set(handles.cb_fixedsou,'Visible','off')
    set(handles.file_soured,'Visible','off')
    set(handles.cb_soured,'Visible','off')
    set(handles.text16,'Visible','off')
end



% TRF - datum files
fils=dir(['../../DATA/GLOB/TRF/DATUM/' '*.txt']);
optdir='';
for a=1:length(fils)
    optdir=[optdir '|' fils(a).name];
end
set(handles.file_datumsta,'String',optdir(2:length(optdir)));
set(handles.file_datumsta,'Value',1);
clear fils optdir

% TRF - session-wise estimation of chosen stations
fils=dir(['../../DATA/GLOB/TRF/REDUCE/' '*.txt']);
optdir='';
for a=1:length(fils)
    optdir=[optdir '|' fils(a).name];
end
set(handles.file_antred,'String',optdir(2:length(optdir)));
set(handles.file_antred,'Value',1);
clear fils optdir

% TRF - station discontinuities
fils=dir(['../../DATA/GLOB/TRF/DISCONT/' '*.mat']);
optdir='';
for a=1:length(fils)
    optdir=[optdir '|' fils(a).name];
end
set(handles.file_discont,'String',optdir(2:length(optdir)));
set(handles.file_discont,'Value',1);
clear fils optdir

% TRF - constant velocity for station with discontinuities
fils=dir(['../../DATA/GLOB/TRF/VELOC/' '*.txt']);
optdir='';
for a=1:length(fils)
    optdir=[optdir '|' fils(a).name];
end
set(handles.file_velconst,'String',optdir(2:length(optdir)));
set(handles.file_velconst,'Value',1);
clear fils optdir

% TRF - velocity ties
fils=dir(['../../DATA/GLOB/TRF/VELOC/TIES/' '*.txt']);
optdir='';
for a=1:length(fils)
    optdir=[optdir '|' fils(a).name];
end
set(handles.file_velties,'String',optdir(2:length(optdir)));
set(handles.file_velties,'Value',1);
clear fils optdir


% CRF - datum files
fils=dir(['../../DATA/GLOB/CRF/DATUM/' '*.txt']);
optdir='';
for a=1:length(fils)
    optdir=[optdir '|' fils(a).name];
end
set(handles.file_datumsou,'String',optdir(2:length(optdir)));
set(handles.file_datumsou,'Value',1);
clear fils optdir

% CRF - sources, which will be fixed
fils=dir(['../../DATA/GLOB/CRF/FIXED_SOURCES/' '*.txt']);
optdir='';
for a=1:length(fils)
    optdir=[optdir '|' fils(a).name];
end
set(handles.file_fixedsou,'String',optdir(2:length(optdir)));
set(handles.file_fixedsou,'Value',1);
clear fils optdir

% CRF - sources, which will be session-wise reduced
fils=dir(['../../DATA/GLOB/CRF/REDUCE/' '*.txt']);
optdir='';
for a=1:length(fils)
    optdir=[optdir '|' fils(a).name];
end
set(handles.file_soured,'String',optdir(2:length(optdir)));
set(handles.file_soured,'Value',1);
clear fils optdir



% UIWAIT makes guiglob_rf wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = guiglob_rf_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% --- Executes on selection change in file_datumsta.
function file_datumsta_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function file_datumsta_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in file_discont.
function file_discont_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function file_discont_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in file_velconst.
function file_velconst_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function file_velconst_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in file_velties.
function file_velties_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function file_velties_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in file_datumsou.
function file_datumsou_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function file_datumsou_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in file_fixedsou.
function file_fixedsou_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function file_fixedsou_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in file_antred.
function file_antred_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function file_antred_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
 

% --- Executes on button press in cb_trf_6p.
function cb_trf_6p_Callback(hObject, eventdata, handles)
if get(handles.cb_trf_6p,'Value') ==1
    set(handles.cb_trf_7p,'Value',0)
elseif get(handles.cb_trf_6p,'Value') ==0
    set(handles.cb_trf_7p,'Value',1)
end

% --- Executes on button press in cb_trf_7p.
function cb_trf_7p_Callback(hObject, eventdata, handles)
if get(handles.cb_trf_7p,'Value') ==1
    set(handles.cb_trf_6p,'Value',0)
elseif get(handles.cb_trf_7p,'Value') ==0
    set(handles.cb_trf_6p,'Value',1)
end



% --- Executes on button press in cb_datumsou.
function cb_datumsou_Callback(hObject, eventdata, handles)
if get(handles.cb_datumsou,'Value') ==1 
    set(handles.text6,'Visible','on')
    set(handles.file_datumsou,'Visible','on')
    set(handles.cb_nnrsou,'Value',0)
else
    set(handles.text6,'Visible','off')
    set(handles.file_datumsou,'Visible','off')
end

% --- Executes on button press in cb_nnrsou.
function cb_nnrsou_Callback(hObject, eventdata, handles)
if get(handles.cb_nnrsou,'Value') ==1 
    set(handles.text6,'Visible','on')
    set(handles.file_datumsou,'Visible','on')
    set(handles.cb_datumsou,'Value',0)
else
    set(handles.text6,'Visible','off')
    set(handles.file_datumsou,'Visible','off')
end


% --- Executes on button press in cb_fixedsou.
function cb_fixedsou_Callback(hObject, eventdata, handles)
if get(handles.cb_fixedsou,'Value') ==1 
    set(handles.text10,'Visible','on')
    set(handles.file_fixedsou,'Visible','on')
else
    set(handles.text10,'Visible','off')
    set(handles.file_fixedsou,'Visible','off')
end


% --- Executes on button press in cb_soured.
function cb_soured_Callback(hObject, eventdata, handles)
if get(handles.cb_soured,'Value') ==1 
    set(handles.text16,'Visible','on')
    set(handles.file_soured,'Visible','on')
else
    set(handles.text16,'Visible','off')
    set(handles.file_soured,'Visible','off')
end



% --- Executes on button press in cb_velconst.
function cb_velconst_Callback(hObject, eventdata, handles)
if get(handles.cb_velconst,'Value') ==1 
    set(handles.text5,'Visible','on')
    set(handles.file_velconst,'Visible','on')
else
    set(handles.text5,'Visible','off')
    set(handles.file_velconst,'Visible','off')
end



% --- Executes on button press in cb_velties.
function cb_velties_Callback(hObject, eventdata, handles)
if get(handles.cb_velties,'Value') ==1 
    set(handles.text13,'Visible','on')
    set(handles.file_velties,'Visible','on')
else
    set(handles.text13,'Visible','off')
    set(handles.file_velties,'Visible','off')
end


% --- Executes on button press in cb_antred.
function cb_antred_Callback(hObject, eventdata, handles)
if get(handles.cb_antred,'Value') ==1 
    set(handles.text15,'Visible','on')
    set(handles.file_antred,'Visible','on')
else
    set(handles.text15,'Visible','off')
    set(handles.file_antred,'Visible','off')
end




% --- Executes on button press in pb_OK.
function pb_OK_Callback(hObject, eventdata, handles)

load('../../OUT/GLOB/paths','paths')
load(['../../OUT/GLOB/parGS'],'parGS');
[g] = globind(parGS);


str=get(handles.file_datumsta,'String');
val=get(handles.file_datumsta,'Value');
paths.datumsta=strcat(str(val,:));

str=get(handles.file_antred,'String');
val=get(handles.file_antred,'Value');
paths.antred=strcat(str(val,:));

str=get(handles.file_discont,'String');
val=get(handles.file_discont,'Value');
paths.discont=strcat(str(val,:));

str=get(handles.file_velconst,'String');
val=get(handles.file_velconst,'Value');
paths.velconst=strcat(str(val,:));

str=get(handles.file_velties,'String');
val=get(handles.file_velties,'Value');
paths.velties=strcat(str(val,:));


str=get(handles.file_datumsou,'String');
val=get(handles.file_datumsou,'Value');
paths.datumsou=strcat(str(val,:));

str=get(handles.file_fixedsou,'String');
val=get(handles.file_fixedsou,'Value');
paths.fixedsou=strcat(str(val,:));

str=get(handles.file_soured,'String');
val=get(handles.file_soured,'Value');
paths.soured=strcat(str(val,:));

if parGS(g.g_coord(1)).id==0 && parGS(g.g_vel(1)).id==0
    paths.datumsta=[''];
end

valred=get(handles.cb_antred,'Value');
if parGS(g.g_coord(1)).id==0 && parGS(g.g_vel(1)).id==0 || valred==0
    paths.antred=[''];
end

valvel=get(handles.cb_velconst,'Value');
if parGS(g.g_coord(1)).id==0 && parGS(g.g_vel(1)).id==0 || valvel==0
    paths.velconst=[''];
end

valvelties=get(handles.cb_velties,'Value');
if parGS(g.g_coord(1)).id==0 && parGS(g.g_vel(1)).id==0 || valvelties==0
    paths.velties=[''];
end

if get(handles.cb_datumsou,'Value') ==0 && get(handles.cb_nnrsou,'Value') ==0
    paths.datumsou=[''];
elseif parGS(g.g_srade(1)).id==0
    paths.datumsou=[''];
end

% station coordinates: NNT/NNR or NNT/NNR/NNS
paths.datumst_scale=[''];
if get(handles.cb_trf_6p,'Value')==1
    paths.datumst_scale=['0']; % NNT/NNR
elseif get(handles.cb_trf_7p,'Value')==1
    paths.datumst_scale=['1']; %NNT/NNR/NNS
end

% source coordinates: NNR+dz or NNR
paths.datumsou_dz=[''];
if get(handles.cb_datumsou,'Value') ==1
    paths.datumsou_dz=['1']; % NNR+dz
elseif get(handles.cb_nnrsou,'Value') ==1
    paths.datumsou_dz=['0']; % NNR
end

if get(handles.cb_fixedsou,'Value') ==0 || parGS(g.g_srade(1)).id==0
    paths.fixedsou=[''];
end

if get(handles.cb_soured,'Value') ==0 || parGS(g.g_srade(1)).id==0
    paths.soured=[''];
end

save('../../OUT/GLOB/paths','paths')

close




% --- Executes on selection change in file_soured.
function file_soured_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function file_soured_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


