function varargout = guiglob_hl(varargin)
% GUIGLOB_HL M-file for guiglob_hl.fig
%      GUIGLOB_HL, by itself, creates a new GUIGLOB_HL or raises the existing
%      singleton*.
%
%      H = GUIGLOB_HL returns the handle to a new GUIGLOB_HL or the handle to
%      the existing singleton*.
%
%      GUIGLOB_HL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIGLOB_HL.M with the given input arguments.
%
%      GUIGLOB_HL('Property','Value',...) creates a new GUIGLOB_HL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guiglob_hl_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guiglob_hl_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guiglob_hl

% Last Modified by GUIDE v2.5 25-Jun-2012 16:55:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guiglob_hl_OpeningFcn, ...
                   'gui_OutputFcn',  @guiglob_hl_OutputFcn, ...
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



% --- Executes just before guiglob_hl is made visible.
function guiglob_hl_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guiglob_hl (see VARARGIN)

% Choose default command line output for guiglob_hl
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guiglob_hl wait for user response (see UIRESUME)
% uiwait(handles.figure1);
path_level='../';
load([path_level 'DATA/GLOB/parGS'],'parGS');


[g] = globind(parGS);

if parGS(g.g_love).id==0
    set(handles.text2,'Visible','off')
    set(handles.text3,'Visible','off')
    set(handles.cb_h2,'Visible','off')
    set(handles.cb_h_diurnal,'Visible','off')
    set(handles.cb_h_diurnal_imag,'Visible','off')
    set(handles.cb_h_long,'Visible','off')
    set(handles.cb_h_long_imag,'Visible','off')
end

if parGS(g.g_shida).id==0
    set(handles.text6,'Visible','off')
    set(handles.text7,'Visible','off')
    set(handles.cb_l2,'Visible','off')
    set(handles.cb_l_diurnal,'Visible','off')
    set(handles.cb_l_diurnal_imag,'Visible','off')
    set(handles.cb_l_long,'Visible','off')
    set(handles.cb_l_long_imag,'Visible','off')
end




% --- Outputs from this function are returned to the command line.
function varargout = guiglob_hl_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;





% --- Executes on button press in cb_h2.
function cb_h2_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox12.
function checkbox12_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox13.
function checkbox13_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox14.
function checkbox14_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox15.
function checkbox15_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox16.
function checkbox16_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox17.
function checkbox17_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox18.
function checkbox18_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox19.
function checkbox19_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox20.
function checkbox20_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox21.
function checkbox21_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox22.
function checkbox22_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox23.
function checkbox23_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox24.
function checkbox24_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox25.
function checkbox25_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox26.
function checkbox26_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox27.
function checkbox27_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox28.
function checkbox28_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox29.
function checkbox29_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox30.
function checkbox30_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox31.
function checkbox31_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox32.
function checkbox32_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox33.
function checkbox33_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox34.
function checkbox34_Callback(hObject, eventdata, handles)




% --- Executes on button press in cb_h_diurnal.
function cb_h_diurnal_Callback(hObject, eventdata, handles)
path_level='../';
load([path_level 'DATA/GLOB/parGS'],'parGS');
[g] = globind(parGS);


if get(handles.cb_h_diurnal,'Value') ==1  || get(handles.cb_l_diurnal,'Value') ==1 || get(handles.cb_h_diurnal_imag,'Value') ==1 || get(handles.cb_l_diurnal_imag,'Value') ==1
    set(handles.checkbox3,'Visible','on')
    set(handles.checkbox4,'Visible','on')
    set(handles.checkbox5,'Visible','on')
    set(handles.checkbox6,'Visible','on')
    set(handles.checkbox7,'Visible','on')
    set(handles.checkbox8,'Visible','on')
    set(handles.checkbox9,'Visible','on')
    set(handles.checkbox10,'Visible','on')
    set(handles.checkbox11,'Visible','on')
    set(handles.checkbox12,'Visible','on')
    set(handles.checkbox13,'Visible','on')
    set(handles.checkbox14,'Visible','on')
    set(handles.checkbox15,'Visible','on')
    set(handles.checkbox16,'Visible','on')
    set(handles.checkbox17,'Visible','on')
    set(handles.checkbox18,'Visible','on')
    set(handles.checkbox19,'Visible','on')
    set(handles.checkbox20,'Visible','on')
    set(handles.checkbox21,'Visible','on')
    set(handles.checkbox22,'Visible','on')
    set(handles.checkbox23,'Visible','on')
    set(handles.checkbox24,'Visible','on')
    set(handles.checkbox25,'Visible','on')
    set(handles.checkbox26,'Visible','on')
    set(handles.checkbox27,'Visible','on')
    set(handles.checkbox28,'Visible','on')
    set(handles.checkbox29,'Visible','on')
    set(handles.checkbox30,'Visible','on')
    set(handles.checkbox31,'Visible','on')
    set(handles.checkbox32,'Visible','on')
    set(handles.checkbox33,'Visible','on')
else
    set(handles.checkbox3,'Visible','off')
    set(handles.checkbox4,'Visible','off')
    set(handles.checkbox5,'Visible','off')
    set(handles.checkbox6,'Visible','off')
    set(handles.checkbox7,'Visible','off')
    set(handles.checkbox8,'Visible','off')
    set(handles.checkbox9,'Visible','off')
    set(handles.checkbox10,'Visible','off')
    set(handles.checkbox11,'Visible','off')
    set(handles.checkbox12,'Visible','off')
    set(handles.checkbox13,'Visible','off')
    set(handles.checkbox14,'Visible','off')
    set(handles.checkbox15,'Visible','off')
    set(handles.checkbox16,'Visible','off')
    set(handles.checkbox17,'Visible','off')
    set(handles.checkbox18,'Visible','off')
    set(handles.checkbox19,'Visible','off')
    set(handles.checkbox20,'Visible','off')
    set(handles.checkbox21,'Visible','off')
    set(handles.checkbox22,'Visible','off')
    set(handles.checkbox23,'Visible','off')
    set(handles.checkbox24,'Visible','off')
    set(handles.checkbox25,'Visible','off')
    set(handles.checkbox26,'Visible','off')
    set(handles.checkbox27,'Visible','off')
    set(handles.checkbox28,'Visible','off')
    set(handles.checkbox29,'Visible','off')
    set(handles.checkbox30,'Visible','off')
    set(handles.checkbox31,'Visible','off')
    set(handles.checkbox32,'Visible','off')
    set(handles.checkbox33,'Visible','off')
end



% --- Executes on button press in cb_l2.
function cb_l2_Callback(hObject, eventdata, handles)


% --- Executes on button press in cb_l_diurnal.
function cb_l_diurnal_Callback(hObject, eventdata, handles)
path_level='../';
load([path_level 'DATA/GLOB/parGS'],'parGS');
[g] = globind(parGS);




if get(handles.cb_l_diurnal,'Value') ==1 || get(handles.cb_h_diurnal,'Value') ==1 || get(handles.cb_h_diurnal_imag,'Value') ==1 || get(handles.cb_l_diurnal_imag,'Value') ==1
    set(handles.checkbox3,'Visible','on')
    set(handles.checkbox4,'Visible','on')
    set(handles.checkbox5,'Visible','on')
    set(handles.checkbox6,'Visible','on')
    set(handles.checkbox7,'Visible','on')
    set(handles.checkbox8,'Visible','on')
    set(handles.checkbox9,'Visible','on')
    set(handles.checkbox10,'Visible','on')
    set(handles.checkbox11,'Visible','on')
    set(handles.checkbox12,'Visible','on')
    set(handles.checkbox13,'Visible','on')
    set(handles.checkbox14,'Visible','on')
    set(handles.checkbox15,'Visible','on')
    set(handles.checkbox16,'Visible','on')
    set(handles.checkbox17,'Visible','on')
    set(handles.checkbox18,'Visible','on')
    set(handles.checkbox19,'Visible','on')
    set(handles.checkbox20,'Visible','on')
    set(handles.checkbox21,'Visible','on')
    set(handles.checkbox22,'Visible','on')
    set(handles.checkbox23,'Visible','on')
    set(handles.checkbox24,'Visible','on')
    set(handles.checkbox25,'Visible','on')
    set(handles.checkbox26,'Visible','on')
    set(handles.checkbox27,'Visible','on')
    set(handles.checkbox28,'Visible','on')
    set(handles.checkbox29,'Visible','on')
    set(handles.checkbox30,'Visible','on')
    set(handles.checkbox31,'Visible','on')
    set(handles.checkbox32,'Visible','on')
    set(handles.checkbox33,'Visible','on')
else
    set(handles.checkbox3,'Visible','off')
    set(handles.checkbox4,'Visible','off')
    set(handles.checkbox5,'Visible','off')
    set(handles.checkbox6,'Visible','off')
    set(handles.checkbox7,'Visible','off')
    set(handles.checkbox8,'Visible','off')
    set(handles.checkbox9,'Visible','off')
    set(handles.checkbox10,'Visible','off')
    set(handles.checkbox11,'Visible','off')
    set(handles.checkbox12,'Visible','off')
    set(handles.checkbox13,'Visible','off')
    set(handles.checkbox14,'Visible','off')
    set(handles.checkbox15,'Visible','off')
    set(handles.checkbox16,'Visible','off')
    set(handles.checkbox17,'Visible','off')
    set(handles.checkbox18,'Visible','off')
    set(handles.checkbox19,'Visible','off')
    set(handles.checkbox20,'Visible','off')
    set(handles.checkbox21,'Visible','off')
    set(handles.checkbox22,'Visible','off')
    set(handles.checkbox23,'Visible','off')
    set(handles.checkbox24,'Visible','off')
    set(handles.checkbox25,'Visible','off')
    set(handles.checkbox26,'Visible','off')
    set(handles.checkbox27,'Visible','off')
    set(handles.checkbox28,'Visible','off')
    set(handles.checkbox29,'Visible','off')
    set(handles.checkbox30,'Visible','off')
    set(handles.checkbox31,'Visible','off')
    set(handles.checkbox32,'Visible','off')
    set(handles.checkbox33,'Visible','off')
end

% --- Executes on button press in cb_h_diurnal_imag.
function cb_h_diurnal_imag_Callback(hObject, eventdata, handles)
% hObject    handle to cb_h_diurnal_imag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_h_diurnal_imag

if get(handles.cb_h_diurnal,'Value') ==1  || get(handles.cb_l_diurnal,'Value') ==1 || get(handles.cb_h_diurnal_imag,'Value') ==1 || get(handles.cb_l_diurnal_imag,'Value') ==1
    set(handles.checkbox3,'Visible','on')
    set(handles.checkbox4,'Visible','on')
    set(handles.checkbox5,'Visible','on')
    set(handles.checkbox6,'Visible','on')
    set(handles.checkbox7,'Visible','on')
    set(handles.checkbox8,'Visible','on')
    set(handles.checkbox9,'Visible','on')
    set(handles.checkbox10,'Visible','on')
    set(handles.checkbox11,'Visible','on')
    set(handles.checkbox12,'Visible','on')
    set(handles.checkbox13,'Visible','on')
    set(handles.checkbox14,'Visible','on')
    set(handles.checkbox15,'Visible','on')
    set(handles.checkbox16,'Visible','on')
    set(handles.checkbox17,'Visible','on')
    set(handles.checkbox18,'Visible','on')
    set(handles.checkbox19,'Visible','on')
    set(handles.checkbox20,'Visible','on')
    set(handles.checkbox21,'Visible','on')
    set(handles.checkbox22,'Visible','on')
    set(handles.checkbox23,'Visible','on')
    set(handles.checkbox24,'Visible','on')
    set(handles.checkbox25,'Visible','on')
    set(handles.checkbox26,'Visible','on')
    set(handles.checkbox27,'Visible','on')
    set(handles.checkbox28,'Visible','on')
    set(handles.checkbox29,'Visible','on')
    set(handles.checkbox30,'Visible','on')
    set(handles.checkbox31,'Visible','on')
    set(handles.checkbox32,'Visible','on')
    set(handles.checkbox33,'Visible','on')
else
    set(handles.checkbox3,'Visible','off')
    set(handles.checkbox4,'Visible','off')
    set(handles.checkbox5,'Visible','off')
    set(handles.checkbox6,'Visible','off')
    set(handles.checkbox7,'Visible','off')
    set(handles.checkbox8,'Visible','off')
    set(handles.checkbox9,'Visible','off')
    set(handles.checkbox10,'Visible','off')
    set(handles.checkbox11,'Visible','off')
    set(handles.checkbox12,'Visible','off')
    set(handles.checkbox13,'Visible','off')
    set(handles.checkbox14,'Visible','off')
    set(handles.checkbox15,'Visible','off')
    set(handles.checkbox16,'Visible','off')
    set(handles.checkbox17,'Visible','off')
    set(handles.checkbox18,'Visible','off')
    set(handles.checkbox19,'Visible','off')
    set(handles.checkbox20,'Visible','off')
    set(handles.checkbox21,'Visible','off')
    set(handles.checkbox22,'Visible','off')
    set(handles.checkbox23,'Visible','off')
    set(handles.checkbox24,'Visible','off')
    set(handles.checkbox25,'Visible','off')
    set(handles.checkbox26,'Visible','off')
    set(handles.checkbox27,'Visible','off')
    set(handles.checkbox28,'Visible','off')
    set(handles.checkbox29,'Visible','off')
    set(handles.checkbox30,'Visible','off')
    set(handles.checkbox31,'Visible','off')
    set(handles.checkbox32,'Visible','off')
    set(handles.checkbox33,'Visible','off')
end
    

% --- Executes on button press in cb_l_diurnal_imag.
function cb_l_diurnal_imag_Callback(hObject, eventdata, handles)
if get(handles.cb_h_diurnal,'Value') ==1  || get(handles.cb_l_diurnal,'Value') ==1 || get(handles.cb_h_diurnal_imag,'Value') ==1 || get(handles.cb_l_diurnal_imag,'Value') ==1
    set(handles.checkbox3,'Visible','on')
    set(handles.checkbox4,'Visible','on')
    set(handles.checkbox5,'Visible','on')
    set(handles.checkbox6,'Visible','on')
    set(handles.checkbox7,'Visible','on')
    set(handles.checkbox8,'Visible','on')
    set(handles.checkbox9,'Visible','on')
    set(handles.checkbox10,'Visible','on')
    set(handles.checkbox11,'Visible','on')
    set(handles.checkbox12,'Visible','on')
    set(handles.checkbox13,'Visible','on')
    set(handles.checkbox14,'Visible','on')
    set(handles.checkbox15,'Visible','on')
    set(handles.checkbox16,'Visible','on')
    set(handles.checkbox17,'Visible','on')
    set(handles.checkbox18,'Visible','on')
    set(handles.checkbox19,'Visible','on')
    set(handles.checkbox20,'Visible','on')
    set(handles.checkbox21,'Visible','on')
    set(handles.checkbox22,'Visible','on')
    set(handles.checkbox23,'Visible','on')
    set(handles.checkbox24,'Visible','on')
    set(handles.checkbox25,'Visible','on')
    set(handles.checkbox26,'Visible','on')
    set(handles.checkbox27,'Visible','on')
    set(handles.checkbox28,'Visible','on')
    set(handles.checkbox29,'Visible','on')
    set(handles.checkbox30,'Visible','on')
    set(handles.checkbox31,'Visible','on')
    set(handles.checkbox32,'Visible','on')
    set(handles.checkbox33,'Visible','on')
else
    set(handles.checkbox3,'Visible','off')
    set(handles.checkbox4,'Visible','off')
    set(handles.checkbox5,'Visible','off')
    set(handles.checkbox6,'Visible','off')
    set(handles.checkbox7,'Visible','off')
    set(handles.checkbox8,'Visible','off')
    set(handles.checkbox9,'Visible','off')
    set(handles.checkbox10,'Visible','off')
    set(handles.checkbox11,'Visible','off')
    set(handles.checkbox12,'Visible','off')
    set(handles.checkbox13,'Visible','off')
    set(handles.checkbox14,'Visible','off')
    set(handles.checkbox15,'Visible','off')
    set(handles.checkbox16,'Visible','off')
    set(handles.checkbox17,'Visible','off')
    set(handles.checkbox18,'Visible','off')
    set(handles.checkbox19,'Visible','off')
    set(handles.checkbox20,'Visible','off')
    set(handles.checkbox21,'Visible','off')
    set(handles.checkbox22,'Visible','off')
    set(handles.checkbox23,'Visible','off')
    set(handles.checkbox24,'Visible','off')
    set(handles.checkbox25,'Visible','off')
    set(handles.checkbox26,'Visible','off')
    set(handles.checkbox27,'Visible','off')
    set(handles.checkbox28,'Visible','off')
    set(handles.checkbox29,'Visible','off')
    set(handles.checkbox30,'Visible','off')
    set(handles.checkbox31,'Visible','off')
    set(handles.checkbox32,'Visible','off')
    set(handles.checkbox33,'Visible','off')
end   
    

% --- Executes on button press in checkbox79.
function checkbox79_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox81.
function checkbox81_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox82.
function checkbox82_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox83.
function checkbox83_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox84.
function checkbox84_Callback(hObject, eventdata, handles)



% --- Executes on button press in cb_h_long.
function cb_h_long_Callback(hObject, eventdata, handles)
if get(handles.cb_h_long,'Value') ==1
    set(handles.checkbox79,'Visible','on')
    set(handles.checkbox81,'Visible','on')
    set(handles.checkbox82,'Visible','on')
    set(handles.checkbox83,'Visible','on')
    set(handles.checkbox84,'Visible','on')
end
if get(handles.cb_h_long,'Value') ==0 && get(handles.cb_h_long_imag,'Value') ==0 && get(handles.cb_l_long,'Value') ==0 && get(handles.cb_l_long_imag,'Value') ==0
    set(handles.checkbox79,'Visible','off')
    set(handles.checkbox81,'Visible','off')
    set(handles.checkbox82,'Visible','off')
    set(handles.checkbox83,'Visible','off')
    set(handles.checkbox84,'Visible','off')
end

% --- Executes on button press in cb_h_long_imag.
function cb_h_long_imag_Callback(hObject, eventdata, handles)
if get(handles.cb_h_long_imag,'Value') ==1
    set(handles.checkbox79,'Visible','on')
    set(handles.checkbox81,'Visible','on')
    set(handles.checkbox82,'Visible','on')
    set(handles.checkbox83,'Visible','on')
    set(handles.checkbox84,'Visible','on')
end
if get(handles.cb_h_long,'Value') ==0 && get(handles.cb_h_long_imag,'Value') ==0 && get(handles.cb_l_long,'Value') ==0 && get(handles.cb_l_long_imag,'Value') ==0
    set(handles.checkbox79,'Visible','off')
    set(handles.checkbox81,'Visible','off')
    set(handles.checkbox82,'Visible','off')
    set(handles.checkbox83,'Visible','off')
    set(handles.checkbox84,'Visible','off')
end


% --- Executes on button press in cb_l_long.
function cb_l_long_Callback(hObject, eventdata, handles)
if get(handles.cb_l_long,'Value') ==1
    set(handles.checkbox79,'Visible','on')
    set(handles.checkbox81,'Visible','on')
    set(handles.checkbox82,'Visible','on')
    set(handles.checkbox83,'Visible','on')
    set(handles.checkbox84,'Visible','on')
end
if get(handles.cb_h_long,'Value') ==0 && get(handles.cb_h_long_imag,'Value') ==0 && get(handles.cb_l_long,'Value') ==0 && get(handles.cb_l_long_imag,'Value') ==0
    set(handles.checkbox79,'Visible','off')
    set(handles.checkbox81,'Visible','off')
    set(handles.checkbox82,'Visible','off')
    set(handles.checkbox83,'Visible','off')
    set(handles.checkbox84,'Visible','off')
end

% --- Executes on button press in cb_l_long_imag.
function cb_l_long_imag_Callback(hObject, eventdata, handles)
if get(handles.cb_l_long_imag,'Value') ==1
    set(handles.checkbox79,'Visible','on')
    set(handles.checkbox81,'Visible','on')
    set(handles.checkbox82,'Visible','on')
    set(handles.checkbox83,'Visible','on')
    set(handles.checkbox84,'Visible','on')
end
if get(handles.cb_h_long,'Value') ==0 && get(handles.cb_h_long_imag,'Value') ==0 && get(handles.cb_l_long,'Value') ==0 && get(handles.cb_l_long_imag,'Value') ==0
    set(handles.checkbox79,'Visible','off')
    set(handles.checkbox81,'Visible','off')
    set(handles.checkbox82,'Visible','off')
    set(handles.checkbox83,'Visible','off')
    set(handles.checkbox84,'Visible','off')
end




% --- Executes on button press in pb_ok.
function pb_ok_Callback(hObject, eventdata, handles)
% hObject    handle to pb_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path_level='../';
load([path_level 'DATA/GLOB/parGS'],'parGS');
[g] = globind(parGS);

parGS_hl.love = [];
    parGS_hl.love.nr=[];
    parGS_hl.love.nr_d=[];
    parGS_hl.love.nr_di=[];
    parGS_hl.love.cond=[];
    parGS_hl.love.nr_long=[];
    parGS_hl.love.nr_longi=[];

parGS_hl.shida =[];
    parGS_hl.shida.nr=[];
    parGS_hl.shida.nr_d=[];
    parGS_hl.shida.nr_di=[];
    parGS_hl.shida.cond=[];
    parGS_hl.shida.nr_long=[];
    parGS_hl.shida.nr_longi=[];
    
parGS_hl.NDFW_cond = [];

if parGS(g.g_love).id==1
    h2=[];
    dih2=[];
    longh2=[];
    
    if get(handles.cb_h2,'Value')==1
        h2=1;
    end


    if get(handles.cb_h_diurnal,'Value') ==1 || get(handles.cb_h_diurnal_imag,'Value') ==1 
        diu_h2(1)=get(handles.checkbox3,'Value');
        diu_h2(2)=get(handles.checkbox4,'Value');
        diu_h2(3)=get(handles.checkbox5,'Value');
        diu_h2(4)=get(handles.checkbox6,'Value');
        diu_h2(5)=get(handles.checkbox7,'Value');
        diu_h2(6)=get(handles.checkbox8,'Value');
        diu_h2(7)=get(handles.checkbox9,'Value');
        diu_h2(8)=get(handles.checkbox10,'Value');
        diu_h2(9)=get(handles.checkbox11,'Value');
        diu_h2(10)=get(handles.checkbox12,'Value');
        diu_h2(11)=get(handles.checkbox13,'Value');
        diu_h2(12)=get(handles.checkbox14,'Value');
        diu_h2(13)=get(handles.checkbox15,'Value');
        diu_h2(14)=get(handles.checkbox16,'Value');
        diu_h2(15)=get(handles.checkbox17,'Value');
        diu_h2(16)=get(handles.checkbox18,'Value');
        diu_h2(17)=get(handles.checkbox19,'Value');
        diu_h2(18)=get(handles.checkbox20,'Value');
        diu_h2(19)=get(handles.checkbox21,'Value');
        diu_h2(20)=get(handles.checkbox22,'Value');
        diu_h2(21)=get(handles.checkbox23,'Value');
        diu_h2(22)=get(handles.checkbox24,'Value');
        diu_h2(23)=get(handles.checkbox25,'Value');
        diu_h2(24)=get(handles.checkbox26,'Value');
        diu_h2(25)=get(handles.checkbox27,'Value');
        diu_h2(26)=get(handles.checkbox28,'Value');
        diu_h2(27)=get(handles.checkbox29,'Value');
        diu_h2(28)=get(handles.checkbox30,'Value');
        diu_h2(29)=get(handles.checkbox31,'Value');
        diu_h2(30)=get(handles.checkbox32,'Value');
        diu_h2(31)=get(handles.checkbox33,'Value');

        [di1,dih2]=find(diu_h2>0);
        if get(handles.cb_h_diurnal,'Value') ==1
            parGS_hl.love.nr_d=dih2 + 5; % only real diurnal tides (5 must be added to the index)
        end
        if get(handles.cb_h_diurnal_imag,'Value') ==1
            parGS_hl.love.nr_di=dih2 + 36; % only imaginary diurnal tides (36 must be added to the index)
        end
    end
    
    if get(handles.cb_h_long,'Value') ==1 || get(handles.cb_h_long_imag,'Value') ==1 
        lon_h2(1)=get(handles.checkbox79,'Value');
        lon_h2(2)=get(handles.checkbox81,'Value');
        lon_h2(3)=get(handles.checkbox82,'Value');
        lon_h2(4)=get(handles.checkbox83,'Value');
        lon_h2(5)=get(handles.checkbox84,'Value');
        
        [lo1,loh2]=find(lon_h2>0);

        if get(handles.cb_h_long,'Value')==1
            parGS_hl.love.nr_long=[loh2 + 67];
        end
        if get(handles.cb_h_long_imag,'Value')==1
            parGS_hl.love.nr_longi=[loh2 + 72];
        end
    end
    
    parGS_hl.love.nr= [h2 parGS_hl.love.nr_d parGS_hl.love.nr_di parGS_hl.love.nr_long parGS_hl.love.nr_longi];  %[1 [1:11 24 26 28 29]+5]; % all love numbers which should be estimated
    
end




if parGS(g.g_shida).id==1
    l2=[];
    dil2=[];

    if get(handles.cb_l2,'Value')==1
        l2=1;
    end


    if get(handles.cb_l_diurnal,'Value') ==1 || get(handles.cb_l_diurnal_imag,'Value') ==1 
        diu_l2(1)=get(handles.checkbox3,'Value');
        diu_l2(2)=get(handles.checkbox4,'Value');
        diu_l2(3)=get(handles.checkbox5,'Value');
        diu_l2(4)=get(handles.checkbox6,'Value');
        diu_l2(5)=get(handles.checkbox7,'Value');
        diu_l2(6)=get(handles.checkbox8,'Value');
        diu_l2(7)=get(handles.checkbox9,'Value');
        diu_l2(8)=get(handles.checkbox10,'Value');
        diu_l2(9)=get(handles.checkbox11,'Value');
        diu_l2(10)=get(handles.checkbox12,'Value');
        diu_l2(11)=get(handles.checkbox13,'Value');
        diu_l2(12)=get(handles.checkbox14,'Value');
        diu_l2(13)=get(handles.checkbox15,'Value');
        diu_l2(14)=get(handles.checkbox16,'Value');
        diu_l2(15)=get(handles.checkbox17,'Value');
        diu_l2(16)=get(handles.checkbox18,'Value');
        diu_l2(17)=get(handles.checkbox19,'Value');
        diu_l2(18)=get(handles.checkbox20,'Value');
        diu_l2(19)=get(handles.checkbox21,'Value');
        diu_l2(20)=get(handles.checkbox22,'Value');
        diu_l2(21)=get(handles.checkbox23,'Value');
        diu_l2(22)=get(handles.checkbox24,'Value');
        diu_l2(23)=get(handles.checkbox25,'Value');
        diu_l2(24)=get(handles.checkbox26,'Value');
        diu_l2(25)=get(handles.checkbox27,'Value');
        diu_l2(26)=get(handles.checkbox28,'Value');
        diu_l2(27)=get(handles.checkbox29,'Value');
        diu_l2(28)=get(handles.checkbox30,'Value');
        diu_l2(29)=get(handles.checkbox31,'Value');
        diu_l2(30)=get(handles.checkbox32,'Value');
        diu_l2(31)=get(handles.checkbox33,'Value');

        [di1,dil2]=find(diu_l2>0);
        if get (handles.cb_l_diurnal,'Value')==1
            parGS_hl.shida.nr_d=dil2 + 7; % only real diurnal tides(7 must be added to the index)
        end
        if get(handles.cb_l_diurnal_imag,'Value') ==1
            parGS_hl.shida.nr_di=dil2 + 38; % only imaginary diurnal tides (38 must be added to the index)
        end
    end
    
    if get(handles.cb_l_long,'Value') ==1 || get(handles.cb_l_long_imag,'Value') ==1 
        lon_l2(1)=get(handles.checkbox79,'Value');
        lon_l2(2)=get(handles.checkbox81,'Value');
        lon_l2(3)=get(handles.checkbox82,'Value');
        lon_l2(4)=get(handles.checkbox83,'Value');
        lon_l2(5)=get(handles.checkbox84,'Value');
        
        [lo1,lol2]=find(lon_l2>0);

        if get(handles.cb_l_long,'Value')==1
            parGS_hl.shida.nr_long=[lol2 + 69];
        end
        if get(handles.cb_l_long_imag,'Value')==1
            parGS_hl.shida.nr_longi=[lol2 + 74];
        end
    end
   
    parGS_hl.shida.nr= [l2 parGS_hl.shida.nr_d parGS_hl.shida.nr_di parGS_hl.shida.nr_long parGS_hl.shida.nr_longi];  %[1 [1:11 24 26 28 29]+5]; % all love numbers which should be estimated
end


path_level='../';
save([path_level 'DATA/GLOB/parGS_hl'],'parGS_hl');

close




