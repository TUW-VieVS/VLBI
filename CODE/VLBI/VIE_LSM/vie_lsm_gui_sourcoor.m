% ************************************************************************
%   Description:
%   function to create a GUI for SESSION SPECIFIC parameterization for
%   estimates of sources' coordinates.
%
%   Reference: 
%
%   Input:	
%       'opt'        structure array     (for info. /DOC/opt.doc)
%
%   Output:
%       'opt'        structure array     (for info. /DOC/opt.doc)
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   30 Aug 2009 by Kamil Teke
%
%   Revision: 
%   06 Dec 2009 by Kamil Teke: header added
%   28 Feb 2017 by Andreas Hellerschmied: opt is stored in handles.opt
% ************************************************************************
function varargout = vie_lsm_gui_sourcoor(varargin)
% VIE_LSM_GUI_SOURCOOR M-file for vie_lsm_gui_sourcoor.fig
%      VIE_LSM_GUI_SOURCOOR, by itself, creates a new VIE_LSM_GUI_SOURCOOR or raises the existing
%      singleton*.
%
%      H = VIE_LSM_GUI_SOURCOOR returns the handle to a new VIE_LSM_GUI_SOURCOOR or the handle to
%      the existing singleton*.
%
%      VIE_LSM_GUI_SOURCOOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIE_LSM_GUI_SOURCOOR.M with the given input arguments.
%
%      VIE_LSM_GUI_SOURCOOR('Property','Value',...) creates a new VIE_LSM_GUI_SOURCOOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vie_lsm_gui_sourcoor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vie_lsm_gui_sourcoor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vie_lsm_gui_sourcoor

% Last Modified by GUIDE v2.5 24-Sep-2009 15:11:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vie_lsm_gui_sourcoor_OpeningFcn, ...
                   'gui_OutputFcn',  @vie_lsm_gui_sourcoor_OutputFcn, ...
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


% --- Executes just before vie_lsm_gui_sourcoor is made visible.
function vie_lsm_gui_sourcoor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vie_lsm_gui_sourcoor (see VARARGIN)

% Choose default command line output for vie_lsm_gui_sourcoor
handles.output = hObject;

handles.opt = varargin{1};


set(handles.cb_opt_pw_sou,'Value',handles.opt.pw_sou)
if get(handles.cb_opt_pw_sou,'Value') == 0
    set(handles.cb_opt_constr_sou,'Visible','off')
end
set(handles.cb_opt_constr_sou,'Value',handles.opt.constr_sou)

D = get(handles.uitable1,'Data');
for isou = 1 : length(handles.opt.source)
    D{isou,1} = handles.opt.source(isou).name;
    D{isou,2} = handles.opt.source(isou).total_obs;
    D{isou,3} = (handles.opt.source(isou).rade_inc == 1);
    D{isou,4} = handles.opt.source(isou).coef_rade;
    D{isou,5} = handles.opt.source(isou).int_rade;
end
set(handles.uitable1,'Data',D);

% UIWAIT makes vie_lsm_gui_sourcoor wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = vie_lsm_gui_sourcoor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in cb_opt_pw_sou.
function cb_opt_pw_sou_Callback(hObject, eventdata, handles)
% hObject    handle to cb_opt_pw_sou (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.opt.pw_sou = get(handles.cb_opt_pw_sou,'Value');
if get(handles.cb_opt_pw_sou,'Value') == 1
    set(handles.cb_opt_constr_sou,'Visible','on')
end
if get(handles.cb_opt_pw_sou,'Value') == 0
    set(handles.cb_opt_constr_sou,'Visible','off')
end
% Hint: get(hObject,'Value') returns toggle state of cb_opt_pw_sou
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in cb_opt_constr_sou.
function cb_opt_constr_sou_Callback(hObject, eventdata, handles)
% hObject    handle to cb_opt_constr_sou (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.opt.constr_sou = get(handles.cb_opt_constr_sou,'Value');
% Hint: get(hObject,'Value') returns toggle state of cb_opt_constr_sou
% Update handles structure
guidata(hObject, handles);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% D = get(handles.uitable1,'Data');
% for isou = 1 : length(handles.opt.source)
%     handles.opt.source(isou).rade_inc = D{isou,3};
%     handles.opt.source(isou).coef_rade = D{isou,4};
%     handles.opt.source(isou).int_rade = D{isou,5};
% end
% % Update handles structure
% guidata(hObject, handles);

function handles = save_uitables(handles)
D = get(handles.uitable1,'Data');
for isou = 1 : length(handles.opt.source)
    handles.opt.source(isou).rade_inc = D{isou,3};
    handles.opt.source(isou).coef_rade = D{isou,4};
    handles.opt.source(isou).int_rade = D{isou,5};
end

% --- Executes on button press in pb_next6.
function pb_next6_Callback(hObject, eventdata, handles)
% hObject    handle to pb_next6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = save_uitables(handles);
close
uiwait(vie_lsm_gui_global(handles.opt));


% --- Executes on button press in pb_back6.
function pb_back6_Callback(hObject, eventdata, handles)
% hObject    handle to pb_back6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = save_uitables(handles);

close
uiwait(vie_lsm_gui_eop(handles.opt));
