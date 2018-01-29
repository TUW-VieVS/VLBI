% ************************************************************************
%   Description:
%   function to create a GUI for SESSION SPECIFIC clock parameterization. 
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
function varargout = vie_lsm_gui_clock(varargin)
% VIE_LSM_GUI_CLOCK M-file for vie_lsm_gui_clock.fig
%      VIE_LSM_GUI_CLOCK, by itself, creates a new VIE_LSM_GUI_CLOCK or raises the existing
%      singleton*.
%
%      H = VIE_LSM_GUI_CLOCK returns the handle to a new VIE_LSM_GUI_CLOCK or the handle to
%      the existing singleton*.
%
%      VIE_LSM_GUI_CLOCK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIE_LSM_GUI_CLOCK.M with the given input arguments.
%
%      VIE_LSM_GUI_CLOCK('Property','Value',...) creates a new VIE_LSM_GUI_CLOCK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vie_lsm_gui_clock_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vie_lsm_gui_clock_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vie_lsm_gui_clock

% Last Modified by GUIDE v2.5 24-Sep-2009 14:55:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vie_lsm_gui_clock_OpeningFcn, ...
                   'gui_OutputFcn',  @vie_lsm_gui_clock_OutputFcn, ...
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


% --- Executes just before vie_lsm_gui_clock is made visible.
function vie_lsm_gui_clock_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vie_lsm_gui_clock (see VARARGIN)

% Choose default command line output for vie_lsm_gui_clock
handles.output = hObject;

handles.opt = varargin{1};

D = get(handles.uitable1,'Data');
for istat = 1 : length(handles.opt.stat)
    D{istat,1}=handles.opt.stat(istat).coef_clk;
    D{istat,2}=handles.opt.stat(istat).int_clk;   
    D{istat,3}=(handles.opt.stat(istat).ref == 1); 
    stat_name(istat,1:8) = handles.opt.stat(istat).name(1:8);
end
set(handles.uitable1,'Data',D);
set(handles.uitable1,'RowName',stat_name);

% UIWAIT makes vie_lsm_gui_clock wait for user response (see UIRESUME)
% uiwait(handles.figure1);

if handles.opt.pw_clk == 0
    set(handles.cb_opt_pw_clk0,'Value',0)
    set(handles.rb_opt_pw_clk1,'Visible','off')
    set(handles.rb_opt_pw_clk2,'Visible','off')
    set(handles.rb_opt_pw_clk3,'Visible','off')
end
if handles.opt.pw_clk == 1
    set(handles.cb_opt_pw_clk0,'Value',1)
    set(handles.rb_opt_pw_clk1,'Value',1)
end
if handles.opt.pw_clk == 2
    set(handles.cb_opt_pw_clk0,'Value',1)
    set(handles.rb_opt_pw_clk2,'Value',1)
end
if handles.opt.pw_clk == 3
    set(handles.cb_opt_pw_clk0,'Value',1)
    set(handles.rb_opt_pw_clk3,'Value',1)
end
set(handles.cb_opt_constr_clk,'Value',handles.opt.constr_clk);

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = vie_lsm_gui_clock_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;
guidata(hObject, handles);


% --- Executes on button press in cb_opt_pw_clk0.
function cb_opt_pw_clk0_Callback(hObject, eventdata, handles)
% hObject    handle to cb_opt_pw_clk0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.cb_opt_pw_clk0,'Value') == 0
    handles.opt.pw_clk = 0;
    set(handles.rb_opt_pw_clk1,'Visible','off')
    set(handles.rb_opt_pw_clk2,'Visible','off')
    set(handles.rb_opt_pw_clk3,'Visible','off')
end
if get(handles.cb_opt_pw_clk0,'Value') == 1
   set(handles.rb_opt_pw_clk1,'Visible','on')
   set(handles.rb_opt_pw_clk2,'Visible','on')
   set(handles.rb_opt_pw_clk3,'Visible','on')
    set(handles.rb_opt_pw_clk1,'Value',0) 
    set(handles.rb_opt_pw_clk2,'Value',0) 
    set(handles.rb_opt_pw_clk3,'Value',0)    
end
% Update handles structure
guidata(hObject, handles);


% Hint: get(hObject,'Value') returns toggle state of cb_opt_pw_clk0

% --- Executes on button press in rb_opt_pw_clk1.
function rb_opt_pw_clk1_Callback(hObject, eventdata, handles)
% hObject    handle to rb_opt_pw_clk1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.rb_opt_pw_clk1,'Value') == 1
    handles.opt.pw_clk = 1;
    set(handles.rb_opt_pw_clk2,'Value',0)
    set(handles.rb_opt_pw_clk3,'Value',0)
end
% Hint: get(hObject,'Value') returns toggle state of rb_opt_pw_clk1
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in rb_opt_pw_clk2.
function rb_opt_pw_clk2_Callback(hObject, eventdata, handles)
% hObject    handle to rb_opt_pw_clk2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.rb_opt_pw_clk2,'Value') == 1
    handles.opt.pw_clk = 2;
    set(handles.rb_opt_pw_clk1,'Value',0)
    set(handles.rb_opt_pw_clk3,'Value',0)
end
% Hint: get(hObject,'Value') returns toggle state of rb_opt_pw_clk2
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in rb_opt_pw_clk3.
function rb_opt_pw_clk3_Callback(hObject, eventdata, handles)
% hObject    handle to rb_opt_pw_clk3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.rb_opt_pw_clk3,'Value') == 1
    handles.opt.pw_clk = 3;
    set(handles.rb_opt_pw_clk1,'Value',0)
    set(handles.rb_opt_pw_clk2,'Value',0)
end
% Hint: get(hObject,'Value') returns toggle state of rb_opt_pw_clk3
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in cb_opt_constr_clk.
function cb_opt_constr_clk_Callback(hObject, eventdata, handles)
% hObject    handle to cb_opt_constr_clk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.opt.constr_clk = get(handles.cb_opt_constr_clk,'Value');
% Hint: get(hObject,'Value') returns toggle state of cb_opt_constr_clk
% Update handles structure
guidata(hObject, handles);

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function handles = save_uitables(handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
D = get(handles.uitable1,'Data');
for istat = 1 : length(handles.opt.stat)
    handles.opt.stat(istat).coef_clk = D{istat,1};
    handles.opt.stat(istat).int_clk = D{istat,2};
    handles.opt.stat(istat).ref = D{istat,3};
end


% --- Executes on button press in pb_next2.
function pb_next2_Callback(hObject, eventdata, handles)
% hObject    handle to pb_next2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = save_uitables(handles);
close
uiwait(vie_lsm_gui_tropo(handles.opt));


% --- Executes on button press in pb_back2.
function pb_back2_Callback(hObject, eventdata, handles)
% hObject    handle to pb_back2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = save_uitables(handles);
close
uiwait(vie_lsm_gui_first(handles.opt));

