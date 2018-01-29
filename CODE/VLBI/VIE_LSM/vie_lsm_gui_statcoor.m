% ************************************************************************
%   Description:
%   function to create a GUI for SESSION SPECIFIC parameterization for
%   estimates of stations' coordinates (TRF)
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
function varargout = vie_lsm_gui_statcoor(varargin)
% VIE_LSM_GUI_STATCOOR M-file for vie_lsm_gui_statcoor.fig
%      VIE_LSM_GUI_STATCOOR, by itself, creates a new VIE_LSM_GUI_STATCOOR or raises the existing
%      singleton*.
%
%      H = VIE_LSM_GUI_STATCOOR returns the handle to a new VIE_LSM_GUI_STATCOOR or the handle to
%      the existing singleton*.
%
%      VIE_LSM_GUI_STATCOOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIE_LSM_GUI_STATCOOR.M with the given input arguments.
%
%      VIE_LSM_GUI_STATCOOR('Property','Value',...) creates a new VIE_LSM_GUI_STATCOOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vie_lsm_gui_statcoor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vie_lsm_gui_statcoor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vie_lsm_gui_statcoor

% Last Modified by GUIDE v2.5 24-Sep-2009 15:11:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vie_lsm_gui_statcoor_OpeningFcn, ...
                   'gui_OutputFcn',  @vie_lsm_gui_statcoor_OutputFcn, ...
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


% --- Executes just before vie_lsm_gui_statcoor is made visible.
function vie_lsm_gui_statcoor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vie_lsm_gui_statcoor (see VARARGIN)

handles.opt = varargin{1};

% Choose default command line output for vie_lsm_gui_statcoor
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


set(handles.cb_opt_stc,'Value',handles.opt.stc)
if get(handles.cb_opt_stc,'Value') == 0
    set(handles.rb_opt_pw_stc0,'Visible','off')
    set(handles.rb_opt_pw_stc1,'Visible','off')
    set(handles.rb_opt_nnt_stc1,'Visible','off')
    set(handles.rb_fix1,'Visible','off')
    set(handles.rb_fix2,'Visible','off')
    set(handles.cb_opt_constr_xyz,'Visible','off')
end
if get(handles.cb_opt_stc,'Value') == 1
    set(handles.rb_opt_pw_stc0,'Visible','on')
    set(handles.rb_opt_pw_stc1,'Visible','on')
    set(handles.rb_opt_nnt_stc1,'Visible','off')
    set(handles.rb_fix1,'Visible','off')
    set(handles.rb_fix2,'Visible','off')
    set(handles.cb_opt_constr_xyz,'Visible','off')
    if handles.opt.pw_stc == 0
        set(handles.rb_opt_pw_stc0,'Value',1)  
        set(handles.rb_opt_nnt_stc1,'Visible','on')
        set(handles.rb_fix1,'Visible','on')
        set(handles.rb_opt_pw_stc1,'Value',0)
        if handles.opt.nnt_stc == 1
           set(handles.rb_opt_nnt_stc1,'Value',1)
           set(handles.rb_fix1,'Value',0)
        end
        if handles.opt.nnt_stc == 0
           set(handles.rb_fix1,'Value',1)
           set(handles.rb_opt_nnt_stc1,'Value',0)
        end
    end
    if handles.opt.pw_stc == 1
        handles.opt.nnt_stc = 0;
        set(handles.rb_opt_pw_stc1,'Value',1) 
        set(handles.rb_opt_pw_stc0,'Value',0)
        set(handles.rb_fix2,'Visible','on')
        set(handles.rb_fix2,'Value',1)
        set(handles.cb_opt_constr_xyz,'Visible','on')
        set(handles.cb_opt_constr_xyz,'Value',handles.opt.constr_xyz)
    end
end

D = get(handles.uitable1,'Data');
for istat = 1 : length(handles.opt.stat)
    stat_name(istat,1:8) = handles.opt.stat(istat).name(1:8);
    D{istat,1}=(handles.opt.stat(istat).nnt_inc == 1);
    D{istat,2}=(handles.opt.stat(istat).nnr_inc == 1);
    D{istat,3}=(handles.opt.stat(istat).nns_inc == 1);
    D{istat,4}=(handles.opt.stat(istat).xyz_inc == 1); 
    D{istat,5}=handles.opt.stat(istat).coef_xyz;
    D{istat,6}=handles.opt.stat(istat).int_xyz;
end
set(handles.uitable1,'Data',D);
set(handles.uitable1,'RowName',stat_name);

if get(handles.rb_opt_nnt_stc1,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[true true true false false false]);
end
if get(handles.rb_fix1,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[false false false true false false]);
end
if get(handles.rb_fix2,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[false false false true true true]);
end

% UIWAIT makes vie_lsm_gui_statcoor wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = vie_lsm_gui_statcoor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in cb_opt_stc.
function cb_opt_stc_Callback(hObject, eventdata, handles)
% hObject    handle to cb_opt_stc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.cb_opt_stc,'Value') == 0
   set(handles.rb_opt_pw_stc0,'Visible','off')
   set(handles.rb_opt_pw_stc1,'Visible','off')
   set(handles.rb_opt_nnt_stc1,'Visible','off')
   set(handles.rb_fix1,'Visible','off')
   set(handles.rb_fix2,'Visible','off')
   set(handles.cb_opt_constr_xyz,'Visible','off')
end
if get(handles.cb_opt_stc,'Value') == 1
   set(handles.rb_opt_pw_stc0,'Visible','on')
   set(handles.rb_opt_pw_stc1,'Visible','on')
   set(handles.rb_opt_nnt_stc1,'Visible','on')
   set(handles.rb_fix1,'Visible','on')
end
handles.opt.stc = get(handles.cb_opt_stc,'Value');

if get(handles.rb_opt_nnt_stc1,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[true true true false false false]);
end
if get(handles.rb_fix1,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[false false false true false false]);
end
if get(handles.rb_fix2,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[false false false true true true]);
end
% Hint: get(hObject,'Value') returns toggle state of cb_opt_stc
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in rb_opt_pw_stc0.
function rb_opt_pw_stc0_Callback(hObject, eventdata, handles)
% hObject    handle to rb_opt_pw_stc0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.rb_opt_pw_stc0,'Value') == 0   
   set(handles.rb_opt_nnt_stc1,'Visible','off')
   set(handles.rb_fix1,'Visible','off')
   set(handles.rb_opt_nnt_stc1,'Value',0)
   set(handles.rb_fix1,'Value',0)
end
if get(handles.rb_opt_pw_stc0,'Value') == 1
   handles.opt.pw_stc = 0;
   set(handles.rb_opt_pw_stc1,'Value',0)
   set(handles.rb_opt_nnt_stc1,'Value',0)
   set(handles.rb_opt_nnt_stc1,'Visible','on')
   set(handles.rb_fix1,'Value',0)
   set(handles.rb_fix1,'Visible','on')
   set(handles.rb_fix2,'Value',0)
   set(handles.rb_fix2,'Visible','off')
   set(handles.cb_opt_constr_xyz,'Visible','off')
end

if get(handles.rb_opt_nnt_stc1,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[true true true false false false]);
end
if get(handles.rb_fix1,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[false false false true false false]);
end
if get(handles.rb_fix2,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[false false false true true true]);
end
% Hint: get(hObject,'Value') returns toggle state of rb_opt_pw_stc0
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in rb_opt_nnt_stc1.
function rb_opt_nnt_stc1_Callback(hObject, eventdata, handles)
% hObject    handle to rb_opt_nnt_stc1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.rb_opt_nnt_stc1,'Value') == 1
    handles.opt.nnt_stc = 1;
    set(handles.rb_fix1,'Value',0)
end
if get(handles.rb_opt_nnt_stc1,'Value') == 0
    handles.opt.nnt_stc = 0;
end

if get(handles.rb_opt_nnt_stc1,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[true true true false false false]);
end
if get(handles.rb_fix1,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[false false false true false false]);
end
if get(handles.rb_fix2,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[false false false true true true]);
end
% Hint: get(hObject,'Value') returns toggle state of rb_opt_nnt_stc1
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in rb_opt_pw_stc1.
function rb_opt_pw_stc1_Callback(hObject, eventdata, handles)
% hObject    handle to rb_opt_pw_stc1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.rb_opt_pw_stc1,'Value') == 0   
   set(handles.rb_fix2,'Visible','off')
   set(handles.cb_opt_constr_xyz,'Visible','off')
end
if get(handles.rb_opt_pw_stc1,'Value') == 1
   handles.opt.pw_stc = 1;
   handles.opt.nnt_stc = 0;
   set(handles.rb_fix2,'Visible','on')
   set(handles.rb_opt_pw_stc0,'Value',0)
   set(handles.rb_fix2,'Value',1)
   set(handles.rb_opt_nnt_stc1,'Visible','off')
   set(handles.rb_fix1,'Visible','off')
   set(handles.cb_opt_constr_xyz,'Visible','on')
end

if get(handles.rb_opt_nnt_stc1,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[true true true false false false]);
end
if get(handles.rb_fix1,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[false false false true false false]);
end
if get(handles.rb_fix2,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[false false false true true true]);
end
% Hint: get(hObject,'Value') returns toggle state of rb_opt_pw_stc1
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in rb_fix1.
function rb_fix1_Callback(hObject, eventdata, handles)
% hObject    handle to rb_fix1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.rb_fix1,'Value') == 1
   handles.opt.nnt_stc = 0;
   set(handles.rb_opt_nnt_stc1,'Value',0)
   set(handles.rb_fix2,'Value',0)
end

if get(handles.rb_opt_nnt_stc1,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[true true true false false false]);
end
if get(handles.rb_fix1,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[false false false true false false]);
end
if get(handles.rb_fix2,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[false false false true true true]);
end
% Hint: get(hObject,'Value') returns toggle state of rb_fix1
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in rb_fix2.
function rb_fix2_Callback(hObject, eventdata, handles)
% hObject    handle to rb_fix2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.rb_fix2,'Value') == 1
   handles.opt.nnt_stc = 0;
   set(handles.rb_fix1,'Value',0)
end

if get(handles.rb_opt_nnt_stc1,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[true true true false false false]);
end
if get(handles.rb_fix1,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[false false false true false false]);
end
if get(handles.rb_fix2,'Value') == 1 
   set(handles.uitable1,'ColumnEditable',[false false false true true true]);
end
% Hint: get(hObject,'Value') returns toggle state of rb_fix2
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in cb_opt_constr_xyz.
function cb_opt_constr_xyz_Callback(hObject, eventdata, handles)
% hObject    handle to cb_opt_constr_xyz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.opt.constr_xyz = get(handles.cb_opt_constr_xyz,'Value');
% Hint: get(hObject,'Value') returns toggle state of cb_opt_constr_xyz
% Update handles structure
guidata(hObject, handles);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function handles = save_uitables(handles)
D = get(handles.uitable1,'Data');
for istat = 1 : length(handles.opt.stat)
    handles.opt.stat(istat).nnt_inc = D{istat,1};
    handles.opt.stat(istat).nnr_inc = D{istat,2};
    handles.opt.stat(istat).nns_inc = D{istat,3};
    handles.opt.stat(istat).xyz_inc = D{istat,4};
    handles.opt.stat(istat).coef_xyz = D{istat,5};
    handles.opt.stat(istat).int_xyz = D{istat,6};
end


% --- Executes on button press in pb_next4.
function pb_next4_Callback(hObject, eventdata, handles)
% hObject    handle to pb_next4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = save_uitables(handles);
close
uiwait(vie_lsm_gui_eop(handles.opt));


% --- Executes on button press in pb_back4.
function pb_back4_Callback(hObject, eventdata, handles)
% hObject    handle to pb_back4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = save_uitables(handles);
close
uiwait(vie_lsm_gui_tropo(handles.opt));
