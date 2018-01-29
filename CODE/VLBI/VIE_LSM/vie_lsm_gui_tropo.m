% ************************************************************************
%   Description:
%   function to create a GUI for SESSION SPECIFIC parameterization for
%   estimates of troposphere (zenith wet delays, troposphere gradients).
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
function varargout = vie_lsm_gui_tropo(varargin)
% VIE_LSM_GUI_TROPO M-file for vie_lsm_gui_tropo.fig
%      VIE_LSM_GUI_TROPO, by itself, creates a new VIE_LSM_GUI_TROPO or raises the existing
%      singleton*.
%
%      H = VIE_LSM_GUI_TROPO returns the handle to a new VIE_LSM_GUI_TROPO or the handle to
%      the existing singleton*.
%
%      VIE_LSM_GUI_TROPO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIE_LSM_GUI_TROPO.M with the given input arguments.
%
%      VIE_LSM_GUI_TROPO('Property','Value',...) creates a new VIE_LSM_GUI_TROPO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vie_lsm_gui_tropo_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vie_lsm_gui_tropo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vie_lsm_gui_tropo

% Last Modified by GUIDE v2.5 18-Feb-2010 14:07:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vie_lsm_gui_tropo_OpeningFcn, ...
                   'gui_OutputFcn',  @vie_lsm_gui_tropo_OutputFcn, ...
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


% --- Executes just before vie_lsm_gui_tropo is made visible.
function vie_lsm_gui_tropo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vie_lsm_gui_tropo (see VARARGIN)

% Choose default command line output for vie_lsm_gui_tropo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vie_lsm_gui_tropo wait for user response (see UIRESUME)
% uiwait(handles.figure1);

handles.opt = varargin{1};

set(handles.cb_opt_constr_zwd,'Value',handles.opt.constr_zwd)
set(handles.cb_opt_constr_rel_ngr,'Value',handles.opt.constr_rel_ngr)
set(handles.cb_opt_constr_rel_egr,'Value',handles.opt.constr_rel_egr)
set(handles.cb_opt_constr_abs_ngr,'Value',handles.opt.constr_abs_ngr)
set(handles.cb_opt_constr_abs_egr,'Value',handles.opt.constr_abs_egr)

D = get(handles.uitable1,'Data');
for istat = 1 : length(handles.opt.stat)
   
    D{istat,1}=handles.opt.stat(istat).coef_zwd;
    D{istat,2}=handles.opt.stat(istat).coef_rel_ngr;
    D{istat,3}=handles.opt.stat(istat).coef_rel_egr;
    
    D{istat,4}=handles.opt.stat(istat).coef_abs_ngr;
    D{istat,5}=handles.opt.stat(istat).coef_abs_egr;
    
    D{istat,6}=handles.opt.stat(istat).int_zwd;
    D{istat,7}=handles.opt.stat(istat).int_ngr;
    D{istat,8}=handles.opt.stat(istat).int_egr;
   
    D{istat,9}=(handles.opt.stat(istat).zwd_inc == 1);
    D{istat,10}=(handles.opt.stat(istat).ngr_inc == 1);
    D{istat,11}=(handles.opt.stat(istat).egr_inc == 1);
   
    stat_name(istat,1:8) = handles.opt.stat(istat).name(1:8);
end
set(handles.uitable1,'Data',D);
set(handles.uitable1,'RowName',stat_name);
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = vie_lsm_gui_tropo_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in cb_opt_constr_zwd.
function cb_opt_constr_zwd_Callback(hObject, eventdata, handles)
% hObject    handle to cb_opt_constr_zwd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.opt.constr_zwd = get(handles.cb_opt_constr_zwd,'Value');
% Hint: get(hObject,'Value') returns toggle state of cb_opt_constr_zwd
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in cb_opt_constr_rel_ngr.
function cb_opt_constr_rel_ngr_Callback(hObject, eventdata, handles)
% hObject    handle to cb_opt_constr_rel_ngr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.opt.constr_rel_ngr = get(handles.cb_opt_constr_rel_ngr,'Value');
% Hint: get(hObject,'Value') returns toggle state of cb_opt_constr_rel_ngr
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in cb_opt_constr_rel_egr.
function cb_opt_constr_rel_egr_Callback(hObject, eventdata, handles)
% hObject    handle to cb_opt_constr_rel_egr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.opt.constr_rel_egr = get(handles.cb_opt_constr_rel_egr,'Value');
% Hint: get(hObject,'Value') returns toggle state of cb_opt_constr_rel_egr
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in cb_opt_constr_abs_ngr.
function cb_opt_constr_abs_ngr_Callback(hObject, eventdata, handles)
% hObject    handle to cb_opt_constr_abs_ngr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.opt.constr_abs_ngr = get(handles.cb_opt_constr_abs_ngr,'Value');
% Hint: get(hObject,'Value') returns toggle state of cb_opt_constr_abs_ngr
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in cb_opt_constr_abs_egr.
function cb_opt_constr_abs_egr_Callback(hObject, eventdata, handles)
% hObject    handle to cb_opt_constr_abs_egr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.opt.constr_abs_egr = get(handles.cb_opt_constr_abs_egr,'Value');
% Hint: get(hObject,'Value') returns toggle state of cb_opt_constr_abs_egr
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
    handles.opt.stat(istat).coef_zwd = D{istat,1};
    handles.opt.stat(istat).coef_rel_ngr = D{istat,2};
    handles.opt.stat(istat).coef_rel_egr = D{istat,3};
    
    handles.opt.stat(istat).coef_abs_ngr = D{istat,4};
    handles.opt.stat(istat).coef_abs_egr = D{istat,5};
    
    handles.opt.stat(istat).int_zwd = D{istat,6};
    handles.opt.stat(istat).int_ngr = D{istat,7};
    handles.opt.stat(istat).int_egr = D{istat,8};
    
    handles.opt.stat(istat).zwd_inc = D{istat,9};
    handles.opt.stat(istat).ngr_inc = D{istat,10};
    handles.opt.stat(istat).egr_inc = D{istat,11};
end


% --- Executes on button press in pb_next3.
function pb_next3_Callback(hObject, eventdata, handles)
% hObject    handle to pb_next3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = save_uitables(handles);
close
uiwait(vie_lsm_gui_statcoor(handles.opt));

% --- Executes on button press in pb_back3.
function pb_back3_Callback(hObject, eventdata, handles)
% hObject    handle to pb_back3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = save_uitables(handles);
close
uiwait(vie_lsm_gui_clock(handles.opt));

