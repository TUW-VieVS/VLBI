% ************************************************************************
%   Description:
%   function to create a GUI for SESSION SPECIFIC EOP parameterization. 
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
function varargout = vie_lsm_gui_eop(varargin)
% VIE_LSM_GUI_EOP M-file for vie_lsm_gui_eop.fig
%      VIE_LSM_GUI_EOP, by itself, creates a new VIE_LSM_GUI_EOP or raises the existing
%      singleton*.
%
%      H = VIE_LSM_GUI_EOP returns the handle to a new VIE_LSM_GUI_EOP or the handle to
%      the existing singleton*.
%
%      VIE_LSM_GUI_EOP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIE_LSM_GUI_EOP.M with the given input arguments.
%
%      VIE_LSM_GUI_EOP('Property','Value',...) creates a new VIE_LSM_GUI_EOP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vie_lsm_gui_eop_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vie_lsm_gui_eop_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vie_lsm_gui_eop

% Last Modified by GUIDE v2.5 24-Sep-2009 15:11:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vie_lsm_gui_eop_OpeningFcn, ...
                   'gui_OutputFcn',  @vie_lsm_gui_eop_OutputFcn, ...
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


% --- Executes just before vie_lsm_gui_eop is made visible.
function vie_lsm_gui_eop_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vie_lsm_gui_eop (see VARARGIN)

% Choose default command line output for vie_lsm_gui_eop
handles.output = hObject;

handles.opt = varargin{1};

% UIWAIT makes vie_lsm_gui_eop wait for user response (see UIRESUME)
% uiwait(handles.figure1);

ED = get(handles.uitable1,'Data');
ED{1,1} = (handles.opt.xpol.model == 1); ED{1,2} = handles.opt.xpol.int; ED{1,3} = (handles.opt.xpol.constrain == 1); ED{1,4} = handles.opt.xpol.coef;
ED{2,1} = (handles.opt.ypol.model == 1); ED{2,2} = handles.opt.ypol.int; ED{2,3} = (handles.opt.ypol.constrain == 1); ED{2,4} = handles.opt.ypol.coef;
ED{3,1} = (handles.opt.dut1.model == 1); ED{3,2} = handles.opt.dut1.int; ED{3,3} = (handles.opt.dut1.constrain == 1); ED{3,4} = handles.opt.dut1.coef;
ED{4,1} = (handles.opt.nutdx.model == 1); ED{4,2} = handles.opt.nutdx.int; ED{4,3} = (handles.opt.nutdx.constrain == 1); ED{4,4} = handles.opt.nutdx.coef;
ED{5,1} = (handles.opt.nutdy.model == 1); ED{5,2} = handles.opt.nutdy.int; ED{5,3} = (handles.opt.nutdy.constrain == 1); ED{5,4} = handles.opt.nutdy.coef;
set(handles.uitable1,'Data',ED);
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = vie_lsm_gui_eop_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% Update handles structure
guidata(hObject, handles);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function handles = save_uitables(handles)
ED = get(handles.uitable1,'Data');
handles.opt.xpol.model = ED{1,1}; handles.opt.xpol.int = ED{1,2}; handles.opt.xpol.constrain = ED{1,3}; handles.opt.xpol.coef = ED{1,4};
handles.opt.ypol.model = ED{2,1}; handles.opt.ypol.int = ED{2,2}; handles.opt.ypol.constrain = ED{2,3}; handles.opt.ypol.coef = ED{2,4};
handles.opt.dut1.model = ED{3,1}; handles.opt.dut1.int = ED{3,2}; handles.opt.dut1.constrain = ED{3,3}; handles.opt.dut1.coef = ED{3,4};
handles.opt.nutdx.model = ED{4,1}; handles.opt.nutdx.int = ED{4,2}; handles.opt.nutdx.constrain = ED{4,3}; handles.opt.nutdx.coef = ED{4,4};
handles.opt.nutdy.model = ED{5,1}; handles.opt.nutdy.int = ED{5,2}; handles.opt.nutdy.constrain = ED{5,3}; handles.opt.nutdy.coef = ED{5,4};


% --- Executes during object deletion, before destroying properties.
function uitable1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pb_next5.
function pb_next5_Callback(hObject, eventdata, handles)
% hObject    handle to pb_next5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = save_uitables(handles);
close
uiwait(vie_lsm_gui_sourcoor(handles.opt));


% --- Executes on button press in pb_back5.
function pb_back5_Callback(hObject, eventdata, handles)
% hObject    handle to pb_back5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = save_uitables(handles);
close
uiwait(vie_lsm_gui_statcoor(handles.opt));
