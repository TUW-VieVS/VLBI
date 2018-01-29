% ************************************************************************
%   Description:
%   function to create a GUI for SESSION SPECIFIC parameterization for 
%   clock options of first solution, determination if main solution will be 
%   carried out or not, outlier test options.
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
function varargout = vie_lsm_gui_first(varargin)
% GUIDETEMPLATE0 M-file for guidetemplate0.fig
%      GUIDETEMPLATE0, by itself, creates a new GUIDETEMPLATE0 or raises the existing
%      singleton*.
%
%      H = GUIDETEMPLATE0 returns the handle to a new GUIDETEMPLATE0 or the handle to
%      the existing singleton*.
%
%      GUIDETEMPLATE0('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIDETEMPLATE0.M with the given input arguments.
%
%      GUIDETEMPLATE0('Property','Value',...) creates a new GUIDETEMPLATE0 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guidetemplate0_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guidetemplate0_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2006 The MathWorks, Inc.

% Edit the above text to modify the response to help guidetemplate0

% Last Modified by GUIDE v2.5 27-Sep-2009 21:59:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guidetemplate0_OpeningFcn, ...
                   'gui_OutputFcn',  @guidetemplate0_OutputFcn, ...
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


% --- Executes just before guidetemplate0 is made visible.
function guidetemplate0_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guidetemplate0 (see VARARGIN)

% Choose default command line output for guidetemplate0
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes guidetemplate0 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

handles.opt = varargin{1};
% % global opt

set(handles.cb_opt_first,'Value',handles.opt.first)
set(handles.cb_opt_treat_breaks,'Value',handles.opt.treat_breaks)
if handles.opt.first == 0
    set(handles.rb_opt_firstclock0,'Visible','off');
    set(handles.rb_opt_firstclock1,'Visible','off');
    set(handles.rb_opt_firstclock2,'Visible','off');
    set(handles.cb_opt_treat_breaks,'Visible','off');
end
if handles.opt.first == 1
    if handles.opt.firstclock == 0
        set(handles.rb_opt_firstclock0,'Value',1);
    end
    if handles.opt.firstclock == 1
        set(handles.rb_opt_firstclock1,'Value',1);
    end
    if handles.opt.firstclock == 2
        set(handles.rb_opt_firstclock2,'Value',1);
    end
end

set(handles.cb_opt_second,'Value',handles.opt.second);
set(handles.cb_opt_basic_outlier,'Value',handles.opt.basic_outlier)
set(handles.cb_opt_simple_outlier,'Value',handles.opt.simple_outlier)
set(handles.et_opt_par_outlier,'String',handles.opt.par_outlier)

for istat = 1 : length(handles.opt.stat)
    stat_name(istat,1:8) = handles.opt.stat(istat).name(1:8);
    if (handles.opt.stat(istat).ref == 1)
        refclk= istat;
    end
end
set(handles.popupmenu1,'String',stat_name);
set(handles.popupmenu1,'Value',refclk);

if handles.opt.second == 0
    set(handles.cb_opt_simple_outlier,'Visible','off')
    set(handles.cb_opt_basic_outlier,'Visible','off')
    set(handles.et_opt_par_outlier,'Visible','off')
    set(handles.coef_text,'Visible','off')
end

% Clock breaks occured stations are written to GUI
if isempty(find([handles.opt.stat.clkbreak],1))
    set(handles.tx_opt_clkbrk,'String','No clock breaks information!')
elseif ~isempty(find([handles.opt.stat.clkbreak],1))
    count = 0;
    for i = 1 : length([handles.opt.stat])
        if ~isempty(find([handles.opt.stat(i).clkbreak],1))
            count = count + 1;
            clk_break_stat(count,1:8) = handles.opt.stat(i).name(1:8);
        end
    end
    set(handles.tx_opt_clkbrk,'String',clk_break_stat)
end

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = guidetemplate0_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in cb_opt_first.
function cb_opt_first_Callback(hObject, eventdata, handles)
% hObject    handle to cb_opt_first (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global opt
handles.opt.first = get(handles.cb_opt_first,'Value');
if handles.opt.first == 1
    set(handles.rb_opt_firstclock0,'Visible','on');
    set(handles.rb_opt_firstclock1,'Visible','on');
    set(handles.rb_opt_firstclock2,'Visible','on');
    set(handles.cb_opt_treat_breaks,'Visible','on');
    if handles.opt.firstclock == 0
        set(handles.rb_opt_firstclock0,'Value',1);
    end
    if handles.opt.firstclock == 1
        set(handles.rb_opt_firstclock1,'Value',1);
    end
    if handles.opt.firstclock == 2
        set(handles.rb_opt_firstclock2,'Value',1);
    end
else
    set(handles.rb_opt_firstclock0,'Visible','off');
    set(handles.rb_opt_firstclock1,'Visible','off');
    set(handles.rb_opt_firstclock2,'Visible','off');
    set(handles.cb_opt_treat_breaks,'Visible','off');
end
% Hint: get(hObject,'Value') returns toggle state of cb_opt_first
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in rb_opt_firstclock0.
function rb_opt_firstclock0_Callback(hObject, eventdata, handles)
% hObject    handle to rb_opt_firstclock0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global opt
if get(handles.rb_opt_firstclock0,'Value') == 1
    handles.opt.firstclock = 0; 
    set(handles.rb_opt_firstclock1,'Value',0);
    set(handles.rb_opt_firstclock2,'Value',0);
end
% Hint: get(hObject,'Value') returns toggle state of rb_opt_firstclock0
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in rb_opt_firstclock1.
function rb_opt_firstclock1_Callback(hObject, eventdata, handles)
% hObject    handle to rb_opt_firstclock1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global opt
if get(handles.rb_opt_firstclock1,'Value') == 1
    handles.opt.firstclock = 1; 
    set(handles.rb_opt_firstclock0,'Value',0);
    set(handles.rb_opt_firstclock2,'Value',0);
end
% Hint: get(hObject,'Value') returns toggle state of rb_opt_firstclock1
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in rb_opt_firstclock2.
function rb_opt_firstclock2_Callback(hObject, eventdata, handles)
% hObject    handle to rb_opt_firstclock2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global opt
if get(handles.rb_opt_firstclock2,'Value') == 1
    handles.opt.firstclock = 2; 
    set(handles.rb_opt_firstclock0,'Value',0);
    set(handles.rb_opt_firstclock1,'Value',0);
end
% Hint: get(hObject,'Value') returns toggle state of rb_opt_firstclock2
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in cb_opt_treat_breaks.
function cb_opt_treat_breaks_Callback(hObject, eventdata, handles)
% hObject    handle to cb_opt_treat_breaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global opt
handles.opt.treat_breaks = get(handles.cb_opt_treat_breaks,'Value');
% Hint: get(hObject,'Value') returns toggle state of cb_opt_treat_breaks
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in cb_opt_second.
function cb_opt_second_Callback(hObject, eventdata, handles)
% hObject    handle to cb_opt_second (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global opt
handles.opt.second = get(handles.cb_opt_second,'Value');
if handles.opt.second == 1
    set(handles.cb_opt_simple_outlier,'Visible','on')
    set(handles.cb_opt_basic_outlier,'Visible','on')
    set(handles.et_opt_par_outlier,'Visible','on')
    set(handles.coef_text,'Visible','on')
end
if handles.opt.second == 0
    set(handles.cb_opt_simple_outlier,'Visible','off')
    set(handles.cb_opt_basic_outlier,'Visible','off')
    set(handles.et_opt_par_outlier,'Visible','off')
    set(handles.coef_text,'Visible','off')
end
% Hint: get(hObject,'Value') returns toggle state of cb_opt_second
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in cb_opt_basic_outlier.
function cb_opt_basic_outlier_Callback(hObject, eventdata, handles)
% hObject    handle to cb_opt_basic_outlier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global opt
handles.opt.basic_outlier = get(handles.cb_opt_basic_outlier,'Value');
if get(handles.cb_opt_basic_outlier,'Value') == 1
    set(handles.cb_opt_simple_outlier,'Value',0)
    handles.opt.simple_outlier = 0;
end
% Hint: get(hObject,'Value') returns toggle state of cb_opt_basic_outlier
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in cb_opt_simple_outlier.
function cb_opt_simple_outlier_Callback(hObject, eventdata, handles)
% hObject    handle to cb_opt_simple_outlier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global opt
handles.opt.simple_outlier = get(handles.cb_opt_simple_outlier,'Value');
if get(handles.cb_opt_simple_outlier,'Value') == 1
    set(handles.cb_opt_basic_outlier,'Value',0)
    handles.opt.basic_outlier = 0;
end
% Hint: get(hObject,'Value') returns toggle state of cb_opt_simple_outlier
% Update handles structure
guidata(hObject, handles);


function et_opt_par_outlier_Callback(hObject, eventdata, handles)
% hObject    handle to et_opt_par_outlier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global opt
handles.opt.par_outlier = get(handles.et_opt_par_outlier,'String');
handles.opt.par_outlier = str2num([handles.opt.par_outlier]);
% Hints: get(hObject,'String') returns contents of et_opt_par_outlier as text
%        str2double(get(hObject,'String')) returns contents of et_opt_par_outlier as a double
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function et_opt_par_outlier_CreateFcn(hObject, eventdata, handles)
% hObject    handle to et_opt_par_outlier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles);



% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global opt
refclk=get(handles.popupmenu1,'Value');
handles.opt.ref_first_clk = refclk;
for istat = 1 : length(handles.opt.stat)
    if istat == refclk
        handles.opt.stat(istat).ref = 1;
    else
        handles.opt.stat(istat).ref = 0;
    end
end
% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles);

% --- Executes on button press in pb_next1.
function pb_next1_Callback(hObject, eventdata, handles)
% hObject    handle to pb_next1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global opt
if handles.opt.second == 1
    close
    uiwait(vie_lsm_gui_clock(handles.opt));
end
if handles.opt.second == 0
    % Save handles.opt to /WORK/opt_tmp.mat
    opt = handles.opt;
    save('opt_tmp.mat', 'opt');
	close
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function tx_opt_clkbrk_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tx_opt_clkbrk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
guidata(hObject, handles);
