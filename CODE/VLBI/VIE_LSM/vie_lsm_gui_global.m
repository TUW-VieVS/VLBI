% ************************************************************************
%   Description:
%   function to create a GUI for SESSION SPECIFIC parameterization for 
%   global estimates (antenna coor., source coor., station velocities, 
%   Love and Shida numbers). 
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
%   05 Oct 2010 by Hana Spicakova: N and b for sinex output
%   12 Nov 2010 by Hana Spicakova: reorganisation of the GUI, option for
%   set up only N, b added (without estimation of the parameters)
%   19 May 2011 by Hana Spicakova: Sinex output added
%   24 Oct 2013 by Hana Krasna: Axis offset added 
%   28 Feb 2017 by Andreas Hellerschmied: opt is stored in handles.opt

% ************************************************************************
function varargout = vie_lsm_gui_global(varargin)
% VIE_LSM_GUI_GLOBAL M-file for vie_lsm_gui_global.fig
%      VIE_LSM_GUI_GLOBAL, by itself, creates a new VIE_LSM_GUI_GLOBAL or raises the existing
%      singleton*.
%
%      H = VIE_LSM_GUI_GLOBAL returns the handle to a new VIE_LSM_GUI_GLOBAL or the handle to
%      the existing singleton*.
%
%      VIE_LSM_GUI_GLOBAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIE_LSM_GUI_GLOBAL.M with the given input arguments.
%
%      VIE_LSM_GUI_GLOBAL('Property','Value',...) creates a new VIE_LSM_GUI_GLOBAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vie_lsm_gui_global_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vie_lsm_gui_global_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vie_lsm_gui_global

% Last Modified by GUIDE v2.5 06-Dec-2018 09:36:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vie_lsm_gui_global_OpeningFcn, ...
                   'gui_OutputFcn',  @vie_lsm_gui_global_OutputFcn, ...
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


% --- Executes just before vie_lsm_gui_global is made visible.
function vie_lsm_gui_global_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vie_lsm_gui_global (see VARARGIN)

% Choose default command line output for vie_lsm_gui_global
handles.output = hObject;

handles.opt = varargin{1};

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vie_lsm_gui_global wait for user response (see UIRESUME)
% uiwait(handles.figure1);


esteop=1; % check if at least one EOP is estimated
if handles.opt.xpol.model==0 && handles.opt.ypol.model==0 && handles.opt.dut1.model==0 && handles.opt.nutdx.model==0 && handles.opt.nutdy.model==0
    esteop=0;
end

% Kamil (09.08.2012)
if handles.opt.outsnx.clk == 0 % Clocks can only be reduced
    set(handles.rb_snx_clk_red,'Value',1); 
end
if handles.opt.outsnx.zwd == 0 % ZWD can be reduced
    set(handles.rb_snx_zwd_red,'Value',1);
    set(handles.rb_snx_zwd_inc,'Value',0); 
elseif handles.opt.outsnx.zwd == 1 % ZWD can be included
     set(handles.rb_snx_zwd_inc,'Value',1);
     set(handles.rb_snx_zwd_red,'Value',0); 
end
if handles.opt.outsnx.tgr == 0 % Trop. gradients can be reduced
    set(handles.rb_snx_tgr_red,'Value',1);
    set(handles.rb_snx_tgr_inc,'Value',0); 
elseif handles.opt.outsnx.tgr == 1 % Trop. gradients can be included
     set(handles.rb_snx_tgr_inc,'Value',1);
     set(handles.rb_snx_tgr_red,'Value',0); 
end
if handles.opt.outsnx.eop == 0  % EOP can be reduced
    set(handles.rb_snx_eop_red,'Value',1);
    set(handles.rb_snx_eop_inc,'Value',0);
elseif handles.opt.outsnx.eop == 1 % EOP can be included
     set(handles.rb_snx_eop_inc,'Value',1);
     set(handles.rb_snx_eop_red,'Value',0); 
end
if handles.opt.outsnx.xyz == 1 % TRF antenna coor. can only be included
    set(handles.rb_snx_xyz_inc,'Value',1); 
end
if handles.opt.outsnx.sou == 1 % CRF source coor. can only be included
    set(handles.rb_snx_sou_inc,'Value',1); 
else
    set(handles.rb_snx_sou_inc,'Value',0);
end
    
set(handles.cb_opt_est_singleses,'Value',handles.opt.est_singleses) %hana
set(handles.cb_opt_global_solve,'Value',handles.opt.global_solve)
set(handles.cb_opt_source,'Value',handles.opt.est_source)
set(handles.cb_opt_vel,'Value',handles.opt.est_vel) %hana
set(handles.ed_refvel,'Value',handles.opt.refvel) %hana
set(handles.cb_opt_ao,'Value',handles.opt.est_AO) %hana
set(handles.cb_ascii_snx,'Value',handles.opt.ascii_snx) %hana

if get(handles.cb_opt_global_solve,'Value') == 0
    set(handles.cb_opt_source,'Visible','off')
    set(handles.cb_opt_vel,'Visible','off')
    set(handles.cb_opt_ao,'Visible','off')
    set(handles.ed_refvel,'Visible','off') %hana
    set(handles.text_refvel,'Visible','off') %hana
    set(handles.warn,'Visible','off') %hana
end
if get(handles.cb_opt_global_solve,'Value') == 1
    set(handles.cb_opt_source,'Visible','on')
    set(handles.cb_opt_vel,'Visible','on')
    set(handles.cb_opt_ao,'Visible','on')
    if get(handles.cb_opt_source,'Value') == 1 & get(handles.cb_opt_global_solve,'Value') == 1
       set(handles.warn,'Visible','on')
    end
    if get(handles.cb_opt_source,'Value') == 0
       set(handles.warn,'Visible','off')
    end
end

%hana
if get(handles.cb_ascii_snx,'Value') == 1
    set(handles.cb_opt_source,'Visible','on')
    set(handles.warn,'Visible','on')
    
    set(handles.rb_snx_clk_red,'Visible','on')
    if handles.opt.pw_zwd==1
        set(handles.rb_snx_zwd_inc,'Visible','on')
        set(handles.rb_snx_zwd_red,'Visible','on')
    end
    if handles.opt.pw_ngr==1 && handles.opt.pw_egr==1
        set(handles.rb_snx_tgr_inc,'Visible','on')
        set(handles.rb_snx_tgr_red,'Visible','on')
    end
    if get(handles.cb_opt_source,'Value') ==1   
        set(handles.rb_snx_sou_inc,'Visible','on')
    end
    if handles.opt.stc==1
        set(handles.rb_snx_xyz_inc,'Visible','on')
    end
    if esteop==1
        set(handles.rb_snx_eop_inc,'Visible','on')
        set(handles.rb_snx_eop_red,'Visible','on')
    end
end

if get(handles.cb_ascii_snx,'Value') == 0
    set(handles.cb_opt_source,'Visible','off')
    set(handles.warn,'Visible','off')
    
    set(handles.rb_snx_clk_red,'Visible','off')
    set(handles.rb_snx_zwd_inc,'Visible','off')
    set(handles.rb_snx_zwd_red,'Visible','off')
    set(handles.rb_snx_tgr_inc,'Visible','off')
    set(handles.rb_snx_tgr_red,'Visible','off')
    set(handles.rb_snx_sou_inc,'Visible','off')
    set(handles.rb_snx_xyz_inc,'Visible','off')
    set(handles.rb_snx_eop_inc,'Visible','off')
    set(handles.rb_snx_eop_red,'Visible','off')
end
% Update handles structure
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = vie_lsm_gui_global_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in cb_opt_global_solve.
function cb_opt_global_solve_Callback(hObject, eventdata, handles)
% hObject    handle to cb_opt_global_solve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.opt.global_solve = get(handles.cb_opt_global_solve,'Value');
handles.opt.ascii_snx = get(handles.cb_ascii_snx,'Value');
if get(handles.cb_opt_global_solve,'Value') == 0  
    if get(handles.cb_ascii_snx,'Value') == 0 %hana
        set(handles.cb_opt_source,'Visible','off')
        set(handles.warn,'Visible','off') %hana
    end
    set(handles.cb_opt_vel,'Visible','off')
    set(handles.cb_opt_ao,'Visible','off')
    set(handles.ed_refvel,'Visible','off') %hana
    set(handles.text_refvel,'Visible','off') %hana
end
if get(handles.cb_opt_global_solve,'Value') == 1
    set(handles.cb_opt_source,'Visible','on')
    set(handles.warn,'Visible','on') %hana
    set(handles.cb_opt_vel,'Visible','on')
    set(handles.cb_opt_ao,'Visible','on')
    if get(handles.cb_opt_vel,'Value') == 0 %hana
        set(handles.ed_refvel,'Visible','off') %hana
        set(handles.text_refvel,'Visible','off') %hana
    else
        set(handles.ed_refvel,'Visible','on') %hana
        set(handles.text_refvel,'Visible','on') %hana
    end
end

% Hint: get(hObject,'Value') returns toggle state of cb_opt_global_solve
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in cb_ascii_snx.
function cb_ascii_snx_Callback(hObject, eventdata, handles)

esteop=1; % check if at least one EOP is estimated
if handles.opt.xpol.model==0 && handles.opt.ypol.model==0 && handles.opt.dut1.model==0 && handles.opt.nutdx.model==0 && handles.opt.nutdy.model==0
    esteop=0;
end

handles.opt.ascii_snx = get(handles.cb_ascii_snx,'Value');
handles.opt.global_solve = get(handles.cb_opt_global_solve,'Value');
if get(handles.cb_ascii_snx,'Value') == 0 & get(handles.cb_opt_global_solve,'Value') == 0
    set(handles.cb_opt_source,'Visible','off')
    set(handles.warn,'Visible','off')
end
if get(handles.cb_ascii_snx,'Value') == 0
    set(handles.rb_snx_clk_red,'Visible','off')
    set(handles.rb_snx_zwd_inc,'Visible','off')
    set(handles.rb_snx_zwd_red,'Visible','off')
    set(handles.rb_snx_tgr_inc,'Visible','off')
    set(handles.rb_snx_tgr_red,'Visible','off')
    set(handles.rb_snx_sou_inc,'Visible','off')
    set(handles.rb_snx_xyz_inc,'Visible','off')
    set(handles.rb_snx_eop_inc,'Visible','off')
    set(handles.rb_snx_eop_red,'Visible','off')
end

if get(handles.cb_ascii_snx,'Value') == 1
    set(handles.cb_opt_source,'Visible','on')
    set(handles.warn,'Visible','on')
    
    set(handles.rb_snx_clk_red,'Visible','on')
    if handles.opt.pw_zwd==1
        set(handles.rb_snx_zwd_inc,'Visible','on')
        set(handles.rb_snx_zwd_red,'Visible','on')
    end
    if handles.opt.pw_ngr==1 && handles.opt.pw_egr==1
        set(handles.rb_snx_tgr_inc,'Visible','on')
        set(handles.rb_snx_tgr_red,'Visible','on')
    end
    if get(handles.cb_opt_source,'Value') ==1   
        set(handles.rb_snx_sou_inc,'Visible','on')
    end
    if handles.opt.stc==1
        set(handles.rb_snx_xyz_inc,'Visible','on')
    end
    if esteop==1
        set(handles.rb_snx_eop_inc,'Visible','on')
        set(handles.rb_snx_eop_red,'Visible','on')
    end
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pb_finish.
function pb_finish_Callback(hObject, eventdata, handles)
  
if  get(handles.rb_snx_zwd_inc,'Value')==1 && handles.opt.pw_zwd==1
    handles.opt.outsnx.zwd = 1;
end
if  get(handles.rb_snx_tgr_inc,'Value')== 1 && handles.opt.pw_ngr==1 && handles.opt.pw_egr==1
    handles.opt.outsnx.tgr = 1;
end
if  get(handles.rb_snx_sou_inc,'Value')==1 && handles.opt.est_source==1
    handles.opt.outsnx.sou = 1;
end

esteop=1; % check if at least one EOP is estimated
if handles.opt.xpol.model==0 && handles.opt.ypol.model==0 && handles.opt.dut1.model==0 && handles.opt.nutdx.model==0 && handles.opt.nutdy.model==0
    esteop=0;
end
if  get(handles.rb_snx_eop_inc,'Value')==1 && esteop==1
    handles.opt.outsnx.eop = 1;
end

if handles.opt.stc==0 % fixed station coordinates
    handles.opt.outsnx.xyz = 0;
end


% Update handles structure
guidata(hObject, handles);

% Save handles.opt to /WORK/opt_tmp.mat
opt = handles.opt;
save('opt_tmp.mat', 'opt');

close


% --- Executes on button press in cb_opt_source.
function cb_opt_source_Callback(hObject, eventdata, handles)

handles.opt.est_source = get(handles.cb_opt_source,'Value');
if get(handles.cb_opt_source,'Value') == 0
    set(handles.rb_snx_sou_inc,'Visible','off')
elseif get(handles.cb_opt_source,'Value') == 1 && get(handles.cb_ascii_snx,'Value') == 1
    set(handles.rb_snx_sou_inc,'Visible','on')
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pb_back7.
function pb_back7_Callback(hObject, eventdata, handles)
close
uiwait(vie_lsm_gui_sourcoor(handles.opt));




% --- Executes on button press in cb_opt_vel.
function cb_opt_vel_Callback(hObject, eventdata, handles)

handles.opt.est_vel = get(handles.cb_opt_vel,'Value');
if get(handles.cb_opt_vel,'Value') == 1
    set(handles.ed_refvel,'Visible','on') %hana
    set(handles.text_refvel,'Visible','on') %hana
end
if get(handles.cb_opt_vel,'Value') == 0
    set(handles.ed_refvel,'Visible','off') %hana
    set(handles.text_refvel,'Visible','off') %hana
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in cb_opt_est_singleses.
function cb_opt_est_singleses_Callback(hObject, eventdata, handles)
handles.opt.est_singleses = get(handles.cb_opt_est_singleses,'Value');
% Update handles structure
guidata(hObject, handles);




function ed_refvel_Callback(hObject, eventdata, handles)
handles.opt.refvel = str2double(get(handles.ed_refvel,'String'));
% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function ed_refvel_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in cb_opt_ao.
function cb_opt_ao_Callback(hObject, eventdata, handles)
handles.opt.est_AO = get(handles.cb_opt_ao,'Value');
% Update handles structure
guidata(hObject, handles);
