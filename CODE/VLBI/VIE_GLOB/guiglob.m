% ************************************************************************
%   Description:
%   GUIGLOB M-file for guiglob.fig
%   1st Graphical User Inteface for vie_glob
%
%   Input:										
%      
%   Output:                stored in VieVS/OUT/GLOB/
%      paths               paths to the data
%      parGS               name of the parameter and if it will be
%                          estimated (1), fixed (0) or reduced (2)
%      maxRMS              max RMS of the sessions, which will be used in
%                          the global adjustment
%
%   External calls: 	
%      globind.m; guiglob_rf.m                					    											
%       
%   Coded for VieVS: 
%   01 Aug 2010 by Hana Spicakova
%
%   Revision: 
%   25 Jul 2011 by Hana Spicakova:   add option to store session-wise N
%                          matrices and b vectors for a backward solution
%                          of reduced parameters.

%%

function varargout = guiglob(varargin)


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guiglob_OpeningFcn, ...
                   'gui_OutputFcn',  @guiglob_OutputFcn, ...
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


% --- Executes just before guiglob is made visible.
function guiglob_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for guiglob
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.path_L2,'String','../../DATA/LEVEL2/'); % apriori path to LEVEL2 data
set(handles.dir_L2,'String','TEST_LEVEL2'); % apriori directory where are the LEVEL2 data
set(handles.dir_out,'String','TEST_OUT'); % apriori subdirectories for output
set(handles.maxRMS,'String','2'); % a posteriori variance of unit weight



% UIWAIT makes guiglob wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = guiglob_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;



function dir_L2_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function dir_L2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dir_out_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function dir_out_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxRMS_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function maxRMS_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function path_L2_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function path_L2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in pb_ok.
function pb_ok_Callback(hObject, eventdata, handles)


% hObject    handle to pb_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

paths.path_in = get(handles.path_L2,'String');
paths.L2 = get(handles.dir_L2,'String');
paths.out = get(handles.dir_out,'String');


maxRMS=str2double(get(handles.maxRMS,'String'));

parGS(1).name = 'pwck';
parGS(2).name = 'rqck';
parGS(3).name = 'zwd';
parGS(4).name = 'ngr';
parGS(5).name = 'egr';
parGS(6).name = 'ant_x';
parGS(7).name = 'ant_y';
parGS(8).name = 'ant_z';
parGS(9).name = 'vx';
parGS(10).name = 'vy';
parGS(11).name = 'vz';
parGS(12).name = 'sra';
parGS(13).name = 'sde';
parGS(14).name = 'xpol';
parGS(15).name = 'ypol';
parGS(16).name = 'dut1';
parGS(17).name = 'dX';
parGS(18).name = 'dY';
parGS(19).name = 'love';
parGS(20).name = 'shida';

[g] = globind(parGS);

est(g.g_clk(1))=0;
est(g.g_clk(2))=0;
est(g.g_zwd)=0;
est(g.g_tgr(1))=0;
est(g.g_tgr(2))=0;
est(g.g_coord(1))=get(handles.rb_acoor_1,'Value');
est(g.g_coord(2))=get(handles.rb_acoor_1,'Value');
est(g.g_coord(3))=get(handles.rb_acoor_1,'Value');
est(g.g_vel(1))=get(handles.rb_avel_1,'Value');
est(g.g_vel(2))=get(handles.rb_avel_1,'Value');
est(g.g_vel(3))=get(handles.rb_avel_1,'Value');
est(g.g_srade(1))=get(handles.rb_srade_1,'Value');
est(g.g_srade(2))=get(handles.rb_srade_1,'Value');
est(g.g_eop(1))=get(handles.rb_xpol_1,'Value');
est(g.g_eop(2))=get(handles.rb_ypol_1,'Value');
est(g.g_eop(3))=get(handles.rb_dut1_1,'Value');
est(g.g_eop(4))=get(handles.rb_dX_1,'Value');
est(g.g_eop(5))=get(handles.rb_dY_1,'Value');
est(g.g_love)=0;
est(g.g_shida)=0;



del(g.g_clk(1))=0;
del(g.g_clk(2))=0;
del(g.g_zwd)=0;
del(g.g_tgr(1))=0;
del(g.g_tgr(2))=0;
del(g.g_coord(1))=get(handles.rb_acoor_0,'Value');
del(g.g_coord(2))=get(handles.rb_acoor_0,'Value');
del(g.g_coord(3))=get(handles.rb_acoor_0,'Value');
del(g.g_vel(1))=get(handles.rb_avel_0,'Value');
del(g.g_vel(2))=get(handles.rb_avel_0,'Value');
del(g.g_vel(3))=get(handles.rb_avel_0,'Value');
del(g.g_srade(1))=get(handles.rb_srade_0,'Value');
del(g.g_srade(2))=get(handles.rb_srade_0,'Value');
del(g.g_eop(1))=get(handles.rb_xpol_0,'Value');
del(g.g_eop(2))=get(handles.rb_ypol_0,'Value');
del(g.g_eop(3))=get(handles.rb_dut1_0,'Value');
del(g.g_eop(4))=get(handles.rb_dX_0,'Value');
del(g.g_eop(5))=get(handles.rb_dY_0,'Value');
del(g.g_love)=get(handles.rb_love_0,'Value');
del(g.g_shida)=get(handles.rb_shida_0,'Value');

red(g.g_clk(1))=get(handles.rb_clock_2,'Value');
red(g.g_clk(2))=get(handles.rb_clock_2,'Value');
red(g.g_zwd)=get(handles.rb_zwd_2,'Value');
red(g.g_tgr(1))=get(handles.rb_tgr_2,'Value');
red(g.g_tgr(2))=get(handles.rb_tgr_2,'Value');
red(g.g_eop(1))=get(handles.rb_xpol_2,'Value');
red(g.g_eop(2))=get(handles.rb_ypol_2,'Value');
red(g.g_eop(3))=get(handles.rb_dut1_2,'Value');
red(g.g_eop(4))=get(handles.rb_dX_2,'Value');
red(g.g_eop(5))=get(handles.rb_dY_2,'Value');

if est(g.g_vel(1))==1
    if est(g.g_coord(1))==0;
        fprintf('\n \n Station coordinates cannot be fixed if station velocities should be estimated!\n The station coordinates will be also estimated! \n\n')
    end
    est(g.g_coord(1))=1;
    est(g.g_coord(2))=1;
    est(g.g_coord(3))=1;
end


for i =1:length(parGS)
    parGS(i).id=[];
end

for i =1:length(parGS)
    if est(i)==1
        parGS(i).id=1;
    elseif del(i)==1
        parGS(i).id=0;
    end
end

for i =[g.g_clk,g.g_zwd,g.g_tgr,g.g_coord,g.g_vel,g.g_srade,g.g_eop]      % change if you add a new parameter
    if red(i)==1
        parGS(i).id=2;
    end
end


valbckwrd=get(handles.cb_bckwrd,'Value');
paths.bckwrdsol=valbckwrd;

%%

save(['../../OUT/GLOB/paths'],'paths');
save(['../../OUT/GLOB/parGS'],'parGS');
save(['../../OUT/GLOB/maxRMS'],'maxRMS');

guiglob_rf
uiwait(guiglob_rf)

close


% --- Executes on button press in cb_bckwrd.
function cb_bckwrd_Callback(hObject, eventdata, handles)
% hObject    handle to cb_bckwrd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_bckwrd


