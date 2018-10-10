function varargout = sched_optimization(varargin)
% SCHED_OPTIMIZATION MATLAB code for sched_optimization.fig
%      SCHED_OPTIMIZATION, by itself, creates a new SCHED_OPTIMIZATION or raises the existing
%      singleton*.
%
%      H = SCHED_OPTIMIZATION returns the handle to a new SCHED_OPTIMIZATION or the handle to
%      the existing singleton*.
%
%      SCHED_OPTIMIZATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCHED_OPTIMIZATION.M with the given input arguments.
%
%      SCHED_OPTIMIZATION('Property','Value',...) creates a new SCHED_OPTIMIZATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sched_optimization_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sched_optimization_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sched_optimization

% Last Modified by GUIDE v2.5 23-Feb-2017 10:55:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sched_optimization_OpeningFcn, ...
                   'gui_OutputFcn',  @sched_optimization_OutputFcn, ...
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


% --- Executes just before sched_optimization is made visible.
function sched_optimization_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sched_optimization (see VARARGIN)
handles.apply = 0;
handles.con = [];

% Choose default command line output for sched_optimization
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sched_optimization wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = sched_optimization_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout{1} = handles.condition_text.String;

delete(handles.figure1);


% --- Executes on selection change in condition1.
function condition1_Callback(hObject, eventdata, handles)
% hObject    handle to condition1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
allcontents = cellstr(get(hObject,'String'));
selected = get(hObject,'Value');
contents = allcontents{selected};

allPossible = {'minimum number of baselines','minimum number of scans','minimum number of different stations','minimum timespan in hours'};
boolSelected = strcmp(allPossible,contents);
    
if strcmp(contents,'remove condition')
    afterContent = {'add new condition',allPossible{:}};
    
    handles.condition1.String = afterContent;
    handles.condition2.String = afterContent;
    handles.condition3.String = afterContent;
    handles.condition4.String = afterContent;
    
    handles.condition2.Enable = 'off';
    handles.text4.Enable = 'off';
    handles.value1.Enable = 'off';
    handles.pushbutton1.Enable = 'off';
elseif strcmp(contents,'add new condition')
    % do nothing
else
    thisContent = {'remove condition',contents};
    handles.condition1.Value = 2;
    afterContent = {'add new condition',allPossible{~boolSelected}};
    handles.condition1.String = thisContent;
    handles.condition2.String = afterContent;
    handles.condition3.String = afterContent;
    handles.condition4.String = afterContent;
    
    
    handles.condition2.Enable = 'on';
    handles.text4.Enable = 'on';
    handles.value1.Enable = 'on';
    handles.pushbutton1.Enable = 'on';
end
writeCondition(hObject,handles);

% Hints: contents = cellstr(get(hObject,'String')) returns condition1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from condition1


% --- Executes during object creation, after setting all properties.
function condition1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to condition1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in condition2.
function condition2_Callback(hObject, eventdata, handles)
% hObject    handle to condition2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
allcontents = cellstr(get(hObject,'String'));
selected = get(hObject,'Value');
contents = allcontents{selected};

allPossible = {'minimum number of baselines','minimum number of scans','minimum number of different stations','minimum timespan in hours'};
boolSelected = strcmp(allPossible,contents);
boolSelected1 = strcmp(allPossible,handles.condition1.String(handles.condition1.Value));

if strcmp(contents,'remove condition')
    handles.condition1.String = {'remove condition',handles.condition1.String{handles.condition1.Value}};
    handles.condition1.Value = 2;
    
    afterContent = {'add new condition',allPossible{~boolSelected1}};
    handles.condition2.String = afterContent;
    handles.condition3.String = afterContent;
    handles.condition4.String = afterContent;
    
    handles.condition3.Enable = 'off';
    handles.text6.Enable = 'off';
    handles.value2.Enable = 'off';
    handles.combination1.Enable = 'off';
elseif strcmp(contents,'add new condition')
    % do nothing
else
    handles.condition1.String = handles.condition1.String(handles.condition1.Value);
    handles.condition1.Value = 1;
    thisContent = {'remove condition',contents};
    handles.condition2.String = thisContent;
    handles.condition2.Value = 2;
    afterContent = {'add new condition',allPossible{~(boolSelected|boolSelected1)}};
    handles.condition3.String = afterContent;
    handles.condition4.String = afterContent;
    
    
    handles.condition3.Enable = 'on';
    handles.text6.Enable = 'on';
    handles.value2.Enable = 'on';
    handles.combination1.Enable = 'on';
end
writeCondition(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns condition2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from condition2


% --- Executes during object creation, after setting all properties.
function condition2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to condition2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in condition3.
function condition3_Callback(hObject, eventdata, handles)
% hObject    handle to condition3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
allcontents = cellstr(get(hObject,'String'));
selected = get(hObject,'Value');
contents = allcontents{selected};

allPossible = {'minimum number of baselines','minimum number of scans','minimum number of different stations','minimum timespan in hours'};
boolSelected = strcmp(allPossible,contents);
boolSelected1 = strcmp(allPossible,handles.condition1.String(handles.condition1.Value));
boolSelected2 = strcmp(allPossible,handles.condition2.String(handles.condition2.Value));


if strcmp(contents,'remove condition')
    handles.condition2.String = {'remove condition',handles.condition2.String{handles.condition2.Value}};
    handles.condition2.Value = 2;
    
    afterContent = {'add new condition',allPossible{~(boolSelected1|boolSelected2)}};
    handles.condition3.String = afterContent;
    handles.condition4.String = afterContent;
    
    handles.condition4.Enable = 'off';
    handles.text7.Enable = 'off';
    handles.value3.Enable = 'off';
    handles.combination2.Enable = 'off';
elseif strcmp(contents,'add new condition')
    % do nothing
else
    handles.condition2.String = handles.condition2.String(handles.condition2.Value);
    handles.condition2.Value = 1;
    thisContent = {'remove condition',contents};
    handles.condition3.String = thisContent;
    handles.condition3.Value = 2;
    afterContent = {'add new condition',allPossible{~(boolSelected|boolSelected1|boolSelected2)}};
    handles.condition4.String = afterContent;
    
    handles.condition4.Enable = 'on';
    handles.text7.Enable = 'on';
    handles.value3.Enable = 'on';
    handles.combination2.Enable = 'on';
end
writeCondition(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns condition3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from condition3


% --- Executes during object creation, after setting all properties.
function condition3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to condition3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in condition4.
function condition4_Callback(hObject, eventdata, handles)
% hObject    handle to condition4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
allcontents = cellstr(get(hObject,'String'));
selected = get(hObject,'Value');
contents = allcontents{selected};

allPossible = {'minimum number of baselines','minimum number of scans','minimum number of different stations','minimum timespan in hours'};
boolSelected = strcmp(allPossible,contents);
boolSelected1 = strcmp(allPossible,handles.condition1.String(handles.condition1.Value));
boolSelected2 = strcmp(allPossible,handles.condition2.String(handles.condition2.Value));
boolSelected3 = strcmp(allPossible,handles.condition3.String(handles.condition3.Value));

if strcmp(contents,'remove condition')
    handles.condition3.String = {'remove condition',handles.condition3.String{handles.condition3.Value}};
    handles.condition3.Value = 2;
    
    afterContent = {'add new condition',allPossible{~(boolSelected1|boolSelected2|boolSelected3)}};
    handles.condition4.String = afterContent;
    
    handles.text8.Enable = 'off';
    handles.value4.Enable = 'off';
    handles.combination3.Enable = 'off';
elseif strcmp(contents,'add new condition')
    % do nothing
else
    handles.condition3.String = handles.condition3.String(handles.condition3.Value);
    handles.condition3.Value = 1;
    thisContent = {'remove condition',contents};
    handles.condition4.String = thisContent;
    handles.condition4.Value = 2;

    handles.text8.Enable = 'on';
    handles.value4.Enable = 'on';
    handles.combination3.Enable = 'on';
end
writeCondition(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns condition4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from condition4


% --- Executes during object creation, after setting all properties.
function condition4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to condition4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function value1_Callback(hObject, eventdata, handles)
% hObject    handle to value1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

v = str2double(get(hObject,'String'));
if isnan(v)
    handles.value1.String = '0';
end
writeCondition(hObject,handles);
% Hints: get(hObject,'String') returns contents of value1 as text
%        str2double(get(hObject,'String')) returns contents of value1 as a double


% --- Executes during object creation, after setting all properties.
function value1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to value1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function value2_Callback(hObject, eventdata, handles)
% hObject    handle to value2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
v = str2double(get(hObject,'String'));
if isnan(v)
    handles.value2.String = '0';
end
writeCondition(hObject,handles);
% Hints: get(hObject,'String') returns contents of value2 as text
%        str2double(get(hObject,'String')) returns contents of value2 as a double


% --- Executes during object creation, after setting all properties.
function value2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to value2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function value3_Callback(hObject, eventdata, handles)
% hObject    handle to value3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
v = str2double(get(hObject,'String'));
if isnan(v)
    handles.value3.String = '0';
end
writeCondition(hObject,handles);
% Hints: get(hObject,'String') returns contents of value3 as text
%        str2double(get(hObject,'String')) returns contents of value3 as a double


% --- Executes during object creation, after setting all properties.
function value3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to value3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function value4_Callback(hObject, eventdata, handles)
% hObject    handle to value4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
v = str2double(get(hObject,'String'));
if isnan(v)
    handles.value4.String = '0';
end
writeCondition(hObject,handles);
% Hints: get(hObject,'String') returns contents of value4 as text
%        str2double(get(hObject,'String')) returns contents of value4 as a double


% --- Executes during object creation, after setting all properties.
function value4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to value4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in combination1.
function combination1_Callback(hObject, eventdata, handles)
% hObject    handle to combination1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
writeCondition(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns combination1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from combination1


% --- Executes during object creation, after setting all properties.
function combination1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to combination1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in combination2.
function combination2_Callback(hObject, eventdata, handles)
% hObject    handle to combination2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
writeCondition(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns combination2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from combination2


% --- Executes during object creation, after setting all properties.
function combination2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to combination2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in combination3.
function combination3_Callback(hObject, eventdata, handles)
% hObject    handle to combination3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
writeCondition(hObject,handles);
% Hints: contents = cellstr(get(hObject,'String')) returns combination3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from combination3


% --- Executes during object creation, after setting all properties.
function combination3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to combination3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [con] = writeCondition(hObject,handles);
allPossible = {'minimum number of baselines','minimum number of scans','minimum number of different stations','minimum timespan in hours'};
translationTable = {'nbl','nscan','nsta','ntime'};

s1 = handles.condition1.String{handles.condition1.Value};
b1 = find(strcmp(s1,allPossible));
n = 0;

if b1>0
    s{1} = translationTable{b1};
    v{1} = handles.value1.String;
    c{1} = '';
    n = 1;
end

s2 = handles.condition2.String{handles.condition2.Value};
b2 = find(strcmp(s2,allPossible));
if b2>0
    s{2} = translationTable{b2};
    v{2} = handles.value2.String;
    c{2} = [' ' handles.combination1.String{handles.combination1.Value} ' '];
    n = 2;
end

s3 = handles.condition3.String{handles.condition3.Value};
b3 = find(strcmp(s3,allPossible));
if b3>0
    s{3} = translationTable{b3};
    v{3} = handles.value3.String;
    c{3} = [' ' handles.combination2.String{handles.combination2.Value} ' '];
    n = 3;
end

s4 = handles.condition4.String{handles.condition4.Value};
b4 = find(strcmp(s4,allPossible));
if b4>0
    s{4} = translationTable{b4};
    v{4} = handles.value4.String;
    c{4} = [' ' handles.combination3.String{handles.combination3.Value} ' '];
    n = 4;
end

str = '';
for i= 1:n
    if i == 1
        str = ['( ' s{i} ' > ' v{i} ' )'];
    else
        str = ['(' str  c{i} '( ' s{i} ' > ' v{i} ' )' ')'];
    end
end

if n>0
    con.type = s;
    con.value = v;
    con.combination = c;
    con.str = str;
    handles.condition = con;
    handles.condition_text.String = str;
    handles.con = con;
else 
    handles.con = [];
end

handles.condition_text.String = str;
handles.pushbutton1.BackgroundColor = [.94 .94 .94];
guidata(hObject, handles);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'),'waiting')
    uiresume(hObject)
else
    delete(hObject)
end
