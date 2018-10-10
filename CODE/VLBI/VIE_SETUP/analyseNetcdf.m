% #########################################################################
% #     analyseNetcdf
% #########################################################################
%
% DESCRIPTION
%   This function opens a GUI to analyse and visualise the content of netCDF data sets.
%
%
% CREATED
%   2016-01-??     Matthias Madzak
%
% REFERENCES
%
%
% COUPLING
% -
%
%
% INPUT
% - vgosFolderPath           - string giving the vgos folder to be analysed (if not given, nothing is loaded yet)
%
%
% OUTPUT
%
% CHANGES:
% - 2016-06-21, A. Hellerschmied: - Fixed several bugs and problems.
%                                 - Treatment of folders without falid data, e.g. "History", added.

function analyseNetcdf(varargin)

fprintf('analyze netcdf started\n');

%% get defaults
fig_pos=[378   211   937   646];
gray=get(0,'DefaultUIControlBackgroundColor');

foldersNoPrint={'Stub', 'CreateTime', 'CreatedBy', 'Program', ...
    'Subroutine','DataOrigin', 'TimeTag', 'TimeTagFile', 'Session',...
    'vgosDB_Version', 'CalcVer'};

%% SET UP GUI (empty)
fig = figure('Position',fig_pos,...
    'Color',gray,...
    'MenuBar','none',...
    'Resize','on',...
    'ResizeFcn',@resize,...
    'toolbar', 'figure',...
    'NumberTitle','off',...
    'Name','Analyse vgosDB (netcdf)',...
    'IntegerHandle','off',...
    'CloseRequestFcn',@cancel,...
    'Visible','off');
fig_invisible=figure('Visible', 'off');
ax_invisible=axes('Parent', fig_invisible);
figure(fig);

% %add Toolbar for Zoom in / out and Pan
% figureToolBar = uimenu('Label','Zoom Functions');
% uimenu(figureToolBar,'Label','Zoom In','Callback','zoom on');
% uimenu(figureToolBar,'Label','Zoom Out','Callback','zoom out');
% uimenu(figureToolBar,'Label','Pan','Callback','pan on');

% tbh = uitoolbar(fig);
%
% % Add a push tool to the toolbar
% a = [.20:.05:0.95];
% img1(:,:,1) = repmat(a,16,1)';
% img1(:,:,2) = repmat(a,16,1);
% img1(:,:,3) = repmat(flipud(a),16,1);
% pth = uipushtool(tbh,'CData',img1,...
%            'TooltipString','My push tool',...
%            'HandleVisibility','off');
% % Add a toggle tool to the toolbar
% img2 = rand(16,16,3);
% tth = uitoggletool(tbh,'CData',img2,'Separator','on',...
%            'TooltipString','Your toggle tool',...
%            'HandleVisibility','off');




% menu bar
f = uimenu('Label','File');
uimenu(f,'Label','Load folder','Callback',@loadNcFolderWithGui);
uimenu(f,'Label','Quit','Callback',@cancel,...
    'Separator','on','Accelerator','Q');


% close button
closebutton = uicontrol('Style','pushbutton',...
    'units', 'pixel',...
    'Position',[fig_pos(3)-70 5 70 20],... % [left bottom width height], lower left to lower left
    'String','Close',...
    'Enable','on',...
    'Callback',@cancel);



% folder edit textbox
% text_vgosFolder=uicontrol('Style', 'text', ...
%     'position', [5 fig_pos(4)-25 100, 20], 'String', 'VGOS folder',...
%     'HorizontalAlignment', 'left');
uipanel_vgosdatabase=uipanel('Title', 'VGOS database','units','pixel',...
    'Position', [5 fig_pos(4)-(5+50) 250 50]);
edit_vgosFolder=uicontrol('Style', 'edit',...
    'units', 'pixel',...
    'Position',[10 fig_pos(4)-40 180 20]);


% reload button
pushbutton_reloadFolder=uicontrol('Style', 'pushbutton',...
    'units', 'pixel',...
    'Position',[10+185 fig_pos(4)-40 50 20],...
    'String', 'Reload',...
    'Callback', @loadVgosFile);

% log textbox
text_log=uicontrol('Style', 'text', 'String', 'Log', 'Position', ...
    [5 100 330 20], 'HorizontalAlignment', 'left');
logTextbox=uicontrol('Style', 'edit',...
    'units', 'pixel',...
    'Position', [5 5 400 100],...
    'max', 3,...
    'min', 1,...
    'horizontalAlignment', 'left','Enable', 'inactive');

% axes (main)
axes_main = axes('Units','pixels', ...
    'parent', fig, ...
    'position', [fig_pos(3)-420 fig_pos(4)-320 400-10 300-10]);


% 'frame' wrapper
uipanel_wrapper=uipanel('Title', 'Wrapper','units','pixel',...
    'Position', [5 fig_pos(4)-(180) 250 120]);
table_wrapper=uitable('Parent', uipanel_wrapper, ...
    'Position', [5 25 240 80],'Data', ...
    {''}, 'ColumnName', {'File'}, ...
    'Columnwidth', {205},...
    'CellSelectionCallback', @selectWrapper);
button_openWrp=uicontrol('Parent',uipanel_wrapper, 'Style','pushbutton',...
    'units', 'pixel',...
    'Position',[5 5 50 20],...
    'String', 'Open',...
    'Callback', @openWrapper);

% 'frame' for vgos data
uipanel_vgosdatacontent=uipanel('Title', 'Database content','units','pixel',...
    'Position', [5 fig_pos(4)-(160+5+120) 250 100]);

text_statlist=uicontrol('Parent', uipanel_vgosdatacontent,...
    'Style', 'text', 'String', 'Stations',...
    'Position', [5 55+15 75 15], ...
    'horizontalAlignment', 'left');
popupmenu_statlist=uicontrol('Parent', uipanel_vgosdatacontent,...
    'Style', 'popupmenu', 'String',' ', 'Value', 1,...
    'Position', [5 40+15 100 15],'Callback', @updateTableVarlist);
text_sourlist=uicontrol('Parent', uipanel_vgosdatacontent,...
    'Style', 'text', 'String', 'Sources',...
    'Position', [3+2+100+5 55+15 75 15], ...
    'horizontalAlignment', 'left');
popupmenu_sourlist=uicontrol('Parent', uipanel_vgosdatacontent,...
    'Style', 'popupmenu', 'String',' ', 'Value', 1,...
    'Position', [3+2+100+5 40+15 100 15]);

text_folderlist=uicontrol('Parent', uipanel_vgosdatacontent,...
    'Style', 'text', 'String', 'Folders',...
    'Position', [5 10+15 50 15], ...
    'horizontalAlignment', 'left');
popupmenu_folderlist=uicontrol('Parent', uipanel_vgosdatacontent,...
    'Style', 'popupmenu', 'String',' ', 'Value', 1,...
    'Position', [5 10 100 15],...
    'Callback', @updatePopupFolders);



% text_subfolderArrow=uicontrol('Style', 'text', 'String', '$\rightarrow$',...
%     'Position', [5+2+100 fig_pos(4)-(40+5+20+15+7+20+15) 15 15]);
text_filelist=uicontrol('Parent', uipanel_vgosdatacontent,...
    'Style', 'text', 'String', 'Files',...
    'Position', [110 10+15 50 15], ...
    'horizontalAlignment', 'left');
popupmenu_filelist=uicontrol('Parent', uipanel_vgosdatacontent,...
    'Style', 'popupmenu', 'String', ' ', ...
    'Value', 1, 'Position', [110 10 100 15],...
    'Callback', @updateTableVarlist);


table_varlist=uitable('Position', [5 5+100+20+90 400 100],'Data', ...
    {'', '', ''}, 'ColumnName', {'Name', 'Type', 'Value'}, ...
    'Columnwidth', {160,90,70}, 'CellSelectionCallback', @selectVariable);
text_varlist=uicontrol('Style', 'text', 'String', 'Variables',...
    'Position', [5 5+100+20+100-1+90 60 20], 'HorizontalAlignment', 'left');
% button for sending a variable to workspace
pushbutton_send2workspace=uicontrol('Style', 'pushbutton',...
    'units', 'pixel',...
    'Position',[5+5+60 5+100+20+100+3+90 180 20],...
    'String', 'Send selected var(s) to workspace',...
    'Callback', @button_sendVar2workspace);

% attributes of variables
table_attlist=uitable('Position', [5 5+100+20 400 90],'Data', ...
    {'', ''}, 'ColumnName', {'Attribute', 'Value'}, ...
    'Columnwidth', {90,260});

%% FILL GUI
% if vgosFolder is given
if nargin>0
    set(edit_vgosFolder,'String', varargin{1});
    loadVgosFile();
end

% Make figure visible and hide handle.
set(fig,'Visible','on')


% --------------- CALLBACK FUNCTIONS --------------------------------------
    function cancel(varargin)
        delete(fig);
    end

    function resize(varargin)
        
        fig_pos=get(fig,'Position');
        try
            set(closebutton, 'Position',[fig_pos(3)-70 5 70 20]);
            set(edit_vgosFolder, 'Position',[10 fig_pos(4)-40 180 20]);
            set(uipanel_wrapper,'Position', [5 fig_pos(4)-(180) 250 120]);
            set(uipanel_vgosdatabase, 'Position',[5 fig_pos(4)-(5+50) 250 50]);
            set(pushbutton_reloadFolder, 'Position',[10+185 fig_pos(4)-40 50 20]);
            set(text_log, 'Position',[5 100 330 20]);
            set(logTextbox, 'Position',[5 5 400 100]);
            set(axes_main, 'Position',[fig_pos(3)-420 fig_pos(4)-320 400-10 300-10]);
            set(uipanel_vgosdatacontent, 'Position',[5 fig_pos(4)-(160+5+120) 250 100]);
            set(text_statlist, 'Position',[5 55+15 75 15]);
            set(popupmenu_statlist, 'Position',[5 40+15 100 15]);
            set(text_sourlist, 'Position',[3+2+100+5 55+15 75 15]);
            set(popupmenu_sourlist, 'Position',[3+2+100+5 40+15 100 15]);
            set(text_folderlist, 'Position',[5 10+15 50 15]);
            set(popupmenu_folderlist, 'Position',[5 10 100 15]);
            set(popupmenu_filelist, 'Position',[110 10 100 15]);
            set(table_varlist, 'Position',[5 5+100+20+90 400 100]);
            set(text_varlist, 'Position',[5 5+100+20+100-1+90 60 20]);
            set(pushbutton_send2workspace, 'Position',[5+5+60 5+100+20+100+3+90 180 20]);
            set(table_attlist, 'Position', [5 5+100+20 400 90]);
        catch
            
        end
    end

    function loadVgosFile(varargin)
        % load data (either programmatically from VieVS GUI or from File ->
        % Load folder
        % get edit textbox entry
        vgosFolder=get(edit_vgosFolder,'String');
        
        if ~strcmp(vgosFolder, '/') && ~strcmp(vgosFolder,'\')
            vgosFolder=[vgosFolder,'/'];
        end
        
        % load data from netcdf files
        if exist(vgosFolder,'dir')
            writeLog('Start loading vgos files');
            % read data
            %             [out_struct, nc_info,errLogical]=read_nc(vgosFolder);
            [out_struct, nc_info]=read_nc(vgosFolder);
            wrprs=dir([vgosFolder, '/*.wrp']);
            
            % save appdata
            setappdata(fig,'out_struct',out_struct);
            
            %             if errLogical==1 % if any error
            %                 writeLog('Error when reading nc files')
            %                 emptyGui(0);
            %             else
            % fill objects (eg popupmenues)
            fillGuiWVgosData(out_struct);
            writeLog('Loading done');
            set(table_wrapper,'Data', {wrprs.name}');
            %             end
        else
            writeLog('Folder not found');
        end
    end

    function button_sendVar2workspace(varargin)
        % This function is exectued when button "send variable to
        % workspace" is clicked. It sends the selcted variable(s) to the
        % matlab main workspace
        
        % get selected indices
        indices_variabTab=getappdata(fig,'indices_variabTab');
        
        if ~isempty(indices_variabTab)
            out_struct=getappdata(fig,'out_struct');
            
            [db,st,so,fo,fi,wr,va] = getGUIselections(varargin);
            
            % for all selected variables
            for kSelVar=1:length(va)
                curSelVar=va{kSelVar};
                curVal=out_struct.(fo).(fi).(curSelVar).val;
                assignin('base', curSelVar, curVal);
                writeLog(sprintf('''%s'' sent to workspace', curSelVar));
            end
        else
            writeLog('No variable selected!');
        end
    end

    function openWrapper(varargin)
        
        % This function is executed when the button to open a wrapper file
        % is clicked
        [db,st,so,fo,fi,wr,va] = getGUIselections(varargin);
        
        % if a wrapper is currently selected
        if ~isempty(wr)
            curWrp=strrep([db,'/',wr{1}],...
                '\','/');
            if ispc
                % try notepad++ (if installed)
                try
                    dosStr=sprintf('start notepad++ %s', curWrp);
                    [status,result] = dos(dosStr);
                catch
                    matlabVersion = ver('MATLAB');
                    if str2double(matlabVersion.Version)<=8
                        dos(['start wordpad ',curWrp]);
                    else
                        % does not work for 7.11 (R2010b):
                        winopen(curWrp) % Open with default editor which is set for *.opt files
                    end
                end
                
            else
                system(['xterm -e ''vi ',curWrp '''']);
            end
            writeLog('Wrapper file opened');
        else
            writeLog('No wrapper selected');
        end
        
    end

    function loadNcFolderWithGui(varargin)
        % File -> Load folder
        vgosDir='../DATA/vgosDB/';
        if ~exist(vgosDir,'dir')
            mkdir(vgosDir);
        end
        out = uipickfiles('FilterSpec', vgosDir);
        if iscell(out)
            out=out';
            % get only foldername % though this removes the path information (if
            % different database to be used -> modify that!)
            
            for iF=1:length(out)
                curSlash=sort([strfind(out{iF},'/'), strfind(out{iF},'\')]);
                out{iF}=out{iF}(curSlash(end)+1:end);
            end
            set(edit_vgosFolder,'String', [vgosDir, num2str(yy2yyyy(out{1}(1:2))),...
                '/', out{1}]);
            loadVgosFile();
            
            %     updateInputFilesBox(hObject, eventdata,handles,out)
        end
        
    end

    function emptyGui(emptyLog,varargin)
        % Empty GUI components (Log is only emptied when emptyLog==1)
        set(text_statlist, 'String', 'stations');
        set(text_sourlist, 'String', 'sources');
        set(popupmenu_statlist, 'Value', 1);
        set(popupmenu_sourlist, 'Value', 1);
        set(popupmenu_folderlist, 'Value', 1);
        set(popupmenu_filelist, 'Value', 1);
        set(popupmenu_statlist, 'String', ' ');
        set(popupmenu_sourlist, 'String', ' ');
        set(popupmenu_folderlist, 'String', ' ');
        set(popupmenu_filelist, 'String', ' ');
        set(table_varlist,'Data', {'', '', ''});
        set(table_wrapper,'Data', {''});
        set(table_attlist,'Data', {'',''});
        if emptyLog==1
            set(logTextbox,'String', '');
        end
    end

    function selectVariable(varargin)
        % This function is exectuted when a variable in the listbox is
        % seltected (then this variable might be plotted)
        % This function is also executed when the variable table is newly
        % filled with data -> indices are then empty
        cla(axes_main);
        % save indices in any case (so that other functions know that
        % nothing is selected
        setappdata(fig,'indices_variabTab',varargin{2}.Indices); % save the slection-indices as appdata
        if ~isempty(varargin{2}.Indices)
            out_struct=getappdata(fig,'out_struct');
            
            [db,st,so,fo,fi,wr,va] = getGUIselections(varargin);
            
            % stat struct contains more entries: get index
            if strcmpi(fo,'stat')
                curVarStruct=out_struct.(fo)(get(popupmenu_statlist,'value')).(fi).(va{1});
            else
                curVarStruct=out_struct.(fo).(fi).(va{1});
            end
            
            % update attribute list
            if isfield(curVarStruct,'attr')
                set(table_attlist,'Data', [{curVarStruct.attr.name}',...
                    {curVarStruct.attr.val}']);
            end
            
            
            
            % try to plot value (first var if multiple selected
            if isfield(curVarStruct,'val')
                if ~ischar(curVarStruct.val)
                    if min(size(curVarStruct.val)) == 1
                        try
                            plot(axes_main, curVarStruct.val);
                        catch
                            writeLog('Cannot plot variable: Unknown problem.');
                        end
                    else
                        writeLog('Cannot plot variable: invalid dimension.');
                    end
                else    
                    writeLog('Cannot plot variable: type=char.');
                end
            end
        end
        
    end

    function selectWrapper(varargin)
        % This function is exectuted when a variable in the wrapper table
        % is seltected
        setappdata(fig,'indices_wrp',varargin{2}.Indices); % save the slection-indices as appdata
    end

    function fillGuiWVgosData(out_struct, varargin)
        % Fill the GUI (popupmenus) with vgos data (eg station/source list,
        % folders...)
        curPopStatStr=cellstr(get(popupmenu_statlist,'String'));
        curPopSourStr=cellstr(get(popupmenu_sourlist,'String'));
        curPopFolderStr=cellstr(get(popupmenu_folderlist,'String'));
        if ~isempty(curPopStatStr)
            curSelStat=curPopStatStr{get(popupmenu_statlist,'Value')};
        end
        if ~isempty(curPopSourStr)
            curSelSour=curPopSourStr{get(popupmenu_sourlist,'Value')};
        end
        if ~isempty(curPopFolderStr)
            curSelFolder=curPopFolderStr{get(popupmenu_folderlist,'Value')};
        end
        
        set(text_statlist,'String', ...
            [num2str(size(out_struct.head.StationList.val',1)), ' stations']);
        set(text_sourlist,'String', ...
            [num2str(size(out_struct.head.SourceList.val',1)), ' sources']);
        set(popupmenu_statlist,'Value', 1);
        set(popupmenu_sourlist,'Value', 1);
        set(popupmenu_statlist,'String', cellstr(out_struct.head.StationList.val'));
        set(popupmenu_sourlist,'String', cellstr(out_struct.head.SourceList.val'));
        
        set(popupmenu_folderlist,'String',fieldnames(out_struct)); % folder list
        
        
        % re-select previously selected station
        if ~isempty(curPopStatStr)
            if sum(strcmpi(cellstr(out_struct.head.StationList.val'),curSelStat))>0
                set(popupmenu_statlist, 'Value', find(strcmpi(cellstr(out_struct.head.StationList.val'),curSelStat)));
            end
        end
        if ~isempty(curPopSourStr)
            if sum(strcmpi(cellstr(out_struct.head.SourceList.val'),curSelSour))>0
                set(popupmenu_sourlist, 'Value', find(strcmpi(cellstr(out_struct.head.SourceList.val'),curSelSour)));
            end
        end
        
        if ~isempty(curPopFolderStr)
            if sum(strcmpi(fieldnames(out_struct),curSelFolder))>0
                set(popupmenu_folderlist, 'Value', find(strcmpi(fieldnames(out_struct),curSelFolder)));
            end
        end
        
        % get subfolder and write to subfolder popupmenu
        updatePopupFolders();
        
        updateTableVarlist();
        
    end

    function updateTableVarlist(varargin)
        % Update the listbox containing the variables (after file is
        % selected in popupmenu)
        cla(axes_main);
        set(table_attlist,'Data',{'',''}); % make variable list empty when variables (might) change
        out_struct=getappdata(fig,'out_struct');
        
        
        [db,st,so,fo,fi,wr,va] = getGUIselections;
        
        if ~isempty(fi)
            % write vars to listbox
            if length(out_struct.(fo))>1
                if strcmpi(fo,'stat')
                    % separate "equal" command to find errors if also other folder
                    % requires an index!
                    statInd=get(popupmenu_statlist,'Value');
                else
                    writeLog('ERROR CHECK HERE');
                end
                curFileContent=out_struct.(fo)(statInd).(fi);
            else
                curFileContent=out_struct.(fo).(fi);
            end
            vars=fieldnames(curFileContent);
            varsGood=vars(~ismember(vars,foldersNoPrint));
            sizeData=cell(length(varsGood),1); % cell for "1x4 double"
            firstValsStr=cell(length(varsGood),1); % cell for 0.4, 0.1, 0.3 ...
            for kVar=1:length(varsGood)
                curFirstValStr='';
                curField=curFileContent.(varsGood{kVar});
                if isfield(curField,'val')
                    sizeCurVar=size(curField.val);
                    curClass=class(curField.val);
                    curString=sprintf('%1.0fx%1.0f %s',...
                        sizeCurVar, curClass);
                    % get the first few entries as user-info-string
                    if isnumeric(curField.val)
                        for kFirstVals=1:min([length(curField.val),3]) % for (up to) first three array entries
                            if curFirstValStr>1
                                curFirstValStr=[curFirstValStr, ', '];
                            end
                            curFirstValStr=[curFirstValStr,num2str(curField.val(kFirstVals))];
                        end
                        % add ... (if more data available)
                        if length(curField.val)>3
                            curFirstValStr=[curFirstValStr,'...'];
                        end
                    elseif ischar(curField.val)
                        curFirstValStr=curField.val(1:min([length(curField.val),3]));
                        curFirstValStr=curFirstValStr(:)'; % make 1xn char array
                        if length(curField.val)>3
                            curFirstValStr=[curFirstValStr,'...'];
                        end
                    else
                        keyboard;
                    end
                else
                    curString='?';
                end
                sizeData{kVar,1}=curString;
                firstValsStr{kVar,1}=curFirstValStr;
            end

            % write data
            set(table_varlist, 'Data', [varsGood, sizeData,firstValsStr])
        else
            set(table_varlist, 'Data', [])
            writeLog(['No valid data in folder: ', fo]);
        end 
    end

    function updatePopupFolders(varargin)
        % Updates the subfolder-popupmenu to the content of the selected
        % folder
        cla(axes_main);
        out_struct=getappdata(fig,'out_struct');
        
        curPopSubfolderStr=cellstr(get(popupmenu_filelist, 'String'));
        if ~isempty(curPopSubfolderStr)
            curSelSubf=curPopSubfolderStr{get(popupmenu_filelist,'Value')};
        end
        
        % (main) folders
        folders=get(popupmenu_folderlist,'String');
        subf=fieldnames(out_struct.(folders{get(popupmenu_folderlist,'Value')}));
        % remove folders which are not needed
        subf=subf(~ismember(subf,foldersNoPrint));
        
        % set subfolder popupmenu
        set(popupmenu_filelist,'Value',1);
        if isempty(subf); subf=' '; end
        set(popupmenu_filelist,'String',subf);
        
        % set to previous value
        if ~isempty(curPopSubfolderStr)
            if sum(strcmpi(subf,curSelSubf))>0
                set(popupmenu_filelist,'Value',find(strcmpi(subf,curSelSubf)));
            end
        end
        
        % update also the file list
        updateTableVarlist();
    end

    function writeLog(inText,varargin)
        % this function writes a log. Usually: "11:40:01 GUI started", if
        % varargin{1}==0, no time is written
        
        writeTime=1;
        
        if nargin>1
            if varargin{1}==0
                writeTime=0;
            end
        end
        if writeTime==1
            t=clock;
            timeStr=sprintf('%02.0f:%02.0f:%02.0f',t(4:6));
        else
            timeStr='';
        end
        
        if isempty(get(logTextbox,'String'))
            newCell=[timeStr, ' ', inText];
        else
            newCell=[{[timeStr, ' ', inText]};...
                get(logTextbox,'String')];
        end
        set(logTextbox,'String', newCell);
    end

% ======================================================================
% OTHER FUNCTIONS
% ======================================================================

    function [db,st,so,fo,fi,wr,va] = getGUIselections(varargin)
        % This function returns the currently selected GUI-options
        %  varargin: nothing required
        %
        % db = (string) VGOS database
        % st = (string) current selected station
        % so = (string) current selected source
        % fo = (string) current selected folder
        % fi = (string) current selected file
        % wr = (cellstr) current selected wrapper,  [] if non is selected
        % va = (cellstr) current selected variable, [] if non is selected
        %
        %
        
        db=get(edit_vgosFolder,'String');
        
        allStations=get(popupmenu_statlist,'String');
        st=allStations{get(popupmenu_statlist,'Value')};
        
        allSources=get(popupmenu_sourlist, 'String');
        so=allSources{get(popupmenu_sourlist,'Value')};
        
        allFolders=get(popupmenu_folderlist,'String');
        fo=allFolders{get(popupmenu_folderlist,'Value')};
        
        allFiles=get(popupmenu_filelist,'String');
        if strcmp(allFiles, ' ')
            fi = [];
        else
            fi=allFiles{get(popupmenu_filelist,'Value')};
        end
        
        allWrappers=get(table_wrapper,'Data');
        indices_wrp=getappdata(fig,'indices_wrp');
        if isempty(indices_wrp)
            wr=[];
        else
            wr=allWrappers(indices_wrp(:,1),1);
        end
        
        allVarnames=get(table_varlist,'Data');
        indices_variabTab=getappdata(fig,'indices_variabTab');
        if isempty(indices_variabTab)
            va=[];
        else
            va=allVarnames(indices_variabTab(:,1),1);
        end
        
    end


end % .m GUI file