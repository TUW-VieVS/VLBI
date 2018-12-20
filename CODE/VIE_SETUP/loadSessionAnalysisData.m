% #########################################################################
% #     loadSessionAnalysisData
% #########################################################################
%
% DESCRITPION
% This function load data for session analysis tool.
% This function loads all data needed for the plotting panel
% "sessionAnalysis" which can be used for e.g. plotting the network of a
% session, the baseline lnegth repeatability,...
%
% need following files: 
% x_      (level3) ... baseline lnegth rep
% antenna (level3) ... map network
% atpa    (level3) ... corrmat
% opt_    (level3) ... opt_.mo
%
% AUTHOR 
%   moved in separate function by A. Girdiuk
%
% INPUT
%   handles      structure from the GUI (also containing data, eg residuals)
%   hObject      ...unused so far.
%
% OUTPUT
%   handles      structure from the GUI (also containing data, eg residuals)
%
% CHANGES
%   2016-08-26, A. Girdiuk: bug-fix
%   2016-09-21, A. Girdiuk: bug-fix
%   2016-10-11, A. Hellerschmied: Fields (soude, soura) added, if missing in loaded x_ file.


function handles=loadSessionAnalysisData(hObject, handles)

% if the "main" load button is clicked on
if strcmpi(get(hObject, 'Tag'), 'pushbutton_plot_sessionAnalysis_load')
    cla(handles.axes_plot_sessionAnalysis);
    set(handles.popupmenu_plot_sessionAnalysis_session, 'Value', 1)
    set(handles.popupmenu_plot_sessionAnalysis_corrMat_secondPar, 'Value', 1)
    set(handles.popupmenu_plot_sessionAnalysis_corrMat_firstPar, 'Value', 1)
    set(handles.popupmenu_plot_sessionAnalysis_session, 'String', ' ')
    set(handles.popupmenu_plot_sessionAnalysis_corrMat_firstPar, 'String', ' ')
    set(handles.popupmenu_plot_sessionAnalysis_corrMat_secondPar, 'String', ' ')
    
    % get all x_ files in folder
    allSubfolders=get(handles.popupmenu_plot_sessionAnalysis_subfolder, 'string');
    curSubfolder=allSubfolders{get(handles.popupmenu_plot_sessionAnalysis_subfolder, 'Value')};
    
    % get all x_ files in folder
    allX_Files=dir(['../DATA/LEVEL3/', curSubfolder, '/x_*.mat']);
    
    if ~isempty(allX_Files)
    
        % enable parts of the interface
        set(handles.popupmenu_plot_sessionAnalysis_session, 'Enable', 'On')
        set(handles.radiobutton_plot_sessionAnalysis_network, 'Enable', 'On')
        set(handles.radiobutton_plot_sessionAnalysis_baselLeRep, 'Enable', 'On')
        if get(handles.radiobutton_plot_sessionAnalysis_baselLeRep, 'Value')
            set(handles.checkbox_plot_sessionAnalysis_baselineNames, 'Enable', 'On')
        end
        set(handles.radiobutton_plot_sessionAnalysis_corMatrix, 'Enable', 'On')
        if get(handles.radiobutton_plot_sessionAnalysis_corMatrix, 'Value')
            set(handles.popupmenu_plot_sessionAnalysis_corrMat_firstPar, 'Enable', 'On')
            set(handles.popupmenu_plot_sessionAnalysis_corrMat_secondPar, 'Enable', 'On')
        end

        % create short names for popupmenu
        sessionnamesShort=strrep({allX_Files.name}, 'x_', '');
        sessionnamesShort=strrep(sessionnamesShort, '.mat', '');
        % save process_list as folder contains
        sessionShortToprocess_list=strrep({allX_Files.name}, 'x_', 'yyyy/');
        handles.data.plot.sessionAnalysis.sessionnamesShort{1}.list=char(strrep(sessionShortToprocess_list, '.mat', ''));

        % create waitbar
        h_waitbar = waitbar(0,'Please wait...');


        % for all res files
        nX_Files=size(allX_Files,1);        
        
        % remove previous data if it exists 
        if isfield(handles.data.plot.sessionAnalysis,'x_files')
            handles.data.plot.sessionAnalysis.x_files(1).x_(1:end)=[];
        end

        for iFile=1:nX_Files

            % load x_ file
            load(['../DATA/LEVEL3/', curSubfolder, '/', allX_Files(iFile).name]);

            % add field "units" which does not exist in older x_ files
            if ~isfield(x_, 'units')
                x_.units=[];
            end
            
            % Check, if all fields in x_ are available.
            % => Init. empty fields, id required:
            if ~isfield(x_, 'soura')
                x_.soura = [];
            end
            if ~isfield(x_, 'soude')
                x_.soude = [];
            end

            % make the same field-order as in the first x_ structure
            % add x_ to handles struct
            if iFile>1
                x_=orderfields(x_, handles.data.plot.sessionAnalysis.x_files(1).x_(iFile-1));
                handles.data.plot.sessionAnalysis.x_files(1).x_(iFile)=x_;
            else
                handles.data.plot.sessionAnalysis.x_files(1).x_=x_;
            end

            % load opt_ file
            load(['../DATA/LEVEL3/', curSubfolder, '/', strrep(allX_Files(iFile).name, 'x_', 'opt_')]);
            
%             if ~isfield(opt_, 'level1OutDir')
%                 level1OutDir=curSubfolder;
%             else
%                 level1OutDir=opt_.level1OutDir;
%             end
            
            % add opt_ to handles struct
            handles.data.plot.sessionAnalysis.optFiles(1).opt(iFile).opt=opt_;
            
            % load antenna file and add to handles struct
%             load(['../DATA/LEVEL1/', level1OutDir, '/', sessionnamesShort{iFile}, '_antenna.mat'])
            load(['../DATA/LEVEL3/', curSubfolder, '/', sessionnamesShort{iFile}, '_antenna.mat'])
            handles.data.plot.sessionAnalysis.antennaFiles(1).antenna(iFile).antenna=antenna;

            % load atpa file and add to handles struct
            load(['../DATA/LEVEL3/', curSubfolder, '/', strrep(allX_Files(iFile).name, 'x_', 'atpa_')]);
            handles.data.plot.sessionAnalysis.atpaFiles(1).atpa(iFile).atpa=atpa_;
            
            % update waitbar
            waitbar(iFile/nX_Files,h_waitbar,sprintf('session %1.0f / %1.0f is loaded.', iFile, nX_Files))
        end % for allres files

        nEntries=size(handles.data.plot.sessionAnalysis.x_files(1).x_,2);
        if nEntries>nX_Files
            handles.data.plot.sessionAnalysis.x_files(1).x_(nEntries:end)=[];
        end

        % close waitbar
        close(h_waitbar);

        % update popupmenus
        % (a) sessions
        set(handles.popupmenu_plot_sessionAnalysis_session, 'String', sessionnamesShort);
    end
    
else % second (or more) load button is clicked
    % get the button (2,3 or 4) which was clicked
    tagOfButton=get(hObject, 'Tag');
    buttonNr=str2double(tagOfButton(end));
    
    set(eval(['handles.popupmenu_plot_sessionAnalysis_session', num2str(buttonNr)]), 'Value', 1)
    set(eval(['handles.popupmenu_plot_sessionAnalysis_session', num2str(buttonNr)]), 'String', ' ')
    
    % get all x_ files in folder
    allSubfolders=get(eval(['handles.popupmenu_plot_sessionAnalysis_subfolder', num2str(buttonNr)]), 'string');
    curSubfolder=allSubfolders{get(eval(['handles.popupmenu_plot_sessionAnalysis_subfolder', num2str(buttonNr)]), 'Value')};
    
    % get all x_ files in folder
    allX_Files=dir(['../DATA/LEVEL3/', curSubfolder, '/x_*.mat']);
    
    % remove previous data if it exists 
    if length(handles.data.plot.sessionAnalysis.x_files)>=buttonNr
        handles.data.plot.sessionAnalysis.x_files(buttonNr).x_(1:end)=[];
    end
    
    if ~isempty(allX_Files)
    
        % enable parts of the interface
        set(eval(['handles.checkbox_plot_sessionAnalysis_add', num2str(buttonNr)]), 'Enable', 'On')
        set(eval(['handles.popupmenu_plot_sessionAnalysis_session', num2str(buttonNr)]), 'Enable', 'On')

        % create short names for popupmenu
        sessionnamesShort=strrep({allX_Files.name}, 'x_', '');
        sessionnamesShort=strrep(sessionnamesShort, '.mat', '');
        
        % save process_list as folder contains
        sessionShortToprocess_list=strrep({allX_Files.name}, 'x_', 'yyyy/');
        handles.data.plot.sessionAnalysis.sessionnamesShort{buttonNr}.list=char(strrep(sessionShortToprocess_list, '.mat', ''));

        % create waitbar
        h_waitbar = waitbar(0,'Please wait...');

        % for all res files
        nX_Files=size(allX_Files,1);
        for iFile=1:nX_Files

            % load x_ file
            load(['../DATA/LEVEL3/', curSubfolder, '/', allX_Files(iFile).name]);

            % add field "units" which does not exist in older x_ files
            if ~isfield(x_, 'units')
                x_.units=[];
            end

            % make the same field-order as in the first x_ structure
            % add x_ to handles struct
            if iFile>1
                x_=orderfields(x_, handles.data.plot.sessionAnalysis.x_files(buttonNr).x_(iFile-1));
                handles.data.plot.sessionAnalysis.x_files(buttonNr).x_(iFile)=x_;
            else
                handles.data.plot.sessionAnalysis.x_files(buttonNr).x_=x_;
            end

            % load opt_ file
            load(['../DATA/LEVEL3/', curSubfolder, '/', ...
                strrep(allX_Files(iFile).name, 'x_', 'opt_')]);
            
%             if ~isfield(opt_, 'level1OutDir')
%                 level1OutDir=curSubfolder;
%             else
%                 level1OutDir=opt_.level1OutDir;
%             end
            
            % add opt_ to handles struct
            handles.data.plot.sessionAnalysis.optFiles(buttonNr).opt(iFile).opt=opt_;
            
            % load antenna file
%             load(['../DATA/LEVEL1/', level1OutDir, '/', sessionnamesShort{iFile}, '_antenna.mat'])
            load(['../DATA/LEVEL3/', curSubfolder, '/', sessionnamesShort{iFile}, '_antenna.mat'])
            handles.data.plot.sessionAnalysis.antennaFiles(buttonNr).antenna(iFile).antenna=antenna;

            % load atpa file
            load(['../DATA/LEVEL3/', curSubfolder, '/atpa_', sessionnamesShort{iFile}, '.mat'])
            handles.data.plot.sessionAnalysis.atpaFiles(buttonNr).atpa(iFile).atpa=atpa_;

            % update waitbar
            waitbar(iFile/nX_Files,h_waitbar,sprintf('session %1.0f / %1.0f is loaded.', iFile, nX_Files))
        end % for allres files

        nEntries=size(handles.data.plot.sessionAnalysis.x_files(buttonNr).x_,2);
        if nEntries>nX_Files
            handles.data.plot.sessionAnalysis.x_files(buttonNr).x_(nEntries:end)=[];
        end

        % close waitbar
        close(h_waitbar);

        % update popupmenus
        % (a) sessions
        set(eval(['handles.popupmenu_plot_sessionAnalysis_session', num2str(buttonNr)]), 'String', sessionnamesShort);
    end

end
