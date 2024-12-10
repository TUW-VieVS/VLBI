% #########################################################################
% #     pushbutton_plot_residuals_writeAmbiguities_Callback
% #########################################################################
%
% DESCRITPION
% Executes on button press in pushbutton_plot_residuals_writeAmbiguities.
%
% AUTHOR 
%   Leo Baldreich, based on
%   pushbutton_plot_residuals_removeOutliers_Callback by Matthias Madzak
%
% INPUT
% hObject    handle to pushbutton_plot_residuals_writeAmbiguities (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% OUTPUT
%
% CHANGES
%

function ext_pushbutton_plot_residuals_writeAmbiguities_Callback(hObject, eventdata, handles)

if size(handles.data.plot.outlierMarksHandle,1)==1
    nSelAmbigs=1;
else
    nSelAmbigs=size(handles.data.plot.outlierMarksHandle(1).XData,2);
end
% if no ambiguity was selected -> write messge box
if nSelAmbigs==0
    msgbox('No value was selected!', 'No value selected', 'warn');
else
%     % (1) ambiguities subfolder
%     allPopupmenuEntriesAMBfolder=get(handles.popupmenu_setInput_ambDir, 'String');
%     if strcmp(get(handles.popupmenu_setInput_ambDir, 'String'), ' ')
%         AMBsubFolder='';
%     else
%         AMBsubFolder=allPopupmenuEntriesAMBfolder{get(handles.popupmenu_setInput_ambDir, 'Value')};
%     end
            % get the AMB filename
            % (2) year
            % firstPlottedMjd=handles.data.plot.res(get(handles.popupmenu_plot_residuals_session, 'Value')).mjd(1);
            % get the year from session name - changed 7.12. S. Boehm
            
            % (3) session name
            chosenSessionInd=get(handles.popupmenu_plot_residuals_session, 'Value');
            allPopupmenuEntriesSessions=get(handles.popupmenu_plot_residuals_session, 'String');
            session = allPopupmenuEntriesSessions{chosenSessionInd};
            
            % ##### Check, if the input dataset file followed the standard naming convention #####
            % => YYMMMDDcc_Nnnn (c...alphabetic character, n...number)
%             flag_std_naming_convention = true;
% 
%             % Total length = 14 char
%             if length(session) ~= 14
%                 flag_std_naming_convention = false; 
%             end
% 
%             if flag_std_naming_convention
%                 % "_" at car. 10
%                 if ~strcmp(session(10), '_')
%                     flag_std_naming_convention = false;
%                 end
%                 % first two characters are numbers:
%                 [~, status_1] = str2num(session(1:2));
%                 % char. 6+7 are numbers:
%                 [~, status_2] = str2num(session(6:7));
%                 % char. 12-14 are numbers:
%                 [~, status_3] = str2num(session(12:14));
%                 if ~(status_1 && status_2 && status_3)
%                     flag_std_naming_convention = false;
%                 end
%             end
% 
%             if flag_std_naming_convention
%                 % ##### Standard naming convention is used for this session: #####
% 
%                 % Get the year from the session name:
%                 if str2double(session(1:2)) > 75
%                     yearStr = ['19', session(1:2)];
%                 else
%                     yearStr = ['20', session(1:2)];
%                 end
% 
%                 
%             else
                % ##### Non-Standard naming convention is used for this session: #####
                % e.g. when using .vso input files
                % => Get yearStr from the opt_ file!

                % #### Load opt_ file from LEVEL 3 (sub-)directory: ####
                % Get sub-dir.:
                allSubfolders   = get(handles.popupmenu_plot_residuals_folder, 'string');
                curSubfolder    = allSubfolders{get(handles.popupmenu_plot_residuals_folder, 'Value')};

                % load opt_ file and get the year:
                load(['../DATA/LEVEL3/', curSubfolder, '/opt_', session, '.mat']);
                yearStr = opt_.data_filepath(end-4:end-1);

                % load _scan:
                load(['../DATA/LEVEL3/', curSubfolder, '/', session, '_scan.mat']);

                % load _parameter:
                load(['../DATA/LEVEL3/', curSubfolder, '/', session, '_parameter.mat']);

%             end
            
            band_letter = parameter.vie_init.vgosDb_observation_parameter(end);
            AMBfolder = ['../DATA/AMB/PU/', yearStr, '/'];
            AMBfilename = [allPopupmenuEntriesSessions{chosenSessionInd}, '_', band_letter, '.AMB'];

            % Check if outlier folder exists => if not => create it!
            if ~exist(AMBfolder, 'dir')
                mkdir(AMBfolder);
            end
            
       % ##### get the values just plotted (station, baseline, all-wise) #####
       % currently only plotting per baseline allows ambig selection
            % ### per station ###
%             if get(handles.radiobutton_plot_residuals_perStat, 'Value')
%                 curValues=sum(handles.data.plot.res(chosenSessionInd).baselineOfObs==get(handles.popupmenu_plot_residuals_station, 'Value'),2);
%                 
%             % ### per baseline ###
%             elseif get(handles.radiobutton_plot_residuals_perBasel, 'Value')
                allBaselines=get(handles.popupmenu_plot_residuals_baseline, 'String');
                curBaseline=allBaselines{get(handles.popupmenu_plot_residuals_baseline, 'Value')};
                if curBaseline == "[all Baselines]"
                    curValues=ones(size(handles.data.plot.res(chosenSessionInd).mjd,1),1);
                else
                    stat1=allBaselines{get(handles.popupmenu_plot_residuals_baseline,'Value')}(1:8);
                    stat2=allBaselines{get(handles.popupmenu_plot_residuals_baseline,'Value')}(10:17);
                    statNr1=~cellfun(@isempty, strfind(handles.data.plot.res(chosenSessionInd).allStatNames, stat1));
                    statNr2=~cellfun(@isempty, strfind(handles.data.plot.res(chosenSessionInd).allStatNames, stat2));
                    valsLogicalsOfFirstStat=handles.data.plot.res(chosenSessionInd).baselineOfObs==find(statNr1);
                    valsLogicalsOfSecondStat=handles.data.plot.res(chosenSessionInd).baselineOfObs==find(statNr2);
                    curValues=sum(valsLogicalsOfFirstStat+valsLogicalsOfSecondStat,2)==2;
                end
%             % ### per source ###
%             elseif get(handles.radiobutton_plot_residuals_perSource, 'Value')
%                 curSource=get(handles.popupmenu_plot_residuals_source, 'Value');
%                 curValues=handles.data.plot.res(chosenSessionInd).source==curSource;
%             
%             % ### all are plotted ###
%             else
%                 curValues=ones(size(handles.data.plot.res(chosenSessionInd).mjd,1),1);
%             end
                
%% ############# GET AMBS ##############
            
            % ##### get epochs (mjd) for ambiguities #####
            
            % Currently plotted values (indices)
            curValuesIndices = find(curValues);
            
            % Selected Data (Ambigs) in residuals plot window:
            ambigEpochsHours = get(handles.data.plot.outlierMarksHandle, 'XData');
%             if(length(outlierEpochsHours)>1)
%                 outlierEpochsHours = outlierEpochsHours{1};
%             end
            if size(ambigEpochsHours,1)>1
                ambigEpochsHours = ambigEpochsHours{1};
            else
                ambigEpochsHours = ambigEpochsHours(1);
            end
                
            
            AmbigValues = get(handles.data.plot.outlierMarksHandle, 'YData');
%             if(length(OutlierValues)>1)
%                 OutlierValues = OutlierValues{1};
%             end
            if(size(AmbigValues,1)>1)
                AmbigValues = AmbigValues{1};
            else
                AmbigValues = AmbigValues(1);
            end

            % Conversion of x-values: "hours from session start" => "MJD": 
            SessionStartTimeMJD =  handles.data.plot.res(chosenSessionInd).mjd(1);
            selectedAmbigEpochsMJD = SessionStartTimeMJD + ambigEpochsHours / 24;

            % Get Baseline indices for selected ambiguities with a 2-dimensional search
            % approach:
            
            % 1.) ### Compare x-values (MJD): ###
            
            % Delete multiple entries in vector "selectedOutlierEpochsMJD":
            ambigEpochsMJD = unique(selectedAmbigEpochsMJD);
            numOfambigEpochsMJD = length(ambigEpochsMJD);
            
            xIndices = zeros(nSelAmbigs, 1);
            index = 1;
            
            for indexEpochs = 1 : numOfambigEpochsMJD
                indicesCurrentEpoch = find(handles.data.plot.res(chosenSessionInd).mjd(curValuesIndices) == ambigEpochsMJD(indexEpochs));
                lengthIndicesCurrentEpoch = length(indicesCurrentEpoch);
                xIndices(index : (index + lengthIndicesCurrentEpoch - 1)) = indicesCurrentEpoch;
                index = index + lengthIndicesCurrentEpoch;
            end
            
            % 2.) ### Compare First/Main Solution Values: ###
            
            % Get all y-values of current selection (baseline, station, all values, sources):
            % Get x-Values (first or main solution):

            % if first solution is chosen
            if get(handles.radiobutton_plot_residuals_firstSolution, 'Value')
                valForXIndices = handles.data.plot.res(chosenSessionInd).firstVal(curValuesIndices);
            % if main solution is chosen
            else
                valForXIndices = handles.data.plot.res(chosenSessionInd).mainVal(curValuesIndices);
            end
            
            
            
            %valForXIndices = handles.data.plot.res(chosenSessionInd).mainVal(curValuesIndices);
            valForXIndices = valForXIndices(xIndices);

            % Get MJD:
            mjdForXIndices = handles.data.plot.res(chosenSessionInd).mjd(curValuesIndices);
            mjdForXIndices = mjdForXIndices(xIndices);

            % Get y Indices:
            yIndices = zeros(length(xIndices), 1);

            for index = 1 : length(AmbigValues)
               temopYIndices = (abs(valForXIndices) == abs(AmbigValues(index)));
               yIndices = temopYIndices + yIndices;
            end

            yIndices = find(yIndices);
            ambigEpochs = mjdForXIndices(yIndices);

            %outlierEpochs=handles.data.plot.res(chosenSessionInd).mjd(curValuesIndices(...
            %    get(handles.data.plot.outlierMarksHandle, 'XData')));
            
            % Get Baselines:
            baselineIndices = handles.data.plot.res(chosenSessionInd).baselineOfObs(curValuesIndices,:);
            baselineForXIndices = baselineIndices(xIndices, :);
            ambigBaselineInd = baselineForXIndices(yIndices, :);
            
            % Get Station names:
            allStatNames=handles.data.plot.res(chosenSessionInd).allStatNames;
            
            % Get Source names:
            sourceIndices = handles.data.plot.res(chosenSessionInd).source(curValuesIndices,:);
            sourceForXIndices = sourceIndices(xIndices, :);
            ambigSourceInd = sourceForXIndices(yIndices, :);
            allSourceNames=handles.data.plot.res(chosenSessionInd).allSourceNames; 
            
%             outlierBaselineInd=handles.data.plot.res(chosenSessionInd).baselineOfObs(curValuesIndices(...
%                 get(handles.data.plot.outlierMarksHandle, 'XData'),:),:);
            
            % Assigne Station names to Station indices:
            ambigBaselines=cell(size(ambigBaselineInd,1),3); 
            ambigBaselines(:,1)=allStatNames(ambigBaselineInd(:,1));
            ambigBaselines(:,2)=allStatNames(ambigBaselineInd(:,2));
            ambigBaselines(:,3)=allSourceNames(ambigSourceInd); 
            
            valForXIndices = valForXIndices(ismember(baselineForXIndices,ambigBaselineInd,"rows"));

    % is user sure to write ambigs?
    ambspace = zeros(nSelAmbigs,1);
    for iAmb=1:size(ambigEpochs,1)
        obs_at_epoch = [scan([scan.mjd] == ambigEpochs(iAmb)).obs];
        obs_at_epoch_at_baseline = obs_at_epoch([obs_at_epoch.i1]==ambigBaselineInd(iAmb,1) & [obs_at_epoch.i2]==ambigBaselineInd(iAmb,2));
        ambspace(iAmb) = obs_at_epoch_at_baseline.ambspace;
    end
    if max(ambspace) ~= min(ambspace)
       msgbox('Selected observations do not all have the same ambiguity spacing! Select a different set of observations.','Error','error');               
       return
    else
        ambspace = ambspace(1)*1e9;
    end
    opts.Interpreter = 'tex';
    prompt = ['Ambiguity spacing is ', num2str(ambspace), ' ns (\approx', num2str(ambspace*30), ' cm). Enter size of multiplicator:'];
    multiplicator = inputdlg({prompt},'Muliplicator',[1 50],{''},opts);
%     choice = questdlg(sprintf('Write out ambiguities for %1.0f selected observations?', ...
%         nSelAmbigs), 'Sure?', ...
%         'Yes','No','Yes');
    % Handle response
%     switch choice
%             case 'Yes'
    if isempty(multiplicator)
        return
    end
            multiplicator = str2double(multiplicator);
            if isnan(multiplicator)
                msgbox('Invalid input','Error','error');               
                return             
            end
            
%% ##### WRITE DATA TO AMB FILE #####
            
            % append data or create new (depending if file exists)
            fid=fopen([AMBfolder, AMBfilename], 'a');    % 'a' is OK for both append or create new for writing
            
            for iAmb=1:size(ambigEpochs,1)
                fprintf(fid, '%8s %8s %18.12f %8s %g\n', ambigBaselines{iAmb,1},...
                    ambigBaselines{iAmb,2}, ambigEpochs(iAmb), ambigBaselines{iAmb,3}, ambspace*multiplicator);
                
                % TEST
                fprintf(1, '%8s %8s %18.12f %8s %g\n', ambigBaselines{iAmb,1},...
                    ambigBaselines{iAmb,2}, ambigEpochs(iAmb), ambigBaselines{iAmb,3}, ambspace*multiplicator);
            end
            fclose(fid);
            
            % sucessful msgbox
            msgbox('Ambiguities(s) sucessfully written to AMB file', 'Done', 'help');               
    end
end
