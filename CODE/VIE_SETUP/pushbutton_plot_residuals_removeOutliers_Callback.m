% #########################################################################
% #     pushbutton_plot_residuals_removeOutliers_Callback
% #########################################################################
%
% DESCRITPION
% Executes on button press in pushbutton_plot_residuals_removeOutliers.
%
% AUTHOR 
%   Matthias Madzak
%
% INPUT
% hObject    handle to pushbutton_plot_residuals_removeOutliers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% OUTPUT
%
% CHANGES
%   2014-06-06, A. Hellerschmied: Adoption to use "hours since session
%       start" as time scale in the plotting residuals window.
%   2015-03-03, A. Hellerschmied: Bug fix - If Outlier folder does not
%       exist, it is created automatically.


function pushbutton_plot_residuals_removeOutliers_Callback(hObject, eventdata, handles)

% % for now - write message box
% msgbox(sprintf('Not working yet!, But this button should remove %1.0f outliers', ...
%     size(get(handles.data.plot.outlierMarksHandle, 'Ydata'), 2)), ':-)', 'warn');

% nSelOutliers=size(get(handles.data.plot.outlierMarksHandle, 'Ydata'), 1);
if size(handles.data.plot.outlierMarksHandle,1)==1
    nSelOutliers=1;
else
    nSelOutliers=size(handles.data.plot.outlierMarksHandle(1).XData,2);
end
% if no outlier was selected -> write messge box
if nSelOutliers==0
    msgbox('No value was selected!', 'No value selected', 'warn');
else
    % (1) outlier subfolder
    allPopupmenuEntriesOUTLIERfolder=get(handles.popupmenu_setInput_outDir, 'String');
    if strcmp(get(handles.popupmenu_setInput_outDir, 'String'), ' ')
        OUTLIERsubFolder='';
    else
        OUTLIERsubFolder=allPopupmenuEntriesOUTLIERfolder{get(handles.popupmenu_setInput_outDir, 'Value')};
    end
    % is user sure to remove outliers?
    choice = questdlg(sprintf('Definitely remove the %1.0f selected outliers?\nSelected (in File->Set Input Files) outlier subfolder: ''%s''', ...
        nSelOutliers, OUTLIERsubFolder), 'Sure?', ...
        'Yes','No','Yes');
    % Handle response
    switch choice
        case 'Yes'
            % get the OUTLIER filename
            % (2) year
            % firstPlottedMjd=handles.data.plot.res(get(handles.popupmenu_plot_residuals_session, 'Value')).mjd(1);
            % get the year from session name - changed 7.12. S. Boehm
            
            % (3) session name
            chosenSessionInd=get(handles.popupmenu_plot_residuals_session, 'Value');
            allPopupmenuEntriesSessions=get(handles.popupmenu_plot_residuals_session, 'String');
            session = allPopupmenuEntriesSessions{chosenSessionInd};
            
            % ##### Check, if the input dataset file followed the standard naming convention #####
            % => YYMMMDDcc_Nnnn (c...alphabetic character, n...number)
            flag_std_naming_convention = true;

            % Total length = 14 char
            if length(session) ~= 14
                flag_std_naming_convention = false; 
            end

            if flag_std_naming_convention
                % "_" at car. 10
                if ~strcmp(session(10), '_')
                    flag_std_naming_convention = false;
                end
                % first two characters are numbers:
                [~, status_1] = str2num(session(1:2));
                % char. 6+7 are numbers:
                [~, status_2] = str2num(session(6:7));
                % char. 12-14 are numbers:
                [~, status_3] = str2num(session(12:14));
                if ~(status_1 && status_2 && status_3)
                    flag_std_naming_convention = false;
                end
            end

            if flag_std_naming_convention
                % ##### Standard naming convention is used for this session: #####

                % Get the year from the session name:
                if str2double(session(1:2)) > 75
                    yearStr = ['19', session(1:2)];
                else
                    yearStr = ['20', session(1:2)];
                end

                
            else
                % ##### Non-Standard naming convention is used for this session: #####
                % e.g. when using .vso input files
                % => Get yearStr from the opt_ file!

                % #### Load opt_ file from LEVEL 3 (sub-)directory: ####
                % Get sub-dir.:
                allSubfolders   = get(handles.popupmenu_plot_residuals_folder, 'string');
                curSubfolder    = allSubfolders{get(handles.popupmenu_plot_residuals_folder, 'Value')};

                % load opt_ file and get the year:
                load(['../DATA/LEVEL3/', curSubfolder, '/opt_', session, '.mat']);
                yearStr = opt_.year;

            end
            
            OUTfolder = ['../DATA/OUTLIER/', OUTLIERsubFolder, '/', yearStr, '/'];
            OUTfilename = [allPopupmenuEntriesSessions{chosenSessionInd}, '.OUT'];

            % Check if outlier folder exists => if not => create it!
            if ~exist(OUTfolder, 'dir')
                mkdir(OUTfolder);
            end
            
       % ##### get the values just plotted (station, baseline, all-wise) #####
            % ### per station ###
            if get(handles.radiobutton_plot_residuals_perStat, 'Value')
                curValues=sum(handles.data.plot.res(chosenSessionInd).baselineOfObs==get(handles.popupmenu_plot_residuals_station, 'Value'),2);
                
            % ### per baseline ###
            elseif get(handles.radiobutton_plot_residuals_perBasel, 'Value')
                allBaselines=get(handles.popupmenu_plot_residuals_baseline, 'String');
                stat1=allBaselines{get(handles.popupmenu_plot_residuals_baseline,'Value')}(1:8);
                stat2=allBaselines{get(handles.popupmenu_plot_residuals_baseline,'Value')}(10:17);
                statNr1=~cellfun(@isempty, strfind(handles.data.plot.res(chosenSessionInd).allStatNames, stat1));
                statNr2=~cellfun(@isempty, strfind(handles.data.plot.res(chosenSessionInd).allStatNames, stat2));
                valsLogicalsOfFirstStat=handles.data.plot.res(chosenSessionInd).baselineOfObs==find(statNr1);
                valsLogicalsOfSecondStat=handles.data.plot.res(chosenSessionInd).baselineOfObs==find(statNr2);
                curValues=sum(valsLogicalsOfFirstStat+valsLogicalsOfSecondStat,2)==2;
                
            % ### per source ###
            elseif get(handles.radiobutton_plot_residuals_perSource, 'Value')
                curSource=get(handles.popupmenu_plot_residuals_source, 'Value');
                curValues=handles.data.plot.res(chosenSessionInd).source==curSource;
            
            % ### all are plotted ###
            else
                curValues=ones(size(handles.data.plot.res(chosenSessionInd).mjd,1),1);
            end
                
%% ############# GET OUTLIERS ##############
            
            % ##### get epochs (mjd) for outliers #####
            
            % Currently plotted values (indices)
            curValuesIndices = find(curValues);
            
            % Selected Data (Outliers) in residuals plot window:
            outlierEpochsHours = get(handles.data.plot.outlierMarksHandle, 'XData');
%             if(length(outlierEpochsHours)>1)
%                 outlierEpochsHours = outlierEpochsHours{1};
%             end
            if size(outlierEpochsHours,1)>1
                outlierEpochsHours = outlierEpochsHours{1};
            else
                outlierEpochsHours = outlierEpochsHours(1);
            end
                
            
            OutlierValues = get(handles.data.plot.outlierMarksHandle, 'YData');
%             if(length(OutlierValues)>1)
%                 OutlierValues = OutlierValues{1};
%             end
            if(size(OutlierValues,1)>1)
                OutlierValues = OutlierValues{1};
            else
                OutlierValues = OutlierValues(1);
            end

            % Conversion of x-values: "hours from session start" => "MJD": 
            SessionStartTimeMJD =  handles.data.plot.res(chosenSessionInd).mjd(1);
            selectedOutlierEpochsMJD = SessionStartTimeMJD + outlierEpochsHours / 24;

            % Get Baseline indices for selected outliers with a 2-dimensional search
            % approach:
            
            % 1.) ### Compare x-values (MJD): ###
            
            % Delete multiple entries in vector "selectedOutlierEpochsMJD":
            outlierEpochsMJD = unique(selectedOutlierEpochsMJD);
            numOfoutlierEpochsMJD = length(outlierEpochsMJD);
            
            xIndices = zeros(nSelOutliers, 1);
            index = 1;
            
            for indexEpochs = 1 : numOfoutlierEpochsMJD
                indicesCurrentEpoch = find(handles.data.plot.res(chosenSessionInd).mjd(curValuesIndices) == outlierEpochsMJD(indexEpochs));
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

            for index = 1 : length(OutlierValues)
               temopYIndices = (abs(valForXIndices) == abs(OutlierValues(index)));
               yIndices = temopYIndices + yIndices;
            end

            yIndices = find(yIndices);
            outlierEpochs = mjdForXIndices(yIndices);

            %outlierEpochs=handles.data.plot.res(chosenSessionInd).mjd(curValuesIndices(...
            %    get(handles.data.plot.outlierMarksHandle, 'XData')));
            
            % Get Baselines:
            baselineIndices = handles.data.plot.res(chosenSessionInd).baselineOfObs(curValuesIndices,:);
            baselineForXIndices = baselineIndices(xIndices, :);
            outlierBaselineInd = baselineForXIndices(yIndices, :);
            
            % Get Station names:
            allStatNames=handles.data.plot.res(chosenSessionInd).allStatNames;
            
            % Get Source names:
            sourceIndices = handles.data.plot.res(chosenSessionInd).source(curValuesIndices,:);
            sourceForXIndices = sourceIndices(xIndices, :);
            outlierSourceInd = sourceForXIndices(yIndices, :);
            allSourceNames=handles.data.plot.res(chosenSessionInd).allSourceNames; 
            
%             outlierBaselineInd=handles.data.plot.res(chosenSessionInd).baselineOfObs(curValuesIndices(...
%                 get(handles.data.plot.outlierMarksHandle, 'XData'),:),:);
            
            % Assigne Station names to Station indices:
            outlierBaselines=cell(size(outlierBaselineInd,1),3); 
            outlierBaselines(:,1)=allStatNames(outlierBaselineInd(:,1));
            outlierBaselines(:,2)=allStatNames(outlierBaselineInd(:,2));
            outlierBaselines(:,3)=allSourceNames(outlierSourceInd); 
           
            
%% ##### WRITE DATA TO OUTLIER FILE #####
            
            % append data or create new (depending if file exists)
            fid=fopen([OUTfolder, OUTfilename], 'a');    % 'a' is OK for both append or create new for writing
            
            for iOutlier=1:size(outlierEpochs,1)
                fprintf(fid, '%8s %8s %18.12f %8s\n', outlierBaselines{iOutlier,1},...
                    outlierBaselines{iOutlier,2}, outlierEpochs(iOutlier), outlierBaselines{iOutlier,3});
                
                % TEST
                fprintf(1, '%8s %8s %18.12f %8s\n', outlierBaselines{iOutlier,1},...
                    outlierBaselines{iOutlier,2}, outlierEpochs(iOutlier), outlierBaselines{iOutlier,3});
            end
            fclose(fid);
            
            % sucessful msgbox
            msgbox('Outlier(s) sucessfully written to OUTLIER file', 'Done', 'help');               
    end
end
