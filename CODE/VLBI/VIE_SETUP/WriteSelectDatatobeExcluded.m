% #########################################################################
% #     WriteSelectDatatobeExcluded
% #########################################################################
%
% DESCRITPION
% This function writes the user-selected data to the OPT file. If no
% OPT file exists for this session, one with default values (no exclusions)
% is created.
%
% AUTHOR 
%       A. Girdiuk
%
% INPUT
%   handles     structure from the GUI (also containing data, eg residuals)
%   DataType    station, source 
%   TimeSpan    time interval to be excluded
%
% OUTPUT
%   writes OPT-file in chosen directory only
%
% COUPLING
% - handles2OptfilePath.m
%
% CHANGES
%   2016-09-21, A. Girdiuk: typos-fix
%   2016-10-11, A. Hellerschmied: Use function handles2OptfilePath.m to get the OPT filename and folder
%   2017-10-11, A. Hellerschmied: More digits used for MJD epochs to remove data
function WriteSelectDatatobeExcluded(handles, DataType, TimeSpan)

% get selected file
[OPTfolder, OPTFileName, selectedOPTdir] = handles2OptfilePath(handles);
wantedOPTfile   = [OPTfolder, OPTFileName];

choice = questdlg(sprintf('Exclude from ''%s'' interval %2.3f - %2.3f \nOPT subfolder is: ''%s''', DataType, TimeSpan, selectedOPTdir), 'Would you like to exclude next data?', 'Yes', 'No', 'Yes');

switch choice
    case 'Yes'
        % check if the directory exists
        if ~exist(OPTfolder, 'dir')
            mkdir(OPTfolder)
        end

        % if OPT file does not exist at all
        if ~exist(wantedOPTfile, 'file')
            writeNewOptFile(wantedOPTfile);
        end
        
        % search "Reference clock in file
        
        % load all content
        fid=fopen(wantedOPTfile, 'r');   
        lineInd=1;
        while ~feof(fid)
            curStr=fgetl(fid);
            optContent{lineInd}=curStr;
            lineInd=lineInd+1;
        end        
        fclose(fid);
        
        if get(handles.radiobutton_plot_residuals_perStat,'Value')
            % station to be excluded
            mask = 'STATIONS TO BE EXCLUDED:';
            indtobeExcluded = find(~cellfun(@isempty,strfind(optContent,mask)));
            indtobeExcluded = indtobeExcluded(1); % second one is supposed to be commented
       
            % number station already exluded:
            numinOPT = str2num(optContent{indtobeExcluded}(25:end));
            optContent(indtobeExcluded+numinOPT+2:length(optContent)+1) = optContent(indtobeExcluded+numinOPT+1:end);
            % put new line to "newly created" line
            optContent{indtobeExcluded}(25:end) = sprintf('%2i',[numinOPT+1]);
            optContent{indtobeExcluded+numinOPT+1}=sprintf('%8s  %12.6f-%12.6f',DataType, TimeSpan);
        end
        
        if get(handles.radiobutton_plot_residuals_perSource,'Value')
            % station to be excluded
            mask = 'SOURCES TO BE EXCLUDED:';
            indtobeExcluded = find(~cellfun(@isempty,strfind(optContent,mask)));
            indtobeExcluded = indtobeExcluded(1); % second one is supposed to be commented
       
            % number station already exluded:
            numinOPT = str2num(optContent{indtobeExcluded}(24:end));
            optContent(indtobeExcluded+numinOPT+2:length(optContent)+1) = optContent(indtobeExcluded+numinOPT+1:end);
            % put new line to "newly created" line
            optContent{indtobeExcluded}(24:end) = sprintf('%2i',[numinOPT+1]);
            optContent{indtobeExcluded+numinOPT+1}=sprintf('%8s  %12.6f-%12.6f',DataType, TimeSpan);
        end
        
        % write the new content to file
        fid=fopen(wantedOPTfile, 'w');
        for iLine=1:length(optContent)
            fprintf(fid, '%s\n', optContent{iLine});
        end
        fclose(fid);
        
        % Status Msg. to CW:
        fprintf(' - Observations excluded via OPT file: %s\n', wantedOPTfile);
end
