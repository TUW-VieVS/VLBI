% #########################################################################
% #     writeClockbreakToOptFile
% #########################################################################
%
% DESCRITPION
% This function writes the user-selected clock-break to the OPT file. If no
% OPT file exists for this session, one with default values (no exclusions)
% is created.
% Works with vievs2_0 interface.
%
% AUTHOR 
%   Matthias Madzak
%
% INPUT
%   handles     structure from the GUI (also containing data, eg residuals)
%   x           x-value which was selected by the user to be a clock break
%
% OUTPUT
%   handles
%
% COUPLING
% - handles2OptfilePath.m
%
% CHANGES
%   2014-06-06, A. Hellerschmied: Adoption to use "hours since session
%       start" as time scale in the plotting residuals window.
%   2014-09-03, H. Krasna: new OPT file is created by calling the
%       writeNewOptFile function
%   2016-09-10, A. Girdiuk: small bug-fix: written OPT-file is not shorter than original
%   2016-10-11, A. Hellerschmied: Use function handles2OptfilePath.m to get the OPT filename and folder
%

function writeClockbreakToOptFile(handles,x)


% ##### Find MJD for picked Clock Break #####

% get values which actually include the chosen station
allStations         = get(handles.popupmenu_plot_residuals_station, 'String');
chosenStation       = allStations{get(handles.popupmenu_plot_residuals_station, 'Value')};

curSessionInd       = get(handles.popupmenu_plot_residuals_session, 'Value');
allBaselInd         = handles.data.plot.res(curSessionInd).baselineOfObs;
allBasel            = handles.data.plot.res(curSessionInd).allStatNames(allBaselInd);
resValsOfCurStatLog = sum(~cellfun(@isempty, strfind(allBasel, chosenStation)),2);

% get all mjds that are actually plotted
allMjd      = handles.data.plot.res(get(handles.popupmenu_plot_residuals_session, 'Value')).mjd;
plottedMjd  = allMjd(logical(resValsOfCurStatLog>=1));

% Calculate MJD from picked value
SessionStartTimeMJD =  handles.data.plot.res(curSessionInd).mjd(1);
mjdOfBreak          = SessionStartTimeMJD + (x / 24);

% Find MJD before and after the picked MJD value
mjdBeforeBreak  = plottedMjd(plottedMjd < mjdOfBreak);
mjdBeforeBreak  = max(mjdBeforeBreak);
mjdAfterBreak   = plottedMjd(plottedMjd > mjdOfBreak);
mjdAfterBreak   = min(mjdAfterBreak);

% Calc. clock break in mid of prev and following mjd
mjdOfBreak = (mjdBeforeBreak + mjdAfterBreak) / 2;

% If MJD before Session-start or after session-end was selected => Error
% and quit funtion!
if isempty(mjdOfBreak)
    errordlg('The selected Clockbreak is not within the valid time period! Clockbreaks has to be within first and last residual value.');
    return; 
end


% ##### Write Clock Break To Opt File #####

% % get OPT file
[OPTfolder, OPTFileName, selectedOPTdir] = handles2OptfilePath(handles);
wantedOPTfile   = [OPTfolder, OPTFileName];

% is user sure to write clock break to OPT file?
choice = questdlg(sprintf('Insert Clock break for station %s at mjd = %1.4f\nOPT subfolder is: ''%s''', chosenStation, mjdOfBreak, selectedOPTdir), 'Insert Break?', 'Yes', 'No', 'Yes');

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

        % read in all data from OPT file
        OPTcontent  = cell(300,1);
        iLine       = 1;
        fid         = fopen(wantedOPTfile);

        while ~feof(fid)
            OPTcontent{iLine} = fgetl(fid);
            iLine = iLine + 1;
        end
        fclose(fid);

        % find the string "CLOCK BREAKS" in OPT file
        lineOfClockBreak=find(~cellfun(@isempty, strfind(OPTcontent(1:iLine-1), 'CLOCK BREAKS: ')));
        % if "CLOCK BREAKS" does not exist, write it to first line
        if isempty(lineOfClockBreak)
            % put everything two lines "down", i.e. to make place
            OPTcontent(2:end)=OPTcontent(1:end-1);
            OPTcontent{1}='CLOCK BREAKS: 1';
            line2WriteBreak=2;
        else
            lineOfClockBreak=lineOfClockBreak(1); % if there are more (eg in the comments!)
            curNumOfBreaks=str2double(OPTcontent{lineOfClockBreak}(15:end));

            % update number of breaks
            OPTcontent{lineOfClockBreak}=[OPTcontent{lineOfClockBreak}(1:14), sprintf('%1.0f', curNumOfBreaks+1)];

            % if there are no breaks yet, put the new break to
            % following line
            if curNumOfBreaks==0
                line2WriteBreak=lineOfClockBreak+1;
            else
                % read current breaks
                existingBreakMjds=zeros(curNumOfBreaks,1);
                for iBreak=1:curNumOfBreaks
                    existingBreakMjds(iBreak)=str2num(OPTcontent{lineOfClockBreak+iBreak}(11:end));
                end

                % find where to put that line
                firstMjdOfBreaksGreaterThanNew=find(existingBreakMjds>mjdOfBreak);
                if isempty(firstMjdOfBreaksGreaterThanNew)
                    line2WriteBreak=lineOfClockBreak+1+size(existingBreakMjds,1);
                else
                    % only take first if more are found
                    firstMjdOfBreaksGreaterThanNew=firstMjdOfBreaksGreaterThanNew(1);

                    % write the new line before the break(s) that
                    % are larger than current one
                    line2WriteBreak=lineOfClockBreak+...
                        firstMjdOfBreaksGreaterThanNew;
                end
            end
        end

        % "shift" down all OPT-files for one 
        OPTcontent(line2WriteBreak+1:end+1)=OPTcontent(line2WriteBreak:end);

        % put new line to "newly created" line
        OPTcontent{line2WriteBreak}=sprintf('%8s  %12.6f', chosenStation, mjdOfBreak);
        
        % i want to write lines that are not empty OR numeric (& is
        % correct)
        nLinesToBeWritten=sum(~cellfun(@isnumeric,OPTcontent));
        % save everything in OPT file
        fid=fopen(wantedOPTfile, 'w');
        for iLine=1:nLinesToBeWritten
            fprintf(fid, '%s\n', OPTcontent{iLine});
        end
        fclose(fid);
        
        % Status Msg. to CW:
        fprintf(' - Clock break added to OPT file: %s\n', wantedOPTfile);

end