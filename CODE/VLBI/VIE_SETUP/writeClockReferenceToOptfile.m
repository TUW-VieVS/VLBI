% 
% DESCRIPTION
% This function writes a (new) reference clock to the OPT file
%
%
% COUPLING
% - handles2OptfilePath.m
%
% CHANGES
% 2014-09-03, H. Krasna: bug corrected; a new OPT file could not be created when the year sub-directory did not exist
% 2016-10-11, A. Hellerschmied: Use function handles2OptfilePath.m to get the OPT filename and folder

function writeClockReferenceToOptfile(handles,newRefclockAntenna)

% get selected file
[OPTfolder, OPTFileName, selectedOPTdir] = handles2OptfilePath(handles);
wantedOPTfile   = [OPTfolder, OPTFileName];
    
choice = questdlg(sprintf('Insert ''%s'' as reference clock?\nOPT subfolder is: ''%s''', newRefclockAntenna, selectedOPTdir), 'Insert clock reference?', 'Yes', 'No', 'Yes');

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
        while ~feof(fid);
            curStr=fgetl(fid);
            optContent{lineInd}=curStr;
            lineInd=lineInd+1;
        end        
        fclose(fid);
        
        indRefClockStat=[];
        for iLine=1:lineInd-1
            if strcmpi('CLOCK REFERENCE:', optContent{iLine}(1:min(length(optContent{iLine}), 16)) )
                indRefClockStat=iLine+1;
                break;
            end
        end
        
        % if there is no clock reference yet - add it at beginning of file
        if isempty(indRefClockStat)
            optContent(3:length(optContent)+2)=optContent(1:end);
            optContent{1}='CLOCK REFERENCE:';
            optContent{2}=newRefclockAntenna;
        else
            optContent{indRefClockStat}='';
            newRefclockAntenna=[newRefclockAntenna, '          '];
            optContent{indRefClockStat}=newRefclockAntenna(1:8);
        end
        
        % write the new content to file
        fid=fopen(wantedOPTfile, 'w');
        for iLine=1:length(optContent)
            fprintf(fid, '%s\n', optContent{iLine});
        end
        fclose(fid);
        
        % Status Msg. to CW:
        fprintf(' - Reference clock updated to %s in OPT file: %s\n', newRefclockAntenna, wantedOPTfile);
end
