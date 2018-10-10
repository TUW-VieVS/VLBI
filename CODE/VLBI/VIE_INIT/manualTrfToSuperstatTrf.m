% This function reads station coordintes from a textfile and puts it values
% to a superstation TRF file.
% created
%   2.8.2012 Matthias Madzak
%   2.2.2016 Younghee Kwak: read additional columns for datum flags and comments 
%                           (11th and 12th columns) like vievsTrf format



function trf=manualTrfToSuperstatTrf(trf,filename)

% open file 
fid=fopen(filename, 'r');

% read data
manualTrfData=textscan(fid, '%8c %f %f %f %f %f %f %f %f %f %f %s', 'commentstyle', '%');
    
fclose(fid);
    
% for all lines
for iLine=1:size(manualTrfData{1},1)
    % try to find corresponding station in superstat TRF
    curStatInSuperstatLog=strcmp({trf.name}, manualTrfData{1}(iLine,1:8)) | ...
        strcmp({trf.name}, strrep(manualTrfData{1}(iLine,1:8), ' ', '_')) | ...
        strcmp({trf.name}, strrep(manualTrfData{1}(iLine,1:8), '_', ' '));
    
    if sum(curStatInSuperstatLog)>1
        fprintf('ERROR: Station %s in manual TRF found in more than one station in superstation file!\n', ...
            manualTrfData{1}(iLine,1:8));
        statIndInSuperstat=find(curStatInSuperstatLog);
    elseif sum(curStatInSuperstatLog)==0
        % write the station to new entry in superstat file
        fprintf('ERROR: Station %s in manual TRF was not found in superstation file!\nStation correction (loading...) information might not be used properly!\n', ...
            manualTrfData{1}(iLine,1:8));
        statIndInSuperstat=size(trf,2)+1;
    else
        statIndInSuperstat=find(curStatInSuperstatLog);
    end
    
    % make a loop (if more stations with that name is found
    for iFoundStat=1:length(statIndInSuperstat)
        curInd=statIndInSuperstat(iFoundStat);
        
        % get break number
        if curInd>size(trf,2) % if we have to add a new field
            breakNr=1;
            trf(curInd).name=manualTrfData{1}(iLine,1:8);
        elseif ~isfield(trf(curInd), 'manualTrf') % for the very first time only!
            breakNr=1;
        elseif ~isfield(trf(curInd).manualTrf, 'break') % if there is no break for "manualTrf" in the current station
                breakNr=1;
        else
            breakNr=size(trf(curInd).manualTrf.break,2)+1;
        end
        
        % write data to superstat trf
        trf(curInd).manualTrf.break(breakNr).x=manualTrfData{2}(iLine);
        trf(curInd).manualTrf.break(breakNr).y=manualTrfData{3}(iLine);
        trf(curInd).manualTrf.break(breakNr).z=manualTrfData{4}(iLine);
        trf(curInd).manualTrf.break(breakNr).vx=manualTrfData{5}(iLine);
        trf(curInd).manualTrf.break(breakNr).vy=manualTrfData{6}(iLine);
        trf(curInd).manualTrf.break(breakNr).vz=manualTrfData{7}(iLine);
        trf(curInd).manualTrf.break(breakNr).epoch=manualTrfData{8}(iLine);
        trf(curInd).manualTrf.break(breakNr).start=manualTrfData{9}(iLine);
        trf(curInd).manualTrf.break(breakNr).end=manualTrfData{10}(iLine);
        trf(curInd).manualTrf.break(breakNr).indatum=manualTrfData{11}(iLine);
   
    end
end