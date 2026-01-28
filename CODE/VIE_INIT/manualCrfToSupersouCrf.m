% This function reads source coordintes from a textfile and puts it values
% to a supersource CRF file.



function crf=manualCrfToSupersouCrf(crf,filename)

% open file 
fid=fopen(filename, 'r');
    % read manualCrfData
    manualCrfData=textscan(fid, '%21c  %11c  %f %f %f %03s %f %f %f %f %f %f %f %f %f %f %8c', ...
    'headerlines', 10, 'delimiter', '\n', 'CommentStyle', '#');
%      manualCrfData=textscan(fid, '%21c  %11c  %f %f %f %03s %f %f %f %f %f %f %f %f %f %f %f', ...
%     'headerlines', 22, 'delimiter', '\n'); % icrf3
fclose(fid);


    
% for all lines
for iLine=1:size(manualCrfData{1},1)
    % try to find corresponding souion in supersou TRF
     curSouInSupersouLog=strcmp({crf.IVSname}, manualCrfData{17}(iLine,1:8));
     %curSouInSupersouLog=strcmp({crf.IERSname}, manualCrfData{2}(iLine,1:8)); %icrf3
    
    if sum(curSouInSupersouLog)>1
        fprintf('ERROR: Source %s in manual CRF found in more than one source in supersource file!\n', ...
        manualCrfData{17}(iLine,1:8));
        souIndInSupersou=find(curSouInSupersouLog);
    elseif sum(curSouInSupersouLog)==0
        % write the souion to new entry in supersou file
        fprintf('ERROR: IVS Source %s in manual CRF was not found in supersource file! Try IERS name', ...
            manualCrfData{17}(iLine,1:8));
        % TRY IERS name
        curSouInSupersouLog=strcmp({crf.IERSname}, manualCrfData{2}(iLine,1:8));
        souIndInSupersou=find(curSouInSupersouLog);
    else
        souIndInSupersou=find(curSouInSupersouLog);
    end
    
    curInd=souIndInSupersou(1);

    RA = deg2rad((manualCrfData{3}(iLine) + manualCrfData{4}(iLine)/60 + manualCrfData{5}(iLine)/3600)*15); %[rad]
    if strncmp('-', manualCrfData{6}(iLine),1)
        DE = deg2rad((str2double(manualCrfData{6}(iLine)) - manualCrfData{7}(iLine)/60 - manualCrfData{8}(iLine)/3600)); %[rad]
    else
        DE = deg2rad((str2double(manualCrfData{6}(iLine)) + manualCrfData{7}(iLine)/60 + manualCrfData{8}(iLine)/3600)); %[rad]
    end

        
    if strcmp('D', manualCrfData{2}(iLine,end))
        crf(curInd).manualCrf.defining = 1;
    else
        crf(curInd).manualCrf.defining = 0;
    end
    crf(curInd).manualCrf.ra = RA; %[rad]
    crf(curInd).manualCrf.de = DE; %[rad]
    crf(curInd).manualCrf.ra_sigma = deg2rad(manualCrfData{9}(iLine)*15/3600); %[rad]
    crf(curInd).manualCrf.de_sigma = deg2rad(manualCrfData{10}(iLine)/3600); %[rad]
    crf(curInd).manualCrf.corr = manualCrfData{11}(iLine); 
    crf(curInd).manualCrf.meanMjd = manualCrfData{12}(iLine); 
    crf(curInd).manualCrf.firstMjd = manualCrfData{13}(iLine); 
    crf(curInd).manualCrf.lastMjd = manualCrfData{14}(iLine); 
    crf(curInd).manualCrf.numberSess = manualCrfData{15}(iLine); 
    crf(curInd).manualCrf.numberObs = manualCrfData{16}(iLine); 

end