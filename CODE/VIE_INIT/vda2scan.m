function [antenna,sources,scan] = vda2scan(pfilepath, psession_name, trffile, crffile, trf, crf)

obs_file_name = [pfilepath psession_name];

band = '1'; %2

% ##### Open file #####
flagBZ2 = 0;
if ispc
    if exist([obs_file_name '.vda.bz2'],'file') == 2
    
        prog7z = 'C:\Program Files\7-Zip\7z.exe';

        [status7z,result] = system(['"' prog7z '" -y x ' '"' obs_file_name '.vda.bz2"' ' -o' '"' pfilepath '"']);

        if status7z ==1
            fprintf('ERROR: Extract contents of .bz2 file: %s %s\n',result, prog7z);
        elseif status7z == 2
            fprintf('%s\n',result);
        end    
    
        flagBZ2 = 1;
    end

elseif isunix
        if exist([obs_file_name '.vda.bz2'],'file') == 2
    
        [statusbunzip2,result] = system("bunzip2 -dkf "+obs_file_name+".vda.bz2");

        if statusbunzip2 == 1
            fprintf('ERROR: Extract contents of .bz2 file: %s %s\n',result);
        elseif statusbunzip2 == 2
            fprintf('%s\n',result);
        end    
    
        flagBZ2 = 1;
        end
end

fid_vda = fopen([obs_file_name '.vda'],'r');

% Loop over all lines
if fid_vda ~= -1
    wholeFile = textscan(fid_vda,'%s','delimiter', '\n', 'whitespace', '');
    wholeFile = wholeFile{1};
else
    error(' couln''t find VDA-file %s',obs_file_name)
end
fclose(fid_vda);


ind=contains(wholeFile,'DATA.1');
DATA1=wholeFile(ind);
nlinesDATA1 = length(DATA1);

ind=contains(wholeFile,'DATA.3');
DATA3=wholeFile(ind);
nlinesDATA3 = length(DATA3);

ind=contains(DATA1,'DATA.1 NUMB_OBS');
currline=DATA1(ind);
seg = split(currline);
nObs=str2double(seg{7});

ind=contains(DATA1,'DATA.1 NUMB_SCA');
currline=DATA1(ind);
seg = split(currline);
nScans=str2double(seg{7});

ind=contains(DATA1,'DATA.1 NUMB_STA');
currline=DATA1(ind);
seg = split(currline);
nSta=str2double(seg{7});

ind=contains(DATA1,'DATA.1 NUM_BAND');
currline=DATA1(ind);
seg = split(currline);
nBand=str2double(seg{7});

valtab=zeros(nObs,14);
 
ind=contains(DATA1,'DATA.1 OBS_TAB');
ourblock=DATA1(ind);
iv=1;
for i=1:size(ourblock,1)
    input_str = ourblock{i};
    seg = split(input_str);
    if char(seg(5)) == '1'
        iscan = str2double(seg(7));
        iobs = str2double(seg(6));
        valtab(iv,1) = iscan;
        valtab(iv,2) = iobs;
    elseif char(seg(5)) == '2'
        istat1 = str2double(seg(7));
        valtab(iv,3) = istat1;
    elseif char(seg(5)) == '3'
        istat2 = str2double(seg(7));
        valtab(iv,4) = istat2;
        iv=iv+1;
    end

end


ind=contains(DATA1,'DATA.1 GR_DELAY');
ourblock=DATA1(ind);
ourblockcell=split(ourblock);
ourblockcell(:,7)=replace(ourblockcell(:,7),'D','e');
bnds=cell2mat(ourblockcell(:,5));
valtab(:,5) = str2double(ourblockcell(bnds==band,7));



ind=contains(DATA1,'DATA.1 GRDELERR');
ourblock=DATA1(ind);
ourblockcell=split(ourblock);
ourblockcell(:,7)=replace(ourblockcell(:,7),'D','e');
bnds=cell2mat(ourblockcell(:,5));
valtab(:,6) = str2double(ourblockcell(bnds==band,7));




valtab(:,7)=0; % computed delay



ind=contains(DATA1,'DATA.1 ION_GDEL');
ourblock=DATA1(ind);
ourblockcell=split(ourblock);
ourblockcell(:,7)=replace(ourblockcell(:,7),'D','e');
bnds=cell2mat(ourblockcell(:,5));
valtab(:,8) = str2double(ourblockcell(bnds==band,7))*10^9; %ns

ind=contains(DATA1,'DATA.1 ION_GERR');
ourblock=DATA1(ind);
ourblockcell=split(ourblock);
ourblockcell(:,7)=replace(ourblockcell(:,7),'D','e');
bnds=cell2mat(ourblockcell(:,5));
valtab(:,9) = str2double(ourblockcell(bnds==band,7))*10^9; %ns

valtab(:,10:11)=zeros(nObs,2); %'q_flag', 'q_flag_ion'


ind=contains(DATA1,'DATA.1 QUALCODE');
ourblock=DATA1(ind);
ourblockcell=split(ourblock);

bnds=cell2mat(ourblockcell(:,6));
valtab(:,12) = str2double(ourblockcell(bnds=='1',7));%band 1 :X
if ~isempty(find(bnds=='2', 1))
    valtab(:,13) = str2double(ourblockcell(bnds=='2',7));%band 2 :S
else
    valtab(:,13)=0;
end


valtab(:,14)=0; % ambiq


ind=contains(DATA1,'DATA.1 MJD_OBS');
ourblock=DATA1(ind);
for i=1:size(ourblock,1)
    input_str = ourblock{i};
    seg = split(input_str);
    id = ismember(valtab(:,1),str2double(seg(3)));
    valtab(id,15)=str2double(seg(7)); %days
    
end

ind=contains(DATA1,'DATA.1 UTC_OBS');
ourblock=DATA1(ind);
for i=1:size(ourblock,1)
    input_str = ourblock{i};
    seg = split(input_str);
    
        id = ismember(valtab(:,1),str2double(seg(3)));
        D=find(char(seg(7))=='D');
        seg{7}(1,D)='e';
        valtab(id,16)=str2double(seg(7)); %s
    
end

ind=contains(DATA1,'DATA.1 SOU_IND');
ourblock=DATA1(ind);
for i=1:size(ourblock,1)
    input_str = ourblock{i};
    seg = split(input_str);
    
        id = ismember(valtab(:,1),str2double(seg(3)));
        valtab(id,17)=str2double(seg(7));    
end

ind=contains(DATA1,'DATA.1 SRCNAMES');
ourblock=DATA1(ind);
for i=1:size(ourblock,1)
    input_str = ourblock{i};
    seg = split(input_str);
    
    src(i).id=str2double(seg(6));
    src(i).name=(seg(7));
end


valtab(:,18)=valtab(:,15)+valtab(:,16)./86400;
[ti(:,1), ti(:,2), ti(:,3), ti(:,4), ti(:,5), ti(:,6)] = mjd2date(valtab(:,18));
[~, ti(:,7)]=dday(ti(:,1),ti(:,2),ti(:,3),ti(:,4),ti(:,5));




ind=contains(wholeFile,'DATA.4');
DATA4=wholeFile(ind);

ind=contains(DATA4,'DATA.4 EFF_FREQ');
ourblock=DATA4(ind);
ourblockcell=split(ourblock);
ourblockcell(:,7)=replace(ourblockcell(:,7),'D','e');
bnds=cell2mat(ourblockcell(:,6));
obser=cell2mat(ourblockcell(:,5)); % 1 - gr.del, 2 - ph.del, 3 - ph.rate
valtab(:,19)=str2double(ourblockcell(bnds==band & obser=='1' ,7)); %Hz



[anttab] = vdaDATA3_2scan(DATA3);


scantab=zeros(nScans,nSta);
for i=1:nScans
    isc= find(valtab(:,1)==i);
    indstat = unique([valtab(isc,3); valtab(isc,4)]);
    scantab(i,indstat)=1;
    
end

% % control
% sum(scantab)


scanantT=zeros(nScans, nSta);
scanantP=zeros(nScans, nSta);
scanantE=zeros(nScans, nSta);
scanantC=zeros(nScans, nSta);


for i=1:nSta
    id=[]; isc=[];
    id=find(anttab(:,2)==i);
    isc=find(scantab(:,i)==1);
    scanantT(isc,i)=anttab(id,3);
    scanantP(isc,i)=anttab(id,4);
    scanantE(isc,i)=anttab(id,6);
    scanantC(isc,i)=anttab(id,10);
end







%% PREALLOCATING
scan(nScans+1)=struct('mjd', [], 'stat', [], 'tim', [], ...
    'nobs', [], 'space', [], 'obs', [], 'iso', []); % +1: not working otherwise - is deleted after loop
space0.source = zeros(3,3);
space0.xp=0; space0.yp=0; space0.era=0; space0.xnut=0; space0.ynut=0;
space0.t2c=zeros(3,3);


subStruct_stat=struct('x', [], 'temp', [], 'pres', [], 'e', [], 'az', ...
    [], 'zd', [], 'zdry', [], 'zwet', [], 'cab', [], 'axkt', [], 'therm', [], ...
    'pantd', [], 'trop', []);

scan(nScans+1).obs=struct('i1', [], 'i2', [], 'obs', [], 'sig', [], 'com', ...
    [], 'delion', [], 'sgdion', [], 'q_flag', [], 'q_flag_ion', [], 'q_code_X', [], 'q_code_S', [], 'eff_freq', []);



% FLAGS
% if parameter.obs_restrictions.suppression_flags
%     'vda2scan.m line 350 - flags from psolve applied: SHOULD BE MOVED TO VIE_LSM'
%     flags_file_name = replace(obs_file_name,'VDA','FLAGS_PSOLVE');
%     fid_flags = fopen([flags_file_name '_edit.txt'],'r');
% 
%     % Loop over all lines
%     if fid_flags ~= -1
%         FLG_PSOLVE = textscan(fid_flags,'%f','delimiter', '\n', 'whitespace', '', 'CommentStyle','#');
%         FLG_PSOLVE = FLG_PSOLVE{1};
%         valtab(FLG_PSOLVE,10) = 1; % 1 - bad Qflag
%         fclose(fid_flags);
% 
%     else
%         fprintf('couln''t find PSOLVE FLAG-file %s_edit.txt',flags_file_name)
%     end
% end



for i=1:nScans   
    iv=[]; currscan=[];
    iv=ismember(valtab(:,1),i);
    numob = sum(iv);
    currscan=valtab(iv,:);
    currtim=ti(iv,:);
    
    scan(i).nobs=numob;
    scan(i).obs_type='q';
    scan(i).iso=currscan(1,17);
    scan(i).mjd=currscan(1,18); %day
    scan(i).tim=currtim(1,:)'; %tim doy
    
    
    
    % loop over stations
    for j=1:nSta
        subStruct_stat(j).temp=scanantT(i,j) - 273.15; % K -> °C
        subStruct_stat(j).pres=scanantP(i,j);
        subStruct_stat(j).e=scanantE(i,j);
        subStruct_stat(j).cab=scanantC(i,j);

        if subStruct_stat(j).temp==0
            subStruct_stat(j).temp=[];
        end
        if subStruct_stat(j).pres==0
            subStruct_stat(j).pres=[];
        end
        if subStruct_stat(j).e==0
            subStruct_stat(j).e=[];
        end
        if subStruct_stat(j).cab==0
            subStruct_stat(j).cab=[];
        end
        
    end
    scan(i).stat=subStruct_stat;

    
    
    
    for j=1:numob
        t=array2table(currscan(:,[3:14 19]),'VariableNames',{'i1','i2','obs','sig','com','delion','sgdion','q_flag','q_flag_ion','q_code_X','q_code_S','amb','eff_freq'});
        scan(i).obs(j) =  table2struct(t(j,:));
        scan(i).obs(j).q_code_X=num2str(scan(i).obs(j).q_code_X);
        scan(i).obs(j).q_code_S=num2str(scan(i).obs(j).q_code_S);

        

        
        % add corrections to the group delay
        % CABLE
        corcab = scan(i).stat(scan(i).obs(j).i2).cab - scan(i).stat(scan(i).obs(j).i1).cab; % [ns]
        if ~isempty(corcab)
             scan(i).obs(j).obs = scan(i).obs(j).obs + corcab*(1e-9);    %[s]  
        end
        % IONO
        if ~isempty(scan(i).obs(j).delion)
            scan(i).obs(j).obs = scan(i).obs(j).obs - scan(i).obs(j).delion*(1e-9); % [s]
            scan(i).obs(j).sig = sqrt(scan(i).obs(j).sig^2 + (scan(i).obs(j).sgdion*(1e-9))^2);
        end
    end
  
    scan(i).space=space0;
    
end  
scan(end)=[];

%%

antenna=vda2antenna(DATA1,nSta,valtab, trf, trffile);
sources=vda2sources(DATA1,valtab, crf, crffile);

if flagBZ2 == 1
    delete([obs_file_name '.vda'])
end

%% check for clock breaks

% ind=contains(DATA4,'DATA.4 NUM_CLBR');
% currline=DATA4(ind);
% seg = split(currline);
% nCB=str2double(seg{7});
% if nCB>0
%     'CLOCK BREAKS IN THE SESSION!!!'
% 
%     ind=contains(DATA4,'DATA.4 STA_CLBR');
%     ourblock=DATA4(ind);
%     ourblockcell=split(ourblock);
%     if nCB==1
%         ourblockcell{7}
%     else
%         ourblockcell{:,7}
%     end
%     
% %     ind=contains(DATA4,'DATA.4 UTC_CLBR');
%     
% end
% 
% 

