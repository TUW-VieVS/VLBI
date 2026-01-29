%
% Computation of effective frequency within a band
%
% Hana Krasna, 2025 Aug 10

function [eFreqGHz,sigma_eFreqGHz,qflag_v] = effective_freq(out_struct,freqband,parameter)


fileoutput=false;
fileoutput_F=false;

if fileoutput
    if ~isfield(out_struct.Observables, 'ChannelInfo_bS')
        fid = fopen(['sessions_missing_ChannelInfo_bS.txt'],'a');
        fprintf(fid,'%s   %s \n', out_struct.head.Session.val, parameter.session_name);
        fclose(fid);
        %return
    end
end


% reference frequency
fileRefv =    ['RefFreq_' freqband];
if ~isfield(out_struct.Observables,fileRefv)
    fprintf('\nReference frequency for %s is missing in the database!\n', freqband);
    if fileoutput
        fid = fopen('sessions_refFreq_missing.txt','a');
        fprintf(fid,'%s   %s   %s \n',freqband, out_struct.head.Session.val, parameter.session_name);
        fclose(fid);
        %return
    end
end
v0 =out_struct.Observables.(fileRefv).RefFreq.val(1); % MHz
if length(v0)>1
    disp('Reference frequency as vector?!?')
end

% ref. channel frequencies: rvi
fileChannInfo = ['ChannelInfo_' freqband];
if ~isfield(out_struct.Observables,fileChannInfo)
    fprintf('\nChannel Info for %s is missing in the database!\n', freqband);
    if fileoutput
        fid = fopen('sessions_ChannelInfo_missing.txt','a');
        fprintf(fid,'%s   %s   %s \n',freqband, out_struct.head.Session.val, parameter.session_name);
        fclose(fid);
        %return
    end
end
ChanFrAtr = {out_struct.Observables.(fileChannInfo).ChannelFreq.attr.name}; % check if channel frequency is same for all observations
idRep = contains(ChanFrAtr,'REPEAT');


if sum(idRep)
    nObs = int64(out_struct.Observables.(fileChannInfo).ChannelFreq.attr(idRep).val);
    rvi = repmat(out_struct.Observables.(fileChannInfo).ChannelFreq.val, 1, nObs);  % MHz   
else
    rvi = out_struct.Observables.(fileChannInfo).ChannelFreq.val;  % MHz
    nObs=size(rvi,2);  
end
nChan=size(rvi,1);

NumChanAtr = {out_struct.Observables.(fileChannInfo).NumChannels.attr.name}; % check if number of channels is same for all observations
idRep = contains(NumChanAtr,'REPEAT');
if sum(idRep)
    NumChanPerObs=repmat(out_struct.Observables.(fileChannInfo).NumChannels.val,nObs,1);
else
    NumChanPerObs=out_struct.Observables.(fileChannInfo).NumChannels.val;
end


%mean sub-band frequency, (sample rate)/4
SamRate = out_struct.Observables.(fileChannInfo).SampleRate.val *1e-6; % MHz

% non-standard values in vgosDB for Sample rate
stats = cellstr(out_struct.head.StationList.val'); % may be station dependent


if isfield (out_struct.head, 'ExpName')
   exper=out_struct.head.ExpName.val';
else
   exper=out_struct.head.Session.val';
end

if fileoutput_F
    fid = fopen(['sampleRate_' num2str(SamRate(1)) '.txt'],'a');
    fprintf(fid,'%s   %s \n', exper, parameter.session_name);
    fclose(fid);
    %return
end



if SamRate(1) == -32.768 % e.g. 08NOV12XA, 18JAN25XC
    if sum(contains(deblank(stats),'VLBA')) > 0 %ug002 USNO-CRF sessions, UG002A/18JAN18XC
        SamRate = 64; %MSamples
    else
        SamRate = 16; % it is a guess
    end
elseif SamRate(1) == -28.672 % ug002 USNO-CRF sessions
    SamRate = 64; %MSamples , guess
elseif SamRate(1) == 0 % E Sessions
    SamRate = 32; %? MSamples (chan.BW 16MHz) 
    %SamRate = 8; %? MSamples (chan.BW 4MHz) 
    % fprintf('Zero sample rate in vgosDB!!!')
    % fid = fopen(['sessions_Samrate0.txt'],'a');
    % fprintf(fid,'%s   %s     %s\n', out_struct.head.Session.val, parameter.session_name, freqband);
    % fclose(fid);
    % return
elseif SamRate(1) == 9.216
    SamRate = 16; %MSamples; guess
elseif ~isempty(find(SamRate==64))    
    if sum(strcmp(deblank(stats),'NYALE13N')) %reduce Nn's channel width from 32 to 8 
        SamRate = 16;
    end
end
% Info> AUSCOPE stations and Ishioka only record USB channels

hBW = SamRate/4;  % MHz

% fileID = fopen('CHANNELBW/channels_SR_-32.768.txt')
% ChanBw = textscan(fileID,'%s  %f %f %s %s');
% fclose(fileID);
% 
% isesBw = strcmpi(ChanBw{4},string(exper));
% if sum(isesBw)>0
%     sesBW = ChanBw{2}(isesBw);
%     hBW = sesBW/2;
% 
%     fid = fopen(['changedChanBW.txt'],'a');
%     fprintf(fid,'%10s %20s %2s %5.3f\n', exper, parameter.session_name, freqband, hBW*2); % channel bandwidth in MHz
%     fclose(fid);
% end


%% WEIGHTS

% NumSamples '# of samples by sideband and channel'
% NumAp '# of AP by sideband and channel'
if isfield (out_struct.Observables.(fileChannInfo), 'NumSamples')
    namparNS = 'NumSamples';
elseif isfield (out_struct.Observables.(fileChannInfo), 'NumAp')
    namparNS = 'NumAp';
        disp('NumAP instead of NumSamples')
        if fileoutput
            fid = fopen(['sessions_NumAp_instead_NumSamples.txt'],'a');
            fprintf(fid,'%s   %s \n', out_struct.head.Session.val, parameter.session_name);
            fclose(fid);
            %return
        end
else    
   disp('Both: NumSamples and NumAp missing!!!')    
        if fileoutput
            fid = fopen(['sessions_NumAp_and_NumSamples_missing.txt'],'a');
            fprintf(fid,'%s   %s     %s\n', out_struct.head.Session.val, parameter.session_name, freqband);
            fclose(fid);
            %return
        end
end

% specify lsb and usb and push frequency to the middle of the channel
% info from S-band
% find LSB and USB according to 0 in S-band
% distinction according to 0 in S-band (!!!assumption!!!: only USB in S-band)
if isfield (out_struct.Observables,'ChannelInfo_bS')
    NS1=out_struct.Observables.ChannelInfo_bS.(namparNS).val(:,:,1);
else
    NS1=out_struct.Observables.(fileChannInfo).(namparNS).val(:,:,1); % S-band missing
end
sNS1 = sum(NS1,2); % sum numSamples over channels
[irowUSB,~]= find(sNS1>0);
[irowLSB,~] = find(sNS1==0);

    % the same number of LSB and USB samples
    % look at the program name
    % the USB/LSB order is actually not needed - the freq. will not be shifted anyway
    if sNS1(1)==sNS1(2)
       prg=out_struct.head.Program.val';
        if contains(prg,'db2vgosDB') % db2vgosDB", then the order is LSB/USB
            irowLSB = 1;
            irowUSB = 2;
        elseif contains(prg,'vgosDbMake') % "vgosDbMake”, then the order is USB/LSB,
            irowLSB = 2;
            irowUSB = 1;
        else
            display('vgosDB generated either with db2vgosDB nor with vgosDbMake!!!')  
        end
    end



NamParNSAtr = {out_struct.Observables.(fileChannInfo).(namparNS).attr.name}; % check if channel frequency is same for all observations
idRep = contains(NamParNSAtr,'REPEAT');

if sum(idRep)
    nNS = int64(out_struct.Observables.(fileChannInfo).(namparNS).attr(idRep).val); %nNS has to equal nObs
    v=out_struct.Observables.(fileChannInfo).(namparNS).val;
    if sum(v,"all")==0
        v(2,:) = v(2,:) + 1; % put the same value (1) for number of samples if 0 are in vgosDB
        irowUSB = 2; % since there is no info in the vgosDB, we define usb channels (ASSUMPTION: USB were observed - maybe it is  wrong?)
        irowLSB = 1;
        disp('Identical NumSamples for USB channels!!!')
            if fileoutput
                fid = fopen(['sessions_identical_NumSamples_USBonly.txt'],'a');
                fprintf(fid,'%s   %s \n', out_struct.head.Session.val, parameter.session_name);
                fclose(fid);
                %return
            end
    end
    Nsam = reshape(repmat(v,nNS,1),2,4,[]);
else
    Nsam = out_struct.Observables.(fileChannInfo).(namparNS).val;  
end

% strange situation which is in the S2 vgosDB (_XG)
if size(Nsam,2) ~= nChan
    Nsam=Nsam./Nsam;
    Nsam(:,1:nChan,:)=[Nsam Nsam(:,1:(nChan-size(Nsam,2)),:)];
end


sNsam = sum(Nsam,1);
sNsam2=squeeze(sNsam)'; % sum (NiLSB + NiUSB)


% channel amplitude
AmpAvail='true';
if isfield (out_struct.Observables.(fileChannInfo), 'ChanAmpPhase')
    namparAmp = 'ChanAmpPhase';
elseif isfield (out_struct.Observables.(fileChannInfo), 'VFRQAM')
        namparAmp = 'VFRQAM';
            disp('VFRQAM instead of ChanAmpPhase')
            if fileoutput
                fid = fopen('sessions_VFRQAM_instead_ChanAmpPhase.txt','a');
                fprintf(fid,'%s   %s \n', out_struct.head.Session.val, parameter.session_name);
                fclose(fid);
                %return
            end
else
    disp('Both: VFRQAM and ChanAmpPhase missing!!! Unit weights are used.')
    AmpAvail='false';
end

if AmpAvail
    if size(out_struct.Observables.(fileChannInfo).(namparAmp).val(1,:,:),3)>1
        ChAmpl = out_struct.Observables.(fileChannInfo).(namparAmp).val(1,:,:); %(0-1)
    else
        ChAmpl = out_struct.Observables.(fileChannInfo).(namparAmp).val; %  S2 vgosDB (_XG)
    end
    sChAmpl = squeeze(ChAmpl)';
end


% ρi = (NiLSB + NiUSB) ∗Ampi    :weights
if AmpAvail
    wi = sNsam2 .* sChAmpl;
else
    wi = squeeze(sum(double(Nsam > 0),1))';
end

tableSB = zeros(nChan,nObs);
for iObs=1:nObs
    NSam1 = Nsam(:,:,iObs);
    
    %idShift = find(NSam1(1,:)~=NSam1(2,:));
    % it happens that the numSample is a little bit different between the
    % LSB and USB. Compare, if it is greater than 10%
    p10 = 0.1.*NSam1(1,:);

    idShift = find(abs(NSam1(1,:)-NSam1(2,:)) > p10);

    if ~isempty(idShift)
        [irowSam,icolSam,~] = find(NSam1(:,idShift)~=0); % check for USB and LSB
        [uicol, ~, ~] = unique(icolSam);
        if length(uicol) ~= length(icolSam) % [4 4 4 4 46 46 46 46; 46 0 0 0 0 0 0 46]
            counts = accumarray(icolSam, 1);
            bothSBsCOL = uicol(counts > 1);
                indcol=[];
                [indcol, ~]=find(ismember(icolSam,bothSBsCOL));
                irowSam(indcol)=[];
                indcol=[];
                [indcol, ~]=find(ismember(idShift,bothSBsCOL));
                idShift(indcol)=[];
        end

        % check over all channels (in case some channel was dropped)
        for ichn = 1:length(irowSam)
            if irowSam(ichn)==irowUSB
                tableSB(idShift(ichn),iObs) = 1;
            elseif irowSam(ichn)==irowLSB
                tableSB(idShift(ichn),iObs) = -1;
                if iObs==1
                    fprintf('LSB only channel identified in %s. Is this expected???\n',freqband)              
                    if fileoutput
                        fid = fopen(['sessions_LSBonly.txt'],'a');
                        fprintf(fid,'%s   %s    %s\n', out_struct.head.Session.val, parameter.session_name,freqband);
                        fclose(fid);
                    end
                end
            end
        end
    else
        if iObs==1
            disp('Channel frequency is not shifted.')
        end
         if fileoutput
            fid = fopen(['sessions_allLSBandUSB.txt'],'a'); % Both LSB and USB observed, so the frequency is not shifted
            fprintf(fid,'%s   %s   %8.0f\n', out_struct.head.Session.val, parameter.session_name, iObs);
            fclose(fid);
         end
    end
end

% if freqband=='bX'
%     tableSB(8,:)=tableSB(8,:).*-1
% end
tabhBW = tableSB.*hBW';

% move reference channel frequencies to the middle of the channels
viAll = rvi+tabhBW; % MHz


fileQualityCode = ['QualityCode_' freqband];
emptyCells = cellfun(@(x) strcmp(x, {''}), deblank(num2cell(out_struct.Observables.(fileQualityCode).QualityCode.val)));
out_struct.Observables.(fileQualityCode).QualityCode.val(emptyCells)=['0'];
QualityCode = out_struct.Observables.(fileQualityCode).QualityCode.val;% Quality Code Flag

eFreq(1:nObs,1) = 0;
qflag_v(1:nObs,1) = 0;
sigma_eFreq(1:nObs,1) = 0;
for iObs= 1 : nObs
 %wi(iObs,:)=1;
 %wi(iObs,[1 8])=2*900000;
    roi = wi(iObs,1:NumChanPerObs(iObs));
    vi = viAll(1:NumChanPerObs(iObs),iObs); %MHz
    vimv0 = vi-v0; %MHz

    s1 = sum(roi);
    s2 =(roi * (vimv0.^2));
    s3 = (roi * vimv0)^2;

    s4 = roi * vimv0;
    s5 = sum(roi'./vi);
    s6 = sum(roi);
    s7 = roi * (vimv0./vi);

    eFreq(iObs,1) = sqrt((s1*s2-s3)/(s4*s5 - s6*s7)); %MHz
    if NumChanPerObs(iObs) == 1
       eFreq(iObs,1) = vi; % effective freq. equals the one channel
    end


    % neef = (sum(roi))^2 / (roi*roi'); %effective sample size
    % stder2 = roi* ((vi-eFreq(iObs,1)).^2) / sum(roi);
    % sigma_eFreq(iObs) = sqrt(stder2 /neef)  / NumChanPerObs(iObs) ; %MHz 
    sigma_eFreq(iObs) = 0 ; %MHz  constant 


    % iono-flag where the effective freq. could not be calculated in the correct way
    if isnan(eFreq(iObs,1)) || eFreq(iObs,1) == 0
        qflag_v(iObs) = -1;
        eFreq(iObs,1) = 0;
    elseif strcmp(QualityCode(iObs),'0')
        qflag_v(iObs) = -1;
    end

end

eFreqGHz = eFreq.*1e-3;
sigma_eFreqGHz= sigma_eFreq.*1e-3;





%% Plots


  % if isfield (out_struct.head, 'ExpName')
  %     exper=out_struct.head.ExpName.val';
  % else
  %     exper=out_struct.head.Session.val';
  % end
% % 
% % % figure(1)
% % % plot(1:nObs,eFreq,'.')
% % % ylabel('[MHz]')
% % % title([ses ' w from NumSamples'])
% % % print('-dpdf' ,'-r500',[ses '_w_NumSamples']);
% % 
% % % figure(2)
% % % plot(1:nObs,eFreq,'.')
% % % ylabel('[MHz]')
% % % title([ses ' w 900000'])
% % % print('-dpdf' ,'-r500',[ses '_w_900000']);
% % 
% %  % figure(3)
% %  % plot(1:nObs,eFreq,'.')
% %  % ylabel('[MHz]')
% %  % title([ses ' w from AP per channel'])
% %  % print('-dpdf' ,'-r500',[ses '_w_APchannel']);
% %  % 
% %  % 
% % 
% fileEF = ['EffFreq_' freqband];
% eFreqDB = [out_struct.ObsDerived.(fileEF).FreqGroupIono.val];
% dEF = eFreq-eFreqDB;
% 
% figure(3)
%  plot(1:nObs,eFreq,'.')
%  ylabel('[MHz]')
%  title([exper])
%  print('-dpdf' ,'-r500',[exper '_' freqband  '_eFreq']);
% 
% 
%   figure(4)
%   plot(1:nObs,eFreqDB,'.')
%   ylabel('[MHz]')
%   title([exper ' vgosDB'])
%   print('-dpdf' ,'-r500',[exper '_' freqband '_eFreq_vgosDB']);
% 
% 
%  figure(5)
%  plot(1:nObs,dEF,'.')
%  ylabel([freqband ' [MHz]'])
%  title([exper ' vievs-vgosDB'])
%  print('-dpdf' ,'-r500',[exper '_' freqband  '_eFreq_diff']);
% 
% 
% f6=figure;
% for i=1:size(rvi,1)
%     plot(rvi(i,:),'.')
%     hold on
% end
% hold off
% ylabel([freqband ' Reference Freq. [MHz]'])
% title([exper ])
%  print(f6,'-dpdf' ,'-r500',[exper '_' freqband  '_refFreq']);


% rvi(rvi==0)=NaN;
%  f7=figure;
% for i=1:size(rvi,1)
%     plot(rvi(i,:),'.')
%     hold on
% end
% errorbar(1:nObs,eFreq,sigma_eFreq,'.')
% 
% hold off
% ylabel([freqband ' Freq. [MHz]'])
% title([exper ])
%  print(f7,'-dpdf' ,'-r500',[exper '_' freqband  '_refFreq']);