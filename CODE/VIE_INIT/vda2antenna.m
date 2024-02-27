function antenna=vda2antenna(DATA1,nSta,valtab, trf, trffile)

chosenTrf=trffile{2};

url_vievswiki_create_superstation = 'http://vievswiki.geo.tuwien.ac.at/doku.php?id=public:vievs_manual:data#create_a_superstation_file';

antenna(nSta)=struct('IDsuper', [], 'in_trf', [], 'name', [],  'x', [], 'y', [], ...
    'z', [], 'vx', [], 'vy', [], 'vz', [], 'epoch', [], 'start', [], ...
    'end', [], 'firstObsMjd', [], 'x_sigma', [], 'y_sigma', [], 'z_sigma', [],...
    'vx_sigma', [], 'vy_sigma', [], 'vz_sigma', [],...
    'thermal', [], 'comments', [], 'domes', [], 'code', [], 'ecc', [], 'gravdef', [], ...
    'ecctype', [], 'axtyp', [], 'offs', [], 'gpt3pres', [],  'gpt3temp', [],  'gpt3e', [], 'gpt3', [], 'noGrad', [], ...
    'cto', [], 'cta', [], 'cnta_dx', [],  ...
    'vmf3', [], 'vmf1', [], 'opl', [], 'numobs', [],  'lastObsMjd', [],  'psd', []);

% for some things no loop is needed
%[antenna.session]=deal(out_struct.head.Session');

ind=contains(DATA1,'DATA.1 SITNAMES');
ourblock=DATA1(ind);
for i=1:size(ourblock,1)
    input_str = ourblock{i};
    seg = split(input_str);
    ourstat(i)=(seg(7));
end

nobsSta=zeros(1,nSta);
ind=contains(DATA1,'DATA.1 NOBS_STA');
ourblock=DATA1(ind);
for i=1:size(ourblock,1)
    input_str = ourblock{i};
    seg = split(input_str);
    nobsSta(i)=str2double(seg(7));
end


                        
for iStat=1:nSta
    % get name of current station
    curName=char(ourstat(iStat));
    curNameLong=[curName, '          '];
    antenna(iStat).name=curNameLong(1:8);
    
    % find current station in trf
    curStatInTrf=strcmp({trf.name},  antenna(iStat).name);
    
    if sum(curStatInTrf)~=1
        error('Station %s not found in the superstation file. Add this station to the superstation file (vievsTRF = backup TRF) by following the steps described at %s\n', curName, url_vievswiki_create_superstation);
    else
        indCurStatInTrf=find(curStatInTrf);
    end
    
    %% write data
    % get ~general infos from superstation file
    antenna(iStat).IDsuper = indCurStatInTrf;
    antenna(iStat).thermal = trf(indCurStatInTrf).antenna_info;
    antenna(iStat).comments = trf(indCurStatInTrf).comments;
    antenna(iStat).domes = trf(indCurStatInTrf).domes;
    antenna(iStat).code = trf(indCurStatInTrf).code; 
    
%     % get first and last mjd
%     curYMDHM=out_struct.stat(iStat).TimeUTC.YMDHM.val;
%     curSecond=out_struct.stat(iStat).TimeUTC.Second.val;
    
    id=[find(valtab(:,3)==iStat); find(valtab(:,4)==iStat)];
    mjdall=(valtab(id,18));
        
    antenna(iStat).firstObsMjd = min(mjdall);
    antenna(iStat).lastObsMjd = max(mjdall);

    
    % Gravitational deformation
    if isfield(trf(indCurStatInTrf), 'gravdef') && ~isempty(trf(indCurStatInTrf).gravdef)      
        % find break
        bnr=find(antenna(iStat).firstObsMjd>=[trf(indCurStatInTrf).gravdef.break.start] & ...
            antenna(iStat).firstObsMjd<=[trf(indCurStatInTrf).gravdef.break.end]);
        
        antenna(iStat).gravdef = trf(indCurStatInTrf).gravdef.break(bnr);
    end
       
    % get trf coordinates from chosen trf
    if isempty(trf(indCurStatInTrf).(chosenTrf))
        if strcmp(chosenTrf, 'vievsTrf') % Station not even in vievsTrf (backup)
            error('Station %s not found in the superstation file (vievsTRF). Add this station to the superstation file by following the steps described at %s\n', antenna(iStat).name, url_vievswiki_create_superstation);
        end
        bnr = [];
    else
        trfToTake=chosenTrf;
        
        % find break
        if ~isempty(trf(indCurStatInTrf).(trfToTake))
            bnr=find(antenna(iStat).firstObsMjd>=[trf(indCurStatInTrf).(trfToTake).break.start] & antenna(iStat).firstObsMjd<=[trf(indCurStatInTrf).(trfToTake).break.end]);
        else
            bnr = [];
        end
    end

    % no break is found for the "trfToTake"
    if isempty(bnr)
%         fprintf('No valid %s coordinates for %s => get vievsTrf coordinates (no NNT/NNR conditions applied!)\n', trffile, curName)
        % so: no valid epoch for official (e.g.. VTRF2008) TRF -> get vievsTrf break
        trfToTake='vievsTrf';
        
        bnr=find(antenna(iStat).firstObsMjd>=[trf(indCurStatInTrf).(trfToTake).break.start] & antenna(iStat).firstObsMjd<=[trf(indCurStatInTrf).(trfToTake).break.end]);

        if isempty(bnr)
            error('Station %s not found in the superstation file (vievsTRF). Add this station to the superstation file by following the steps described at %s\n', antenna(iStat).name, url_vievswiki_create_superstation);
        end
    end
    
    curBreak=trf(indCurStatInTrf).(trfToTake).break(bnr);
    
    % set in_trf
    % if user has chosen vievsTrf
    if strcmp(trffile, 'vievsTrf')
        antenna(iStat).in_trf=trf(indCurStatInTrf).(trffile).break(bnr).indatum;
    else
        if strcmp(trfToTake, 'vievsTrf')
            antenna(iStat).in_trf=0;
        else
            antenna(iStat).in_trf=1;
        end
    end
        
    
    antenna(iStat).x=curBreak.x;
    antenna(iStat).y=curBreak.y;
    antenna(iStat).z=curBreak.z;
    if isfield(curBreak, 'vx')
        antenna(iStat).vx=curBreak.vx;
        antenna(iStat).vy=curBreak.vy;
        antenna(iStat).vz=curBreak.vz;
        antenna(iStat).epoch=curBreak.epoch;
        antenna(iStat).start=curBreak.start;
        antenna(iStat).end=curBreak.end;
    else
        antenna(iStat).vx=0;
        antenna(iStat).vy=0;
        antenna(iStat).vz=0;
        antenna(iStat).epoch=0;
        antenna(iStat).start=0;
        antenna(iStat).end=99999;
    end
    
    if isfield(curBreak, 'x_sigma')
        antenna(iStat).x_sigma=curBreak.x_sigma;
        antenna(iStat).y_sigma=curBreak.y_sigma;
        antenna(iStat).z_sigma=curBreak.z_sigma;
        antenna(iStat).vx_sigma=curBreak.vx_sigma;
        antenna(iStat).vy_sigma=curBreak.vy_sigma;
        antenna(iStat).vz_sigma=curBreak.vz_sigma;
    else
        antenna(iStat).x_sigma=[];
        antenna(iStat).y_sigma=[];
        antenna(iStat).z_sigma=[];
        antenna(iStat).vx_sigma=[];
        antenna(iStat).vy_sigma=[];
        antenna(iStat).vz_sigma=[];
    end
    
    % numobs
    antenna(iStat).numobs=nobsSta(iStat);
    
    if ~isempty(trf(indCurStatInTrf).antenna_info)
        antenna(iStat).axtyp=trf(indCurStatInTrf).antenna_info.mount(4:end);
        antenna(iStat).offs=trf(indCurStatInTrf).antenna_info.axis_offset;
        antenna(iStat).axtyp=trf(indCurStatInTrf).antenna_info.mount(4:end);
        antenna(iStat).axtyp=trf(indCurStatInTrf).antenna_info.mount(4:end);
        antenna(iStat).axtyp=trf(indCurStatInTrf).antenna_info.mount(4:end);
    end
    
    %% ECC/ECCTYPE
    antenna(iStat).ecc=[0 0 0];
    antenna(iStat).ecctype='NEU';
    if ~isempty(trf(indCurStatInTrf).ecc)
        % get start and end of all ecc-breaks
        nEccBreaks=size(trf(indCurStatInTrf).ecc.break,2);
        eccBreaks=zeros(nEccBreaks,2); % one line for one break; two cols: start/end
        
        for iEccBreak=1:nEccBreaks
            eccBreaks(iEccBreak,1)=cal2jd(str2double(trf(indCurStatInTrf).ecc.break(iEccBreak).starting(1:4)), ...
                str2double(trf(indCurStatInTrf).ecc.break(iEccBreak).starting(6:7)), ...
                str2double(trf(indCurStatInTrf).ecc.break(iEccBreak).starting(9:10)))-2400000.5 +...
                str2double(trf(indCurStatInTrf).ecc.break(iEccBreak).starting(12:13))/24+...
                str2double(trf(indCurStatInTrf).ecc.break(iEccBreak).starting(15:16))/60/24;    % start
                
            eccBreaks(iEccBreak,2)=cal2jd(str2double(trf(indCurStatInTrf).ecc.break(iEccBreak).ending(1:4)), ...
                str2double(trf(indCurStatInTrf).ecc.break(iEccBreak).ending(6:7)), ...
                str2double(trf(indCurStatInTrf).ecc.break(iEccBreak).ending(9:10)))-2400000.5 +...
                str2double(trf(indCurStatInTrf).ecc.break(iEccBreak).ending(12:13))/24+...
                str2double(trf(indCurStatInTrf).ecc.break(iEccBreak).ending(15:16))/60/24; % end
        end
        
        % get correct break
        eccBnr=find(antenna(iStat).firstObsMjd>=eccBreaks(:,1) & ...
            antenna(iStat).firstObsMjd<=eccBreaks(:,2));
        
        if isempty(eccBnr)
            fprintf('No valid ECC break found for station %s\n', curName)
        else
            antenna(iStat).ecc=[trf(indCurStatInTrf).ecc.break(eccBnr).FCE, ...
                trf(indCurStatInTrf).ecc.break(eccBnr).SCE, ...
                trf(indCurStatInTrf).ecc.break(eccBnr).TCE];
            antenna(iStat).ecctype=trf(indCurStatInTrf).ecc.break(eccBnr).type_e;
        end
    
    end
    
    % ##### Met. data: #####     

    % Init. flags:
    antenna(iStat).gpt3pres     = 0;
    antenna(iStat).gpt3temp     = 0;
    antenna(iStat).gpt3e        = 0;
    
    % Init. emptx fields => They will be filled in vie_mod (calc_met_data.m)
    antenna(iStat).gpt3.p      = [];
    antenna(iStat).gpt3.T      = [];    
    antenna(iStat).gpt3.dT     = [];
    antenna(iStat).gpt3.Tm     = [];
    antenna(iStat).gpt3.e      = [];
    antenna(iStat).gpt3.ah     = [];
    antenna(iStat).gpt3.aw     = [];
    antenna(iStat).gpt3.lambda = [];
    antenna(iStat).gpt3.undu   = [];
    antenna(iStat).gpt3.Gn_h   = [];
    antenna(iStat).gpt3.Ge_h   = [];
    antenna(iStat).gpt3.Gn_w   = [];
    antenna(iStat).gpt3.Ge_w   = [];
end

