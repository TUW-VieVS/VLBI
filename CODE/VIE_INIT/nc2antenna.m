% ************************************************************************
%     nc2antenna
% ************************************************************************
% Description:
%   This function creates the VieVS antenna struct from netCDF files
%
% Input:	
%   Both variables from the netCDF files - out_struct and nc_info.
%   - wrapper_data: Struct, containing all information loaded from the vgosDB wrapper file
%
% Output:
%   The scan structure array.
% 
% External calls: 	
%   
%       
% Coded for VieVS: 
%   Jul 2012 by Matthias Madzak
%
% Revision: 
%   yyyy-mm-dd,FIRSTNAME SECONDNAME:
%   2015-11-11, M. Madzak: New Goddard structure adopted
%   2016-07-06, A. Hellerschmied: - Header added
%                                 - Evaluation of GPT and GPT2 moved to vie_mod
%   2017-01-23, D. Landskron: Initilization of GPT2 changed to GPT3, GPT removed
%   2017-02-09, D. Landskron: Preallocation extended
%   2017-02-22, A. Hellerschmied: antenna.psd initialized
%   2018-07-06, D. Landskron: vm1 renamed to vmf1 and VMF3 added to the troposphere models 
%   2019-01-15, M. Mikschi: added gravitational deformation
%   2020-12-06, J. Boehm: datum information from manualTRF was enabled
%
% ************************************************************************


function antenna=nc2antenna(out_struct, trf, chosenTrf, wrapper_data)

url_vievswiki_create_superstation = 'http://vievswiki.geo.tuwien.ac.at/doku.php?id=public:vievs_manual:data#create_a_superstation_file';

% number of stations in session
nStat=size(out_struct.head.StationList.val,2);

antenna(nStat)=struct('IDsuper', [], 'in_trf', [], 'name', [],  'x', [], 'y', [], ...
    'z', [], 'vx', [], 'vy', [], 'vz', [], 'epoch', [], 'start', [], ...
    'end', [], 'firstObsMjd', [], 'x_sigma', [], 'y_sigma', [], 'z_sigma', [],...
    'vx_sigma', [], 'vy_sigma', [], 'vz_sigma', [],...
    'thermal', [], 'comments', [], 'domes', [], 'code', [], 'ecc', [], 'gravdef', [], ...
    'ecctype', [], 'axtyp', [], 'offs', [], 'gpt3pres', [],  'gpt3temp', [],  'gpt3e', [], 'gpt3', [], 'noGrad', [], ...
    'cto', [], 'cta', [], 'cnta_dx', [],  ...
    'vmf3', [], 'vmf1', [], 'opl', [], 'numobs', [],  'lastObsMjd', [],  'psd', []);

% for some things no loop is needed
%[antenna.session]=deal(out_struct.head.Session');
                        
                        
for iStat=1:nStat
    % get name of current station
    curName=out_struct.head.StationList.val(:,iStat)';
    curNameLong=[curName, '          '];
    antenna(iStat).name=curNameLong(1:8);
    
    % find current station in trf
    curStatInTrf=strcmp({trf.name}, curName);
    
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
    
    % get first and last mjd
    curYMDHM=out_struct.stat(iStat).TimeUTC.YMDHM.val;
    curSecond=out_struct.stat(iStat).TimeUTC.Second.val;
    
    antenna(iStat).firstObsMjd=cal2jd(double(curYMDHM(1,1)), double(curYMDHM(2,1)), double(curYMDHM(3,1)))-2400000.5 +double(curYMDHM(4,1))/24+double(curYMDHM(5,1))/60/24+curSecond(1)/60/60/24;
    antenna(iStat).lastObsMjd=cal2jd(double(curYMDHM(1,end)), double(curYMDHM(2,end)), double(curYMDHM(3,end)))-2400000.5 +double(curYMDHM(4,end))/24+double(curYMDHM(5,end))/60/24+curSecond(end)/60/60/24;
    
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
        fprintf('No valid %s coordinates for %s => get vievsTrf coordinates (no NNT/NNR conditions applied!)\n', chosenTrf, curName)
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
    if strcmp(chosenTrf, 'manualTrf')  % Johannes 6 Dec 2020
        antenna(iStat).in_trf=trf(indCurStatInTrf).(chosenTrf).break(bnr).indatum; % Johannes 6 Dec 2020
    else % Johannes 6 Dec 2020
        if strcmp(chosenTrf, 'vievsTrf')
            antenna(iStat).in_trf=trf(indCurStatInTrf).(chosenTrf).break(bnr).indatum;
        else
            if strcmp(trfToTake, 'vievsTrf')
                antenna(iStat).in_trf=0;
            else
                antenna(iStat).in_trf=1;
            end
        end
    end    % Johannes 6 Dec 2020
    
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
    nc_filename = get_nc_filename('ObsCrossRef', wrapper_data.Observation.CrossReference.files, 1);
    nobstemp3=sum(sum(out_struct.CrossReference.(nc_filename).Obs2Baseline.val==iStat)); % 808
    
       
    antenna(iStat).numobs=nobstemp3;
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

