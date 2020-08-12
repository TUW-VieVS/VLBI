% This function creates the sources struct from the supersource file as
% well as from netCDF data.
%
% Input
%  - wrapper_data: Struct, containing all information loaded from the vgosDB wrapper file
%
% LOG
%  11.11.2015 M. Madzak  New Goddard structure adopted
%   07 Jul 2016 by D. Mayer: included NNR for only defining sources
%   15 Jul 2016 by D. Mayer: Added IVS name
%   19 Jul 2016 by D. Mayer: Added firstObsMjd and lastObsMjd
%   26 Sep 2016 by H. Krasna: changes related to updated supersource file -crf.commonName deleted. We assumed that source names in the vgosDB are the IVS names - crosschecking between IVS and IERS names is not desirable.
%   13 Nov 2017 by J. Gruber: def. sources update

function sources=nc2sources(out_struct, crf, chosenCrf, wrapper_data)

% get number of sources
nSources=size(out_struct.head.SourceList.val,2);

sources(nSources)=struct('name', [],'IVSname', [], 'IERSname', [], 'ICRFdes', [], 'ra2000', [], 'de2000', [], 'ra_sigma', [], 'de_sigma', [], 'corr', [], 'in_crf', [], 'numobs', [], 'flag_defining', [], 'firstObsMjd', [], 'lastObsMjd', []);

% % for speed improvement: delete all sources in crf which do not occur
% sourceNamesCellStr=cellstr(out_struct.head.SourceName');
% neededSourcesInCrf=ismember({crf.IERSname}, sourceNamesCellStr) | ...
%     ismember({crf.commonName}, sourceNamesCellStr) | ...
%     ismember({crf.designation}, sourceNamesCellStr);

nc_filename = get_nc_filename('ObsCrossRef', wrapper_data.Observation.CrossReference.files, 1);
obs2scan = out_struct.CrossReference.(nc_filename).Obs2Scan.val;


% for all sources
for iSource=1:nSources
    curName=out_struct.head.SourceList.val(:,iSource)';
    
    % find index in crf
    indSourceInCrf=strcmp(deblank({crf.IVSname}),deblank(curName));
    if sum(indSourceInCrf) == 0
        indSourceInCrf=strcmp(deblank({crf.IERSname}),deblank(curName));
        if sum(indSourceInCrf) == 0
        warning('ERROR: source %s not found in supersource file!\nAdd source to supersource file!', curName)
        fprintf('Run: check_sources_in_vgosDB_file or check_sources_in_NGS_file function\n')
        end
    end
    
    % get index
    indSourceInCrf=find(indSourceInCrf);
    
    % write name
    sources(iSource).name=curName;
    sources(iSource).IERSname=crf(indSourceInCrf).IERSname;
    sources(iSource).IVSname = crf(indSourceInCrf).IVSname;
    sources(iSource).ICRFdes=crf(indSourceInCrf).designation(6:end);
    
    % number of observations
    %     out_struct.Observables.Source.Source % source names
    nc_filename = get_nc_filename('SourceCrossRef', wrapper_data.Session.CrossReference.files, 1);
    scansUsingCurSource=out_struct.CrossReference.(nc_filename).Scan2Source.val==iSource; % logicals
    
    
    numobs1=sum(ismember(obs2scan,find(scansUsingCurSource)));
    
    sources(iSource).numobs=numobs1;
    % in_crf
    
    % if there are coordinates for the chosen source in chosen CRF
    if isempty(crf(indSourceInCrf).(chosenCrf))
        sources(iSource).in_crf=0;
        crfToTake='vievsCrf';
    else
        % if vievsCrf was chosen -> take in_crf from there
        if strcmp(chosenCrf, 'vievsCrf')
            sources(iSource).in_crf=crf(indSourceInCrf).vievsCrf.in_crf;
        else
            sources(iSource).in_crf=1;
        end
        crfToTake=chosenCrf;
    end
    %David - always add information about defining sources
    %from the ICRF2
    if isfield(crf(indSourceInCrf),'icrf3sx')
        if ~isempty(crf(indSourceInCrf).icrf3sx)
            sources(iSource).flag_defining=crf(indSourceInCrf).icrf3sx.defining;
        else
            sources(iSource).flag_defining=0;
        end
    else
        sources(iSource).flag_defining=0;
    end
    
    % RA/DE
    sources(iSource).ra2000=      crf(indSourceInCrf).(crfToTake).ra;
    sources(iSource).de2000=      crf(indSourceInCrf).(crfToTake).de;
    
    % sigmas (RA/DE)
    if isfield(crf(indSourceInCrf).(crfToTake), 'ra_sigma')
        sources(iSource).ra_sigma=crf(indSourceInCrf).(crfToTake).ra_sigma;
        sources(iSource).de_sigma=crf(indSourceInCrf).(crfToTake).de_sigma;
    else
        sources(iSource).ra_sigma=[];
        sources(iSource).de_sigma=[];
    end
    if isfield(crf(indSourceInCrf).(crfToTake), 'corr')
        sources(iSource).corr=    crf(indSourceInCrf).(crfToTake).corr;
    else
        sources(iSource).corr=[];
    end
    
    
    %out_struct.CrossReference.SourceCrossRef.Scan2Source.val % source ID per scan
    %out_struct.Scan.TimeUTC.YMDHM.val % time of each scan
    nc_filename = get_nc_filename('TimeUTC', wrapper_data.Scan.Scan.files, 1);
    
    YMDHM_tmp = out_struct.Scan.(nc_filename).YMDHM.val;
    mjd_of_scans = cal2jd(double(YMDHM_tmp(1,:)), double(YMDHM_tmp(2,:)),double(YMDHM_tmp(3,:))) - 2400000.5+double(YMDHM_tmp(4,:))./24 +double(out_struct.Scan.(nc_filename).YMDHM.val(5,:))./24/60 + double(out_struct.Scan.(nc_filename).Second.val')./24/60/60;

    sources(iSource).firstObsMjd = min(mjd_of_scans(scansUsingCurSource));
    sources(iSource).lastObsMjd = max(mjd_of_scans(scansUsingCurSource));
end
