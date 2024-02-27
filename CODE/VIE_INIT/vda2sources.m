% This function creates the sources struct from the supersource file as
% well as from VDA data.

function sources=vda2sources(DATA1,valtab, crf, crffile)

chosenCrf=crffile{2};
ind=contains(DATA1,'DATA.1 NUMB_SOU');
currline=DATA1(ind);
seg = split(currline);
nSources=str2double(seg{7});


ind=contains(DATA1,'DATA.1 SRCNAMES');
ourblock=DATA1(ind);
for i=1:size(ourblock,1)
    input_str = ourblock{i};
    seg = split(input_str);
    oursou(i)=(seg(7));
end


sources(nSources)=struct('name', [],'IVSname', [], 'IERSname', [], 'ICRFdes', [], 'J2000nameS', [], 'ra2000', [], 'de2000', [], 'ra_sigma', [], 'de_sigma', [], 'corr', [], 'in_crf', [], 'numobs', [], 'flag_defining', [], 'firstObsMjd', [], 'lastObsMjd', []);

% for all sources
for iSource=1:nSources
    curName=oursou(iSource);
    curNameLong=[char(curName), '          '];

    indSourceInCrfIERS=strcmp(deblank({crf.IERSname}),deblank(curName));
    if sum(indSourceInCrfIERS) == 0
        indSourceInCrfIERS=strcmp(deblank({crf.IVSname}),deblank(curName));
    end
    curName = crf(indSourceInCrfIERS).IVSname;



    % find index in crf
    indSourceInCrf=strcmp(deblank({crf.IVSname}),deblank(curName));
    if sum(indSourceInCrf) == 0
        warning('ERROR: source %s not found in supersource file!\nAdd source to supersource file!', curName)
        fprintf('Run: check_sources_in_vgosDB_file or check_sources_in_NGS_file function\n')
    end
    
    % get index
    indSourceInCrf=find(indSourceInCrf);
    
    % write name
    sources(iSource).name=curNameLong(1:8);
    sources(iSource).IERSname=crf(indSourceInCrf).IERSname;
    sources(iSource).IVSname = crf(indSourceInCrf).IVSname;
    sources(iSource).ICRFdes=crf(indSourceInCrf).designation(6:end);
    
    id=valtab(:,17)==iSource;
    mjdall=(valtab(id,18));
    numobs1=length(mjdall);
        
    sources(iSource).firstObsMjd = min(mjdall);
    sources(iSource).lastObsMjd = max(mjdall);
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
    

end


% Create a sub-structure in "sources" for quasars sources:
q = sources;
clear sources
sources.q   = q;
sources.s 	= [];

