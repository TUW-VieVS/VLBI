% ************************************************************************
%   Description:
%       This function creates for each parameter, which will be estimated
%       an indices-pair (ar): number of column, where the parameter is in the
%       old N matrix in this session X number of column, where the
%       parameter will be in the new global N matrix.
%
%
%   Input:
%      parGS               name of the parameters which will be
%                          estimated (1), fixed (0) or reduced (2),
%                          and the position of columns in the old N matrix
%      antenna             information about a priori position of the antenna (glob1.an)
%      refantbr            information about the discontinuities at the
%                          stations, with respective coordinates and velocities
%                          from the external file 'VLBI-DISCONT.txt'
%      refnamec            names of antenna, also multiple if there are
%                          breaks in the position, that will be considered
%      qnames              names of sources in this session
%      qrefname            names of all sources, which will be estimated in
%                          the global adjustment
%      ise                 number of this session
%      lnc                 number of stations for estimating coordinates
%      lnv                 number of stations for estimating velocities
%      lns                 number of sources
%
%
%   Output: 
%      ar                  actual (in this session) X reference (in the new global N matrix) indices
%      actp                actual number of est. parameters in this session
%      parGS               positions of the parameters in the new N matrix
%
%   External calls: 	
%      globind               					    											
%       
%   Coded for VieVS: 
%   01 Aug 2010 by Hana Spicakova
%
%   Revision: 
%   16 Sep 2011 by Hana Spicakova   bug by estimation of EOP with different
%   resolution was fixed
%  20 Jun 2012  by Hana Krásná: INTERNAL version
%               added: Love and Shida numbers, FCN per from SET, velocity of sources, acceleration of SSB, seasonal
%               variations in station positions (annual, semi-annual), pole
%               tide Love and Shida numbers
%  25 Sep 2012  by Hana Krásná: relativistic parameter gamma in INTERNAL
%               version
%  04 Oct 2012  by Hana Krásná  small bug fixed which occured if there
%               were no observations in the first interval of station
%               coordinates
%   04 Oct 2013 by Hana Krasna: estimation of antenna axis offset added
%   06 Dec 2013 by Hana Krasna: APL regression coefficient added
%   24 Aug 2015 by Sigrid Boehm: tidal ERP coefficients added to INTERNAL
%   10 Apr 2017 by David Mayer: added errormessage 
% ************************************************************************




function [ar,actp,parGS] = par_newcol(parGS,parGS_hl,antenna,refantbr,refnamec,qnames,qrefname,aostname,rgstname,ise,lnc,lnv,lns,lnao,...
                                               llove,lshida,lFCNset,stseasname,lse,nvsou,laccSSB,lhpole,llpole,lgamma,lrg,...
                                               ltidpm,ltidut);



[g] = globind(parGS);

actp=0; ar=[];
mjd=antenna.firstscan_mjd;

for j = 1: length(refantbr)
    refname(j,:) = refantbr(j,:).name;
end


ncoord=3*lnc;
nveloc=3*lnv;
nsou=2*lns;

% Antenna coordinates
if parGS(g.g_coord(1)).id==1
    % put station coordinates into the right columns
    clear arx ary arz
    arx =[]; ary =[]; arz =[];
    l=1;
    for k = 1:length(antenna.x)
        aname=antenna.name(k,:);
        p=strcmp(cellstr(aname),cellstr(refnamec));
        if sum(p)~=0
            [ref_idm,xx]=find(p==1);
            if length(ref_idm)>1
                p=strcmp(cellstr(aname),cellstr(refname));
                [ref_id,xx]=find(p==1);
                for i=1:length(refantbr(ref_id).break)-1
                    gg = find(refantbr(ref_id).break(i)<mjd && mjd<refantbr(ref_id).break(i+1));
                    if isempty(gg); gg=0; end
                    interv(i)=gg;
                end 
                
                % find only intervals with observations
                [xx1,xx2,xx1]=find(refantbr(ref_id).interv > 0); % HANA 2012-10-04
                id=find(interv(xx2),1);
                ref_idm=ref_idm(id); % index of the column      
                clear id xx2
            
            end

            arx(l,:)=[parGS(g.g_coord(1)).oldcol(k),ref_idm];       %actual_reference indices of station
            ary(l,:)=[parGS(g.g_coord(2)).oldcol(k),ref_idm+lnc];   %actual_reference indices of station
            arz(l,:)=[parGS(g.g_coord(3)).oldcol(k),ref_idm+2*lnc]; %actual_reference indices of station
            l=l+1;
        end
    end
%     if ~exist('arx','var')
%         error('Non of the antennas in session are estimated (probably reduced)');
%     end
    ar=vertcat(arx,ary,arz); % actual X reference indices
    actp=size(ar); actp=actp(1); % actp: actual number of est. parameters in this session
end


% Antenna velocities
if parGS(g.g_vel(1)).id==1
    clear arvx arvy arvz ref_id ref_idm
    arvx =[]; arvy=[]; arvz=[];
    % put station coordinates into the correct columns
    l=1;
    for k = 1:length(antenna.x)
        anames=antenna.name;
        aname=anames(k,:);
        p=strcmp(cellstr(aname),cellstr(refnamec));
        if sum(p)~=0
            [ref_idm,xx]=find(p==1);
            if length(ref_idm)>1
                p=strcmp(cellstr(aname),cellstr(refname));
                [ref_id,xx]=find(p==1);
                for i=1:length(refantbr(ref_id).break)-1
                    gg = find(refantbr(ref_id).break(i)<mjd && mjd<refantbr(ref_id).break(i+1));
                    if isempty(gg); gg=0; end
                    interv(i)=gg;
                end 
                % find only intervals with observations
                 [xx1,xx2,xx1]=find(refantbr(ref_id).interv > 0); % HANA 2012-10-04
                 id=find(interv(xx2),1);
                ref_idm=ref_idm(id); % index of the column
                clear id xx2
            end
            arvx(l,:)=[parGS(g.g_vel(1)).oldcol(k),ncoord+ref_idm];      %actual_reference indices of station
            arvy(l,:)=[parGS(g.g_vel(2)).oldcol(k),ncoord+ref_idm+lnv];   %actual_reference indices of station
            arvz(l,:)=[parGS(g.g_vel(3)).oldcol(k),ncoord+ref_idm+2*lnv]; %actual_reference indices of station
            l=l+1;
        end
    end
    ar=vertcat(ar,arvx,arvy,arvz); % actual X reference indices
    actp=size(ar); actp=actp(1); % actp: actual number of est. parameters in this session
end


% source coordinates
if parGS(g.g_srade(1)).id==1
    clear arra arde ref_idm
    l=1;
    % put source coordinates into the right columns
    for iso = 1:size(qnames,1)
        qname=qnames(iso,:);
        fq=strcmp(cellstr(qname),cellstr(qrefname)); % find the source
        if sum(fq)~=0
            [ref_idm]=find(fq==1);
            clear fq
            arra(l,:)=[parGS(g.g_srade(1)).oldcol(iso),ncoord+nveloc+ref_idm];       %actual_reference indices of source RA
            arde(l,:)=[parGS(g.g_srade(2)).oldcol(iso),ncoord+nveloc+ref_idm+lns];   %actual_reference indices of source De
            l=l+1;
        end
    end
    ar=vertcat(ar,arra,arde); % actual X reference indices
    actp=size(ar); actp=actp(1); % actp: actual number of est. parameters in this session
end

%% INTERNAL
% source velocities
if parGS(g.g_svrade(1)).id==1
    clear arra arde ref_idm
    l=1;
    % put source vel into the right columns
    for iso = 1:size(qnames,1)
        qname=qnames(iso,:);
        fq=strcmp(cellstr(qname),cellstr(qrefname)); % find the source
        if sum(fq)~=0
            [ref_idm]=find(fq==1);
            clear fq
            arvra(l,:)=[parGS(g.g_svrade(1)).oldcol(iso),ncoord+nveloc+nsou+ref_idm];       %actual_reference indices of source vel RA
            arvde(l,:)=[parGS(g.g_svrade(2)).oldcol(iso),ncoord+nveloc+nsou+ref_idm+lns];   %actual_reference indices of source vel De
            l=l+1;
        end
    end
    ar=vertcat(ar,arvra,arvde); % actual X reference indices
    actp=size(ar); actp=actp(1); % actp: actual number of est. parameters in this session
end
%%

globp = ncoord+nveloc+nsou+nvsou;   % number of st. coord + veloc + source coordinates + veloc in the global N-matrix 

% EOP
for ieop=g.g_eop
    if parGS(ieop).id==1
        oldcol = parGS(ieop).oldcol;
        mjd    = parGS(ieop).mjd;
        mjdstep= parGS(ieop).mjdstep;
        newcol = parGS(ieop).newcol;
        for i=1:length(oldcol)
            [c,ia]=intersect(mjdstep,mjd(i));
            areop(i,:)=[oldcol(i),globp+ia];
            newcol(ise,i)=globp+ia;
        end
        ar(actp+1 : actp+length(oldcol),:)=areop;
        actp =size(ar,1);
        globp=globp+length(mjdstep);
        parGS(ieop).newcol=newcol;
        clear mjdstep mjd oldcol areop newcol ia
    end
end

% Antenna Offsets
if parGS(g.g_ao).id==1
    l=1;
    arao=[];
    % put AO into the right columns
    for k = 1:length(antenna.x)
        aname=antenna.name(k,:);
        fq=strcmp(cellstr(aname),cellstr(aostname)); % find the antenna
        if sum(fq)~=0
            [ref_idm]=find(fq==1);
            clear fq
            arao(l,:)=[parGS(g.g_ao).oldcol(k),globp+ref_idm];   
            l=l+1;
        end
    end
    if ~isempty(arao)    
        ar(actp+1:actp+size(arao,1), :)=arao;
        actp=size(ar,1);
        parGS(g.g_ao).newcol=[globp+1:globp+lnao];
        globp=globp+lnao;        
    end
end



% seasonal variations in station positions R
if parGS(g.g_stseaspos_Ar(1)).id==1
    idwave=find(cell2mat(parGS(g.g_stseaspos_Ar(1)).spec)==1);
    lwave=length(idwave);
    nstses=size(stseasname,1);
    
    clear arAcr arAsr ref_idm allsesm_c allsesm_s
    l1=0;
    l2=0;
    % put antennas and amplitudes into the right columns
    for k = 1:length(antenna.x)
        aname=antenna.name(k,:);
        fq=strcmp(cellstr(aname),cellstr(stseasname)); % find the antenna
        if sum(fq)~=0
            [ref_idm]=find(fq==1);
            clear fq
            allsesm_c=cell2mat(parGS(g.g_stseaspos_Ar(1)).oldcol(k));
            allsesm_s=cell2mat(parGS(g.g_stseaspos_Ar(2)).oldcol(k));
            for i=1:lwave
                l1=l1+1;
                arAcr(l1,:)=[allsesm_c(idwave(i)),globp+(ref_idm-1)*2*lwave+i];       %actual_reference indices 
            end
            for i=1:lwave
                l2=l2+1;
                arAsr(l2,:)=[allsesm_s(idwave(i)),globp+(ref_idm-1)*2*lwave+i+lwave];   %actual_reference indices 
            end
        end
    end
    if l1>0
        ar=vertcat(ar,arAcr,arAsr); % actual X reference indices
        actp=size(ar,1); % actp: actual number of est. parameters in this session
    end

    if ise==lse
        for ian=1:nstses
            for iw=1:lwave
                stc(ian).newcol(iw)=globp+(ian-1)*2*lwave+iw;
                sts(ian).newcol(iw)=globp+(ian-1)*2*lwave+iw+lwave;
            end
        end
        parGS(g.g_stseaspos_Ar(1)).newcol={stc.newcol};
        parGS(g.g_stseaspos_Ar(2)).newcol={sts.newcol};
    end
    globp=globp+2*lwave*nstses;

end

% seasonal variations in station positions E
if parGS(g.g_stseaspos_Ae(1)).id==1
    idwave=find(cell2mat(parGS(g.g_stseaspos_Ae(1)).spec)==1);
    lwave=length(idwave);
    nstses=size(stseasname,1);
    
    clear arAce arAse ref_idm allsesm_c allsesm_s stc sts
    l1=0;
    l2=0;
    % put antennas and amplitudes into the right columns
    for k = 1:length(antenna.x)
        aname=antenna.name(k,:);
        fq=strcmp(cellstr(aname),cellstr(stseasname)); % find the antenna
        if sum(fq)~=0
            [ref_idm]=find(fq==1);
            clear fq
            allsesm_c=cell2mat(parGS(g.g_stseaspos_Ae(1)).oldcol(k));
            allsesm_s=cell2mat(parGS(g.g_stseaspos_Ae(2)).oldcol(k));
            for i=1:lwave
                l1=l1+1;
                arAce(l1,:)=[allsesm_c(idwave(i)),globp+(ref_idm-1)*2*lwave+i];       %actual_reference indices 
            end
            for i=1:lwave
                l2=l2+1;
                arAse(l2,:)=[allsesm_s(idwave(i)),globp+(ref_idm-1)*2*lwave+i+lwave];   %actual_reference indices 
            end
        end
    end
    if l1>0
        ar=vertcat(ar,arAce,arAse); % actual X reference indices
        actp=size(ar,1); % actp: actual number of est. parameters in this session
    end
    
    if ise==lse
        for ian=1:nstses
            for iw=1:lwave
                stc(ian).newcol(iw)=globp+(ian-1)*2*lwave+iw;
                sts(ian).newcol(iw)=globp+(ian-1)*2*lwave+iw+lwave;
            end
        end
        parGS(g.g_stseaspos_Ae(1)).newcol={stc.newcol};
        parGS(g.g_stseaspos_Ae(2)).newcol={sts.newcol};
    end
    globp=globp+2*lwave*nstses;
end


% seasonal variations in station positions N
if parGS(g.g_stseaspos_An(1)).id==1
    idwave=find(cell2mat(parGS(g.g_stseaspos_An(1)).spec)==1);
    lwave=length(idwave);
    nstses=size(stseasname,1);
    
    clear arAcn arAsn ref_idm allsesm_c allsesm_s stc sts
    l1=0;
    l2=0;
    % put antennas and amplitudes into the right columns
    for k = 1:length(antenna.x)
        aname=antenna.name(k,:);
        fq=strcmp(cellstr(aname),cellstr(stseasname)); % find the antenna
        if sum(fq)~=0
            [ref_idm]=find(fq==1);
            clear fq
            allsesm_c=cell2mat(parGS(g.g_stseaspos_An(1)).oldcol(k));
            allsesm_s=cell2mat(parGS(g.g_stseaspos_An(2)).oldcol(k));
            for i=1:lwave
                l1=l1+1;
                arAcn(l1,:)=[allsesm_c(idwave(i)),globp+(ref_idm-1)*2*lwave+i];       %actual_reference indices 
                stc(ref_idm).newcol(i)=globp+(ref_idm-1)*2*lwave+i;
            end
            for i=1:lwave
                l2=l2+1;
                arAsn(l2,:)=[allsesm_s(idwave(i)),globp+(ref_idm-1)*2*lwave+i+lwave];   %actual_reference indices 
                sts(ref_idm).newcol(i)=globp+(ref_idm-1)*2*lwave+i+lwave;
            end
        end
    end
    if l1>0
        ar=vertcat(ar,arAcn,arAsn); % actual X reference indices
        actp=size(ar,1); % actp: actual number of est. parameters in this session
    end
    
    if ise==lse
        for ian=1:nstses
            for iw=1:lwave
                stc(ian).newcol(iw)=globp+(ian-1)*2*lwave+iw;
                sts(ian).newcol(iw)=globp+(ian-1)*2*lwave+iw+lwave;
            end
        end
        parGS(g.g_stseaspos_An(1)).newcol={stc.newcol};
        parGS(g.g_stseaspos_An(2)).newcol={sts.newcol};
    end

    
    globp=globp+2*lwave*nstses;
end

if parGS(g.g_hpole).id==1
    ar(actp+1,:)=[parGS(g.g_hpole).oldcol,globp+1];
    actp=size(ar,1);
    globp=globp+lhpole;
    parGS(g.g_hpole).newcol=globp(end);
end

if parGS(g.g_lpole).id==1
    ar(actp+1,:)=[parGS(g.g_lpole).oldcol,globp+1];
    actp=size(ar,1);
    globp=globp+llpole;
    parGS(g.g_lpole).newcol=globp(end);
end


% APL RgC
if parGS(g.g_aplrg).id==1
    l=1;
    arrg=[];
    % put APL RgC into the right columns
    for k = 1:length(antenna.x)
        aname=antenna.name(k,:);
        fq=strcmp(cellstr(aname),cellstr(rgstname)); % find the antenna
        if sum(fq)~=0
            [ref_idm]=find(fq==1);
            clear fq
            arrg(l,:)=[parGS(g.g_aplrg).oldcol(k),globp+ref_idm];   
            l=l+1;
        end
    end
    if ~isempty(arrg)    
        ar(actp+1:actp+size(arrg,1), :)=arrg;
        actp=size(ar,1);
        parGS(g.g_aplrg).newcol=[globp+1:globp+lrg];
        globp=globp+lrg;        
    end
end



%% INTERNAL

% Love
if parGS(g.g_love).id==1
     for i=1:llove
         ilo=parGS_hl.love.nr(i); % the indices of the Love numbers that will be estimated according to paramGS.m
         ardlove(i,:)= [parGS(g.g_love).oldcol(ilo),globp+i];
         newcol_h(i)=globp+i;
     end
    ar(actp+1:actp+llove,:)=ardlove;
    actp=size(ar,1);
    globp=globp+llove;
    parGS(g.g_love).newcol=newcol_h;
end

% Shida
if parGS(g.g_shida).id==1
     for i=1:lshida
         ish=parGS_hl.shida.nr(i); % the indices of the Shida numbers that will be estimated according to paramGS.m
         ardshida(i,:)= [parGS(g.g_shida).oldcol(ish),globp+i];
         newcol_l(i)=globp+i;
     end
    ar(actp+1:actp+lshida,:)=ardshida;
    actp=size(ar,1);
    globp=globp+lshida;
    parGS(g.g_shida).newcol=newcol_l;
end

% FCN period from solid Earth tides
if parGS(g.g_FCNset).id==1
    ar(actp+1,:)=[parGS(g.g_FCNset).oldcol,globp+1];
    actp=size(ar,1);
    globp=globp+lFCNset;
    parGS(g.g_FCNset).newcol=globp(end);
end


% acceleration of SSB
if parGS(g.g_accSSB).id==1
    oldcol=parGS(g.g_accSSB).oldcol;
    for i=1:length(oldcol)
        aracc(i,:)=[oldcol(i),globp+i];
        newcol(i)=globp+i;
    end
    ar(actp+1:actp+length(oldcol),:)=aracc;
    actp=size(ar,1);
    globp=globp+laccSSB;
    parGS(g.g_accSSB).newcol=newcol;
    clear newcol
end



if parGS(g.g_gamma).id==1
    ar(actp+1,:)=[parGS(g.g_gamma).oldcol,globp+1];
    actp=size(ar,1);
    globp=globp+lgamma;
    parGS(g.g_gamma).newcol=globp(end);
end




% tidal ERP variations
if parGS(g.g_tidpm).id==1
    ar(actp+1:actp+ltidpm , :) = [parGS(g.g_tidpm).oldcol', (globp+1 : globp+ltidpm)' ];
    actp=size(ar,1);
    parGS(g.g_tidpm).newcol= globp +1 :globp + ltidpm;
    globp=globp+ltidpm;
    
end    
    
if parGS(g.g_tidut).id==1
    ar(actp+1:actp+ltidut, :) = [parGS(g.g_tidut).oldcol' ,(globp+1 : globp+ltidut)'];
    actp=size(ar,1);
    parGS(g.g_tidut).newcol= globp +1 :globp + ltidut;
    globp=globp+ltidut;
end  

%%