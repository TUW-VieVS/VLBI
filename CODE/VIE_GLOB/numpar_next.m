% ************************************************************************
%   Description:
%   This function finds the time steps, for which the EOP will be estimated.
%
%
%   Input:										
%      parGS               name of the parameter and if it will be
%                          estimated (1), fixed (0) or reduced (2),
%      lse                 number of sessions
%      ses                 names of the sessions
%      path                path to the LEVEL2 data
%      *'_par_glob.mat'    glob2; LEVEL2 data

%   Output:                
%      lmjd_eop            number of estimated EOPs parametres (from all
%                          sessions)
%      parGS               if EOP will be estimated, time steps for them are added
%
%   External calls: 	
%      globind               					    											
%       
%   Coded for VieVS: 
%   01 Aug 2010 by Hana Spicakova
%
%   Revision: 21 Jun 2012 by Hana Krásná: INTERNAL version with Love and Shida numbers
%               SET, FCN period from SET, velocity of sources, acceleration of SSB, seasonal
%               variations in station positions (annual, semi-annual), pole
%               tide Love and Shida numbers
%  25 Sep 2012  by Hana Krásná: relativistic parameter gamma in INTERNAL
%               version
%  24 Aug 2015  by Sigrid Boehm: tidal ERP coefficients added to INTERNAL
% ************************************************************************


function [lmjd_eop,llove,lshida,lFCNset,laccSSB,lstseaspos,lhpole,llpole,lgamma,...
          ltidpm,ltidut,parGS] = numpar_next(parGS,parGS_hl,lse,ses,path)

[g] = globind(parGS);

for ieop=g.g_eop
    parGS(ieop).mjdi=[];
end
    
for ise=1:lse
    load ([path ses{ise} '_par_glob.mat']);
    
    x_=glob2.x;
    
    iseeop(1).mjd = x_.xpol.mjd;
    iseeop(2).mjd = x_.ypol.mjd;
    iseeop(3).mjd = x_.dut1.mjd;
    iseeop(4).mjd = x_.nutdx.mjd;
    iseeop(5).mjd = x_.nutdy.mjd;
    
    for ieop=g.g_eop
        if parGS(ieop).id==1; 
            parGS(ieop).mjdi(length(parGS(ieop).mjdi)+1 : length(parGS(ieop).mjdi)+length(iseeop(ieop-(g.g_eop(1)-1)).mjd))=iseeop(ieop-(g.g_eop(1)-1)).mjd;
        end
    end
    clear glob2
end


lmjd_eop=0; % number of estimated EOPs parametres (from all sessions)

for ieop=g.g_eop
    if parGS(ieop).id==1
        mjdstep = unique(parGS(ieop).mjdi);
        lmjd_eop=lmjd_eop+length(mjdstep);
        parGS(ieop).mjdstep=mjdstep;
    end
end

lstseaspos_R=0;
if parGS(g.g_stseaspos_Ar(1)).id==1
    n=cell2mat(parGS(g.g_stseaspos_Ar(1)).spec);
    lstseaspos_R=length(find(n==1)) * 2; % cos and sin amplitude
    % has to be multiplicated with the number of station in the main
    % program
end

lstseaspos_E=0;
if parGS(g.g_stseaspos_Ae(1)).id==1
    n=cell2mat(parGS(g.g_stseaspos_Ae(1)).spec);
    lstseaspos_E=length(find(n==1)) * 2; % cos and sin amplitude
    % has to be multiplicated with the number of station in the main
    % program
end


lstseaspos_N=0;
if parGS(g.g_stseaspos_An(1)).id==1
    n=cell2mat(parGS(g.g_stseaspos_An(1)).spec);
    lstseaspos_N=length(find(n==1)) * 2; % cos and sin amplitude
    % has to be multiplicated with the number of station in the main
    % program
end

lstseaspos = lstseaspos_R + lstseaspos_E + lstseaspos_N;


lhpole=0;
if parGS(g.g_hpole).id==1
    lhpole=1;
end

llpole=0;
if parGS(g.g_lpole).id==1
    llpole=1;
end

% tidal ERP variations
ltidpm=0;
if parGS(g.g_tidpm).id==1
    ltidpm = length(parGS(g.g_tidpm).spectid)*2 + length(parGS(g.g_tidpm).spectidret)*2;
end

ltidut=0;
if parGS(g.g_tidut).id==1
    ltidut = length(parGS(g.g_tidut).spectid)*2;
end

%% INTERNAL
llove=0;
if parGS(g.g_love).id==1
     llove=length(parGS_hl.love.nr);
end

lshida=0;
if parGS(g.g_shida).id==1
    lshida=length(parGS_hl.shida.nr);
end

lFCNset=0;
if parGS(g.g_FCNset).id==1
    lFCNset=1;
end

laccSSB=0;
if parGS(g.g_accSSB).id==1
    laccSSB=3;
end



lgamma=0;
if parGS(g.g_gamma).id==1
    lgamma=1;
end
%%