% ************************************************************************
%   Description:
%   The variable g stores the name of the parameters and its order within
%   subgroups.
%
%   Input:										
%      parGS               name of the parameter and if it will be
%                          estimated (1), fixed (0) or reduced (2)
%   Output:                
%      g                   name of subroups and order of the parameters
%
%   External calls: 	
%      -               					    											
%       
%   Coded for VieVS: 
%   01 Aug 2010 by Hana Spicakova
%
%   Revision: 
%   20 Jun 2012 by Hana Krásná
%               Love and Shida numbers (solid Earth tide), FCN, velocity of sources, acceleration of SSB, seasonal
%               variations in station positions (annual, semi-annual), pole
%               tide Love and Shida numbers
%  25 Sep 2012  by Hana Krásná: relativistic parameter gamma
%   04 Oct 2013 by Hana Krasna: antenna axis offset added
%   06 Dec 2013 by Hana Krasna: apl regression coefficient added
%   25 Aug 2015 by Sigrid Boehm: tidal ERP coefficients added 
%   21 Oct 2015 by Hana Krasna: seasonal variations in station positions (annual, semi-annual), pole
%               tide Love and Shida numbers, apl regression coefficients added
% ************************************************************************

function [g] = globind(parGS)

for i=1:length(parGS)
    if strcmp(parGS(i).name,'pwck')==1
        g_pwck=i;
    
    elseif strcmp(parGS(i).name,'zwd')==1
        g_zwd=i;
    elseif strcmp(parGS(i).name,'ngr')==1
        g_ngr=i;
    
    elseif strcmp(parGS(i).name,'ant_x')==1
        g_coordx=i;
    
    elseif strcmp(parGS(i).name,'vx')==1
        g_velx=i;
    
    elseif strcmp(parGS(i).name,'sra')==1
        g_sra=i;
        
    elseif strcmp(parGS(i).name,'xpol')==1
        g_xpol=i;
        
    elseif strcmp(parGS(i).name,'ao')==1
        g_ao=i;

    elseif strcmp(parGS(i).name,'stseaspos_Acr')==1
        g_stseaspos_Acr=i;
    elseif strcmp(parGS(i).name,'stseaspos_Ace')==1
        g_stseaspos_Ace=i;
    elseif strcmp(parGS(i).name,'stseaspos_Acn')==1
        g_stseaspos_Acn=i;
   
    elseif strcmp(parGS(i).name,'hpole')==1
        g_hpole=i;
    elseif strcmp(parGS(i).name,'lpole')==1
        g_lpole=i;
   
    elseif strcmp(parGS(i).name,'aplrg')==1
        g_aplrg=i;
    elseif strcmp(parGS(i).name,'tidpm')==1
        g_tidpm=i;
    elseif strcmp(parGS(i).name,'tidut')==1
        g_tidut=i;   
        
    elseif strcmp(parGS(i).name,'love')==1
        g_love=i;
    elseif strcmp(parGS(i).name,'shida')==1
        g_shida=i;
    elseif strcmp(parGS(i).name,'FCNset')==1
        g_FCNset=i;
    
    elseif strcmp(parGS(i).name,'accSSB')==1
        g_accSSB=i;
    
    elseif strcmp(parGS(i).name,'svra')==1
        g_svra=i;

    elseif strcmp(parGS(i).name,'gamma')==1
        g_gamma=i;

    elseif strcmp(parGS(i).name,'bdco')==1
        g_bdco=i;

    end
end

g.g_clk   = [g_pwck,g_pwck+1];
g.g_zwd   = [g_zwd];
g.g_tgr   = [g_ngr, g_ngr+1];
g.g_coord = [g_coordx,g_coordx+1,g_coordx+2];
g.g_vel   = [g_velx, g_velx+1, g_velx+2];
g.g_srade = [g_sra, g_sra+1];
g.g_eop   = [g_xpol:1:g_xpol+4];
g.g_ao   = [g_ao];
g.g_stseaspos_Ar=[g_stseaspos_Acr, g_stseaspos_Acr+3];
g.g_stseaspos_Ae=[g_stseaspos_Ace, g_stseaspos_Ace+3];
g.g_stseaspos_An=[g_stseaspos_Acn, g_stseaspos_Acn+3];
g.g_hpole = [g_hpole];
g.g_lpole = [g_lpole];
g.g_aplrg = [g_aplrg];
g.g_tidpm = [g_tidpm];
g.g_tidut = [g_tidut];
g.g_love  = [g_love];
g.g_shida = [g_shida];
g.g_FCNset = [g_FCNset];
g.g_accSSB= [g_accSSB];
g.g_svrade= [g_svra, g_svra+1];
g.g_gamma = [g_gamma];
g.g_bdco = [g_bdco];



