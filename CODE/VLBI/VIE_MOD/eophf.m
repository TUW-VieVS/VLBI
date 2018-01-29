% ************************************************************************
%   Description:
%   Applies high frequency eop variations for
%       - diurnal and subdiurnal ocean tides
%       - Earth's triaxiality
%   Decisions which models are used is done in VIE_SETUP (mod_qu)
% 
%   Reference: 
%   IERS Conventions 2000, Chapter 5 (as of 4 Nov 2009)
%
%   Input:										
%      MJD      (n,1)      TT date mod. jul. date (n,1)
%      parameter           VieVS structure array includes settings
% 
%   Output:
%      CORX     (n,1)      corrections for X-coordinate [rad]
%      CORY     (n,1)      corrections for Y-corrdinate [rad]
%      CORUT    (n,1)      corrections for dUT1 [s]
% 
%   External calls: 	
%   as2rad, eop_eanes, pm_libration, ut_libration 					    											
%       
%   Coded for VieVS: 
%   04 Nov 2009 by Lucia Plank
%
%   Revision:
%   26 May 2010 by Lucia Plank: adapted for new eop_eanes.m
%   13 Aug 2010 by Lucia Plank: libration added
%   07 Feb 2012 by Lucia Plank: Combi_IGG_Bonn model added
% ************************************************************************
function [CORX,CORY,CORUT,parameter] = eophf(TT,parameter)

CORUT = TT*0;
CORX  = TT*0;
CORY  = TT*0;
parameter.vie_mod.eophf =0;
ocf = parameter.vie_mod.eopoc;
% diurnal and semidiurnal ocean tide correction
switch ocf
    case 'none'
        disp('no high frequency EOP model applied')
    case 'Combi_IGG_Bonn'
        disp ('EOP: Combi_IGG_Bonn')
        parameter.vie_mod.eophf =1;
        [ocx,ocy,ocut] = eophf_combi(TT);         %[as,sec]
        
        CORUT = CORUT + ocut;
        CORX  = CORX + as2rad(ocx);
        CORY  = CORY + as2rad(ocy);
        
    otherwise
        parameter.vie_mod.eophf =1;
        [ocx,ocy,ocut] = eop_eanes(TT,ocf);         %[as,sec]
        %[ocx2,ocy2,ocut2] = eop_ocean(TT,ocf);         %[as,sec]

        CORUT = CORUT + ocut;
        CORX  = CORX + as2rad(ocx);
        CORY  = CORY + as2rad(ocy);
end

% correction due to the Earth's triaxiality
% px, py (recommended by the IERS conventions 2003, chap. 5.4.2) 
if parameter.vie_mod.lib_pm == 1
    parameter.vie_mod.eophf =1;
    if strcmp(ocf,'boehm08.dat')||strcmp(ocf,'boehm09.dat')
        disp('THE EFFECT OF EARTHS TRIAXIALITY ON POLAR MOTION IS ALREADY INCLUDED IN YOUR')
        disp('SUBDAILY OCEAN TIDES MODEL')
        disp('pm_gravi is not applied !!!')
        parameter.vie_mod.lib_pm = 0;
    else
        [trix,triy] = pm_libration(TT);   % [rad]
        CORX = CORX + trix;
        CORY = CORY + triy;
    end
end
if parameter.vie_mod.lib_ut == 1
    parameter.vie_mod.eophf =1;
    if strcmp(ocf,'boehm08.dat')||strcmp(ocf,'boehm09.dat')
        disp('THE EFFECT OF EARTHS TRIAXIALITY ON UT1 IS ALREADY INCLUDED IN YOUR')
        disp('SUBDAILY OCEAN TIDES MODEL')
        disp('spinlibV is not applied !!!')
        parameter.vie_mod.lib_ut = 0;
    else
        [triut] = ut_libration(TT);   % [sec]
        CORUT = CORUT + triut;
    end
end

