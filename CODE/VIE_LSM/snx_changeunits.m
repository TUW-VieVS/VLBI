% ************************************************************************
%   Description:
%   function which change the units in the N matrix and b vector
%   In VieVS we are computing in cm and mas. Sinex requires m, mas, ms, rad
%   ...
%
%   Reference: 
%
%   Input:	
%   N_sinex          N matrix (units of VieVS)
%   b_sinex          b vector (units of VieVS)
%   col_sinex        number of columns
%   outsnx           1/0
%
%   Output:
%   N_sinex          N matrix (units of Sinex)
%   b_sinex          b vector (units of Sinex)
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   01 May 2011 by Hana Spicakova
%
%	Revision:
%   2025-02-02 by H.Wolf: added orbital elements to be written in the sinex file
%
% ************************************************************************ 

function [N_sinex, b_sinex] = snx_changeunits(N_sinex,b_sinex,col_sinex,outsnx)  

    c_zwd=[col_sinex.zwd.col];
    c_tgr=[[col_sinex.ngr.col] [col_sinex.egr.col]];
    c_sou=[col_sinex.ra col_sinex.de];
    c_eop4=[col_sinex.xp col_sinex.yp col_sinex.dX col_sinex.dY];
    c_dut1=[col_sinex.dut1];
    c_KepEle_a=[col_sinex.orb_sma];
    c_KepEle_5=[col_sinex.orb_ecc col_sinex.orb_inc col_sinex.orb_raan col_sinex.orb_argp col_sinex.orb_argl];
    c_SRP = [col_sinex.srp_D0 col_sinex.srp_Y0 col_sinex.srp_B0 col_sinex.srp_DC col_sinex.srp_YC col_sinex.srp_BC col_sinex.srp_DS col_sinex.srp_YS col_sinex.srp_BS ]; %nm/s2
    
    c_xyz = [];
    if isfield(col_sinex, 'coorx') && isfield(col_sinex, 'coory') && isfield(col_sinex, 'coorz')
        c_xyz = [col_sinex.coorx col_sinex.coory col_sinex.coorz];
    elseif isfield(col_sinex, 'coorx_sat') && isfield(col_sinex, 'coory_sat') && isfield(col_sinex, 'coorz_sat')
        c_xyz = [col_sinex.coorx_sat col_sinex.coory_sat col_sinex.coorz_sat];
        % check if also stations from quasar coordinates are estimated
        if isfield(col_sinex, 'coorx_qu') && isfield(col_sinex, 'coory_qu') && isfield(col_sinex, 'coorz_qu')
            c_xyz = [c_xyz col_sinex.coorx_qu col_sinex.coory_qu col_sinex.coorz_qu];
        end
    elseif isfield(col_sinex, 'coorx_qu') && isfield(col_sinex, 'coory_qu') && isfield(col_sinex, 'coorz_qu')
        c_xyz = [col_sinex.coorx_qu col_sinex.coory_qu col_sinex.coorz_qu];
    end


    % zwd
    if outsnx.zwd==1
        N_sinex(c_zwd,c_zwd) = N_sinex(c_zwd,c_zwd).*10000; %1/cm^2 --> 1/m^2
        if outsnx.tgr==1
            N_sinex(c_zwd,c_tgr) = N_sinex(c_zwd,c_tgr).*10000; %1/cm^2 --> 1/m^2
        end
        if outsnx.sou==1
            N_sinex(c_zwd,c_sou) = N_sinex(c_zwd,c_sou).*(100*1000*3600*180/pi); %1/mas.cm --> 1/rad.m
        end
        if outsnx.xyz==1
            N_sinex(c_zwd,c_xyz) = N_sinex(c_zwd,c_xyz).*10000; %1/cm^2 --> 1/m^2
        end
        if outsnx.eop==1
            N_sinex(c_zwd,c_eop4) = N_sinex(c_zwd,c_eop4).*100; %1/cm.mas --> 1/m.mas
            N_sinex(c_zwd,c_dut1) = N_sinex(c_zwd,c_dut1).*100.*15; %1/cm.mas --> 1/m.ms
        end
        if outsnx.orb ==1
            N_sinex(c_zwd,c_KepEle_a) = N_sinex(c_zwd,c_KepEle_a).*10000; %1/cm^2 --> 1/m^2
			N_sinex(c_zwd,c_KepEle_5) = N_sinex(c_zwd,c_KepEle_5).*100; %1/cm --> 1/m or 1/mas.cm --> 1/mas.m
            N_sinex(c_zwd,c_SRP) = N_sinex(c_zwd,c_SRP).*100; %1/(nm/s2*cm) --> 1/(nm/s2 *m)
        end
        b_sinex(c_zwd) = b_sinex(c_zwd).*100; %1/cm --> 1/m
    end

    % tgr
    if outsnx.tgr==1
        N_sinex(c_tgr,c_tgr) = N_sinex(c_tgr,c_tgr).*10000; %1/cm^2 --> 1/m^2
        if outsnx.zwd==1
           N_sinex(c_tgr,c_zwd) = N_sinex(c_tgr,c_zwd).*10000; %1/cm^2 --> 1/m^2
        end
        if outsnx.sou==1
           N_sinex(c_tgr,c_sou) = N_sinex(c_tgr,c_sou).*(100*1000*3600*180/pi); %1/cm.mas --> 1/m.rad
        end
        if outsnx.xyz==1
           N_sinex(c_tgr,c_xyz) = N_sinex(c_tgr,c_xyz).*10000; %1/cm^2 --> 1/m^2
        end
        if outsnx.eop==1
            N_sinex(c_tgr,c_eop4) = N_sinex(c_tgr,c_eop4).*100; %1/cm.mas --> 1/m.mas
            N_sinex(c_tgr,c_dut1) = N_sinex(c_tgr,c_dut1).*100.*15; %1/cm.mas --> 1/m.ms
        end
        if outsnx.orb ==1
            N_sinex(c_tgr,c_KepEle_a) = N_sinex(c_tgr,c_KepEle_a).*10000; %1/cm^2 --> 1/m^2
			N_sinex(c_tgr,c_KepEle_5) = N_sinex(c_tgr,c_KepEle_5).*100; %1/cm --> 1/m or 1/mas.cm --> 1/mas.m
            N_sinex(c_tgr,c_SRP) = N_sinex(c_tgr,c_SRP).*100; %1/(nm/s2*cm) --> 1/(nm/s2 *m)
        end
        b_sinex(c_tgr) = b_sinex(c_tgr).*100; %1/cm --> 1/m
    end

    % sou
    if outsnx.sou==1
        N_sinex(c_sou,c_sou) = N_sinex(c_sou,c_sou).*(1000*3600*180/pi)^2; %1/mas^2 --> 1/rad^2
        if outsnx.zwd==1
           N_sinex(c_sou,c_zwd) = N_sinex(c_sou,c_zwd).*(100*1000*3600*180/pi); %1/cm.mas --> 1/m.rad
        end
        if outsnx.tgr==1
           N_sinex(c_sou,c_tgr) = N_sinex(c_sou,c_tgr).*(100*1000*3600*180/pi); %1/cm.mas --> 1/m.rad
        end
        if outsnx.xyz==1
           N_sinex(c_sou,c_xyz) = N_sinex(c_sou,c_xyz).*(100*1000*3600*180/pi); %1/cm.mas --> 1/m.rad
        end
        if outsnx.eop==1
           N_sinex(c_sou,c_eop4) = N_sinex(c_sou,c_eop4).*(1000*3600*180/pi); %1/mas^2 --> 1/mas.rad
           N_sinex(c_sou,c_dut1) = N_sinex(c_sou,c_dut1).*(1000*3600*180/pi*15); %1/mas^2 --> 1/ms.rad
        end
        if outsnx.orb ==1
            N_sinex(c_sou,c_KepEle_a) = N_sinex(c_sou,c_KepEle_a).*100*(1000*3600*180/pi); %1/cm.mas --> 1/m.rad
			N_sinex(c_sou,c_KepEle_5) = N_sinex(c_sou,c_KepEle_5).*(1000*3600*180/pi); %1/mas --> 1/rad or 1/mas^2 --> 1/mas.rad
            N_sinex(c_sou,c_SRP) = N_sinex(c_sou,c_SRP).*(1000*3600*180/pi); %1/(nm/s2*mas) --> 1/(nm/s2 *rad)
        end
        b_sinex(c_sou) = b_sinex(c_sou).*(1000*3600*180/pi); %1/mas --> 1/rad
    end
    
    if outsnx.xyz==1
        N_sinex(c_xyz,c_xyz) = N_sinex(c_xyz,c_xyz).*10000; %1/cm^2 --> 1/m^2
        if outsnx.zwd==1
           N_sinex(c_xyz,c_zwd) = N_sinex(c_xyz,c_zwd).*10000; %1/cm^2 --> 1/m^2
        end
        if outsnx.tgr==1
           N_sinex(c_xyz,c_tgr) = N_sinex(c_xyz,c_tgr).*10000; %1/cm^2 --> 1/m^2
        end
        if outsnx.sou==1
            N_sinex(c_xyz,c_sou) = N_sinex(c_xyz,c_sou).*(100*1000*3600*180/pi); %1/mas.cm --> 1/rad.m
        end
        if outsnx.eop==1
            N_sinex(c_xyz,c_eop4) = N_sinex(c_xyz,c_eop4).*100; %1/cm.mas --> 1/m.mas
            N_sinex(c_xyz,c_dut1) = N_sinex(c_xyz,c_dut1).*100.*15; %1/cm.mas --> 1/m.ms
        end 
        if outsnx.orb ==1
            N_sinex(c_xyz,c_KepEle_a) = N_sinex(c_xyz,c_KepEle_a).*10000; %1/cm^2 --> 1/m^2
			N_sinex(c_xyz,c_KepEle_5) = N_sinex(c_xyz,c_KepEle_5).*100; %1/cm --> 1/m or 1/mas.cm --> 1/mas.m
            N_sinex(c_xyz,c_SRP) = N_sinex(c_xyz,c_SRP).*100; %1/(nm/s2*cm) --> 1/(nm/s2 *m)
        end
        b_sinex(c_xyz) = b_sinex(c_xyz).*100; %1/cm --> 1/m
    end
        
    if outsnx.eop==1
        N_sinex(c_dut1,c_dut1) = N_sinex(c_dut1,c_dut1).*225; %1/mas^2 --> 1/ms^2
        N_sinex(c_dut1,c_eop4) = N_sinex(c_dut1,c_eop4).*15; %1/mas^2 --> 1/mas.ms
        N_sinex(c_eop4,c_dut1) = N_sinex(c_eop4,c_dut1).*15; %1/mas^2 --> 1/mas.ms
        if outsnx.zwd==1
            N_sinex(c_eop4,c_zwd) = N_sinex(c_eop4,c_zwd).*100; %1/cm.mas --> 1/m.mas
            N_sinex(c_dut1,c_zwd) = N_sinex(c_dut1,c_zwd).*100.*15; %1/cm.mas --> 1/m.ms
        end
        if outsnx.tgr==1
            N_sinex(c_eop4,c_tgr) = N_sinex(c_eop4,c_tgr).*100; %1/cm.mas --> 1/m.mas
            N_sinex(c_dut1,c_tgr) = N_sinex(c_dut1,c_tgr).*100.*15; %1/cm.mas --> 1/m.ms
        end
        if outsnx.sou==1
            N_sinex(c_eop4,c_sou) = N_sinex(c_eop4,c_sou).*(1000*3600*180/pi); %1/mas^2 --> 1/mas.rad
            N_sinex(c_dut1,c_sou) = N_sinex(c_dut1,c_sou).*(1000*3600*180/pi*15); %1/mas^2 --> 1/ms.rad
        end
        if outsnx.xyz==1
            N_sinex(c_eop4,c_xyz) = N_sinex(c_eop4,c_xyz).*100; %1/cm.mas --> 1/m.mas
            N_sinex(c_dut1,c_xyz) = N_sinex(c_dut1,c_xyz).*100.*15; %1/cm.mas --> 1/m.ms
        end
        if outsnx.orb ==1
            N_sinex(c_eop4,c_KepEle_a) = N_sinex(c_eop4,c_KepEle_a).*100; %1/cm.mas --> 1/m.mas
            N_sinex(c_dut1,c_KepEle_a) = N_sinex(c_dut1,c_KepEle_a).*100*15; %1/cm.mas --> 1/m.ms
			N_sinex(c_dut1,c_KepEle_5) = N_sinex(c_dut1,c_KepEle_5).*15; %1/mas --> 1/ms or 1/mas^2 --> 1/mas.ms
            N_sinex(c_dut1,c_SRP) = N_sinex(c_dut1,c_SRP).*15; %1/(nm/s2*mas) --> 1/(nm/s2 *ms)
        end
        b_sinex(c_dut1) = b_sinex(c_dut1).*15; %1/mas --> 1/ms
    end

    if outsnx.orb==1
        %KepEle1 -> a -> cm zu m
        N_sinex(c_KepEle_a,c_KepEle_a) = N_sinex(c_KepEle_a,c_KepEle_a).*10000; %1/cm^2 --> 1/m^2 
		N_sinex(c_KepEle_a,c_KepEle_5) = N_sinex(c_KepEle_a,c_KepEle_5).*100; %1/([].cm) --> 1/([].m) or 1/(mas.cm) --> 1/(mas.m)
        N_sinex(c_KepEle_5,c_KepEle_a) = N_sinex(c_KepEle_5,c_KepEle_a).*100; %1/cm --> 1/m or 1/mas.cm --> 1/mas.m ,

        N_sinex(c_KepEle_a,c_SRP) = N_sinex(c_KepEle_a,c_SRP).*100; %1/(nm/s2*cm) --> 1/(nm/s2 *m)
        N_sinex(c_SRP,c_KepEle_a) = N_sinex(c_SRP,c_KepEle_a).*100; %1/(nm/s2*cm) --> 1/(nm/s2 *m)

	    if outsnx.zwd==1
            N_sinex(c_KepEle_a,c_zwd) = N_sinex(c_KepEle_a,c_zwd).*10000; %1/cm^2 --> 1/m^2
            N_sinex(c_KepEle_5,c_zwd) = N_sinex(c_KepEle_5,c_zwd).*100; %1/cm --> 1/m or 1/mas.cm --> 1/mas.m
            N_sinex(c_SRP,c_zwd) = N_sinex(c_SRP,c_zwd).*100; %1/(cm nm/s2) --> 1/(m nm/s2)
        end
        if outsnx.tgr==1
            N_sinex(c_KepEle_a,c_tgr) = N_sinex(c_KepEle_a,c_tgr).*10000; %1/cm^2 --> 1/m^2
            N_sinex(c_KepEle_5,c_tgr) = N_sinex(c_KepEle_5,c_tgr).*100; %1/cm --> 1/m or 1/mas.cm --> 1/mas.m
            N_sinex(c_SRP,c_tgr) = N_sinex(c_SRP,c_tgr).*100; %1/(cm nm/s2) --> 1/(m nm/s2)
        end
        if outsnx.sou==1
            N_sinex(c_KepEle_a,c_sou) = N_sinex(c_KepEle_a,c_sou).*(1000*3600*180/pi)*100; %1/cm.mas --> 1/m.rad
       		N_sinex(c_KepEle_5,c_sou) = N_sinex(c_KepEle_5,c_sou).*(1000*3600*180/pi); %1/mas --> 1/rad or 1/mas^2 --> 1/mas.rad
            N_sinex(c_SRP,c_sou) = N_sinex(c_SRP,c_sou).*(1000*3600*180/pi); %1/(nm/s2*mas) --> 1/(nm/s2 *rad)
        end
        if outsnx.xyz==1
            N_sinex(c_KepEle_a,c_xyz) = N_sinex(c_KepEle_a,c_xyz).*10000; %1/cm^2 --> 1/m^2
	    	N_sinex(c_KepEle_5,c_xyz) = N_sinex(c_KepEle_5,c_xyz).*100; %1/cm --> 1/m or 1/mas.cm --> 1/mas.m
            N_sinex(c_SRP,c_xyz) = N_sinex(c_SRP,c_xyz).*100; %1/(cm nm/s2) --> 1/(m nm/s2)
        end
        if outsnx.eop==1
            N_sinex(c_KepEle_a,c_eop4) = N_sinex(c_KepEle_a,c_eop4).*100; %1/cm.mas --> 1/m.mas
            N_sinex(c_KepEle_a,c_dut1) = N_sinex(c_KepEle_a,c_dut1).*100*15; %1/cm.mas --> 1/m.ms
            N_sinex(c_KepEle_5,c_dut1) = N_sinex(c_KepEle_5,c_dut1).*15; %1/mas --> 1/ms or 1/mas^2 --> 1/mas.ms
            N_sinex(c_SRP,c_dut1) = N_sinex(c_SRP,c_dut1).*15; %1/(mas nm/s2) --> 1/(ms nm/s2) 
        end 
        b_sinex(c_KepEle_a) = b_sinex(c_KepEle_a).*100; %1/cm --> 1/m
    end
end    
