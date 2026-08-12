% ************************************************************************
%   Description:
%   function which finds columns for reduction and columns which will be written
%   into the SINEX file
%
%   Reference: 
%
%   Input:										
%       x_           information about the columns
%       outsnx       1/0
%
%   Output:
%       col_est      vector with the parameters going into SINEX file
%       col_red      vector with the parameters which will be reduced from
%                    the N_sinex
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   01 May 2011 by Hana Spicakova
%
%   Changes: 
%   15 May 2025 by H. Wolf: changes for different possibilities of estimation of station  coordinates (all, quasar only, ..) 
%
% ************************************************************************ 


function [col_red, col_est] = snx_split(x_,outsnx, parameter)


    % clock parameters will be always reduced, in this version they cannot
    % be written into SINEX file
    col_red = [[x_.pwclk.col] [x_.rqclk.col] [x_.bdclko.col]];
    
    % station coordinates cannot be reduced in SINEX file in this version
    if parameter.lsmopt.stc_all == 1
        col_est = [[x_.coorx.col] [x_.coory.col] [x_.coorz.col]];
    elseif parameter.lsmopt.stc_sat == 1
        col_est = [[x_.coorx_sat.col] [x_.coory_sat.col] [x_.coorz_sat.col]];
    elseif parameter.lsmopt.stc_qu == 1
        col_est = [[x_.coorx_qu.col] [x_.coory_qu.col] [x_.coorz_qu.col]];
    elseif parameter.lsmopt.stc_qs
        col_est = [[x_.coorx_sat.col] [x_.coory_sat.col] [x_.coorz_sat.col] [x_.coorx_qu.col] [x_.coory_qu.col] [x_.coorz_qu.col]];
    end      
    
    if outsnx.zwd == 0
        col_red = [col_red [x_.zwd.col]];
    else
        col_est = [col_est [x_.zwd.col]];
    end
    
    if outsnx.tgr == 0
        col_red = [col_red [x_.ngr.col] [x_.egr.col]];
    else
        col_est = [col_est [x_.ngr.col] [x_.egr.col]];
    end
    
    if outsnx.sou ==1
        col_est = [col_est  [x_.col_soura] [x_.col_soude]];
    end
    
    if outsnx.eop == 0
        col_red = [col_red [x_.xpol.col] [x_.ypol.col] [x_.dut1.col] [x_.nutdx.col] [x_.nutdy.col]];
    else
        col_est = [col_est [x_.xpol.col] [x_.ypol.col] [x_.dut1.col] [x_.nutdx.col] [x_.nutdy.col]];
    end
    
    if outsnx.orb==0
       col_red = [col_red [x_.ORB.sma.col] [x_.ORB.ecc.col] [x_.ORB.inc.col] [x_.ORB.raan.col] [x_.ORB.argp.col] [x_.ORB.argl.col] [x_.SRP.D0.col] [x_.SRP.Y0.col] [x_.SRP.B0.col] [x_.SRP.DC.col] [x_.SRP.YC.col] [x_.SRP.BC.col] [x_.SRP.DS.col] [x_.SRP.YS.col] [x_.SRP.BS.col] ];
    else
       col_est = [col_est [x_.ORB.sma.col] [x_.ORB.ecc.col] [x_.ORB.inc.col] [x_.ORB.raan.col] [x_.ORB.argp.col] [x_.ORB.argl.col] [x_.SRP.D0.col] [x_.SRP.Y0.col] [x_.SRP.B0.col] [x_.SRP.DC.col] [x_.SRP.YC.col] [x_.SRP.BC.col] [x_.SRP.DS.col] [x_.SRP.YS.col] [x_.SRP.BS.col] ];
    end

end
