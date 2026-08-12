% ************************************************************************
%   Description:
%   function to copute the partial derivatives of the time delay (for
%   satellites) w.r.t. the Earth Orientation Parameters
%
%   Reference: 
%
%   Input:										
%       'crfScPos'        (1,3)           position of satellite in CRS
%       'trsStation1'     (3,1)           position of station 1 in TRS
%       'trsStation2'     (3,1)           position of station 2 in TRS 
%       'crsStation1'     (1,3)           position of station 1 in CRS
%       'crsStation2'     (1,3)           position of station 2 in CRS
%       'dQdxp'           (3,3)           partial derivative of t2c w.r.t. pole x
%       'dQdyp'           (3,3)           partial derivative of t2c w.r.t. pole y
%       'dQdut'           (3,3)           partial derivative of t2c w.r.t. dut1 
%       'dQddX'           (3,3)           partial derivative of t2c w.r.t. celestial X 
%       'dQddY'           (3,3)           partial derivative of t2c w.r.t. celestial X
%       'v2'              (1,3)           velocity of station 2 
%   Output:
%       'pdeop'          (5,1)            partial derivatives of near-field time delay w.r.t. EOPs  %[sec/rad]          
%      
%   External calls: 	
%   
%   Coded for VieVS: 
%   07 Aug 2025 by Helene Wolf
%
%   Revision: 
%
% ************************************************************************

function [pdeop] = nfd_deop(crfScPos, trsStation1, trsStation2, crsStation1, crsStation2, dQdxp, dQdyp, dQdut, dQddX, dQddY, v2)
    global c omega  

    % vector station-spacecraft
    L1  = crfScPos' - crsStation1'; 
    L2  = crfScPos' - crsStation2';
    nL1 = norm(L1);
    nL2 = norm(L2);

    %  partial deriv. of station velocity w.r.t. EOP:
    dv2dut = dQdut *[-omega*trsStation2(2);omega*trsStation2(1);0]; 
    dv2dxp = dQdxp *[-omega*trsStation2(2);omega*trsStation2(1);0];
    dv2dyp = dQdyp *[-omega*trsStation2(2);omega*trsStation2(1);0]; 
    dv2ddX = dQddX *[-omega*trsStation2(2);omega*trsStation2(1);0]; 
    dv2ddY = dQddY *[-omega*trsStation2(2);omega*trsStation2(1);0]; 

    % partial deriv. of station position w.r.t. EOP:
    dw1dut = dQdut * -trsStation1';
    dw2dut = dQdut * -trsStation2'; 
    dw1dxp = dQdxp * -trsStation1';
    dw2dxp = dQdxp * -trsStation2';
    dw1dyp = dQdyp * -trsStation1';
    dw2dyp = dQdyp * -trsStation2';
    dw1ddX = dQddX * -trsStation1';
    dw2ddX = dQddX * -trsStation2'; 
    dw1ddY = dQddY * -trsStation1';
    dw2ddY = dQddY * -trsStation2'; 

    nL1_dut = L1/norm(L1) * dw1dut ; 
    nL2_dut = L2/norm(L2) * dw2dut ; 
    nL1_dxp = L1/norm(L1) * dw1dxp ; 
    nL2_dxp = L2/norm(L2) * dw2dxp ; 
    nL1_dyp = L1/norm(L1) * dw1dyp ; 
    nL2_dyp = L2/norm(L2) * dw2dyp ; 
    nL1_ddX = L1/norm(L1) * dw1ddX ; 
    nL2_ddX = L2/norm(L2) * dw2ddX ; 
    nL1_ddY = L1/norm(L1) * dw1ddY ; 
    nL2_ddY = L2/norm(L2) * dw2ddY ;

    %xpol
    A_dxp = nL2_dxp - nL1_dxp; 
    B_dxp = dw2dxp' * v2 + L2 * dv2dxp; 
    C1_dxp = (nL1_dxp * L2 * v2 + nL1 * dw2dxp' * v2 + nL1 * L2 * dv2dxp) / nL2; 
    C2_dxp = (nL1 * L2 * v2 * nL2_dxp)/ nL2^2; 

    %ypol
    A_dyp = nL2_dyp - nL1_dyp; 
    B_dyp = dw2dyp' * v2 + L2 * dv2dyp; 
    C1_dyp = (nL1_dyp * L2 * v2 + nL1 * dw2dyp' * v2 + nL1 * L2 * dv2dyp) / nL2; 
    C2_dyp = (nL1 * L2 * v2 * nL2_dyp)/ nL2^2; 

    %ut1
    A_dut = nL2_dut - nL1_dut; 
    B_dut = dw2dut' * v2 + L2 * dv2dut; 
    C1_dut = (nL1_dut * L2 * v2 + nL1 * dw2dut' * v2 + nL1 * L2 * dv2dut) / nL2; 
    C2_dut = (nL1 * L2 * v2 * nL2_dut)/ nL2^2; 

    %dX
    A_ddX = nL2_ddX - nL1_ddX; 
    B_ddX = dw2ddX' * v2 + L2 * dv2ddX;  
    C1_ddX = (nL1_ddX * L2 * v2 + nL1 * dw2ddX' * v2 + nL1 * L2 * dv2ddX) / nL2; 
    C2_ddX = (nL1 * L2 * v2 * nL2_ddX)/ nL2^2; 
    
    %dY
    A_ddY = nL2_ddY - nL1_ddY; 
    B_ddY = dw2ddY' * v2 + L2 * dv2ddY; 
    C1_ddY = (nL1_ddY * L2 * v2 + nL1 * dw2ddY' * v2 + nL1 * L2 * dv2ddY) / nL2; 
    C2_ddY = (nL1 * L2 * v2 * nL2_ddY)/ nL2^2; 
 
    pxp = (1/c) * A_dxp - 1/c^2*(B_dxp - C1_dxp + C2_dxp); %[sec/rad]
    pyp = (1/c) * A_dyp - 1/c^2*(B_dyp - C1_dyp + C2_dyp); %[sec/rad]
    put = (1/c) * A_dut - 1/c^2*(B_dut - C1_dut + C2_dut); %[sec/rad]
    pdX = (1/c) * A_ddX - 1/c^2*(B_ddX - C1_ddX + C2_ddX); %[sec/rad]
    pdY = (1/c) * A_ddY - 1/c^2*(B_ddY - C1_ddY + C2_ddY); %[sec/rad]

    pdeop = [pxp, pyp, put, pdX, pdY];
end