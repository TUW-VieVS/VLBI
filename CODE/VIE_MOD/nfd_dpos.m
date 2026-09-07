% ************************************************************************
%   Description:
%   function to compute the partial derivatives of the near-field delay
%   model with respect to the satellite position in satellite fixed frame
%   and with respect to the station positions
%
%   Reference:
%
%   Input:
%       'crsStation1'             (3,1)           station coordinates of station 1 (CRS)
%       'crsStation2'             (3,1)           station coordinates of station 2 (CRS)
%       't2c'                     (3,1)           transformation matrix 
%       'crfScPos'                (3,1)           satellite position (CRS)
%       'crfScVel'                (3,1)           satellite velocity (CRS)
%       'v2'                      (3,1)           velocity of station 2 (CRS)
%  
%
%   Output:
%       'psat_rsw'              (3,1)           partial derivative of delay wrt satellite position in RSW-frame
%       'psat_gcrf'             (3,1)           partial derivative of delay wrt satellite position in GCRF-frame
%       'psat_ntw'              (3,1)           partial derivative of delay wrt satellite position in NTW-frame
%       'psat_trf'              (3,3)           partial derivative of delay wrt satellite position in TRF-frame
%       'ps1'                   (3,1)           partial derivative of delay wrt to coord of station1 in TRF
%       'ps2'                   (3,1)           partial derivative of delay wrt to coord of station2 in TRF
%
%   External calls:
%
%   Coded for VieVS:
%   10 Aug 2026 by H. Wolf 
%
%   Revision:
%
% ************************************************************************
function [psat_gcrf, psat_rsw, psat_ntw, psat_trf, ps1, ps2] = nfd_dpos(crsStation1, crsStation2, t2c, crfScPos, crfScVel, v2)
    global c

    L1  = crfScPos' - crsStation1';
    L2  = crfScPos' - crsStation2';

    nL1 = norm(L1);
    nL2 = norm(L2);

    % analytical PD in GCRF:
    dudws_part1 = (crfScPos - crsStation2) ./ nL2   -  (crfScPos - crsStation1) ./ nL1; 
    dudws_part2 = v2 - (crfScPos - crsStation1)./nL1 .* 1/nL2 .* (L2*v2)  +  nL1 .* (crfScPos - crsStation2) * 1/nL2^3 * (L2*v2)  - (nL1 * 1/nL2) .* v2; 
    psat_gcrf = dudws_part1./ c - dudws_part2./ c^2; % [sec/m] % PD of delay time du w.r.t. satellite pos. in GCRF [sec/m], (3x1 vector)
    psat_gcrf = psat_gcrf * c; % * 100 / 100; % Unit conversion: [sec/m] => [cm/cm] = []; => estimates will be in [cm]

    % Rotation of PD to RSW system:
    [~, ~, transmatRSW] = rv2rsw(crfScPos, crfScVel);
    psat_rsw = transmatRSW*psat_gcrf;

    % Rotation of PD to NTW system:
    [~, ~, transmatNTW] = rv2ntw(crfScPos, crfScVel);
    psat_ntw = transmatNTW*psat_gcrf;
    
    % Rotation of PD to TRF system:
    psat_trf = t2c' * psat_gcrf;

    %partial derivative of du0 w.r.t. station coordinates (in TRF!)
    ps1 = +t2c'*(L1'/norm(L1));  % partial derivative wrt to station 1, unit: [] => estimates will be in [cm]
    ps2 = -t2c'*(L2'/norm(L2));  % partial derivative wrt to station 2, unit: [] => estimates will be in [cm]
    
    % % Nummerical derivation of du0 w.r.t. ws1 (du0/dws1)
    s = 0.5;
    crfScPosp = crfScPos + [s; 0; 0]; % plus s m in ws1
    L1p     = crfScPosp' - crsStation1';
    L2p     = crfScPosp' - crsStation2';
    crfScPosm = crfScPos - [s; 0; 0]; % minus s m in ws1
    L1m     = crfScPosm' - crsStation1';
    L2m     = crfScPosm' - crsStation2';
    dudws1_num     = ( ((norm(L2p) - norm(L1p))) - ((norm(L2m) - norm(L1m))) ) / (2*s);

    crfScPosp = crfScPos + [0; s; 0]; % plus s m in ws2
    L1p     = crfScPosp' - crsStation1';
    L2p     = crfScPosp' - crsStation2';
    crfScPosm = crfScPos - [0; s; 0]; % minus s m in ws2
    L1m     = crfScPosm' - crsStation1';
    L2m     = crfScPosm' - crsStation2';
    dudws2_num     = ( ((norm(L2p) - norm(L1p))) - ((norm(L2m) - norm(L1m))) ) / (2*s);

    crfScPosp = crfScPos + [0; 0; s]; % plus s m in ws3
    L1p     = crfScPosp' - crsStation1';
    L2p     = crfScPosp' - crsStation2';
    crfScPosm = crfScPos - [0; 0; s]; % minus s m in ws3
    L1m     = crfScPosm' - crsStation1';
    L2m     = crfScPosm' - crsStation2';
    dudws3_num     = ( ((norm(L2p) - norm(L1p))) - ((norm(L2m) - norm(L1m))) ) / (2*s);

    npdSatPosGCRF = [dudws1_num; dudws2_num; dudws3_num];

    % check if analytical and numerical solutions are equal
    gcrDfDiffPd = psat_gcrf - npdSatPosGCRF; 
    if gcrDfDiffPd(1,1)> 0.00001 || gcrDfDiffPd(2,1)> 0.00001 || gcrDfDiffPd(3,1)> 0.00001
       disp('Analytical and numerical derivatives of tau wrt satellite position are not equal!')
    end
end