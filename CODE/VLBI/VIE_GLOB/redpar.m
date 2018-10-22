% ************************************************************************
%   Description:
%   Find position of all parameters in this session in the old N matrix,
%   which will be reduced
%
%
%   Input:
%      parGS               name of the parameter and if it will be
%                          estimated (1), fixed (0) or reduced (2),
%      antfo_old           position of the reduced antenna in the old N matrix
%      reducsou_old        position of the reduced source in the old N matrix
%
%   Output:     
%      redpos_all          position of all parameters, which will be reduced
%      redcoor             position of antennas, which will be reduced
%      redsou              position of sources, which will be reduced
%                          (positions antenna and sources are needed for backward solution)
%
%   Coded for VieVS: 
%   21 Jul 2011 by Hana Spicakova
%
%   Revision: 
%%



function [redpos_all redcoor redsou] = redpar(parGS,antfo_old,reducsou_old)
[g] = globind(parGS);

redpos=[]; redcoor=[]; redsou=[];

% position of the parameters, which will be reduced
ipar = find([parGS.id]==2); 
redpos = [parGS(ipar).oldcol];

% stations with few observations, which will be session-wise reduced
if parGS(g.g_coord(1)).id==1
    redcoorx=parGS(g.g_coord(1)).oldcol(antfo_old);
    redcoory=parGS(g.g_coord(2)).oldcol(antfo_old);
    redcoorz=parGS(g.g_coord(3)).oldcol(antfo_old);
    redcoor=[redcoorx, redcoory, redcoorz];
end

% sources, which will be session-wise reduced (i.e., special handling sources)
if parGS(g.g_srade(1)).id==1
    redsoura=parGS(g.g_srade(1)).oldcol(reducsou_old);
    redsoude=parGS(g.g_srade(2)).oldcol(reducsou_old);
    redsou=[redsoura, redsoude];
end

redpos_all = [redpos redcoor redsou];


