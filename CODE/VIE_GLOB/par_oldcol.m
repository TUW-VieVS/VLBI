% ************************************************************************
%   Description:
%   This function saves numbers of the column/row for the parameters in 
%   the old N matrix (+ mjd for EOP). (i.e. LEVEL2 N matrix)
%  
%
%
%   Input:
%      x_                  glob2.x (LEVEL2 data)
%      parGS               name of the parameter and if it will be
%                          estimated (1), fixed (0) or reduced (2),
%
%   Output:                
%      parGS               nambers of the old columns,
%                          for EOP it stores also mjd
%
%   External calls: 	
%      globind               					    											
%       
%   Coded for VieVS: 
%   01 Aug 2010 by Hana Spicakova
%
%   Revision: 
%   21 Jul 2011 by Hana Spicakova:  save mjd also for zwd and troposheric
%                                  gradients - for backward solution
%   21 Jun 2012 by Hana Krásná: INTERNAL version with Love and Shida numbers
%               SET, FCN period from SET, velocity of sources, acceleration of SSB, seasonal
%               variations in station positions (annual, semi-annual), pole
%               tide Love and Shida numbers
%  25 Sep 2012  by Hana Krásná: relativistic parameter gamma in INTERNAL
%               version
%   04 Oct 2013 by Hana Krasna: antenna axis offset added
%   06 Dec 2013 by Hana Krasna: APL regression coefficient added
%   24 Aug 2015 by Sigrid Boehm: tidal ERP coefficients added to INTERNAL
% ************************************************************************

function [parGS]=par_oldcol(x_,parGS)

[g] = globind(parGS);

parGS(g.g_clk(1)).oldcol = [x_.pwclk.col];
parGS(g.g_clk(2)).oldcol = [x_.rqclk.col];

parGS(g.g_zwd).oldcol = [x_.zwd.col];
parGS(g.g_zwd).mjd = [x_.zwd.mjd];

parGS(g.g_tgr(1)).oldcol = [x_.ngr.col];
parGS(g.g_tgr(2)).oldcol = [x_.egr.col];
parGS(g.g_tgr(1)).mjd = [x_.ngr.mjd];
parGS(g.g_tgr(2)).mjd = [x_.egr.mjd];


parGS(g.g_coord(1)).oldcol = [x_.coorx.col];
parGS(g.g_coord(2)).oldcol = [x_.coory.col];
parGS(g.g_coord(3)).oldcol = [x_.coorz.col];

parGS(g.g_vel(1)).oldcol = x_.col_vx;
parGS(g.g_vel(2)).oldcol = x_.col_vy;
parGS(g.g_vel(3)).oldcol = x_.col_vz;

parGS(g.g_srade(1)).oldcol = x_.col_soura;
parGS(g.g_srade(2)).oldcol = x_.col_soude;

parGS(g.g_eop(1)).oldcol = [x_.xpol.col];
parGS(g.g_eop(1)).mjd = x_.xpol.mjd;

parGS(g.g_eop(2)).oldcol = [x_.ypol.col];
parGS(g.g_eop(2)).mjd = x_.ypol.mjd;

parGS(g.g_eop(3)).oldcol = [x_.dut1.col];
parGS(g.g_eop(3)).mjd = x_.dut1.mjd;

parGS(g.g_eop(4)).oldcol = [x_.nutdx.col];
parGS(g.g_eop(4)).mjd = x_.nutdx.mjd;

parGS(g.g_eop(5)).oldcol = [x_.nutdy.col];
parGS(g.g_eop(5)).mjd = x_.nutdy.mjd;

parGS(g.g_ao).oldcol = x_.col_AO;

parGS(g.g_love).oldcol = x_.col_love;
parGS(g.g_shida).oldcol =x_.col_shida;
parGS(g.g_FCNset).oldcol =x_.col_FCNset;

parGS(g.g_accSSB).oldcol = x_.col_accSSB;

parGS(g.g_svrade(1)).oldcol = x_.col_souvra;
parGS(g.g_svrade(2)).oldcol = x_.col_souvde;

parGS(g.g_stseaspos_Ar(1)).oldcol = {x_.col_ssp_Acr.col};
parGS(g.g_stseaspos_Ar(2)).oldcol = {x_.col_ssp_Asr.col};
parGS(g.g_stseaspos_Ae(1)).oldcol = {x_.col_ssp_Ace.col};
parGS(g.g_stseaspos_Ae(2)).oldcol = {x_.col_ssp_Ase.col};
parGS(g.g_stseaspos_An(1)).oldcol = {x_.col_ssp_Acn.col};
parGS(g.g_stseaspos_An(2)).oldcol = {x_.col_ssp_Asn.col};

parGS(g.g_hpole).oldcol = x_.col_hpole;
parGS(g.g_lpole).oldcol = x_.col_lpole;

parGS(g.g_gamma).oldcol = x_.col_gamma;

parGS(g.g_aplrg).oldcol = x_.col_rg;

% tidal ERP variations
if parGS(g.g_tidpm).id==1 || parGS(g.g_tidut).id==1
    parGS(g.g_tidpm).oldcol = [x_.col_tidap(parGS(g.g_tidpm).spectid) ...
                               x_.col_tidbp(parGS(g.g_tidpm).spectid) ...
                               x_.col_tidam(parGS(g.g_tidpm).spectidret)...
                               x_.col_tidbm(parGS(g.g_tidpm).spectidret)];
                           
    parGS(g.g_tidut).oldcol = [x_.col_tidutc(parGS(g.g_tidut).spectid) ...
                               x_.col_tiduts(parGS(g.g_tidut).spectid)];
end

if isfield(x_,'bdclko')
    parGS(g.g_bdco).oldcol = [x_.bdclko.col];
else
    parGS(g.g_bdco).oldcol = [];
end
