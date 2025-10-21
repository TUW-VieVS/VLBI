% ************************************************************************
%   Description:
%   function to set undefined fields in case a old parameter file is loaded
%
%
%   Input: 
%           parameter struct
%  
%
%   Output:
%           parameter struct
%
%
%   Coded for VieVS:
%   3 June 2025 by H. Wolf 
%
%   Revision:
%   2025-Sep-25, P. Urban: ambiguity parameters added
%
% ************************************************************************

function [parameter] = checkParameterFile(parameter)

    if ~isfield(parameter.lsmopt, 'stc_sat')
        parameter.lsmopt.stc_sat=0;
    end
    if ~isfield(parameter.lsmopt, 'stc_qu')
        parameter.lsmopt.stc_qu=0;
    end
    if ~isfield(parameter.lsmopt, 'stc_all')
        parameter.lsmopt.stc_all=1;
    end
    if ~isfield(parameter.lsmopt, 'stc_qs')
        parameter.lsmopt.stc_qs=0;
    end
    if ~isfield(parameter.lsmopt, 'nnt_stc')
        parameter.lsmopt.nnt_stc=0;
    end
    if ~isfield(parameter.lsmopt, 'nnr_stc')
        parameter.lsmopt.nnr_stc=0;
    end
    if ~isfield(parameter.lsmopt, 'nns_stc')
        parameter.lsmopt.nns_stc=0;
    end
    if ~isfield(parameter.lsmopt, 'stc_qs_snx_sat')
        parameter.lsmopt.stc_qs_snx_sat=0;
    end
    if ~isfield(parameter.lsmopt, 'stc_qs_snx_qu')
        parameter.lsmopt.stc_qs_snx_qu=0;
    end
    if (parameter.lsmopt.nnt_stc==1 || parameter.lsmopt.nnr_stc==1 || parameter.lsmopt.nns_stc==1) && ~isfield(parameter.lsmopt, 'addDatumCd')
        parameter.lsmopt.addDatumCd=1;
    elseif ~isfield(parameter.lsmopt, 'addDatumCd')
        parameter.lsmopt.addDatumCd=0;
    end
    if ~isfield(parameter.lsmopt.outsnx, 'orb')
        parameter.lsmopt.outsnx.orb=0;
    end
    if ~isfield(parameter.lsmopt,'SatPos')
        parameter.lsmopt.SatPos.pw_sat=0;
        parameter.lsmopt.SatPos.sat_pos_int=0;
        parameter.lsmopt.SatPos.constr_sat=0;
        parameter.lsmopt.SatPos.sat_pos_coef=0;
        parameter.lsmopt.SatPos.sat_pos_coef=0;
        parameter.lsmopt.SatPos.sat_pos_est_ref_frame='';

    end
    if ~isfield(parameter.lsmopt, 'KepEle')
        parameter.lsmopt.KepEle.estKepEle=0;
        parameter.lsmopt.KepEle.estKepEle_FRP=0;
    end
    if ~isfield(parameter.lsmopt, 'trf_excldatum')
        parameter.lsmopt.trf_excldatum = 0;
        parameter.lsmopt.trf_excldatum_file = ' ';
    end
    if ~isfield(parameter.lsmopt, 'remove_sources_from_list')
        parameter.lsmopt.remove_sources_from_list = 0;
    end
    if ~isfield(parameter.lsmopt, 'datum')
        parameter.lsmopt.datum = 'trf';
    end
    if ~isfield(parameter.lsmopt, 'pw_sou_select')
        parameter.lsmopt.pw_sou_select = 'notcat';
    end
    if ~isfield(parameter.lsmopt, 'est_sourceNNR_selection')
        parameter.lsmopt.est_sourceNNR_selection = 0;
    end
    if ~isfield(parameter.lsmopt, 'est_sourceNNR_defining')
        parameter.lsmopt.est_sourceNNR_defining = 1;
    end
    if parameter.lsmopt.est_sourceNNR_defining == 0 && ( ~isfield(parameter.lsmopt, 'est_sourceNNR_selection') || parameter.lsmopt.est_sourceNNR_selection == 0)
        parameter.lsmopt.est_sourceNNR_defining = 1;
    end
    if ~isfield(parameter.lsmopt, 'bdco_auto')
        parameter.lsmopt.bdco_auto = 1;
    end

    % to ensure compatibility with parameter files before the vievs update (10/2020)
    % where the estimation of scale offset was added. Can be removed in the future!
    if ~isfield(parameter.lsmopt, 'est_scale') 
        parameter.lsmopt.est_scale=0;
    end
    
    % to ensure compatibility with parameter files before the vievs update (10/2020)
    % where the bas-dep clock offsets were added. Can be removed in the future!
    if ~isfield(parameter.lsmopt, 'est_bdco') 
       parameter.lsmopt.est_bdco=0; 
    end

    if ~isfield(parameter.vie_mod, 'cntrl')
        parameter.vie_mod.cntrl = 0;
    end
    if ~isfield(parameter.vie_mod, 'cntrlm')
        parameter.vie_mod.cntrlm = ' ';
    end
    if ~isfield(parameter, 'createSkyPlots')
        parameter.createSkyPlots = 0;
    end
    
    % ambiguity related parameters 
    if ~isfield(parameter.vie_init, 'res_compute')
        parameter.vie_init.res_compute = 0;
    end
    if ~isfield(parameter.lsmopt, 'res_compute')
        parameter.lsmopt.res_compute = 0;
    end
    if ~isfield(parameter.vie_init, 'res_apply')
        parameter.vie_init.res_apply = 0;
    end
    if ~isfield(parameter.lsmopt, 'res_apply')
        parameter.lsmopt.res_apply = 0;
    end
    if ~isfield(parameter.vie_init, 'amb')
        parameter.vie_init.amb='observation_database_amb';
    end
 
end
