% #########################################################################
% #     saveParamFile
% #########################################################################
%
% DESCRITPION
% This function saves the structure "parameter", which contains all
% required parameters to run the VieVS batch mode, to a *.mat file.
%
% AUTHOR 
%   Matthias Madzak
%
% INPUT
%   handles      structure from the GUI (also containing data, eg residuals)
%   hObject      ...unused so far.
%   fullOutName  filename of saved parameter file (including filepath!)
%
% OUTPUT
%
% CHANGES
%   2014-11-05, A. Hellerschmied: Exctact code from "vievs2_3.m" file and
%     create this file
%   2014-11-05, A. Hellerschmied: Parameter "use_opt_files" added.
%   2014-12-09, A. Hellerschmied: Add field "use_opt_files" to struct "parameter.vie_init" 
%   2015-01-13, M. Madzak: contact details are saved for SINEX output
%   2015-02-03, A. Hellerschmied: Modification related to GUI update and
%       bug fix
%   2015-07-14, A. Hellerschmied: New mean pole model added (IERS 2015 model, "cmp2015") 
%   2015-08-21, S. Boehm: Parameter parameter.lsmopt.est_erpsp added (tidal ERP variations).
%   2015-10-21, H. Krasna: seasonal variations in station positions (annual, semi-annual), pole
%               tide Love and Shida numbers, APL regression coefficient,
%               GIA uplift added 
%   2016-06-11, H. Krasna: FCN period from nutation and amplitudes added
%   2016-07-07, A. Hellerschmied: Added checkbox "checkbox_estimation_leastSquares_sources_ICRF2_def"
%   2016-08-03, A. Hellerschmied: Added selection of SP3 file via "edit_models_sc_sp3_file" (parameter.vie_init.sc_orbit_file_path_name, parameter.vie_init.sc_orbit_file_type)
%   2016-08-31, A. Girdiuk: New ephemerides is added
%   2016-12-15, H. Krasna: Elevation dependent noise (P matrix) added
%   2017-01-16, H. Krasna: "checkbox_run_globalPram_stseaspos", and "checkbox_run_globalPram_hlpole" added
%   2017/01/19, D. Landskron: choice of temperature source transferred from "Models - Troposphere" to "Models - Station models", and bug fixed
%   2017/01/23, D. Landskron: general shape of Troposphere changed
%   2017/02/19, H. Krasna: checkbox_run_globalPram_APLrg,
%               checkbox_parameters_statCorr_APLrg, checkbox_parameters_statCorr_GIA added
%   2017-02-23, A. Hellerschmied: "handles.edit_run_add_noise" added
%   2017-03-01, A. Hellerschmied: "handles.checkbox_run_runOptions_manuallyFindBreaks" removed
%   2017-03-02, H. Krasna: possibility to apply for antenna axis offset altitude correction in vie_mod
%   2017-03-20, A. Hellerschmied: Changes required for new EOP data loading routine
%   2017-09-13, D. Landskron: 'tropSource' shifted into 'vie_init' 
%   2018-01-11, D. Landskron: external troposphere modeling removed
%   2018-07-06, D. Landskron: VMF3 added to the troposphere models 
%   2018-11-29, D. Landskron: structure of observation restrictions standardized
%   2018-03-06, D. Landskron: suffix checkbox added to the sinex files

function saveParamFile(hObject, handles, fullOutName)

% preallocate
auto_save_parameterfileparameter=struct('vie_init', [], 'vie_mod', [], 'lsmopt', []);

% ===========================
% OPT files
% ===========================

% OPT file flag (apply opt file options?):
parameter.opt.use_opt_files=get(handles.checkbox_setInput_useOptFiles, 'Value');

% OPT file directory:
allOptDirs=get(handles.popupmenu_setInput_optDir, 'String');
if ~isempty(allOptDirs)
    parameter.opt.opt_file_dir = allOptDirs{get(handles.popupmenu_setInput_optDir, 'Value')};
else
    parameter.opt.opt_file_dir = '';
end


% ===========================
% Outlier files
% ===========================

% outlier file directory:
allOutDirs = get(handles.popupmenu_setInput_outDir, 'String');
if strcmp(allOutDirs, ' ')
    parameter.outlier.out_file_dir = '';
else
    parameter.outlier.out_file_dir = allOutDirs{get(handles.popupmenu_setInput_outDir, 'Value')};
end

% remove outlier option:
parameter.outlier.flag_remove_outlier = get(handles.checkbox_setInput_eliminOutliers, 'Value');

% ===========================
% Observation restrictions
% ===========================

% Set cut-off elevation:
minElevInputNum = str2double(get(handles.edit_parameter_obsRestr_cutOff, 'String'));
if isnan(minElevInputNum)
    parameter.obs_restrictions.cut_off_elev = 0; % [rad]
else
    parameter.obs_restrictions.cut_off_elev = minElevInputNum * pi/180; % [rad]
end

% Quality code limit:
qLimitInputNum = str2double(get(handles.edit_parameter_obsRestr_qualityCode, 'String'));
if isnan(qLimitInputNum)
    parameter.obs_restrictions.Qlim=0;
else
    parameter.obs_restrictions.Qlim=qLimitInputNum;
end


% ===========================
% VIE_INIT
% ===========================

% opt directory
allOptDirs=get(handles.popupmenu_setInput_optDir, 'String');
if ~isempty(allOptDirs)
    parameter.vie_init.diropt=allOptDirs{get(handles.popupmenu_setInput_optDir, 'Value')};
end

% use OPT file option
parameter.vie_init.use_opt_files=get(handles.checkbox_setInput_useOptFiles, 'Value');

% outlier directory
allOutDirs=get(handles.popupmenu_setInput_outDir, 'String');
if strcmp(allOutDirs, ' ')
    parameter.vie_init.dirout='';
else
    parameter.vie_init.dirout=allOutDirs{get(handles.popupmenu_setInput_outDir, 'Value')};
end

% remove outlier option
parameter.vie_init.rm_outlier=get(handles.checkbox_setInput_eliminOutliers, 'Value');

% ambiguity correction (on or off)
parameter.vie_init.ambiguity_correction=get(handles.checkbox_ambiguity_correction,'Value');

% ionospheric correction (on or off)
parameter.vie_init.iono_correction=get(handles.checkbox_ionosphere_corrected,'Value');

% observation parameter
parameter.vie_init.vgosDb_observation_parameter=get(handles.edit_vgosdb_observation ,'String');

% institutes
parameter.vie_init.vgosDb_institute=get(handles.edit_vgosdb_wrapper_institute ,'String');

% wrapper version
parameter.vie_init.vgosDb_wrapper_version=get(handles.edit_wrapper_version_number ,'String');

% trf
% if superstations file is ticked
if get(handles.radiobutton_parameters_refFrames_superstationTRF, 'Value')
    allPopupmenuEntries=get(handles.popupmenu_parameters_refFrames_superstationTRF, 'String');
    % save filename
    parameter.vie_init.trf{1}=handles.data.superstationFile;
    % and field
    parameter.vie_init.trf{2}=allPopupmenuEntries{get(handles.popupmenu_parameters_refFrames_superstationTRF, 'Value')};  % only the struct-fieldname is saved to parameter struct
else % manual txt trf file was chosen
    allPopupmenuEntries=get(handles.popupmenu_parameters_refFrames_otherTRF, 'String');
    parameter.vie_init.trf{1}=['../TRF/', allPopupmenuEntries{get(handles.popupmenu_parameters_refFrames_otherTRF, 'Value')}];
    parameter.vie_init.trf{2}=' ';
end

% crf
% if supersource file is ticked
if get(handles.radiobutton_parameters_refFrames_supersourceCRF, 'Value')
    allPopupmenuEntries=get(handles.popupmenu_parameters_refFrames_supersourceCRF, 'String');
    % save filename
    parameter.vie_init.crf{1}=handles.data.supersourceFile;
    % and field
    parameter.vie_init.crf{2}=allPopupmenuEntries{get(handles.popupmenu_parameters_refFrames_supersourceCRF, 'Value')};
else
    allPopupmenuEntries=get(handles.popupmenu_parameters_refFrames_otherCRF, 'String');
    parameter.vie_init.crf{1}=['../CRF/', allPopupmenuEntries{get(handles.popupmenu_parameters_refFrames_otherCRF, 'Value')}];
    parameter.vie_init.crf{2}=' ';
end

% station info file (not variable yet)
% parameter.vie_init.sta_info='../TRF/STAT_INFO';

jetangleInput=get(handles.edit_parameter_obsRestr_jetang,'String');
if strcmp(jetangleInput,'none')
    parameter.vie_init.ex_jet=0;
else
    parameter.vie_init.ex_jet=1;
    parameter.vie_init.jetang=str2double(jetangleInput);
end

% tropospheric correction
if get(handles.radiobutton_parameters_troposphere_indModeling, 'Value')
    parameter.vie_init.tropSource.name='indModeling';
elseif get(handles.radiobutton_parameters_troposphere_raytr, 'Value')
    parameter.vie_init.tropSource.name='raytr';
end

% ionospheric correction
if get(handles.radiobutton_parameters_iono_fromNGS, 'Value')
    parameter.vie_init.iono='observation_database';
else % external file was chosen
    parameter.vie_init.iono='ext';
    
    % get all iono folders
    allIonoFolders=get(handles.popupmenu_parameters_iono_ext, 'String');
    if iscell(allIonoFolders)
        parameter.vie_init.ionoFolder=allIonoFolders{get(handles.popupmenu_parameters_iono_ext, 'Value')};
    else
        parameter.vie_init.ionoFolder=allIonoFolders(get(handles.popupmenu_parameters_iono_ext, 'Value'));
    end
end

% load orbit data file
sp3_file_name_path_str = get(handles.edit_models_sc_sp3_file,'String');
if ~isempty(sp3_file_name_path_str)
    if ~strcmp(sp3_file_name_path_str, ' ')
        parameter.vie_init.sc_orbit_file_path_name  = sp3_file_name_path_str;
        parameter.vie_init.sc_orbit_file_type       = 'sp3';
    else
        parameter.vie_init.sc_orbit_file_path_name  = '';
        parameter.vie_init.sc_orbit_file_type       = '';
    end
else
    parameter.vie_init.sc_orbit_file_path_name  = '';
    parameter.vie_init.sc_orbit_file_type       = '';
end


% VIE_MOD
% =========================

% ephemerides
if     get(handles.radiobutton_parameters_eph_jpl405, 'Value')
    parameter.vie_mod.eph='jpl_405';
elseif get(handles.radiobutton_parameters_eph_jpl421, 'Value')
    parameter.vie_mod.eph='jpl_421';
elseif get(handles.radiobutton_parameters_eph_jpl430, 'Value')
    parameter.vie_mod.eph='jpl_430';
end

% EOP file
if get(handles.radiobutton_parameters_eop_aPriori_05C04, 'Value')
    parameter.vie_mod.EOPfile = 'C04_08_1962_now.txt';
elseif get(handles.radiobutton_parameters_eop_aPriori_08C04, 'Value')
    parameter.vie_mod.EOPfile = 'C04_14_1962_now.txt';
elseif get(handles.radiobutton_parameters_eop_aPriori_finals, 'Value')
    % parameter.vie_mod.EOPfile='eop_finals2000A.txt';
    parameter.vie_mod.EOPfile = 'finals_all_IAU2000.txt';
else % other trf was chosen
    allEopFiles=get(handles.popupmenu_parameters_eop_aPriori_other, 'String');
    if ischar(allEopFiles)
        parameter.vie_mod.EOPfile=allEopFiles;
    else % we have a cell of different entries
        parameter.vie_mod.EOPfile=allEopFiles{get(handles.popupmenu_parameters_eop_aPriori_other, 'Value')};
    end
end

% interpolation - linear?
parameter.vie_mod.linear=get(handles.radiobutton_parameters_eop_interp_lin, 'Value');
if get(handles.radiobutton_parameters_eop_interp_lin, 'Value')==1
 parameter.vie_mod.linear48h=get(handles.checkbox_parameters_eop_interp_lin48h, 'Value');
else
 parameter.vie_mod.linear48h=0;
end

% tidal UT
parameter.vie_mod.tidalUT=get(handles.checkbox_parameters_eop_tidalUtVariations, 'Value');
parameter.vie_mod.tidalUT35=get(handles.rb_eop_inter_UT1R, 'Value');

% ocean tides
    % parameter.vie_mod.eopoc='interpf (Conventions)';
allHfFiles=get(handles.popupmenu_parameters_eop_oceanTideModel, 'String');
if ischar(allHfFiles)
    parameter.vie_mod.eopoc=allHfFiles;
else % we have a cell with more entries than one
    parameter.vie_mod.eopoc=allHfFiles{get(handles.popupmenu_parameters_eop_oceanTideModel, 'Value')};
end
% if ocean tides were unticked -> use model "none"
if get(handles.checkbox_parameters_eop_inclHf_oceanTidesOther, 'value')==0
    parameter.vie_mod.eopoc='none';
end


% libration
parameter.vie_mod.lib_pm=get(handles.checkbox_parameters_eop_inclHf_LibrationXpYp, 'value');
parameter.vie_mod.lib_ut=get(handles.checkbox_parameters_eop_inclHf_LibrationUt1, 'value');

% nutation offsets at all
parameter.vie_mod.dXdY=get(handles.checkbox_parameters_eop_models_inclAPrioriNutOffs, 'Value');

% precession/nutation model
if get(handles.radiobutton_parameters_eop_precNutModel_2000a, 'value')
    parameter.vie_mod.nutmod='IAU_2000A';
else % IAU 2006/2000A
    parameter.vie_mod.nutmod='IAU_2006/2000A';
end

% station corrections
parameter.vie_mod.cts=get(handles.checkbox_parameters_statCorr_solidEarthTides, 'Value');
parameter.vie_mod.cto=get(handles.checkbox_parameters_statCorr_tidalOceanLoad, 'Value');
parameter.vie_mod.cta=get(handles.checkbox_parameters_statCorr_tidalAtmoLoad, 'Value');
parameter.vie_mod.cnta=get(handles.checkbox_parameters_statCorr_nonTidalAtmoLoad, 'Value');
parameter.vie_mod.ctp=get(handles.checkbox_parameters_statCorr_poleTides, 'Value');
parameter.vie_mod.ctop=get(handles.checkbox_parameters_statCorr_oceanPoleTides, 'Value');
parameter.vie_mod.chl=get(handles.checkbox_parameters_statCorr_hydroLoading, 'Value');
parameter.vie_mod.therm=get(handles.checkbox_parameters_statCorr_thermalDef, 'Value');
parameter.vie_mod.gravDef=get(handles.checkbox_parameters_statCorr_gravitationalDef, 'Value');
parameter.vie_mod.crg = get(handles.checkbox_parameters_statCorr_APLrg, 'Value');
if parameter.vie_mod.crg == 1
    parameter.vie_mod.cnta = 0; parameter.vie_mod.cta = 0;
end
parameter.vie_mod.gia = get(handles.checkbox_parameters_statCorr_GIA, 'Value');


% Mean pole model:
if get(handles.radiobutton_parameters_statCorr_poleModel_lin, 'Value') % linear model (IERS 2003)
    parameter.vie_mod.ctpm='linear';
elseif get(handles.radiobutton_parameters_statCorr_poleModel_cub, 'Value') % cubic model (IERS 2010)
    parameter.vie_mod.ctpm='cubic';
elseif get(handles.radiobutton_parameters_statCorr_poleModel_iers2015, 'Value') % IERS 2015 model
    parameter.vie_mod.ctpm='cmp2015';
end

% temperature source for antenna deformation
if get(handles.checkbox_parameters_statCorr_temp_fromInSitu, 'Value')
    parameter.vie_init.tp='in situ';
elseif get(handles.checkbox_parameters_statCorr_temp_GPT3, 'Value')
    parameter.vie_init.tp='gpt3';
else
    error('There is something wrong here!')
end

% loading models
allTidOceanLoadFiles=get(handles.popupmenu_parameters_statCorr_tidalOceanLoad, 'String');
allTidAtmoLoadFiles=get(handles.popupmenu_parameters_statCorr_tidalAtmoOceanLoad, 'String');
allNTidAtmoLoadFiles=get(handles.popupmenu_parameters_statCorr_nonTidalAtmoOceanLoad, 'String');
allHydroLoadFiles=get(handles.popupmenu_parameters_statCorr_hydroLoading, 'String');
allAPLrgLoadFiles=get(handles.popupmenu_parameters_statCorr_APLrg, 'String');
allGIALoadFiles=get(handles.popupmenu_parameters_statCorr_GIA, 'String');



if strcmp(allTidOceanLoadFiles, ' ')
    parameter.vie_mod.ocm=' ';
else
    parameter.vie_mod.ocm=allTidOceanLoadFiles{get(handles.popupmenu_parameters_statCorr_tidalOceanLoad, 'Value')};
end
if strcmp(allTidAtmoLoadFiles, ' ')
    parameter.vie_mod.ctam=' ';
else
    parameter.vie_mod.ctam=allTidAtmoLoadFiles{get(handles.popupmenu_parameters_statCorr_tidalAtmoOceanLoad, 'Value')};
end
if strcmp(allNTidAtmoLoadFiles, ' ')
    parameter.vie_mod.cntam=' ';
else
    parameter.vie_mod.cntam=allNTidAtmoLoadFiles{get(handles.popupmenu_parameters_statCorr_nonTidalAtmoOceanLoad, 'Value')};
end
if strcmp(allHydroLoadFiles, ' ')
    parameter.vie_mod.chlm=' ';
else
    parameter.vie_mod.chlm=allHydroLoadFiles{get(handles.popupmenu_parameters_statCorr_hydroLoading, 'Value')};
end
if strcmp(allAPLrgLoadFiles, ' ')
    parameter.vie_mod.crgm=' ';
else
    parameter.vie_mod.crgm=allAPLrgLoadFiles{get(handles.popupmenu_parameters_statCorr_APLrg, 'Value')};
end
if strcmp(allGIALoadFiles, ' ')
    parameter.vie_mod.giam=' ';
else
    parameter.vie_mod.giam=allGIALoadFiles{get(handles.popupmenu_parameters_statCorr_GIA, 'Value')};
end



% zenith delay
if get(handles.radiobutton_parameters_troposphere_zhd_fromInSitu, 'Value')
    parameter.vie_init.zhd='in situ';
elseif get(handles.radiobutton_parameters_troposphere_zhd_no, 'Value')
    parameter.vie_init.zhd='no';
elseif get(handles.radiobutton_parameters_troposphere_zhd_VMF3, 'Value')
    parameter.vie_init.zhd='vmf3';
elseif get(handles.radiobutton_parameters_troposphere_zhd_VMF1, 'Value')
    parameter.vie_init.zhd='vmf1';
elseif get(handles.radiobutton_parameters_troposphere_zhd_GPT3, 'Value')
    parameter.vie_init.zhd='gpt3';
end
if get(handles.radiobutton_parameters_troposphere_zwd_no, 'Value')
    parameter.vie_init.zwd='no';
elseif get(handles.radiobutton_parameters_troposphere_zwd_fromInSitu, 'Value')
    parameter.vie_init.zwd='in situ';
elseif get(handles.radiobutton_parameters_troposphere_zwd_VMF3, 'Value')
    parameter.vie_init.zwd='vmf3';
elseif get(handles.radiobutton_parameters_troposphere_zwd_VMF1, 'Value')
    parameter.vie_init.zwd='vmf1';
elseif get(handles.radiobutton_parameters_troposphere_zwd_GPT3, 'Value')
    parameter.vie_init.zwd='gpt3';
end

% gradients
if get(handles.radiobutton_parameters_troposphere_gradients_h_no, 'Value')
    parameter.vie_mod.apgm_h='no';
elseif get(handles.radiobutton_parameters_troposphere_gradients_h_GRAD, 'Value')
    parameter.vie_mod.apgm_h='grad';
elseif get(handles.radiobutton_parameters_troposphere_gradients_h_GPT3, 'Value')
    parameter.vie_mod.apgm_h='gpt3';
elseif get(handles.radiobutton_parameters_troposphere_gradients_h_DAO, 'Value')
    parameter.vie_mod.apgm_h='dao';
end
if get(handles.radiobutton_parameters_troposphere_gradients_w_no, 'Value')
    parameter.vie_mod.apgm_w='no';
elseif get(handles.radiobutton_parameters_troposphere_gradients_w_GRAD, 'Value')
    parameter.vie_mod.apgm_w='grad';
elseif get(handles.radiobutton_parameters_troposphere_gradients_w_GPT3, 'Value')
    parameter.vie_mod.apgm_w='gpt3';
end

% mapping function
if get(handles.radiobutton_parameters_troposphere_mfh_VMF3, 'value')
    parameter.vie_mod.mfh='vmf3';
elseif get(handles.radiobutton_parameters_troposphere_mfh_VMF1, 'value')
    parameter.vie_mod.mfh='vmf1';
elseif get(handles.radiobutton_parameters_troposphere_mfh_GPT3, 'value')
    parameter.vie_mod.mfh='gpt3';
end
if get(handles.radiobutton_parameters_troposphere_mfw_VMF3, 'value')
    parameter.vie_mod.mfw='vmf3';
elseif get(handles.radiobutton_parameters_troposphere_mfw_VMF1, 'value')
    parameter.vie_mod.mfw='vmf1';
elseif get(handles.radiobutton_parameters_troposphere_mfw_GPT3, 'value')
    parameter.vie_mod.mfw='gpt3';
end

% Source Structure
allSSCatalogFiles=get(handles.popupmenu_parameters_ss_catalog, 'String');
if strcmp(allSSCatalogFiles, ' ')
    parameter.vie_mod.sou_cat=' ';
else
    parameter.vie_mod.sou_cat=allSSCatalogFiles{get(handles.popupmenu_parameters_ss_catalog, 'Value')};
end
parameter.vie_mod.ssou=get(handles.checkbox_parameters_ss_applyss,'Value');
parameter.vie_mod.write_jet=get(handles.checkbox_parameters_ss_writejet, 'Value');


% VIE_LSM
% =========================



% use OPT file option
parameter.lsmopt.use_opt_files=get(handles.checkbox_setInput_useOptFiles, 'Value');

% LEVEL 1 output sub-directory
if get(handles.checkbox_run_outDirs_diffSubs, 'Value')
    % different subs were chosen
    parameter.lsmopt.level1OutDir=get(handles.edit_run_outDirs_level1, 'String');
else
    % One sub-directory
    parameter.lsmopt.level1OutDir=get(handles.edit_run_outDirs_oneSub, 'String');
end

% opt and outlier directory (copied)
parameter.lsmopt.diropt=parameter.vie_init.diropt;
parameter.lsmopt.dirout=parameter.vie_init.dirout;

% first / main solution
parameter.lsmopt.first=get(handles.checkbox_run_runFirstSolution, 'Value');
parameter.lsmopt.second=get(handles.checkbox_run_runMainSolution, 'Value');

% baseline dependent weighting
parameter.lsmopt.bsldep=get(handles.checkbox_run_bslDepWeights, 'Value');

% elevation dependent noise
parameter.lsmopt.eldep_noise=get(handles.checkbox_run_ElevDepNoise, 'Value');

% add constant noise [cm]
parameter.lsmopt.add_const_noise_cm = str2double(get(handles.edit_run_add_noise, 'String'));


% clock treatment of first solution
if get(handles.radiobutton_run_oneOffsPerClock, 'Value')
    parameter.lsmopt.firstclock=0;
elseif get(handles.radiobutton_run_oneOffsAndRatePerClock, 'Value')
    parameter.lsmopt.firstclock=1;
else % offset+rate+quadratic
    parameter.lsmopt.firstclock=2;
end

% treat clock breaks
parameter.lsmopt.treat_breaks=get(handles.checkbox_estimation_leastSquares_clocks_useClockBreaks, 'value');

% estimation of clocks
if get(handles.checkbox_estimation_leastSquares_clocks, 'Value')==0
    parameter.lsmopt.pw_clk=0;
else
    if get(handles.radiobutton_estimation_leastSquares_pwlClock, 'Value')
        parameter.lsmopt.pw_clk=1;
    elseif get(handles.radiobutton_estimation_leastSquares_pwlAndRateClock, 'value')
        parameter.lsmopt.pw_clk=2;
    else % offset+rate+quadratic
        parameter.lsmopt.pw_clk=3;
    end
end

% clock constraints
parameter.lsmopt.constr_clk=get(handles.checkbox_estimation_leastSquares_clocksRelConstr, 'Value');
clockConstrCoeffNum=str2double(get(handles.edit_estimation_leastSquares_clockConstr, 'String'));
if isnan(clockConstrCoeffNum)
    parameter.lsmopt.coef_clk=0.5;
else
    parameter.lsmopt.coef_clk=clockConstrCoeffNum;
end

% clock interval
clockIntNum=str2double(get(handles.edit_estimation_leastSquares_clockInterval, 'String'));
if isnan(clockIntNum)
    parameter.lsmopt.int_clk=60;
else
    parameter.lsmopt.int_clk=clockIntNum;
end

parameter.lsmopt.ref=0; % set 1st clock as reference (could be chagned in OPT file)
parameter.lsmopt.ref_first_clk=1; %...

% baseline dependent clock offset
parameter.lsmopt.est_bdco = get(handles.checkbox_estimation_leastSquares_clocksBasDepOffset, 'Value');
bdco_mino=str2double(get(handles.edit_estimation_leastSquares_clockBasDepO_minNobs, 'String'));
if isnan(bdco_mino)
    parameter.lsmopt.bdco_minobs = 5; % min obs at the baseline
else
    parameter.lsmopt.bdco_minobs = bdco_mino;
end
parameter.lsmopt.bdco_fromOPT = get(handles.radiobutton_estimation_leastSquares_basdepClockoff_OPT, 'Value'); % 1 - read from OPT file, 0 - create the list in vie_lsm
parameter.lsmopt.bdco_nrlist = [];    

% estimation of zwd
parameter.lsmopt.pw_zwd=get(handles.checkbox_estimation_leastSquares_tropo_zwd, 'Value');

% zwd contraint
parameter.lsmopt.constr_zwd=get(handles.checkbox_estimation_leastSquares_tropo_zwdRelConstr, 'Value');
zwdConstrNum=str2double(get(handles.edit_estimation_leastSquares_tropo_zwdRelConstr, 'String'));
if isnan(zwdConstrNum)
    parameter.lsmopt.coef_zwd=0.7;
else
    parameter.lsmopt.coef_zwd=zwdConstrNum;
end

% zwd interval
zwdIntervalNum=str2double(get(handles.edit_estimation_leastSquares_tropo_zwdInterval, 'String'));
if isnan(zwdIntervalNum)
    parameter.lsmopt.int_zwd=30;
else
    parameter.lsmopt.int_zwd=zwdIntervalNum;
end

% estimation of north gradients
parameter.lsmopt.pw_ngr=get(handles.checkbox_estimation_leastSquares_tropo_ngr, 'Value');

% constraints of north gradients
parameter.lsmopt.constr_rel_ngr=get(handles.checkbox_estimation_leastSquares_tropo_ngrRelConstr, 'Value');
parameter.lsmopt.constr_abs_ngr=get(handles.checkbox_estimation_leastSquares_tropo_ngrAbsConstr, 'Value');
ngrRelConstrNum=str2double(get(handles.edit_estimation_leastSquares_tropo_ngrRelConstr, 'String'));
if isnan(ngrRelConstrNum)
    parameter.lsmopt.coef_rel_ngr=0.7;
else
    parameter.lsmopt.coef_rel_ngr=ngrRelConstrNum;
end
ngrAbsConstrNum=str2double(get(handles.edit_estimation_leastSquares_tropo_ngrAbsConstr, 'String'));
if isnan(ngrAbsConstrNum)
    parameter.lsmopt.coef_abs_ngr=0.7;
else
    parameter.lsmopt.coef_abs_ngr=ngrAbsConstrNum;
end

% interval of north gradients
ngrIntervalNum=str2double(get(handles.edit_estimation_leastSquares_tropo_ngrInterval, 'String'));
if isnan(ngrIntervalNum)
    parameter.lsmopt.int_ngr=360;
else
    parameter.lsmopt.int_ngr=ngrIntervalNum;
end

% estimation of east gradients
parameter.lsmopt.pw_egr=get(handles.checkbox_estimation_leastSquares_tropo_egr, 'Value');

% constraints of east gradients
parameter.lsmopt.constr_rel_egr=get(handles.checkbox_estimation_leastSquares_tropo_egrRelConstr, 'Value');
parameter.lsmopt.constr_abs_egr=get(handles.checkbox_estimation_leastSquares_tropo_egrAbsConstr, 'Value');
egrRelConstrNum=str2double(get(handles.edit_estimation_leastSquares_tropo_egrRelConstr, 'String'));
if isnan(egrRelConstrNum)
    parameter.lsmopt.coef_rel_egr=0.7;
else
    parameter.lsmopt.coef_rel_egr=egrRelConstrNum;
end
egrAbsConstrNum=str2double(get(handles.edit_estimation_leastSquares_tropo_egrAbsConstr, 'String'));
if isnan(egrAbsConstrNum)
    parameter.lsmopt.coef_abs_egr=0.7;
else
    parameter.lsmopt.coef_abs_egr=egrAbsConstrNum;
end

% interval of east gradients
egrIntervalNum=str2double(get(handles.edit_estimation_leastSquares_tropo_egrInterval, 'String'));
if isnan(egrIntervalNum)
    parameter.lsmopt.int_egr=360;
else
    parameter.lsmopt.int_egr=egrIntervalNum;
end

% station coordinates estimation
parameter.lsmopt.stc=get(handles.checkbox_estimation_leastSquares_coordinates_estimate, 'Value');
parameter.lsmopt.pw_stc=0; % is changed in sessionwise parameterization
parameter.lsmopt.nnt_stc=get(handles.checkbox_estimation_leastSquares_coordinates_NNT, 'Value');
parameter.lsmopt.nnr_stc=get(handles.checkbox_estimation_leastSquares_coordinates_NNR, 'Value');
parameter.lsmopt.sca_stc=get(handles.checkbox_estimation_leastSquares_coordinates_NNS, 'Value');
parameter.lsmopt.datum='trf';
if get(handles.radiobutton_estimation_leastSquares_coordinates_datum_all, 'Value')
    parameter.lsmopt.datum='all';
else
    parameter.lsmopt.datum='trf';
end
parameter.lsmopt.constr_xyz=1;
parameter.lsmopt.coef_xyz=10;
parameter.lsmopt.int_xyz=360;


% estimation of sources
parameter.lsmopt.pw_sou                 = get(handles.checkbox_estimation_leastSquares_sources_est, 'Value');
parameter.lsmopt.est_sourceNNR          = get(handles.checkbox_estimation_leastSquares_sources_NNR, 'Value');
parameter.lsmopt.est_sourceNNR_defining = get(handles.checkbox_estimation_leastSquares_sources_ICRF2_def, 'Value');

% interval of sources
sourceIntervalNum=str2double(get(handles.edit_estimation_leastSquares_sources_interval, 'String'));
if isnan(sourceIntervalNum)
    parameter.lsmopt.sour_int_rade=1440;
else
    parameter.lsmopt.sour_int_rade=sourceIntervalNum;
end

% constraints of sources
parameter.lsmopt.constr_sou=get(handles.checkbox_estimation_leastSquares_sources_constr, 'Value');
sourceConstrCoefNum=str2double(get(handles.edit_estimation_leastSquares_sources_constr, 'String'));
if isnan(sourceConstrCoefNum)
    parameter.lsmopt.sour_coef_rade=1.0e-4;
else
    parameter.lsmopt.sour_coef_rade=sourceConstrCoefNum;
end

%special parameters for source estimation

parameter.lsmopt.UseSourceAbsConstrNNR = get(handles.checkbox_estimation_leastSquares_sources_abs_constr, 'Value');
if parameter.lsmopt.UseSourceAbsConstrNNR
    sourceAbsConstrNNR=str2double(get(handles.edit_estimation_leastSquares_sources_abs_constr, 'String'));
    if isnan(sourceAbsConstrNNR)
        warning('Inappropriat input at Estimation/Leaste squares/Source coordinates ')
    else
        parameter.lsmopt.sourceAbsConstrNNR=sourceAbsConstrNNR;
    end
else
    parameter.lsmopt.sourceAbsConstrNNR = 1;
end

parameter.lsmopt.use_min_num_obs_per_est_source = get(handles.checkbox_estimation_leastSquares_sources_obs_per_source, 'Value');
if parameter.lsmopt.use_min_num_obs_per_est_source
    parameter.lsmopt.min_num_obs_per_est_source=str2double(get(handles.edit_estimation_leastSquares_sources_obs_per_source, 'String'));
    if isnan(parameter.lsmopt.min_num_obs_per_est_source)
        warning('Inappropriat input at Estimation/Leaste squares/Source coordinates')
    end
else
    parameter.lsmopt.min_num_obs_per_est_source = 0;
end





% prepare N and b for global solution
parameter.lsmopt.global_solve=get(handles.checkbox_run_prepareGlobParam, 'Value');
parameter.lsmopt.est_vel=get(handles.checkbox_run_globalPram_statVel, 'Value');
parameter.lsmopt.est_source=get(handles.checkbox_run_globalPram_source, 'Value');
refEpochNum=str2double(get(handles.edit_run_globalPram_statVelEpoch, 'String'));
if isnan(refEpochNum)
    parameter.lsmopt.refvel=2000;
else
    parameter.lsmopt.refvel=refEpochNum;
end
if get(handles.checkbox_run_globalPram_axisOffset, 'Value')
    parameter.lsmopt.est_AO=1;
else
    parameter.lsmopt.est_AO=0;
end

 % seasonal variations of station positions
if get(handles.checkbox_run_globalPram_stseaspos, 'Value')
    parameter.lsmopt.est_stsespos =1;
else
    parameter.lsmopt.est_stsespos =0;
end


% Tidal ERP variations
parameter.lsmopt.est_erpsp = 0;
if get(handles.checkbox_run_globalPram_tidERPvar, 'Value')
    parameter.lsmopt.est_erpsp = 1;
end

% Pole tide Love and Shida number
if get(handles.checkbox_run_globalPram_hlpole, 'Value')
    parameter.lsmopt.est_hpole =1;
    parameter.lsmopt.est_lpole =1;
else
    parameter.lsmopt.est_hpole =0;
    parameter.lsmopt.est_lpole =0;
end


% APL regression coefficients
parameter.lsmopt.est_rg =0;
if get(handles.checkbox_run_globalPram_APLrg, 'Value')
    parameter.lsmopt.est_rg =1;   
end



% sinex output
parameter.lsmopt.ascii_snx=get(handles.checkbox_run_sinex_write, 'Value');
parameter.lsmopt.outsnx.clk=get(handles.radiobutton_run_sinex_clockParam_incl, 'Value');
if get(handles.checkbox_estimation_leastSquares_tropo_zwd, 'Value')
    parameter.lsmopt.outsnx.zwd=get(handles.radiobutton_run_sinex_zwd_incl, 'Value');
else
    parameter.lsmopt.outsnx.zwd=0;
end
if get(handles.checkbox_estimation_leastSquares_tropo_ngr, 'Value') && ...
        get(handles.checkbox_estimation_leastSquares_tropo_egr, 'Value')
    parameter.lsmopt.outsnx.tgr=get(handles.radiobutton_run_sinex_tropoParam_incl, 'Value');
else
    parameter.lsmopt.outsnx.tgr=0;
end
if get(handles.checkbox_run_sinex_sources, 'Value')
    parameter.lsmopt.outsnx.sou=get(handles.radiobutton_run_sinex_sources_incl, 'Value');
else
    parameter.lsmopt.outsnx.sou=0;
end
if get(handles.checkbox_estimation_leastSquares_coordinates_estimate, 'Value')
    parameter.lsmopt.outsnx.xyz=get(handles.radiobutton_run_sinex_stationCoords_incl, 'Value');
else
    parameter.lsmopt.outsnx.xyz=0;
end
if (get(handles.checkbox_estimation_leastSquares_eop_xpEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_ypEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_dut1Est, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_nutdxEst, 'Value')+...
        get(handles.checkbox_estimation_leastSquares_eop_nutdyEst, 'Value'))==0
    parameter.lsmopt.outsnx.eop=0;
else
    parameter.lsmopt.outsnx.eop=get(handles.radiobutton_run_sinex_eop_incl, 'Value');
end

parameter.lsmopt.outsnx.bdco = 0; % 0 = reduce bas-dep clk offset in sinex

parameter.lsmopt.addSnxSource=get(handles.checkbox_run_sinex_sources, 'Value');

if get(handles.checkbox_run_sinex_changeAnalystsName, 'Value')
    parameter.lsmopt.outsnx.changeAnalystsName=1;
    parameter.lsmopt.outsnx.firstname=get(handles.edit_run_sinex_firstname, 'String');
    parameter.lsmopt.outsnx.lastname=get(handles.edit_run_sinex_lastname, 'String');
    parameter.lsmopt.outsnx.email=get(handles.edit_run_sinex_email, 'String');
else
    parameter.lsmopt.outsnx.changeAnalystsName=0;
    parameter.lsmopt.outsnx.firstname='';
    parameter.lsmopt.outsnx.lastname='';
    parameter.lsmopt.outsnx.email='';
end

if get(handles.checkbox_run_sinex_addSuffix, 'Value')
    parameter.lsmopt.outsnx.addSuffix=1;
    parameter.lsmopt.outsnx.suffix=get(handles.edit_run_sinex_suffix, 'String');
else
    parameter.lsmopt.outsnx.addSuffix=0;
    parameter.lsmopt.outsnx.suffix='';
end

% EOP estimation
% xpol
parameter.lsmopt.xpol.model=get(handles.checkbox_estimation_leastSquares_eop_xpEst, 'Value');
parameter.lsmopt.xpol.constrain=get(handles.checkbox_estimation_leastSquares_eop_constr_xp, 'Value');
xpolIntNum=str2double(get(handles.edit_estimation_leastSquares_eop_interval_xp, 'String'));
if isnan(xpolIntNum)
    parameter.lsmopt.xpol.int=1440;
else
    parameter.lsmopt.xpol.int=xpolIntNum;
end
xpolConstrNum=str2double(get(handles.edit_estimation_leastSquares_constr_xp, 'String'));
if isnan(xpolConstrNum)
    parameter.lsmopt.xpol.coef=1.0e-4;
else
    parameter.lsmopt.xpol.coef=xpolConstrNum;
end

% ypol
parameter.lsmopt.ypol.model=get(handles.checkbox_estimation_leastSquares_eop_ypEst, 'Value');
parameter.lsmopt.ypol.constrain=get(handles.checkbox_estimation_leastSquares_eop_constr_yp, 'Value');
ypolIntNum=str2double(get(handles.edit_estimation_leastSquares_eop_interval_yp, 'String'));
if isnan(ypolIntNum)
    parameter.lsmopt.ypol.int=1440;
else
    parameter.lsmopt.ypol.int=xpolIntNum;
end
ypolConstrNum=str2double(get(handles.edit_estimation_leastSquares_constr_yp, 'String'));
if isnan(ypolConstrNum)
    parameter.lsmopt.ypol.coef=1.0e-4;
else
    parameter.lsmopt.ypol.coef=ypolConstrNum;
end

% dut1
parameter.lsmopt.dut1.model=get(handles.checkbox_estimation_leastSquares_eop_dut1Est, 'Value');
parameter.lsmopt.dut1.constrain=get(handles.checkbox_estimation_leastSquares_eop_constr_dut1, 'Value');
dut1IntNum=str2double(get(handles.edit_estimation_leastSquares_eop_interval_dut1, 'String'));
if isnan(dut1IntNum)
    parameter.lsmopt.dut1.int=1440;
else
    parameter.lsmopt.dut1.int=dut1IntNum;
end
dut1ConstrNum=str2double(get(handles.edit_estimation_leastSquares_constr_dut1, 'String'));
if isnan(dut1ConstrNum)
    parameter.lsmopt.dut1.coef=1.0e-4;
else
    parameter.lsmopt.dut1.coef=dut1ConstrNum;
end
% nutdx
parameter.lsmopt.nutdx.model=get(handles.checkbox_estimation_leastSquares_eop_nutdxEst, 'Value');
parameter.lsmopt.nutdx.constrain=get(handles.checkbox_estimation_leastSquares_eop_constr_nutdx, 'Value');
nutdxIntNum=str2double(get(handles.edit_estimation_leastSquares_eop_interval_nutdx, 'String'));
if isnan(nutdxIntNum)
    parameter.lsmopt.nutdx.int=1440;
else
    parameter.lsmopt.nutdx.int=nutdxIntNum;
end
nutdxConstrNum=str2double(get(handles.edit_estimation_leastSquares_constr_nutdx, 'String'));
if isnan(nutdxConstrNum)
    parameter.lsmopt.nutdx.coef=1.0e-4;
else
    parameter.lsmopt.nutdx.coef=nutdxConstrNum;
end

% nutdy
parameter.lsmopt.nutdy.model=get(handles.checkbox_estimation_leastSquares_eop_nutdyEst, 'Value');
parameter.lsmopt.nutdy.constrain=get(handles.checkbox_estimation_leastSquares_eop_constr_nutdy, 'Value');
nutdyIntNum=str2double(get(handles.edit_estimation_leastSquares_eop_interval_nutdy, 'String'));
if isnan(nutdyIntNum)
    parameter.lsmopt.nutdy.int=1440;
else
    parameter.lsmopt.nutdy.int=nutdyIntNum;
end
nutdyConstrNum=str2double(get(handles.edit_estimation_leastSquares_constr_nutdy, 'String'));
if isnan(nutdyConstrNum)
    parameter.lsmopt.nutdy.coef=1.0e-4;
else
    parameter.lsmopt.nutdy.coef=nutdyConstrNum;
end

% outlier handling
% parameter
outParamNum=str2double(get(handles.edit_run_outlierTestC, 'String'));
if isnan(outParamNum)
    parameter.lsmopt.par_outlier=5;
else
    parameter.lsmopt.par_outlier=outParamNum;
end
parameter.lsmopt.simple_outlier=get(handles.checkbox_run_simpleOutlierTest, 'Value');
parameter.lsmopt.basic_outlier=get(handles.checkbox_run_normalOutlierTest, 'Value');

% define if parameters should be estimated
parameter.lsmopt.est_singleses=get(handles.checkbox_run_estParameters, 'Value');

% allow station/sessionwise parameterization
parameter.lsmopt.control_gui_vie_lsm=get(handles.checkbox_run_allowStationwise, 'Value');



% save the parameter file to file (before: create path if not exists)
allSlashFindings=sort([strfind(fullOutName, '\'), strfind(fullOutName, '/')]); % find \ and /
if ~isempty(allSlashFindings)
    outPath=fullOutName(1:allSlashFindings(end));
    if ~exist(outPath, 'dir')
        mkdir(outPath)
    end
end


% do you want to apply antenna axis offset altitude correction a priori in
% vie_mod?
parameter.vie_mod.aoaltcorr = 0;

% estimation of single session parameters in VIE_LSM
parameter.lsmopt.est_scale = 0;


% for which parameters which are not in GUI yet should the NEQ be prepared for the estimation in VIE_GLOB? (1) / (0)?
parameter.lsmopt.est_love = 0;
parameter.lsmopt.est_shida = 0;
parameter.lsmopt.est_FCNset = 0; % FCN period from Solid Earth tides
parameter.lsmopt.est_accSSB = 0; % acceleration of the Solar System Barycentre
parameter.lsmopt.est_source_velo = 0;
parameter.lsmopt.est_gamma = 0;
parameter.lsmopt.est_FCNnut = 0; % FCN period from nutation
parameter.lsmopt.est_FCNnutAmpl = 0; % FCN amplitudes from nutation


save(fullOutName, 'parameter');