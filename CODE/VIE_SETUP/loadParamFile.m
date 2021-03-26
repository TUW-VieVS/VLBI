% #########################################################################
% #     loadParamFile
% #########################################################################
%
% DESCRITPION
% This function loads the structure "parameter", which contains all
% required parameters to run the VieVS bath mode, from a *.mat file.
%
% AUTHOR
%   Matthias Madzak
%
% INPUT
%   handles      structure from the GUI (also containing data, eg residuals)
%   hObject      ...unused so far.
%   fullOutName  filename of loaded parameter file (including filepath!)
%
% OUTPUT
%
% CHANGES
%   2014-11-05, A. Hellerschmied: Exctact code from "vievs2_3.m" file and
%     create this m-file
%   2014-11-05, A. Hellerschmied: Parameter "use_opt_files" added.
%   2014-12-03, M. Madzak & Anastasiia: Bug-fixes
%   2014-12-09, A. Hellerschmied: Add field "use_opt_files" to struct "parameter.vie_init"
%   2015-01-08, M. Madzak & Anastasiia: Minor bug-fix
%   2015-01-14, M. Madzak: Minor bug fix
%   2015-02-03, A. Hellerschmied: Modification related to GUI update and
%       bug fix
%   2015-07-14, A. Hellerschmied: New mean pole model added (IERS 2015 model, "cmp2015")
%   2016-02-29, D. Landskron: the path of external tropo files is now also loaded
%   2016-07-07, A. Hellerschmied: Added checkbox "checkbox_estimation_leastSquares_sources_ICRF2_def"
%   2016-07-11, D. Mayer: Added edit field "edit_estimation_leastSquares_sources_abs_constr"
%   2016-08-03, A. Hellerschmied: Added sedit field "edit_models_sc_sp3_file"
%   2016-08-31, A. Girdiuk: New ephemerides is added
%   2016-12-20, A. Hellerschmied: "checkbox_run_bslDepWeights" and "checkbox_run_ElevDepNoise" added
%   2017-01-16, H. Krasna: "checkbox_run_globalPram_axisOffset",
%       "checkbox_run_globalPram_stseaspos",
%       "checkbox_run_globalPram_tidERPvar", and
%       "checkbox_run_globalPram_hlpole" added
%   2017/01/17, D. Landskron: choice of temperature source transferred from "Models - Troposphere" to "Models - Station models", and bug fixed
%   2017/01/23, D. Landskron: general shape of Troposphere changed
%   2017/02/07, D. Mayer: added exception handling to troposphere pannel
%       --> old parameter files can still be used
%   2017/02/19, H. Krasna: checkbox_run_globalPram_APLrg, popupmenu_parameters_statCorr_APLrg,
%                          checkbox_parameters_statCorr_APLrg, checkbox_parameters_statCorr_GIA,
%                          popupmenu_parameters_statCorr_GIA added
%   2017-02-23, A. Hellerschmied: "handles.edit_run_add_noise" added
%   2017-03-01, A. Hellerschmied: "handles.checkbox_run_runOptions_manuallyFindBreaks" removed
%   2017-03-10, D. Landskron: bug corrected with the removing of manually find clock breaks
%   2017-03-20, A. Hellerschmied: Changes required for new EOP data loading routine
%   2017-09-13, D. Landskron: 'tropSource' shifted into 'vie_init' 
%   2018-01-11, D. Landskron: external troposphere modeling removed
%   2018-07-06, D. Landskron: VMF3 added to the troposphere models 
%   2018-11-29, D. Landskron: structure of observation restrictions standardized
%   2018-03-06, D. Landskron: suffix checkbox added to the sinex files
%
function loadParamFile(hObject, handles, fullFileName)

% load parameter file
load(fullFileName);

% VIE_INIT
% ================
% set opt directory if there is a directory like this
if sum(strcmp(get(handles.popupmenu_setInput_optDir, 'String'), parameter.vie_init.diropt))==1
    set(handles.popupmenu_setInput_optDir, 'Value', find(strcmp(get(handles.popupmenu_setInput_optDir, 'String'), parameter.vie_init.diropt)));
else
    msgbox('OPT directory not found!', 'Warning', 'warn');
end

% set use-OPT-file option

% "Older" parameter files may not contain this option => set to "1" and
% print warning msg!
if ~isfield(parameter.vie_init, 'use_opt_files')
    parameter.vie_init.use_opt_files = 1;
    parameter.lsmopt.use_opt_files = 1;   % NOTE: This is a Vie_LSM parameter!!
    msgbox('Loaded parameter file (from older VieVS version) does not contain the "use OPT file" option!', 'Warning', 'warn');
end
set(handles.checkbox_setInput_useOptFiles, 'Value', parameter.vie_init.use_opt_files);
% Enable/Disable popupmenu to select the OPT file
if parameter.vie_init.use_opt_files == 0
    set(handles.popupmenu_setInput_optDir, 'Enable', 'off')
else
    set(handles.popupmenu_setInput_optDir, 'Enable', 'on')
end


% set outlier directory if there is a directory like this
if sum(strcmp(get(handles.popupmenu_setInput_outDir, 'String'), parameter.vie_init.dirout))==1
    set(handles.popupmenu_setInput_outDir, 'Value', find(strcmp(get(handles.popupmenu_setInput_outDir, 'String'), parameter.vie_init.dirout)));
else
    % take empty when saved parameter is '' (only ' ' is found)
    if strcmp(parameter.vie_init.dirout, '')
        set(handles.popupmenu_setInput_outDir, 'Value', find(strcmp(get(handles.popupmenu_setInput_outDir, 'String'), ' ')));
    else
        msgbox('Outlier directory not found!', 'Warning', 'warn');
    end
end

% set remove outlier option
set(handles.checkbox_setInput_eliminOutliers, 'Value', parameter.vie_init.rm_outlier)

% set TRF
% if superstations (.mat) file was taken
if strcmp(parameter.vie_init.trf{1}(end-3:end), '.mat')
    % if file exists
    if exist(parameter.vie_init.trf{1}, 'file')
        % load file and set popupmenu
        handles=loadSuperstationFile(hObject, handles, parameter.vie_init.trf{1});
        
        % select chosen trf (field)
        logVTRF2008Found=~cellfun(@isempty, strfind(get(handles.popupmenu_parameters_refFrames_superstationTRF, 'String'), parameter.vie_init.trf{2}));
        if sum(logVTRF2008Found)>0
            set(handles.popupmenu_parameters_refFrames_superstationTRF, 'Value', find(logVTRF2008Found));
        end
    else
        msgbox('TRF file does not exist!\nBe sure to select a proper one!', 'TRF file not found', 'warn');
    end
    set(handles.radiobutton_parameters_refFrames_superstationTRF, 'Value', 1)
    set(handles.pushbutton_parameters_refFrames_superstationTRF_chose, 'Enable', 'On')
    set(handles.pushbutton_parameters_refFrames_superstationTRF_create, 'Enable', 'On')
    set(handles.popupmenu_parameters_refFrames_superstationTRF, 'Enable', 'On')
    set(handles.popupmenu_parameters_refFrames_otherTRF, 'Enable', 'Off')
else % txt trf file was taken
    % try to find it in "other-popupmenu"
    logOtherTrf=~cellfun(@isempty, strfind(get(handles.popupmenu_parameters_refFrames_otherTRF, 'String'), parameter.vie_init.trf{1}(8:end)));
    if sum(logOtherTrf)>0
        set(handles.popupmenu_parameters_refFrames_otherTRF, 'Value', find(logOtherTrf));
    else
        msgbox('Other (.txt) TRF was not found.\nProbably the file does not exist anymore. Be sure to select a proper one!', 'TRF file not found', 'warn');
    end
    set(handles.radiobutton_parameters_refFrames_otherTRF, 'Value', 1)
    set(handles.pushbutton_parameters_refFrames_superstationTRF_chose, 'Enable', 'Off')
    set(handles.pushbutton_parameters_refFrames_superstationTRF_create, 'Enable', 'Off')
    set(handles.popupmenu_parameters_refFrames_superstationTRF, 'Enable', 'Off')
    set(handles.popupmenu_parameters_refFrames_otherTRF, 'Enable', 'On')
    
end

% for debug:
% parameter.vie_init.crf{1}
% ans =
% ../CRF/supersource.mat

% parameter.vie_init.crf{2}
% ans =
% icrf2

% set CRF
% if supersource file was chosen
if strcmpi(parameter.vie_init.crf{1}(end-2:end), 'mat')
    if exist(parameter.vie_init.crf{1}, 'file')
        % load file and set popupmenu
        handles=loadSupersourceFile(hObject, handles, parameter.vie_init.crf{1});
        
        % select chosen crf (field)
        %logCRFFound=~cellfun(@isempty, strfind(get(handles.popupmenu_parameters_refFrames_supersourceCRF, 'String'), parameter.vie_init.crf{2}));
        logCRFFound=strcmp(get(handles.popupmenu_parameters_refFrames_supersourceCRF, 'String'), parameter.vie_init.crf{2});
        if sum(logCRFFound)>0
            set(handles.popupmenu_parameters_refFrames_supersourceCRF, 'Value', find(logCRFFound));
        end
    else
        msgbox('CRF file does not exist!\nBe sure to select a proper one!', 'CRF file not found', 'warn');
    end
    set(handles.radiobutton_parameters_refFrames_supersourceCRF, 'Value', 1)
    set(handles.pushbutton_parameters_refFrames_superstationCRF_chose, 'Enable', 'On')
    set(handles.pushbutton_parameters_refFrames_createSupersourceFile, 'Enable', 'On')
    set(handles.popupmenu_parameters_refFrames_supersourceCRF, 'Enable', 'On')
    set(handles.popupmenu_parameters_refFrames_otherCRF, 'Enable', 'Off')
else
    % manual CRF was chosen
    logOtherCrf=~cellfun(@isempty, strfind(get(handles.popupmenu_parameters_refFrames_otherCRF, 'String'), parameter.vie_init.crf{1}(8:end)));
    if sum(logOtherCrf)>0
        set(handles.popupmenu_parameters_refFrames_otherCRF, 'Value', find(logOtherCrf));
    else
        msgbox('Other (.txt) CRF was not found.\nProbably the file does not exist anymore. Be sure to select a proper one!', 'CRF file not found', 'warn');
    end
    set(handles.radiobutton_parameters_refFrames_otherCRF, 'Value', 1)
    set(handles.pushbutton_parameters_refFrames_superstationCRF_chose, 'Enable', 'Off')
    set(handles.pushbutton_parameters_refFrames_createSupersourceFile, 'Enable', 'Off')
    set(handles.popupmenu_parameters_refFrames_supersourceCRF, 'Enable', 'Off')
    set(handles.popupmenu_parameters_refFrames_otherCRF, 'Enable', 'On')
end


% set min elevation and quality max
set(handles.edit_parameter_obsRestr_cutOff, 'String', parameter.obs_restrictions.cut_off_elev*180/pi)
set(handles.edit_parameter_obsRestr_qualityCode, 'String', parameter.obs_restrictions.Qlim)

% set tropodelay model
try
    switch parameter.vie_init.tropSource.name
        
        case 'indModeling'   % individual modeling
            set(handles.radiobutton_parameters_troposphere_indModeling, 'Value', 1)
            set(handles.radiobutton_parameters_troposphere_raytr, 'Value', 0)
            
            % set en/disable
            set(handles.radiobutton_parameters_troposphere_zhd_no, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_zhd_fromInSitu, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_zhd_VMF3, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_zhd_VMF1, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_zhd_GPT3, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_zwd_no, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_zwd_fromInSitu, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_zwd_VMF3, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_zwd_VMF1, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_zwd_GPT3, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_mfh_VMF3, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_mfh_VMF1, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_mfh_GPT3, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_mfw_VMF1, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_mfw_GPT3, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_gradients_h_no, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_gradients_h_GRAD, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_gradients_h_GPT3, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_gradients_h_DAO, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_gradients_w_no, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_gradients_w_GRAD, 'Enable', 'on')
            set(handles.radiobutton_parameters_troposphere_gradients_w_GPT3, 'Enable', 'on')
            
        case 'raytr'
            set(handles.radiobutton_parameters_troposphere_indModeling, 'Value', 0)
            set(handles.radiobutton_parameters_troposphere_raytr, 'Value', 1)
            
            % set en/disable
            set(handles.radiobutton_parameters_troposphere_zhd_no, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_zhd_fromInSitu, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_zhd_VMF3, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_zhd_VMF1, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_zhd_GPT3, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_zwd_no, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_zwd_fromInSitu, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_zwd_VMF3, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_zwd_VMF1, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_zwd_GPT3, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_mfh_VMF3, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_mfh_VMF1, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_mfh_GPT3, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_mfw_VMF1, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_mfw_GPT3, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_gradients_h_no, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_gradients_h_GRAD, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_gradients_h_GPT3, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_gradients_h_DAO, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_gradients_w_no, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_gradients_w_GRAD, 'Enable', 'off')
            set(handles.radiobutton_parameters_troposphere_gradients_w_GPT3, 'Enable', 'off')
            
    end
catch
    
    % set en/disable
    set(handles.radiobutton_parameters_troposphere_zhd_no, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_zhd_fromInSitu, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_zhd_VMF3, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_zhd_VMF1, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_zhd_GPT3, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_zwd_no, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_zwd_fromInSitu, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_zwd_VMF3, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_zwd_VMF1, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_zwd_GPT3, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_mfh_VMF3, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_mfh_VMF1, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_mfh_GPT3, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_mfw_VMF1, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_mfw_GPT3, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_gradients_h_no, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_gradients_h_GRAD, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_gradients_h_GPT3, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_gradients_h_DAO, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_gradients_w_no, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_gradients_w_GRAD, 'Enable', 'on')
    set(handles.radiobutton_parameters_troposphere_gradients_w_GPT3, 'Enable', 'on')
    
    set(handles.radiobutton_parameters_troposphere_indModeling, 'Value', 1)
    set(handles.radiobutton_parameters_troposphere_raytr, 'Value', 0)
            
    warning('Setting of troposphere delay model failed --> individual models used');
    
end

% set iono delay source:
if isfield(parameter.vie_init, 'iono') % older parameter files do not have that field
    switch parameter.vie_init.iono
        case 'ngs'
            set(handles.radiobutton_parameters_iono_fromNGS, 'Value', 1)
            set(handles.popupmenu_parameters_iono_ext, 'Enable', 'off')
            set(handles.pushbutton_parameters_iono_create, 'Enable', 'off')
        case 'observation_database'
            set(handles.radiobutton_parameters_iono_fromNGS, 'Value', 1)
            set(handles.popupmenu_parameters_iono_ext, 'Enable', 'off')
            set(handles.pushbutton_parameters_iono_create, 'Enable', 'off')
        otherwise % external iono delay
            set(handles.radiobutton_parameters_iono_ext, 'Value', 1)
            set(handles.popupmenu_parameters_iono_ext, 'Enable', 'on')
            set(handles.pushbutton_parameters_iono_create, 'Enable', 'on')
            % if the external iono file-folder is found: set this as selected
            if sum(strcmp(get(handles.popupmenu_parameters_iono_ext, 'String'), parameter.vie_init.ionoFolder))==1
                set(handles.popupmenu_parameters_iono_ext, 'Value', find(strcmp(get(handles.popupmenu_parameters_iono_ext, 'String'), parameter.vie_init.ionoFolder)));
            else
                msgbox(sprintf('Folder of external ionosperic file not found!\nBe sure to select a proper one!'), 'Warning', 'warn');
            end
            
    end
else
    set(handles.radiobutton_parameters_iono_fromNGS, 'Value', 1)
end

% Sp3 filename and path
if isfield(parameter.vie_init, 'sc_orbit_file_path_name')
    set(handles.edit_models_sc_sp3_file, 'String', parameter.vie_init.sc_orbit_file_path_name)
end



% VIE_MOD
% ================
% ephemerides
switch parameter.vie_mod.eph
    case 'jpl_430'
        set(handles.radiobutton_parameters_eph_jpl430, 'Value', 1)
    case 'jpl_421'
        set(handles.radiobutton_parameters_eph_jpl421, 'Value', 1)
    case 'jpl_405'
        set(handles.radiobutton_parameters_eph_jpl405, 'Value', 1)
end

% EOP file
switch parameter.vie_mod.EOPfile
    case {'C04 14', 'C04_14_1962_now.txt'}
        set(handles.radiobutton_parameters_eop_aPriori_08C04, 'Value', 1)
        set(handles.popupmenu_parameters_eop_aPriori_other, 'Enable', 'off')
    case {'C04 08', 'C04_08_1962_now.txt'}
        set(handles.radiobutton_parameters_eop_aPriori_05C04, 'Value', 1)
        set(handles.popupmenu_parameters_eop_aPriori_other, 'Enable', 'off')
    case {'eop_finals2000A.txt', 'finals_all_IAU2000.txt'}
        set(handles.radiobutton_parameters_eop_aPriori_finals, 'Value', 1)
        set(handles.popupmenu_parameters_eop_aPriori_other, 'Enable', 'off')
    case 'C04 05'
        warning('EOP C04 05 time series not available any more. C04 14 is used instead.');
        set(handles.radiobutton_parameters_eop_aPriori_08C04, 'Value', 1)
        set(handles.popupmenu_parameters_eop_aPriori_other, 'Enable', 'off')
    otherwise % external EOP file
        set(handles.radiobutton_parameters_eop_aPriori_other, 'Value', 1)
        set(handles.popupmenu_parameters_eop_aPriori_other, 'Enable', 'on')
        % try to find the proper textfile in the popupmenu
        allEntriesInPopup=get(handles.popupmenu_parameters_eop_aPriori_other,'String');
        foundInd=find(strcmpi(allEntriesInPopup, parameter.vie_mod.EOPfile));
        if isempty(foundInd)
            msgbox('Chosen special EOP file was not found in the EOP directory!','EOP file not found');
        else
            set(handles.popupmenu_parameters_eop_aPriori_other, 'Value', foundInd);
        end
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




% a priori nutation offsets
set(handles.checkbox_parameters_eop_models_inclAPrioriNutOffs, 'Value', parameter.vie_mod.dXdY)

% interpolation
set(handles.radiobutton_parameters_eop_interp_lin, 'Value', parameter.vie_mod.linear)
if isfield(parameter.vie_mod,'linear48h')
    set(handles.checkbox_parameters_eop_interp_lin48h, 'Value', parameter.vie_mod.linear48h)
else
    set(handles.checkbox_parameters_eop_interp_lin48h, 'Value', 0)
end
set(handles.radiobutton_parameters_eop_interp_lag, 'Value', ~parameter.vie_mod.linear)

% tidal UT
set(handles.checkbox_parameters_eop_tidalUtVariations, 'Value', parameter.vie_mod.tidalUT)
if parameter.vie_mod.tidalUT
    set(handles.rb_eop_inter_UT1R, 'Enable', 'on')
    set(handles.rb_eop_inter_UT1S, 'Enable', 'on')
else
    set(handles.rb_eop_inter_UT1R, 'Enable', 'off')
    set(handles.rb_eop_inter_UT1S, 'Enable', 'off')
end
set(handles.rb_eop_inter_UT1R, 'Value', parameter.vie_mod.tidalUT35)
if parameter.vie_mod.tidalUT35==0
    set(handles.rb_eop_inter_UT1S,'Value',1);
end
% eopoc
% if model == 'none' -> untick 'ocean tides'
if strcmp(parameter.vie_mod.eopoc, 'none')
    set(handles.checkbox_parameters_eop_inclHf_oceanTidesOther, 'Value', 0);
    set(handles.popupmenu_parameters_eop_oceanTideModel, 'Enable', 'off');
else
    % try to find entry in popupmenu
    indChosenEopoc=find(strcmp(get(handles.popupmenu_parameters_eop_oceanTideModel, 'String'), parameter.vie_mod.eopoc));
    if ~isempty(indChosenEopoc)
        set(handles.popupmenu_parameters_eop_oceanTideModel, 'Value', indChosenEopoc)
    else
        msgbox('Chosen Eopoc file not found. Be sure to select a proper one.', 'Manual eopoc file', 'warn');
    end
end

% libration
set(handles.checkbox_parameters_eop_inclHf_LibrationXpYp, 'Value', parameter.vie_mod.lib_pm)
set(handles.checkbox_parameters_eop_inclHf_LibrationUt1, 'Value', parameter.vie_mod.lib_ut)

% precession/nutation model
if strcmp(parameter.vie_mod.nutmod, 'IAU_2006/2000A')
    set(handles.radiobutton_parameters_eop_precNutModel_20062000A, 'value', 1)
else
    set(handles.radiobutton_parameters_eop_precNutModel_2000a, 'value', 1)
end

% station corrections
% solid earth tides
set(handles.checkbox_parameters_statCorr_solidEarthTides, 'Value', parameter.vie_mod.cts)

% ocean tides
set(handles.checkbox_parameters_statCorr_tidalOceanLoad, 'value', parameter.vie_mod.cto)
if parameter.vie_mod.cto==1
    set(handles.popupmenu_parameters_statCorr_tidalOceanLoad, 'Enable', 'on')
else
    set(handles.popupmenu_parameters_statCorr_tidalOceanLoad, 'Enable', 'off')
end
% find chosen model in current popupmenu
indOfModelFound=find(strcmp(get(handles.popupmenu_parameters_statCorr_tidalOceanLoad, 'string'), parameter.vie_mod.ocm));
if ~isempty(indOfModelFound)
    set(handles.popupmenu_parameters_statCorr_tidalOceanLoad, 'Value', indOfModelFound)
else
    msgbox(sprintf('Tidal ocean loading model file not found.\nBe sure to select a proper one!'), 'Warning!', 'warn');
end

% tidal atmosphere loading
set(handles.checkbox_parameters_statCorr_tidalAtmoLoad, 'Value', parameter.vie_mod.cta)
if parameter.vie_mod.cta==1
    set(handles.popupmenu_parameters_statCorr_tidalAtmoOceanLoad, 'Enable', 'on')
else
    set(handles.popupmenu_parameters_statCorr_tidalAtmoOceanLoad, 'Enable', 'off')
end
% find chosen model in current popupmenu
indOfModelFound=find(strcmp(get(handles.popupmenu_parameters_statCorr_tidalAtmoOceanLoad, 'string'), parameter.vie_mod.ctam));
if ~isempty(indOfModelFound)
    set(handles.popupmenu_parameters_statCorr_tidalAtmoOceanLoad, 'Value', indOfModelFound)
else
    msgbox(sprintf('Tidal atmospheric loading model file not found.\nBe sure to select a proper one!'), 'Warning!', 'warn');
end

% non tidal atmosphere loading
set(handles.checkbox_parameters_statCorr_nonTidalAtmoLoad, 'Value', parameter.vie_mod.cnta)
if parameter.vie_mod.cnta==1
    set(handles.popupmenu_parameters_statCorr_nonTidalAtmoOceanLoad, 'Enable', 'on')
else
    set(handles.popupmenu_parameters_statCorr_nonTidalAtmoOceanLoad, 'Enable', 'off')
end
% find chosen model in current popupmenu
indOfModelFound=find(strcmp(get(handles.popupmenu_parameters_statCorr_nonTidalAtmoOceanLoad, 'string'), parameter.vie_mod.cntam));
if ~isempty(indOfModelFound)
    set(handles.popupmenu_parameters_statCorr_nonTidalAtmoOceanLoad, 'Value', indOfModelFound)
else
    msgbox(sprintf('Non-tidal atmosphere loading model file not found.\nBe sure to select a proper one!'), 'Warning!', 'warn');
end

% pole tides
set(handles.checkbox_parameters_statCorr_poleTides, 'Value', parameter.vie_mod.ctp)
set(handles.checkbox_parameters_statCorr_oceanPoleTides, 'Value', parameter.vie_mod.ctop)

% mean pole model
if strcmp(parameter.vie_mod.ctpm, 'cubic')
    set(handles.radiobutton_parameters_statCorr_poleModel_cub, 'Value', 1)
elseif strcmp(parameter.vie_mod.ctpm, 'linear')
    set(handles.radiobutton_parameters_statCorr_poleModel_lin, 'Value', 1)
elseif strcmp(parameter.vie_mod.ctpm, 'cmp2015')
    set(handles.radiobutton_parameters_statCorr_poleModel_iers2015, 'Value', 1)
end
% set en/disable
if parameter.vie_mod.ctp+parameter.vie_mod.ctop==0
    set(handles.radiobutton_parameters_statCorr_poleModel_lin, 'Enable', 'off')
    set(handles.radiobutton_parameters_statCorr_poleModel_cub, 'Enable', 'off')
    set(handles.radiobutton_parameters_statCorr_poleModel_iers2015, 'Enable', 'off')
else
    set(handles.radiobutton_parameters_statCorr_poleModel_lin, 'Enable', 'on')
    set(handles.radiobutton_parameters_statCorr_poleModel_cub, 'Enable', 'on')
    set(handles.radiobutton_parameters_statCorr_poleModel_iers2015, 'Enable', 'on')
end

% hydrology loading
if isfield(parameter.vie_mod, 'chl') % this questin could be deleted when we are sure that no old parameter files (without this field) are used!!
    set(handles.checkbox_parameters_statCorr_hydroLoading, 'Value', parameter.vie_mod.chl)
    if parameter.vie_mod.chl==1
        set(handles.popupmenu_parameters_statCorr_hydroLoading, 'Enable', 'on')
    else
        set(handles.popupmenu_parameters_statCorr_hydroLoading, 'Enable', 'off')
    end
    % find chosen model in current popupmenu
    indOfModelFound=find(strcmp(get(handles.popupmenu_parameters_statCorr_hydroLoading, 'string'), parameter.vie_mod.chlm));
    if ~isempty(indOfModelFound)
        set(handles.popupmenu_parameters_statCorr_hydroLoading, 'Value', indOfModelFound)
    else
        msgbox(sprintf('Hydrology loading model file not found.\nBe sure to select a proper one!'), 'Warning!', 'warn');
    end
else
    % if this model does not exist in parameter file (old parameter file)
    % set this model to 0 (=not use)
    set(handles.checkbox_parameters_statCorr_hydroLoading, 'Value', 0)
    % and set the fodler to chose to disable
    set(handles.popupmenu_parameters_statCorr_hydroLoading, 'Enable', 'off')
end

% thermal deformation
set(handles.checkbox_parameters_statCorr_thermalDef, 'Value', parameter.vie_mod.therm)
if strcmp(parameter.vie_init.tp, 'in situ')
    set(handles.checkbox_parameters_statCorr_temp_fromInSitu, 'Value', 1)
elseif strcmp(parameter.vie_init.tp, 'gpt3')
    set(handles.checkbox_parameters_statCorr_temp_GPT3, 'Value', 1)
end
% set en/disable
if parameter.vie_mod.therm==0
    set(handles.checkbox_parameters_statCorr_temp_fromInSitu, 'Enable', 'off')
    set(handles.checkbox_parameters_statCorr_temp_GPT3, 'Enable', 'off')
else
    set(handles.checkbox_parameters_statCorr_temp_fromInSitu, 'Enable', 'on')
    set(handles.checkbox_parameters_statCorr_temp_GPT3, 'Enable', 'on')
end

% gravitational deformation
if isfield(parameter.vie_mod, 'gravDef')
    set(handles.checkbox_parameters_statCorr_gravitationalDef, 'Value', parameter.vie_mod.gravDef)
else
    set(handles.checkbox_parameters_statCorr_gravitationalDef, 'Value', 0)
end

% APL regression coefficients
if isfield(parameter.vie_mod, 'crg') % this questin could be deleted when we are sure that no old parameter files (without this field) are used!!
    set(handles.checkbox_parameters_statCorr_APLrg, 'Value', parameter.vie_mod.crg)
    if parameter.vie_mod.crg==1
        set(handles.popupmenu_parameters_statCorr_APLrg, 'Enable', 'on')
    else
        set(handles.popupmenu_parameters_statCorr_APLrg, 'Enable', 'off')
    end
    % find chosen model in current popupmenu
    indOfModelFound=find(strcmp(get(handles.popupmenu_parameters_statCorr_APLrg, 'string'), parameter.vie_mod.crgm));
    if ~isempty(indOfModelFound)
        set(handles.popupmenu_parameters_statCorr_APLrg, 'Value', indOfModelFound)
    else
        msgbox(sprintf('APL regression coeff. file not found.\nBe sure to select a proper one!'), 'Warning!', 'warn');
    end
else
    % if this model does not exist in parameter file (old parameter file)
    % set this model to 0 (=not use)
    set(handles.checkbox_parameters_statCorr_APLrg, 'Value', 0)
    % and set the fodler to chose to disable
    set(handles.popupmenu_parameters_statCorr_APLrg, 'Enable', 'off')
end



% GIA uplift rates
if isfield(parameter.vie_mod, 'gia') % this questin could be deleted when we are sure that no old parameter files (without this field) are used!!
    set(handles.checkbox_parameters_statCorr_GIA, 'Value', parameter.vie_mod.gia)
    if parameter.vie_mod.gia==1
        set(handles.popupmenu_parameters_statCorr_GIA, 'Enable', 'on')
    else
        set(handles.popupmenu_parameters_statCorr_GIA, 'Enable', 'off')
    end
    % find chosen model in current popupmenu
    indOfModelFound=find(strcmp(get(handles.popupmenu_parameters_statCorr_GIA, 'string'), parameter.vie_mod.giam));
    if ~isempty(indOfModelFound)
        set(handles.popupmenu_parameters_statCorr_GIA, 'Value', indOfModelFound)
    else
        msgbox(sprintf('GIA uplift rates file not found.\nBe sure to select a proper one!'), 'Warning!', 'warn');
    end
else
    % if this model does not exist in parameter file (old parameter file)
    % set this model to 0 (=not use)
    set(handles.checkbox_parameters_statCorr_GIA, 'Value', 0)
    % and set the fodler to chose to disable
    set(handles.popupmenu_parameters_statCorr_GIA, 'Enable', 'off')
end




% zenith delay
try
    switch parameter.vie_init.zhd
        case 'no'
            set(handles.radiobutton_parameters_troposphere_zhd_no, 'Value', 1)
        case 'in situ'
            set(handles.radiobutton_parameters_troposphere_zhd_fromInSitu, 'Value', 1)
        case 'vmf3'
            set(handles.radiobutton_parameters_troposphere_zhd_VMF3, 'Value', 1)
        case 'vmf1'
            set(handles.radiobutton_parameters_troposphere_zhd_VMF1, 'Value', 1)
        case 'gpt3'
            set(handles.radiobutton_parameters_troposphere_zhd_GPT3, 'Value', 1)
    end
catch 
    set(handles.radiobutton_parameters_troposphere_zhd_fromInSitu, 'Value', 1)
    warning('Setting of ZHD failed --> in situ is set');
end
try
    switch parameter.vie_init.zwd
        case 'no'
            set(handles.radiobutton_parameters_troposphere_zwd_no, 'Value', 1)
        case 'in situ'
            set(handles.radiobutton_parameters_troposphere_zwd_fromInSitu, 'Value', 1)
        case 'vmf3'
            set(handles.radiobutton_parameters_troposphere_zwd_VMF3, 'Value', 1)
        case 'vmf1'
            set(handles.radiobutton_parameters_troposphere_zwd_VMF1, 'Value', 1)
        case 'gpt3'
            set(handles.radiobutton_parameters_troposphere_zwd_GPT3, 'Value', 1)
    end
catch 
    set(handles.radiobutton_parameters_troposphere_zwd_no, 'Value', 1)
    warning('Setting of ZWD failed --> no is set');
end

%  mapping function
try
    switch parameter.vie_mod.mfh
        case 'vmf3'
            set(handles.radiobutton_parameters_troposphere_mfh_VMF3, 'Value', 1)
        case 'vmf1'
            set(handles.radiobutton_parameters_troposphere_mfh_VMF1, 'Value', 1)
        case 'gpt3'
            set(handles.radiobutton_parameters_troposphere_mfh_GPT3, 'Value', 1)
    end
catch
    set(handles.radiobutton_parameters_troposphere_mfh_VMF1, 'Value', 1)
    warning('Setting of hydrostatic mapping function failed --> VMF1 is set');
end
try
    switch parameter.vie_mod.mfw
        case 'vmf3'
            set(handles.radiobutton_parameters_troposphere_mfw_VMF3, 'Value', 1)
        case 'vmf1'
            set(handles.radiobutton_parameters_troposphere_mfw_VMF1, 'Value', 1)
        case 'gpt3'
            set(handles.radiobutton_parameters_troposphere_mfw_GPT3, 'Value', 1)
    end
catch
    set(handles.radiobutton_parameters_troposphere_mfw_VMF1, 'Value', 1)
    warning('Setting of wet mapping function failed --> VMF1 is set');
end

% gradients
try
    switch parameter.vie_mod.apgm_h
        case 'no'
            set(handles.radiobutton_parameters_troposphere_gradients_h_no, 'Value', 1)
        case 'grad'
            set(handles.radiobutton_parameters_troposphere_gradients_h_GRAD, 'Value', 1)
        case 'gpt3'
            set(handles.radiobutton_parameters_troposphere_gradients_h_GPT3, 'Value', 1)
        case 'dao'
            set(handles.radiobutton_parameters_troposphere_gradients_h_DAO, 'Value', 1)
    end
catch
    set(handles.radiobutton_parameters_troposphere_gradients_h_no, 'Value', 1)
    warning('Setting of hydrostatic gradients failed --> no gradients are set');
end
try
    switch parameter.vie_mod.apgm_w
        case 'no'
            set(handles.radiobutton_parameters_troposphere_gradients_w_no, 'Value', 1)
        case 'grad'
            set(handles.radiobutton_parameters_troposphere_gradients_w_GRAD, 'Value', 1)
        case 'gpt3'
            set(handles.radiobutton_parameters_troposphere_gradients_w_GPT3, 'Value', 1)
    end
catch
    set(handles.radiobutton_parameters_troposphere_gradients_w_no, 'Value', 1)
    warning('Setting of wet gradients failed --> no gradients are set');
end



% source structure
% find chosen model in current popupmenu
indOfModelFound=find(strcmp(get(handles.popupmenu_parameters_ss_catalog, 'string'), parameter.vie_mod.sou_cat));
if ~isempty(indOfModelFound)
    set(handles.popupmenu_parameters_ss_catalog, 'Value', indOfModelFound)
else
    msgbox(sprintf('Source structure catalog file not found.\nBe sure to select a proper one!'), 'Warning!', 'warn');
end
set(handles.checkbox_parameters_ss_applyss, 'Value', parameter.vie_mod.ssou)
set(handles.checkbox_parameters_ss_writejet, 'Value', parameter.vie_mod.write_jet)


% VIE_LSM
% ================

% set output directory
% is not saved - so cannot be loaded! set(handles.edit_run_outDirs_oneSub, 'String', parameter.lsmopt.dirout)

% set first solution (0|1) and main solution (0|1)
set(handles.checkbox_run_runFirstSolution, 'Value', parameter.lsmopt.first)
if parameter.lsmopt.first==0
    % set disable
    set(handles.radiobutton_run_oneOffsPerClock, 'enable', 'off')
    set(handles.radiobutton_run_oneOffsAndRatePerClock, 'enable', 'off')
    set(handles.radiobutton_run_oneOffsAndRateAndQuPerClock, 'enable', 'off')
else
    % set enable
    set(handles.radiobutton_run_oneOffsPerClock, 'enable', 'on')
    set(handles.radiobutton_run_oneOffsAndRatePerClock, 'enable', 'on')
    set(handles.radiobutton_run_oneOffsAndRateAndQuPerClock, 'enable', 'on')
end

set(handles.checkbox_run_runMainSolution, 'Value', parameter.lsmopt.second)
if parameter.lsmopt.second==0
    % set disable
    set(handles.checkbox_run_simpleOutlierTest, 'Enable', 'off')
    set(handles.checkbox_run_normalOutlierTest, 'enable', 'off')
    set(handles.text_run_runOptions_c, 'enable', 'off')
    set(handles.edit_run_outlierTestC, 'Enable', 'off')
else
    % set enable
    set(handles.checkbox_run_simpleOutlierTest, 'Enable', 'on')
    set(handles.checkbox_run_normalOutlierTest, 'enable', 'on')
    if (parameter.lsmopt.simple_outlier+parameter.lsmopt.basic_outlier) >= 1
        set(handles.text_run_runOptions_c, 'enable', 'on')
        set(handles.edit_run_outlierTestC, 'Enable', 'on')
    else
        set(handles.text_run_runOptions_c, 'enable', 'off')
        set(handles.edit_run_outlierTestC, 'Enable', 'off')
    end
end

% set outlier test options
set(handles.checkbox_run_simpleOutlierTest, 'Value', parameter.lsmopt.simple_outlier)
set(handles.checkbox_run_normalOutlierTest, 'Value', parameter.lsmopt.basic_outlier)
set(handles.edit_run_outlierTestC, 'String', num2str(parameter.lsmopt.par_outlier))

% baseline dependent weighting
if isfield(parameter.lsmopt, 'bsldep')
    set(handles.checkbox_run_bslDepWeights, 'Value', parameter.lsmopt.bsldep)
else
    set(handles.checkbox_run_bslDepWeights, 'Value', 0)
end

% elevation dependent noise
if isfield(parameter.lsmopt, 'eldep_noise')
    set(handles.checkbox_run_ElevDepNoise, 'Value', parameter.lsmopt.eldep_noise)
else
    set(handles.checkbox_run_ElevDepNoise, 'Value', 0)
end

% add constant noise [cm]
if isfield(parameter.lsmopt, 'add_const_noise_cm')
    set(handles.edit_run_add_noise, 'String', num2str(parameter.lsmopt.add_const_noise_cm))
else
    set(handles.edit_run_add_noise, 'String', 1) % default: 1 cm
end


% set clock treatment of first solution
if parameter.lsmopt.firstclock==2
    set(handles.radiobutton_run_oneOffsAndRateAndQuPerClock, 'Value', 1)
elseif parameter.lsmopt.firstclock==1
    set(handles.radiobutton_run_oneOffsAndRatePerClock, 'Value', 1)
elseif parameter.lsmopt.firstclock==0
    set(handles.radiobutton_run_oneOffsPerClock, 'Value', 1)
end

% set clock treatment of main solution
set(handles.edit_estimation_leastSquares_clockInterval, 'String', parameter.lsmopt.int_clk)
set(handles.text308, 'String', sprintf('after %s minutes', num2str(parameter.lsmopt.int_clk)))
set(handles.checkbox_estimation_leastSquares_clocksRelConstr, 'Value', parameter.lsmopt.constr_clk)
set(handles.edit_estimation_leastSquares_clockConstr, 'String', num2str(parameter.lsmopt.coef_clk))
if parameter.lsmopt.pw_clk==0
    set(handles.checkbox_estimation_leastSquares_clocks, 'Value', 0)
    
    % set disable
    set(handles.radiobutton_estimation_leastSquares_pwlClock, 'Enable', 'off')
    set(handles.radiobutton_estimation_leastSquares_pwlAndRateClock, 'Enable', 'off')
    set(handles.radiobutton_estimation_leastSquares_pwlRateAndQuClock, 'Enable', 'off')
    set(handles.text_estimation_leastSquares_clockInt, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_clockInterval, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_clocksRelConstr, 'Enable', 'off')
    set(handles.text_estimation_leastSquares_clockConstr, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_clockConstr, 'Enable', 'off')
    set(handles.text308, 'Enable', 'off')
    snxState='off';
else % clock should be estimated
    % set estimation model (pwl/+rate/+quadratic)
    switch parameter.lsmopt.pw_clk
        case 1
            set(handles.radiobutton_estimation_leastSquares_pwlClock, 'Value', 1)
        case 2
            set(handles.radiobutton_estimation_leastSquares_pwlAndRateClock, 'Value', 1)
        case 3
            set(handles.radiobutton_estimation_leastSquares_pwlRateAndQuClock, 'Value', 1)
    end
    % set enable
    set(handles.radiobutton_estimation_leastSquares_pwlClock, 'Enable', 'on')
    set(handles.radiobutton_estimation_leastSquares_pwlAndRateClock, 'Enable', 'on')
    set(handles.radiobutton_estimation_leastSquares_pwlRateAndQuClock, 'Enable', 'on')
    set(handles.text_estimation_leastSquares_clockInt, 'Enable', 'on')
    set(handles.edit_estimation_leastSquares_clockInterval, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_clocksRelConstr, 'Enable', 'on')
    if get(handles.checkbox_estimation_leastSquares_clocksRelConstr, 'Value')==1
        set(handles.text_estimation_leastSquares_clockConstr, 'Enable', 'on')
        set(handles.edit_estimation_leastSquares_clockConstr, 'Enable', 'on')
        set(handles.text308, 'Enable', 'on')
    else
        set(handles.text_estimation_leastSquares_clockConstr, 'Enable', 'off')
        set(handles.edit_estimation_leastSquares_clockConstr, 'Enable', 'off')
        set(handles.text308, 'Enable', 'off')
    end
    if parameter.lsmopt.pw_clk==1
        set(handles.radiobutton_estimation_leastSquares_pwlClock, 'Value', 1)
    elseif parameter.lsmopt.pw_clk==2
        set(handles.radiobutton_estimation_leastSquares_pwlAndRateClock, 'Value', 1)
    elseif parameter.lsmopt.pw_clk==3
        set(handles.radiobutton_estimation_leastSquares_pwlRateAndQuClock, 'Value', 1)
    end
    if parameter.lsmopt.ascii_snx==1
        snxState='on';
    else
        snxState='off';
    end
end

% set baseline-dependent clock offsets
if isfield(parameter.lsmopt, 'est_bdco')
    if parameter.lsmopt.est_bdco==1
        set(handles.checkbox_estimation_leastSquares_clocksBasDepOffset, 'Value', 1);
        set(handles.radiobutton_estimation_leastSquares_basdepClockoff_OPT, 'Enable', 'on')
        set(handles.radiobutton_estimation_leastSquares_basdepClockoff_automatic, 'Enable', 'on')
        set(handles.edit_estimation_leastSquares_clockBasDepO_minNobs, 'Enable', 'on');
        set(handles.text403, 'Enable', 'on');
        set(handles.text404, 'Enable', 'on');
        if parameter.lsmopt.bdco_fromOPT
            set(handles.radiobutton_estimation_leastSquares_basdepClockoff_OPT, 'Value',1);
        else
            set(handles.radiobutton_estimation_leastSquares_basdepClockoff_OPT, 'Value',0);
            set(handles.edit_estimation_leastSquares_clockBasDepO_minNobs, 'String', num2str(parameter.lsmopt.bdco_minobs));
        end
    else
        set(handles.checkbox_estimation_leastSquares_clocksBasDepOffset, 'Value', 0);
    end
else
    set(handles.checkbox_estimation_leastSquares_clocksBasDepOffset, 'Value', 0);
    set(handles.radiobutton_estimation_leastSquares_basdepClockoff_OPT, 'Enable', 'off')
    set(handles.radiobutton_estimation_leastSquares_basdepClockoff_automatic, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_clockBasDepO_minNobs, 'Enable', 'off');
    set(handles.text403, 'Enable', 'off');
    set(handles.text404, 'Enable', 'off');
end


set(handles.text_run_sinex_clockParam, 'Enable', snxState)
set(handles.radiobutton_run_sinex_clockParam_incl, 'Enable', 'off') % is never written to SINEX
set(handles.radiobutton_run_sinex_clockParam_excl, 'Enable', snxState)

set(handles.checkbox_run_sinex_sources, 'Enable', snxState)
set(handles.checkbox_run_sinex_changeAnalystsName, 'Enable', snxState)
set(handles.text300, 'Enable', snxState)
set(handles.checkbox_run_sinex_addSuffix, 'Enable', snxState)

% estimation of ZWDs
% input interval of zwd
set(handles.edit_estimation_leastSquares_tropo_zwdInterval, 'String', parameter.lsmopt.int_zwd)
set(handles.text309, 'String', sprintf('after %s minutes', num2str(parameter.lsmopt.int_zwd)))
% rel constr (0|1)
set(handles.checkbox_estimation_leastSquares_tropo_zwdRelConstr, 'Value', parameter.lsmopt.constr_zwd)
% set constraint
set(handles.edit_estimation_leastSquares_tropo_zwdRelConstr, 'String', num2str(parameter.lsmopt.coef_zwd))
if parameter.lsmopt.pw_zwd==0 % if zwd should be estimtated
    set(handles.checkbox_estimation_leastSquares_tropo_zwd, 'Value', 0)
    
    % set disable
    set(handles.text_estimation_leastSquares_tropo_zwdInt, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_tropo_zwdInterval, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_tropo_zwdRelConstr, 'Enable', 'off')
    set(handles.text_estimation_leastSquares_tropo_zwdRelConstr, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_tropo_zwdRelConstr, 'Enable', 'off')
else
    % set enable
    set(handles.text_estimation_leastSquares_tropo_zwdInt, 'Enable', 'on')
    set(handles.edit_estimation_leastSquares_tropo_zwdInterval, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_tropo_zwdRelConstr, 'Enable', 'on')
    if get(handles.checkbox_estimation_leastSquares_tropo_zwdRelConstr, 'Value')==1
        set(handles.text_estimation_leastSquares_tropo_zwdRelConstr, 'Enable', 'on')
        set(handles.edit_estimation_leastSquares_tropo_zwdRelConstr, 'Enable', 'on')
    end
end

% north gradients
set(handles.checkbox_estimation_leastSquares_tropo_ngr, 'Value', parameter.lsmopt.pw_ngr)
% NGR interval
set(handles.edit_estimation_leastSquares_tropo_ngrInterval, 'String', num2str(parameter.lsmopt.int_ngr))
set(handles.text310, 'String', sprintf('after %s minutes', num2str(parameter.lsmopt.int_ngr)))
% NGR rel. constraint (1|0)
set(handles.checkbox_estimation_leastSquares_tropo_ngrRelConstr, 'Value', parameter.lsmopt.constr_rel_ngr)
% NGR rel constraint
set(handles.edit_estimation_leastSquares_tropo_ngrRelConstr, 'String', num2str(parameter.lsmopt.coef_rel_ngr))
% NGR abs. constraint (1|0)
set(handles.checkbox_estimation_leastSquares_tropo_ngrAbsConstr, 'Value', parameter.lsmopt.constr_abs_ngr)
% NGR abs. constraint
set(handles.edit_estimation_leastSquares_tropo_ngrAbsConstr, 'String', num2str(parameter.lsmopt.coef_abs_ngr))

% set (1|0) for estimation of north gradients
if parameter.lsmopt.pw_ngr==0
    % set disable
    set(handles.text_estimation_leastSquares_tropo_ngrInt, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_tropo_ngrInterval, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_tropo_ngrRelConstr, 'Enable', 'off')
    set(handles.text_estimation_leastSquares_tropo_ngrRelConstr, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_tropo_ngrRelConstr, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_tropo_ngrAbsConstr, 'Enable', 'off')
    set(handles.text_estimation_leastSquares_tropo_ngrAbsConstr, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_tropo_ngrAbsConstr, 'Enable', 'off')
else
    % set enable
    set(handles.text_estimation_leastSquares_tropo_ngrInt, 'Enable', 'on')
    set(handles.edit_estimation_leastSquares_tropo_ngrInterval, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_tropo_ngrRelConstr, 'Enable', 'on')
    if get(handles.checkbox_estimation_leastSquares_tropo_ngrRelConstr, 'Value')==1
        set(handles.text_estimation_leastSquares_tropo_ngrRelConstr, 'Enable', 'on')
        set(handles.edit_estimation_leastSquares_tropo_ngrRelConstr, 'Enable', 'on')
    else
        set(handles.text_estimation_leastSquares_tropo_ngrRelConstr, 'Enable', 'off')
        set(handles.edit_estimation_leastSquares_tropo_ngrRelConstr, 'Enable', 'off')
    end
    set(handles.checkbox_estimation_leastSquares_tropo_ngrAbsConstr, 'Enable', 'on')
    if get(handles.checkbox_estimation_leastSquares_tropo_ngrAbsConstr, 'Value')==1
        set(handles.text_estimation_leastSquares_tropo_ngrAbsConstr, 'Enable', 'on')
        set(handles.edit_estimation_leastSquares_tropo_ngrAbsConstr, 'Enable', 'on')
    else
        set(handles.text_estimation_leastSquares_tropo_ngrAbsConstr, 'Enable', 'off')
        set(handles.edit_estimation_leastSquares_tropo_ngrAbsConstr, 'Enable', 'off')
    end
end

% east gradients
set(handles.checkbox_estimation_leastSquares_tropo_egr, 'Value', parameter.lsmopt.pw_egr)
% EGR interval
set(handles.edit_estimation_leastSquares_tropo_egrInterval, 'String', num2str(parameter.lsmopt.int_egr))
set(handles.text313, 'String', sprintf('after %s minutes', num2str(parameter.lsmopt.int_egr)))
% EGR rel. constraint (1|0)
set(handles.checkbox_estimation_leastSquares_tropo_egrRelConstr, 'Value', parameter.lsmopt.constr_rel_egr)
% EGR rel constraint
set(handles.edit_estimation_leastSquares_tropo_egrRelConstr, 'String', num2str(parameter.lsmopt.coef_rel_egr))
% EGR abs. constraint (1|0)
set(handles.checkbox_estimation_leastSquares_tropo_egrAbsConstr, 'Value', parameter.lsmopt.constr_abs_egr)
% EGR abs. constraint
set(handles.edit_estimation_leastSquares_tropo_egrAbsConstr, 'String', num2str(parameter.lsmopt.coef_abs_egr))

% set (1|0) for estimation of east gradients
if parameter.lsmopt.pw_egr==0
    % set disable
    set(handles.text_estimation_leastSquares_tropo_egrInt, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_tropo_egrInterval, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_tropo_egrRelConstr, 'Enable', 'off')
    set(handles.text_estimation_leastSquares_tropo_egrRelConstr, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_tropo_egrRelConstr, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_tropo_egrAbsConstr, 'Enable', 'off')
    set(handles.text_estimation_leastSquares_tropo_egrAbsConstr, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_tropo_egrAbsConstr, 'Enable', 'off')
else
    % set enable
    set(handles.text_estimation_leastSquares_tropo_egrInt, 'Enable', 'on')
    set(handles.edit_estimation_leastSquares_tropo_egrInterval, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_tropo_egrRelConstr, 'Enable', 'on')
    if get(handles.checkbox_estimation_leastSquares_tropo_egrRelConstr, 'Value')==1
        set(handles.text_estimation_leastSquares_tropo_egrRelConstr, 'Enable', 'on')
        set(handles.edit_estimation_leastSquares_tropo_egrRelConstr, 'Enable', 'on')
    else
        set(handles.text_estimation_leastSquares_tropo_egrRelConstr, 'Enable', 'off')
        set(handles.edit_estimation_leastSquares_tropo_egrRelConstr, 'Enable', 'off')
    end
    set(handles.checkbox_estimation_leastSquares_tropo_egrAbsConstr, 'Enable', 'on')
    if get(handles.checkbox_estimation_leastSquares_tropo_egrAbsConstr, 'Value')==1
        set(handles.text_estimation_leastSquares_tropo_egrAbsConstr, 'Enable', 'on')
        set(handles.edit_estimation_leastSquares_tropo_egrAbsConstr, 'Enable', 'on')
    else
        set(handles.text_estimation_leastSquares_tropo_egrAbsConstr, 'Enable', 'off')
        set(handles.edit_estimation_leastSquares_tropo_egrAbsConstr, 'Enable', 'off')
    end
end

% estimate station coordinates
set(handles.checkbox_estimation_leastSquares_coordinates_estimate, 'Value', parameter.lsmopt.stc)
if parameter.lsmopt.stc==1
    set(handles.checkbox_estimation_leastSquares_coordinates_NNT, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_coordinates_NNR, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_coordinates_NNS, 'Enable', 'on')
    if parameter.lsmopt.ascii_snx==1
        snxState='on';
    else
        snxState='off';
    end
else
    set(handles.checkbox_estimation_leastSquares_coordinates_NNT, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_coordinates_NNR, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_coordinates_NNS, 'Enable', 'off')
    snxState='off';
end
set(handles.checkbox_estimation_leastSquares_coordinates_NNT, 'Value', parameter.lsmopt.nnt_stc)
set(handles.checkbox_estimation_leastSquares_coordinates_NNR, 'Value', parameter.lsmopt.nnr_stc)
set(handles.checkbox_estimation_leastSquares_coordinates_NNS, 'Value', parameter.lsmopt.sca_stc)
if ~isfield(parameter.lsmopt, 'datum')
    parameter.lsmopt.datum = 'trf';
    msgbox('Loaded parameter file (from older VieVS version) does not contain the "datum (trf/all)" option!', 'Warning', 'warn');
end
if strcmp(parameter.lsmopt.datum,'trf')
    set(handles.radiobutton_estimation_leastSquares_coordinates_datum_trf, 'Value', 1)
elseif strcmp(parameter.lsmopt.datum,'all')
    set(handles.radiobutton_estimation_leastSquares_coordinates_datum_all, 'Value', 1)
else
    msgbox('"parameter.lsmopt.datum" does not contain a valid option (trf/all)!', 'Warning', 'warn');
end
    
set(handles.radiobutton_run_sinex_stationCoords_incl, 'Enable', snxState)
set(handles.text_run_sinex_stationCoords, 'Enable', snxState)

% sources
set(handles.checkbox_estimation_leastSquares_sources_est, 'Value', parameter.lsmopt.pw_sou)
set(handles.checkbox_estimation_leastSquares_sources_constr, 'Value', parameter.lsmopt.constr_sou)
set(handles.edit_estimation_leastSquares_sources_constr, 'String', num2str(parameter.lsmopt.sour_coef_rade))
set(handles.edit_estimation_leastSquares_sources_interval, 'String', num2str(parameter.lsmopt.sour_int_rade))
set(handles.text321, 'String', sprintf('after %s minutes', num2str(parameter.lsmopt.sour_int_rade)))

set(handles.checkbox_estimation_leastSquares_sources_NNR, 'Value', parameter.lsmopt.est_sourceNNR)
if parameter.lsmopt.est_sourceNNR
    set(handles.checkbox_estimation_leastSquares_sources_ICRF2_def, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_sources_abs_constr, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_sources_obs_per_source, 'Enable', 'on')
else
    set(handles.checkbox_estimation_leastSquares_sources_ICRF2_def, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_sources_abs_constr, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_sources_obs_per_source, 'Enable', 'off')
end

if isfield(parameter.lsmopt, 'est_sourceNNR_defining')
    set(handles.checkbox_estimation_leastSquares_sources_ICRF2_def, 'Value', parameter.lsmopt.est_sourceNNR_defining)
else
    set(handles.checkbox_estimation_leastSquares_sources_ICRF2_def, 'Value', 0)
end

if isfield(parameter.lsmopt, 'sourceAbsConstrNNR')
    set(handles.edit_estimation_leastSquares_sources_abs_constr, 'String', num2str(parameter.lsmopt.sourceAbsConstrNNR))
else
    set(handles.edit_estimation_leastSquares_sources_abs_constr, 'String', '1')
end

if isfield(parameter.lsmopt, 'UseSourceAbsConstrNNR')
    if parameter.lsmopt.UseSourceAbsConstrNNR
        set(handles.checkbox_estimation_leastSquares_sources_abs_constr, 'Value', 1)
        set(handles.edit_estimation_leastSquares_sources_abs_constr, 'Enable', 'on')
    else
        set(handles.checkbox_estimation_leastSquares_sources_abs_constr, 'Value', 0)
        set(handles.edit_estimation_leastSquares_sources_abs_constr, 'Enable', 'off')
    end
else
    set(handles.checkbox_estimation_leastSquares_sources_abs_constr, 'Value', 0)
    set(handles.edit_estimation_leastSquares_sources_abs_constr, 'Enable', 'off')
end


if isfield(parameter.lsmopt, 'min_num_obs_per_est_source')
    set(handles.edit_estimation_leastSquares_sources_obs_per_source, 'String', num2str(parameter.lsmopt.min_num_obs_per_est_source))
else
    set(handles.edit_estimation_leastSquares_sources_obs_per_source, 'String', '5')
end

if isfield(parameter.lsmopt, 'use_min_num_obs_per_est_source')
    if parameter.lsmopt.use_min_num_obs_per_est_source
        set(handles.checkbox_estimation_leastSquares_sources_obs_per_source, 'Value', 1)
        set(handles.edit_estimation_leastSquares_sources_obs_per_source, 'Enable', 'on')
    else
        set(handles.checkbox_estimation_leastSquares_sources_obs_per_source, 'Value', 0)
        set(handles.edit_estimation_leastSquares_sources_obs_per_source, 'Enable', 'off')
    end

else
    set(handles.checkbox_estimation_leastSquares_sources_obs_per_source, 'Value', 1)
    set(handles.edit_estimation_leastSquares_sources_obs_per_source, 'Enable', 'on')
end



if parameter.lsmopt.pw_sou==0 || parameter.lsmopt.est_sourceNNR==1
    % set disable
    set(handles.edit_estimation_leastSquares_sources_interval, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_sources_constr, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_sources_constr, 'Enable', 'off')
else
    set(handles.edit_estimation_leastSquares_sources_interval, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_sources_constr, 'Enable', 'on')
    if get(handles.checkbox_estimation_leastSquares_sources_constr, 'value')==1
        set(handles.edit_estimation_leastSquares_sources_constr, 'Enable', 'on')
    else
        set(handles.edit_estimation_leastSquares_sources_constr, 'Enable', 'off')
    end
end

% prepare parameters for global solution
set(handles.checkbox_run_prepareGlobParam, 'Value', parameter.lsmopt.global_solve)
set(handles.checkbox_run_globalPram_statVel, 'Value', parameter.lsmopt.est_vel)
set(handles.checkbox_run_globalPram_source, 'Value', parameter.lsmopt.est_source)
set(handles.edit_run_globalPram_statVelEpoch, 'String', num2str(parameter.lsmopt.refvel))
if parameter.lsmopt.global_solve==0
    set(handles.checkbox_run_globalPram_source, 'Enable', 'off')
    set(handles.text_run_runOptions_attentionSources, 'Enable', 'off')
    set(handles.checkbox_run_globalPram_statVel, 'Enable', 'off')
    set(handles.edit_run_globalPram_statVelEpoch, 'Enable', 'off')
    set(handles.text_run_runOptions_refEpochInYrs, 'Enable', 'off')
    set(handles.checkbox_run_globalPram_axisOffset, 'Enable', 'off')
    set(handles.checkbox_run_globalPram_stseaspos, 'Enable','off')
    set(handles.checkbox_run_globalPram_tidERPvar,'Enable','off')
    set(handles.checkbox_run_globalPram_hlpole, 'Enable','off')
    set(handles.checkbox_run_globalPram_APLrg, 'Enable','off')
else
    set(handles.checkbox_run_globalPram_source, 'Enable', 'on')
    set(handles.text_run_runOptions_attentionSources, 'Enable', 'on')
    set(handles.checkbox_run_globalPram_statVel, 'Enable', 'on')
    if parameter.lsmopt.est_vel==1
        set(handles.edit_run_globalPram_statVelEpoch, 'Enable', 'on')
        set(handles.text_run_runOptions_refEpochInYrs, 'Enable', 'on')
    else
        set(handles.edit_run_globalPram_statVelEpoch, 'Enable', 'off')
        set(handles.text_run_runOptions_refEpochInYrs, 'Enable', 'off')
    end
    set(handles.checkbox_run_globalPram_axisOffset, 'Enable', 'on')
    set(handles.checkbox_run_globalPram_stseaspos, 'Enable','on')
    set(handles.checkbox_run_globalPram_tidERPvar,'Enable','on')
    set(handles.checkbox_run_globalPram_hlpole, 'Enable','on')
    set(handles.checkbox_run_globalPram_APLrg, 'Enable','on')
end
set(handles.checkbox_run_globalPram_axisOffset, 'Value', parameter.lsmopt.est_AO)
set(handles.checkbox_run_globalPram_stseaspos, 'Value', parameter.lsmopt.est_stsespos)
set(handles.checkbox_run_globalPram_tidERPvar, 'Value', parameter.lsmopt.est_erpsp)
set(handles.checkbox_run_globalPram_hlpole, 'Value', parameter.lsmopt.est_hpole)
set(handles.checkbox_run_globalPram_APLrg, 'Value', parameter.lsmopt.est_rg)


% SINEX output
set(handles.checkbox_run_sinex_write, 'Value', parameter.lsmopt.ascii_snx)
set(handles.radiobutton_run_sinex_clockParam_incl, 'Value', parameter.lsmopt.outsnx.clk)
set(handles.radiobutton_run_sinex_zwd_incl, 'Value', parameter.lsmopt.outsnx.zwd)
set(handles.radiobutton_run_sinex_tropoParam_incl, 'Value', parameter.lsmopt.outsnx.tgr)
set(handles.radiobutton_run_sinex_sources_incl, 'Value', 1)
set(handles.radiobutton_run_sinex_stationCoords_incl, 'Value', parameter.lsmopt.outsnx.xyz)
set(handles.radiobutton_run_sinex_eop_incl, 'Value', parameter.lsmopt.outsnx.eop)
set(handles.checkbox_run_sinex_sources, 'Value', parameter.lsmopt.addSnxSource)
set(handles.edit_run_sinex_firstname, 'String', parameter.lsmopt.outsnx.firstname)
set(handles.edit_run_sinex_lastname, 'String', parameter.lsmopt.outsnx.lastname)
set(handles.edit_run_sinex_email, 'String', parameter.lsmopt.outsnx.email)
set(handles.edit_run_sinex_firstname,'Enable', 'on');
set(handles.edit_run_sinex_lastname,'Enable', 'on');
set(handles.edit_run_sinex_email,'Enable', 'on');
if isfield(parameter.lsmopt.outsnx, 'suffix')
    set(handles.edit_run_sinex_suffix, 'String', parameter.lsmopt.outsnx.suffix)
    set(handles.edit_run_sinex_suffix,'Enable', 'on');
end

try
    set(handles.checkbox_run_sinex_changeAnalystsName, 'Value', parameter.lsmopt.outsnx.changeAnalystsName)
catch
    set(handles.checkbox_run_sinex_changeAnalystsName, 'Value', 0)
    set(handles.edit_run_sinex_firstname,'Enable', 'off');
    set(handles.edit_run_sinex_lastname,'Enable', 'off');
    set(handles.edit_run_sinex_email,'Enable', 'off');
    warning('Information about analyst name for sinex not in parameter file');
end
try 
    set(handles.checkbox_run_sinex_addSuffix, 'Value', parameter.lsmopt.outsnx.addSuffix)
catch
    set(handles.checkbox_run_sinex_addSuffix, 'Value', 0)
    set(handles.edit_run_sinex_suffix,'Enable', 'off');
    warning('Information about adding a suffix to the sinex file not in parameter file');
end

if parameter.lsmopt.addSnxSource && parameter.lsmopt.ascii_snx
    set(handles.text_run_sinex_sources, 'Enable', 'on')
    set(handles.radiobutton_run_sinex_sources_incl, 'Enable', 'on')
end

set(handles.radiobutton_run_sinex_clockParam_excl, 'Value', ~parameter.lsmopt.outsnx.clk)
set(handles.radiobutton_run_sinex_zwd_excl, 'Value', ~parameter.lsmopt.outsnx.zwd)
set(handles.radiobutton_run_sinex_tropoParam_excl, 'Value', ~parameter.lsmopt.outsnx.tgr)
set(handles.radiobutton_run_sinex_eop_excl, 'Value', ~parameter.lsmopt.outsnx.eop)

% sinex zwd enabling
if (parameter.lsmopt.pw_zwd == 0) || ...
        parameter.lsmopt.ascii_snx == 0
    newState='off';
else
    newState='on';
end
set(handles.text_run_sinex_zwd, 'Enable', newState)
set(handles.radiobutton_run_sinex_zwd_incl, 'Enable', newState)
set(handles.radiobutton_run_sinex_zwd_excl, 'Enable', newState)

% sinex gradients enabling
if (parameter.lsmopt.pw_ngr+parameter.lsmopt.pw_egr == 0) || ...
        parameter.lsmopt.ascii_snx == 0
    newState='off';
else
    newState='on';
end
set(handles.text_run_sinex_tropoGradients, 'Enable', newState)
set(handles.radiobutton_run_sinex_tropoParam_incl, 'Enable', newState)
set(handles.radiobutton_run_sinex_tropoParam_excl, 'Enable', newState)

% EOP
% xpol
set(handles.checkbox_estimation_leastSquares_eop_xpEst, 'Value', parameter.lsmopt.xpol.model)
set(handles.edit_estimation_leastSquares_eop_interval_xp, 'String', num2str(parameter.lsmopt.xpol.int))
set(handles.text316, 'String', sprintf('after %s minutes', num2str(parameter.lsmopt.xpol.int)))
set(handles.checkbox_estimation_leastSquares_eop_constr_xp, 'value', parameter.lsmopt.xpol.constrain)
set(handles.edit_estimation_leastSquares_constr_xp, 'String', num2str(parameter.lsmopt.xpol.coef))
if parameter.lsmopt.xpol.model==0
    set(handles.edit_estimation_leastSquares_eop_interval_xp, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_eop_constr_xp, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_constr_xp, 'Enable', 'off')
else
    set(handles.edit_estimation_leastSquares_eop_interval_xp, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_eop_constr_xp, 'Enable', 'on')
    if parameter.lsmopt.xpol.constrain==1
        set(handles.edit_estimation_leastSquares_constr_xp, 'Enable', 'on')
    else
        set(handles.edit_estimation_leastSquares_constr_xp, 'Enable', 'off')
    end
end
% ypol
set(handles.checkbox_estimation_leastSquares_eop_ypEst, 'Value', parameter.lsmopt.ypol.model)
set(handles.edit_estimation_leastSquares_eop_interval_yp, 'String', num2str(parameter.lsmopt.ypol.int))
set(handles.text317, 'String', sprintf('after %s minutes', num2str(parameter.lsmopt.ypol.int)))
set(handles.checkbox_estimation_leastSquares_eop_constr_yp, 'value', parameter.lsmopt.ypol.constrain)
set(handles.edit_estimation_leastSquares_constr_yp, 'String', num2str(parameter.lsmopt.ypol.coef))
if parameter.lsmopt.ypol.model==0
    set(handles.edit_estimation_leastSquares_eop_interval_yp, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_eop_constr_yp, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_constr_yp, 'Enable', 'off')
else
    set(handles.edit_estimation_leastSquares_eop_interval_yp, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_eop_constr_yp, 'Enable', 'on')
    if parameter.lsmopt.ypol.constrain==1
        set(handles.edit_estimation_leastSquares_constr_yp, 'Enable', 'on')
    else
        set(handles.edit_estimation_leastSquares_constr_yp, 'Enable', 'off')
    end
end
% dut1
set(handles.checkbox_estimation_leastSquares_eop_dut1Est, 'Value', parameter.lsmopt.dut1.model)
set(handles.edit_estimation_leastSquares_eop_interval_dut1, 'String', num2str(parameter.lsmopt.dut1.int))
set(handles.text318, 'String', sprintf('after %s minutes', num2str(parameter.lsmopt.dut1.int)))
set(handles.checkbox_estimation_leastSquares_eop_constr_dut1, 'value', parameter.lsmopt.dut1.constrain)
set(handles.edit_estimation_leastSquares_constr_dut1, 'String', num2str(parameter.lsmopt.dut1.coef))
if parameter.lsmopt.dut1.model==0
    set(handles.edit_estimation_leastSquares_eop_interval_dut1, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_eop_constr_dut1, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_constr_dut1, 'Enable', 'off')
else
    set(handles.edit_estimation_leastSquares_eop_interval_dut1, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_eop_constr_dut1, 'Enable', 'on')
    if parameter.lsmopt.dut1.constrain==1
        set(handles.edit_estimation_leastSquares_constr_dut1, 'Enable', 'on')
    else
        set(handles.edit_estimation_leastSquares_constr_dut1, 'Enable', 'off')
    end
end
% nutdx
set(handles.checkbox_estimation_leastSquares_eop_nutdxEst, 'Value', parameter.lsmopt.nutdx.model)
set(handles.edit_estimation_leastSquares_eop_interval_nutdx, 'String', num2str(parameter.lsmopt.nutdx.int))
set(handles.text319, 'String', sprintf('after %s minutes', num2str(parameter.lsmopt.nutdx.int)))
set(handles.checkbox_estimation_leastSquares_eop_constr_nutdx, 'value', parameter.lsmopt.nutdx.constrain)
set(handles.edit_estimation_leastSquares_constr_nutdx, 'String', num2str(parameter.lsmopt.nutdx.coef))
if parameter.lsmopt.nutdx.model==0
    set(handles.edit_estimation_leastSquares_eop_interval_nutdx, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_eop_constr_nutdx, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_constr_nutdx, 'Enable', 'off')
else
    set(handles.edit_estimation_leastSquares_eop_interval_nutdx, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_eop_constr_nutdx, 'Enable', 'on')
    if parameter.lsmopt.nutdx.constrain==1
        set(handles.edit_estimation_leastSquares_constr_nutdx, 'Enable', 'on')
    else
        set(handles.edit_estimation_leastSquares_constr_nutdx, 'Enable', 'off')
    end
end
% nutdy
set(handles.checkbox_estimation_leastSquares_eop_nutdyEst, 'Value', parameter.lsmopt.nutdy.model)
set(handles.edit_estimation_leastSquares_eop_interval_nutdy, 'String', num2str(parameter.lsmopt.nutdy.int))
set(handles.text320, 'String', sprintf('after %s minutes', num2str(parameter.lsmopt.nutdy.int)))
set(handles.checkbox_estimation_leastSquares_eop_constr_nutdy, 'value', parameter.lsmopt.nutdy.constrain)
set(handles.edit_estimation_leastSquares_constr_nutdy, 'String', num2str(parameter.lsmopt.nutdy.coef))
if parameter.lsmopt.nutdy.model==0
    set(handles.edit_estimation_leastSquares_eop_interval_nutdy, 'Enable', 'off')
    set(handles.checkbox_estimation_leastSquares_eop_constr_nutdy, 'Enable', 'off')
    set(handles.edit_estimation_leastSquares_constr_nutdy, 'Enable', 'off')
else
    set(handles.edit_estimation_leastSquares_eop_interval_nutdy, 'Enable', 'on')
    set(handles.checkbox_estimation_leastSquares_eop_constr_nutdy, 'Enable', 'on')
    if parameter.lsmopt.nutdy.constrain==1
        set(handles.edit_estimation_leastSquares_constr_nutdy, 'Enable', 'on')
    else
        set(handles.edit_estimation_leastSquares_constr_nutdy, 'Enable', 'off')
    end
end

% if neither eop is estimtated -> disable options in sinex output
if (parameter.lsmopt.xpol.model+...
        parameter.lsmopt.ypol.model+...
        parameter.lsmopt.dut1.model+...
        parameter.lsmopt.nutdx.model+...
        parameter.lsmopt.nutdy.model)==0 || ...
        parameter.lsmopt.ascii_snx==0
    set(handles.radiobutton_run_sinex_eop_incl, 'Enable', 'off');
    set(handles.radiobutton_run_sinex_eop_excl, 'Enable', 'off');
    set(handles.text_run_sinex_eop, 'Enable', 'off')
else
    set(handles.radiobutton_run_sinex_eop_incl, 'Enable', 'on');
    set(handles.radiobutton_run_sinex_eop_excl, 'Enable', 'on');
    set(handles.text_run_sinex_eop, 'Enable', 'on')
end

% define if parameters should be estimated
set(handles.checkbox_run_estParameters, 'Value', parameter.lsmopt.est_singleses);

% set session-wise parameterization option
set(handles.checkbox_run_allowStationwise, 'Value', parameter.lsmopt.control_gui_vie_lsm)

% set use of clock breaks and manually finding breaks
set(handles.checkbox_estimation_leastSquares_clocks_useClockBreaks, 'Value', parameter.lsmopt.treat_breaks)
