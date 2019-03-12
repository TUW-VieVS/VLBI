% creates a protocol about chosen models and parameterization
% the information is taken from struct SESSION_parameter.mat in LEVEL0,
% which has been created by VIE_SETUP

% created for VieVS by Hana Spicakova
% 19 Jan 2012

% updates:
% 27 March 2012 by Hana Spicakova: hydrology loading added
% 15 Nov 2012 by Hana Kr�sn�: distinction between pressure and temperature
% values
% 19 Dec 2012 by Hana Kr�sn�: changes according to changes of constrain's representation in
% the GUI 2.1
% 18 March 2013 by Hana Kr�sn�: change according to new supersource file
% 25 Jan 2017 by Daniel Landskron: output changed due to troposphere model change
% 13 Sep 2017 by Daniel Landskron: 'tropSource' shifted into 'vie_init' 
% 11 May 2018 by Daniel Landskron: remaining variable 'ext' removed
% 29 Nov 2018 by Daniel Landskron: structure of observation restrictions standardized
%%

function create_input_protocol(parameter, outfile)
    
    if nargin<2
        outfile='input_protocol.txt';
    end
    fid=fopen(outfile, 'w');
    clear input_protocol

    curDate=clock;          % current date and time
    
    fprintf(fid, 'created [yr mon day h min sec]:');
    fprintf(fid, '  %4.0f %2.0f %2.0f  %2.0f %2.0f %2.0f \n', curDate);
    
    fprintf(fid, '\n\nselected input files and models\n');
    fprintf(fid, '___________________________________________\n\n');
    
    fprintf(fid, '\n OPT directory: %s ', parameter.vie_init.diropt );
    if iscell(parameter.vie_init.dirout)
        fprintf(fid, '\n OUTLIER directory: %s ', parameter.vie_init.dirout{1} );        
    else
        fprintf(fid, '\n OUTLIER directory: %s ', parameter.vie_init.dirout );
    end
    fprintf(fid, '\n remove outliers: %1.0f ', parameter.vie_init.rm_outlier );
    
    fprintf(fid, '\n\n a priori TRF file: %s ', parameter.vie_init.trf{1});
    fprintf(fid, '\n TRF field (for main station file): %s ', parameter.vie_init.trf{2});
    %fprintf(fid, '\n\n a priori CRF: %s ', parameter.vie_init.crf(8:end));
    fprintf(fid, '\n\n a priori CRF: %s ', parameter.vie_init.crf{1}); % change according to supersource file
    fprintf(fid, '\n CRF field (for main station file): %s ', parameter.vie_init.crf{2}); % change according to supersource file
    
    cutoff = parameter.obs_restrictions.cut_off_elev/pi*180; % rad --> deg
    fprintf(fid, '\n\n cut-off elevation angle: %2.0f degree', cutoff);
    fprintf(fid, '\n quality code limit: %1.0f', parameter.obs_restrictions.Qlim);

    fprintf(fid, '\n info about temperature: %s', parameter.vie_init.tp);
    fprintf(fid, '\n info about ionosphere: %s', parameter.vie_init.iono);
    
    
    
    fprintf(fid, '\n\n ephemerides: %s', parameter.vie_mod.eph);
    
    % EOP
    if parameter.vie_mod.linear==1; eopIntrpl=' linear';
    else eopIntrpl=' lagrange'; end
    
    if parameter.vie_mod.dXdY==0; eop_dXdY = ' - but nutation offsets dX, dY were set to zero!';
    else eop_dXdY = ''; end
    
    fprintf(fid, '\n\n EOP file: %s %s', parameter.vie_mod.EOPfile, eop_dXdY);
    fprintf(fid, '\n ocean tides: %s', parameter.vie_mod.eopoc);
    fprintf(fid, '\n libration in xpol, ypol: %1.0f', parameter.vie_mod.lib_pm);
    fprintf(fid, '\n libration in UT1: %1.0f', parameter.vie_mod.lib_ut);
    fprintf(fid, '\n EOP interpolation: %s', eopIntrpl);
    fprintf(fid, '\n tidal UT variations: %1.0f', parameter.vie_mod.tidalUT);
    
    fprintf(fid, '\n\n precession/nutation model: %s', parameter.vie_mod.nutmod);
    

    % stations corrections
    fprintf(fid, '\n\n solid Earth tides: %1.0f ', parameter.vie_mod.cts );
    fprintf(fid, '\n tidal ocean loading: %1.0f ', parameter.vie_mod.cto );
    fprintf(fid, ', %s ', parameter.vie_mod.ocm );
    fprintf(fid, '\n tidal atmosphere loading: %1.0f ', parameter.vie_mod.cta );
    fprintf(fid, ', %s ', parameter.vie_mod.ctam );
    fprintf(fid, '\n non-tidal atmosphere loading: %1.0f ', parameter.vie_mod.cnta );
    fprintf(fid, ', %s ', parameter.vie_mod.cntam );
    fprintf(fid, '\n pole tide: %1.0f ', parameter.vie_mod.ctp );
    if parameter.vie_mod.ctp;  fprintf(fid, ', mean pole model: %s ', parameter.vie_mod.ctpm ); end
    fprintf(fid, '\n ocean pole tide: %1.0f ', parameter.vie_mod.ctop );
    if parameter.vie_mod.ctop;  fprintf(fid, ', mean pole model: %s ', parameter.vie_mod.ctpm ); end
    fprintf(fid, '\n hydrology loading: %1.0f ', parameter.vie_mod.chl );
    fprintf(fid, ', %s ', parameter.vie_mod.chlm );

    
    fprintf(fid, '\n\n antenna thermal deformation: %1.0f', parameter.vie_mod.therm);
    fprintf(fid, '\n hydr. a priori troposphere gradients: %s ', parameter.vie_mod.apgm_h);
    fprintf(fid, '\n wet a priori troposphere gradients: %s ', parameter.vie_mod.apgm_w);
    fprintf(fid, '\n hydrostatic mapping functions: %s ', parameter.vie_mod.mfh);
    fprintf(fid, '\n wet mapping functions: %s ', parameter.vie_mod.mfw);

%% parametrization

    fprintf(fid, '\n\n\n\nselected parametrization\n');
    fprintf(fid, '___________________________________________\n');

    
   % OUTLIERS 
    outl=0; coef=0; outltest='';
    if parameter.lsmopt.simple_outlier || parameter.lsmopt.basic_outlier
        outl = 1;
        coef = parameter.lsmopt.par_outlier;
        if parameter.lsmopt.simple_outlier
             outltest = 'simple';
        else if parameter.lsmopt.basic_outlier
                outltest = 'basic';
            end
        end
    end
   
    fprintf(fid, '\n\n outlier detection: %1.0f', outl);
    if outl == 1
        fprintf(fid, '\n (outliers will be stored but not removed in this run)');
        fprintf(fid, '\n %s outlier test with a coefficient of %2.1f', outltest, coef);
    end
    
    % CLOCK
    if parameter.lsmopt.pw_clk==3
        clk=' piece-wise linear offsets + one rate + one squared rate per station';
    elseif parameter.lsmopt.pw_clk==2
        clk=' piece-wise linear offsets + one rate per station';
    elseif parameter.lsmopt.pw_clk==1
        clk=' piece-wise linear offsets per station';
    else %parameter.lsmopt.pw_clk==0
        clk=' not estimated'; 
    end
    clk_con=[''];
    if parameter.lsmopt.pw_clk~=0
        if parameter.lsmopt.constr_clk
            clk_con=[num2str(parameter.lsmopt.int_clk) ' min interval and constraints ' num2str(parameter.lsmopt.coef_clk) ' cm after ' num2str(parameter.lsmopt.int_clk) ' min' ];
        else
            clk_con=['no constraints'];
        end
    end
    
    fprintf(fid, '\n\n clocks');
    fprintf(fid, '\n parametrization: %s', clk);
    fprintf(fid, '\n constraints: %s', clk_con);

    
    % ZWD
    if parameter.lsmopt.pw_zwd==1
        zwd=[num2str(parameter.lsmopt.int_zwd) ' min offsets'];
        zwd_con=['relative constr. ' num2str(parameter.lsmopt.coef_zwd) ' cm after '  num2str(parameter.lsmopt.int_zwd) ' min' ];
    else zwd='not estimated'; zwd_con=[''];
    end
    
    
    fprintf(fid, '\n\n zenith wet delay');
    fprintf(fid, '\n parametrization: %s', zwd);
    fprintf(fid, '\n constraints: %s', zwd_con);
    
    % NGR
    if parameter.lsmopt.pw_ngr==1 
        if parameter.lsmopt.constr_rel_ngr==1 && parameter.lsmopt.constr_abs_ngr==0
            ngr=[' ' num2str(parameter.lsmopt.int_ngr) ' min offsets'];
            ngr_con=['relative constr. ' num2str(parameter.lsmopt.coef_rel_ngr) ' cm after ' num2str(parameter.lsmopt.int_ngr) ' min'];
        elseif parameter.lsmopt.constr_rel_ngr==1 && parameter.lsmopt.constr_abs_ngr==1
            ngr=[' ' num2str(parameter.lsmopt.int_ngr) ' min offsets'];
            ngr_con=['relative constr. ' num2str(parameter.lsmopt.coef_rel_ngr) 'cm after ' num2str(parameter.lsmopt.int_ngr) ' min and absolute constr. ' num2str(parameter.lsmopt.coef_abs_ngr) ' mm'];
        elseif parameter.lsmopt.constr_rel_ngr==0 && parameter.lsmopt.constr_abs_ngr==1
            ngr=[' ' num2str(parameter.lsmopt.int_ngr) ' min offsets'];
            ngr_con=['absolute constr. ' num2str(parameter.lsmopt.coef_abs_ngr) ' mm'];
        elseif parameter.lsmopt.constr_rel_ngr==0;
            ngr=[' ' num2str(parameter.lsmopt.int_ngr) ' min offsets'];
            ngr_con=['no constraints'];
        end
    else ngr='not estimated'; ngr_con=[''];
    end
    
    fprintf(fid, '\n\n north tropospheric gradients');
    fprintf(fid, '\n parametrization: %s', ngr);
    fprintf(fid, '\n constraints: %s', ngr_con);
    
    % EGR
    if parameter.lsmopt.pw_egr==1 
        if parameter.lsmopt.constr_rel_egr==1 && parameter.lsmopt.constr_abs_egr==0
            egr=[' ' num2str(parameter.lsmopt.int_egr) ' min offsets'];
            egr_con=['relative constr. ' num2str(parameter.lsmopt.coef_rel_egr) ' cm after ' num2str(parameter.lsmopt.int_egr) ' min'];
        elseif parameter.lsmopt.constr_rel_egr==1 && parameter.lsmopt.constr_abs_egr==1
            egr=[' ' num2str(parameter.lsmopt.int_egr) ' min offsets'];
            egr_con=['relative constr. ' num2str(parameter.lsmopt.coef_rel_egr) ' cm after ' num2str(parameter.lsmopt.int_egr) ' min and absolute constr. ' num2str(parameter.lsmopt.coef_abs_egr) ' mm'];
        elseif parameter.lsmopt.constr_rel_egr==0 && parameter.lsmopt.constr_abs_egr==1
            egr=[' ' num2str(parameter.lsmopt.int_egr) ' min offsets'];
            egr_con=['absolute constr. ' num2str(parameter.lsmopt.coef_abs_egr) ' mm'];
        elseif parameter.lsmopt.constr_rel_egr==0;
            egr=[' ' num2str(parameter.lsmopt.int_egr) ' min offsets'];
            egr_con=['no constraints'];
        end
    else egr='not estimated'; egr_con=[''];
    end
    
    fprintf(fid, '\n\n east tropospheric gradients');
    fprintf(fid, '\n parametrization: %s', egr);
    fprintf(fid, '\n constraints: %s', egr_con);
   
    
    % STATIONS
    if parameter.lsmopt.stc==1
        stat=['one offset per session'];
        if isfield(parameter.lsmopt,'nnt_stc'); stat_con = ['NNT']; end
        if isfield(parameter.lsmopt,'nnr_stc'); stat_con = [stat_con ' NNR']; end
        if isfield(parameter.lsmopt,'nns_stc'); stat_con = [stat_con ' NNS']; end
    else stat=['not estimated'];
        stat_con=[''];
    end
    
    fprintf(fid, '\n\n stations');
    fprintf(fid, '\n parametrization: %s', [stat ',' parameter.lsmopt.datum]);
    fprintf(fid, '\n conditions: %s', stat_con);
    
    % EOP
    eop(1) = parameter.lsmopt.xpol.model;
    eop(2) = parameter.lsmopt.ypol.model;
    eop(3) = parameter.lsmopt.dut1.model;
    eop(4) = parameter.lsmopt.nutdx.model;
    eop(5) = parameter.lsmopt.nutdy.model;
    
    if eop(1)==1
        eop_int(1) = parameter.lsmopt.xpol.int;
        if parameter.lsmopt.xpol.constrain
            eop_con(1) = parameter.lsmopt.xpol.coef;
        else eop_con(1) = 0;
        end
    else eop_int(1)=0; eop_con(1)=0;
    end
    
    
    if eop(2)==1;
        eop_int(2) = parameter.lsmopt.ypol.int;
        if parameter.lsmopt.ypol.constrain
            eop_con(2) = parameter.lsmopt.ypol.coef;
        else eop_con(2) = 0;
        end
    else eop_int(2)=0; eop_con(2)=0;
    end
    
    if eop(3)==1;
        eop_int(3) = parameter.lsmopt.dut1.int;
        if parameter.lsmopt.dut1.constrain
            eop_con(3) = parameter.lsmopt.dut1.coef;
        else eop_con(3) = 0;
        end
    else eop_int(3)=0; eop_con(3)=0;
    end
    
    if eop(4)==1;
        eop_int(4) = parameter.lsmopt.nutdx.int;
        if parameter.lsmopt.nutdx.constrain
                eop_con(4) = parameter.lsmopt.nutdx.coef;
        else eop_con(4) = 0;
        end
    else eop_int(4)=0; eop_con(4)=0;
    end
    
    if eop(5)==1;
        eop_int(5) = parameter.lsmopt.nutdy.int;
        if parameter.lsmopt.nutdy.constrain
                eop_con(5) = parameter.lsmopt.nutdy.coef;
        else eop_con(5) = 0;
        end
    else eop_int(5)=0; eop_con(5)=0;
    end
    
    fprintf(fid, '\n\n EOP:            x pole - y pole -   dut1 -     dX -     dY');
    fprintf(fid, '\n estimated:      %6.0f - %6.0f - %6.0f - %6.0f - %6.0f', eop);
    fprintf(fid, '\n intervals:      %6.0f - %6.0f - %6.0f - %6.0f - %6.0f min', eop_int);
    fprintf(fid, '\n constraints:    %3.7g - %3.7g - %3.7g - %3.7g - %3.7g mas/day', eop_con);
    
    
    % SOURCES
    if parameter.lsmopt.pw_sou==1
        sou1=[' as default only one source is estimated - you can change it in the session-wise GUIs'];
        sou = [' ' num2str(parameter.lsmopt.sour_int_rade) ' min offsets'];
        if parameter.lsmopt.constr_sou ==1
            sou_con = ['relative constr. ' num2str(parameter.lsmopt.sour_coef_rade) ' mas/day'];
        else
            sou_con=[''];
        end
    else
        sou=[' not estimated'];
        sou1=['']; sou_con=[''];
    end
    
    fprintf(fid, '\n\n sources');
    fprintf(fid, '\n parametrization: %s', sou);
    fprintf(fid, '\n constraints: %s', sou_con);
    fprintf(fid, '\n%s', sou1);
    
    fprintf(fid, '\n\n\n');

    fclose(fid);


