% ************************************************************************
%   Description:
%   This function saves all estimates into a matlab structure 'globsol'
%
%   Input:										
%      globsol              a structure with estimates and all relavant
%                           information about the global adjustment
%      ses_time             order of sessions according to time
%      paths                paths to directories
%
%   Output:                
%      'glob_results_*.txt' TXT file with estimates in
%                           VieVS/OUT/GLOB/_ESTIMATES
%
%
%
%   External calls: 	
%      -              					    											
%       
%   Coded for VieVS: 
%   03 Aug 2010 by Hana Spicakova
%
%   Revision: 
%   05 Oct 2010 by Hana Spicakova: writes datum definition for sources
%   12 Oct 2010 by Hana Spicakova: writes datum definition for stations
%   20 Jun 2012 by Hana Krásná: 
%               added: Love and Shida numbers, FCN per from SET, velocity of sources, acceleration of SSB, seasonal
%               variations in station positions (annual, semi-annual), pole
%               tide Love and Shida numbers
%   25 Sep 2012 by Hana Krásná: relativistic parameter gamma
%   24 Oct 2012 by Hana Krásná: change units of RA
%   04 Oct 2013 by Hana Krasna: antenna axis offset added
%   06 Dec 2013 by Hana Krasna: APL regression coefficient added
%   ?? Apr 2017 by David Mayer: Minor changes (??)

%%


function createTXT(globsol,ses_time,badses,badses_mo,paths)

if exist([paths.path_out '_ESTIMATES/' paths.out])~=7
   mkdir([paths.path_out '_ESTIMATES/' paths.out])
end

fid=fopen([paths.path_out '_ESTIMATES/' paths.out '/glob_results_' paths.L2 '.txt'],'wt');

fprintf(fid,' %% Parameterisation \n');
fprintf(fid,'DIRIN:  %s \n',paths.L2);
fprintf(fid,'DIROUT: %s \n',paths.out);
fprintf(fid,'TRF datum: %s \n',paths.datumsta);
fprintf(fid,'TRF reduced ant: %s \n',paths.antred);
fprintf(fid,'TRF constant vel: %s \n',paths.velconst);
fprintf(fid,'TRF velocity ties: %s \n',paths.velties);
fprintf(fid,'TRF discontinuity: %s \n',paths.discont);
fprintf(fid,'CRF datum: %s \n',paths.datumsou);
fprintf(fid,'CRF fixed sou: %s \n',paths.fixedsou);
fprintf(fid,'CRF reduced sou: %s \n',paths.soured);

fprintf(fid,' %% Sessions in the solution \n');
for i=1:size(ses_time,2)
    fprintf(fid,'%s   \n',ses_time{i});
end

fprintf(fid,'\n %% Number of sessions in the global adjustment\n');
lse=size(globsol.sessions,2);
fprintf(fid,'%1.0f\n',lse);

fprintf(fid,'\n %% Maximal RMS of the sessions in the solution\n');
fprintf(fid,'\n %% applicable only if you run single-session solution in VIE_LSM\n');
maxRMS=globsol.maxRMS;
fprintf(fid,'%3.2f\n',maxRMS);

fprintf(fid,'\n %% Sessions which were excluded from the solution (RMS > %3.2f) \n',maxRMS);
for i=1:size(badses,2)
	fprintf(fid,'%s      %10.2f\n',badses{i},badses_mo(i)');
end


fprintf(fid,'\n\n %% Corrections to station coordinates and formal errors in [cm]\n');
fprintf(fid,'\n %% station         dx         dy         dz           mx         my         mz       epoch    start     end  \n');

if globsol.antenna.id_dxyz==1
    [refname_coord,ind]=sortrows(globsol.antenna.refname_coord);
    a = char(refname_coord);
    dxyz     = globsol.antenna.dxyz(ind,:);
    sigma_xyz= globsol.antenna.sigma_xyz(ind,:);
    eb       = globsol.antenna.epoch_break(ind,:);
    for i=1:size(globsol.antenna.dxyz,1)
        fprintf(fid,'%c%c%c%c%c%c%c%c    %9.3f  %9.3f  %9.3f    %9.3f  %9.3f  %9.3f      %6.0f   %6.0f  %6.0f\n' ,a(i,:), dxyz(i,:), sigma_xyz(i,:), eb(i,:));
    end
    clear a
end
fprintf(fid,'\n\n %% Corrections to station velocities and formal errors in [cm/y]\n');
fprintf(fid,'\n %% station        dvx        dvy        dvz          mvx        mvy        mvz       epoch    start     end  \n');

if globsol.antenna.id_dvxyz==1
    [refname_vel,ind]=sortrows(globsol.antenna.refname_coord);
    a = char(refname_vel);
    dvxyz     = globsol.antenna.dvxyz(ind,:);
    sigma_vxyz= globsol.antenna.sigma_vxyz(ind,:);
    eb        = globsol.antenna.epoch_break(ind,:);
    for i=1:size(globsol.antenna.dvxyz,1)
        fprintf(fid,'%c%c%c%c%c%c%c%c    %9.3f  %9.3f  %9.3f    %9.3f  %9.3f  %9.3f      %6.0f   %6.0f  %6.0f\n' ,a(i,:), dvxyz(i,:), sigma_vxyz(i,:), eb(i,:));
    end
    clear a
end

if size(globsol.antenna.dxyz,1)>0 || size(globsol.antenna.dvxyz,1)>0
    fprintf(fid,'\n %% datum definition: %s\n',globsol.antenna.datumdef);
    datum=sortrows(globsol.antenna.datum);
    for i=1:size(globsol.antenna.datum,1)
        fprintf(fid,'%c%c%c%c%c%c%c%c \n', datum(i,:));
    end
end

clear ind

fprintf(fid,'\n\n____________________________________________________________________\n');
fprintf(fid,'\n\n\n %% Corrections to source coordinates and formal errors: RA in [ms], De in [mas]\n');
fprintf(fid,'\n %% source            RA          De            mRA         mDe \n');

if globsol.source.id==1
    [refname_s,ind]=sortrows(globsol.source.refname.IVS);
    a = char(refname_s);
    drade = globsol.source.drade(ind,:);
    sigma_rade = globsol.source.sigma_rade(ind,:);
    for i=1:size(globsol.source.drade,1)
        fprintf(fid,'%c%c%c%c%c%c%c%c     %10.4f  %10.4f     %10.4f  %10.4f \n' ,a(i,:), drade(i,:), sigma_rade(i,:));
    end

    
    fprintf(fid,'\n\n____________________________________________________________________\n');
    fprintf(fid,'\n\n\n %% Source velocities and formal errors: RA in [ms/y] and De in [mas/y]\n');
    fprintf(fid,'\n %% source            vRA          vDe            mvRA         mvDe \n');
    if globsol.source.id_vel==1
        dvrade = globsol.source.dvrade(ind,:);
        sigma_vrade = globsol.source.sigma_vrade(ind,:);
        for i=1:size(globsol.source.drade,1)
            fprintf(fid,'%c%c%c%c%c%c%c%c     %10.4f  %10.4f     %10.4f  %10.4f \n' ,a(i,:), dvrade(i,:), sigma_vrade(i,:));
        end
    end
    clear a
    
    if size(globsol.source.datum,1)>0 
        fprintf(fid,'\n %% sources included in NNR (+ dz) condition\n');
        nnr=sortrows(globsol.source.datum);
        for i=1:size(globsol.source.datum,1)
            fprintf(fid,'%c%c%c%c%c%c%c%c \n', nnr(i,:));
        end
    end

    if size(globsol.source.fixed_apr,1)>0 
        fprintf(fid,'\n %% sources fixed to apriori coordinates\n');
        fso=sortrows(globsol.source.fixed_apr);
        for i=1:size(globsol.source.fixed_apr,1)
            fprintf(fid,'%c%c%c%c%c%c%c%c \n', fso(i,:));
        end
    end
    
    fprintf(fid,'\n %% datum definition: %s\n', globsol.source.datumdef);
    
end

fprintf(fid,'\n\n____________________________________________________________________\n');
fprintf(fid,'\n\n\n %% EOPs were:  0 = fixed,  1 = estimated,  2 = reduced\n');


fprintf(fid,'\n\n %% x pole: mjd of estimates, corrections in [mas], formal errors in [mas]\n');
fprintf(fid,' %1d\n',globsol.xpol.id);
fprintf(fid,' %7.4f       %7.4f    %7.4f\n',[globsol.xpol.mjd globsol.xpol.val globsol.xpol.sigma]');

fprintf(fid,'\n %% y pole: mjd of estimates, corrections in [mas], formal errors in [mas]\n');
fprintf(fid,' %1d\n',globsol.ypol.id);
fprintf(fid,' %7.4f       %7.4f    %7.4f\n',[globsol.ypol.mjd globsol.ypol.val globsol.ypol.sigma]');

fprintf(fid,'\n %% ut1: mjd of estimates, corrections in [ms], formal errors in [ms]\n');
fprintf(fid,' %1d\n',globsol.dut1.id);
fprintf(fid,' %7.4f       %7.4f    %7.4f\n',[globsol.dut1.mjd globsol.dut1.val globsol.dut1.sigma]');

fprintf(fid,'\n %% nutation X: mjd of estimates, corrections in [mas], formal errors in [mas]\n');
fprintf(fid,' %1d\n',globsol.dX.id);
fprintf(fid,' %7.4f       %7.4f    %7.4f\n',[globsol.dX.mjd globsol.dX.val globsol.dX.sigma]');

fprintf(fid,'\n %% nutation Y: mjd of estimates, corrections in [mas], formal errors in [mas]\n');
fprintf(fid,' %1d\n',globsol.dY.id);
fprintf(fid,' %7.4f       %7.4f    %7.4f\n',[globsol.dY.mjd globsol.dY.val globsol.dY.sigma]');

fprintf(fid,'\n\n____________________________________________________________________\n');
fprintf(fid,'\n\n\n %% Corrections to axis offset in [cm], formal error in [cm]\n');
fprintf(fid,'\n %% station     dAO            mAO\n');

    [aoname,ind]=sortrows(globsol.antenna.aoname);
    a = char(aoname);
    ao     = globsol.antenna.ao(ind,:);
    sigma_ao= globsol.antenna.sigma_ao(ind,:);
    for i=1:length(globsol.antenna.ao)
        fprintf(fid,'%c%c%c%c%c%c%c%c    %9.3f     %9.3f   \n' ,a(i,:), ao(i,:), sigma_ao(i,:));
    end
    clear a

fprintf(fid,'\n\n %% Amplitudes of seasonal variations in station position \n');
fprintf(fid,' %% period in [solar days], amplitudes in REN [cm], formal errors in [cm], phases in REN [deg], formal errors in [deg]\n');

if globsol.stseaspos.id==1;
    Pwav = globsol.stseaspos.periods_solardays;
    nwav = length(Pwav);
    for i=1:length(globsol.stseaspos.results)
        fprintf(fid,' %c%c%c%c%c%c%c%c \n',globsol.stseaspos.results(i).aname);
        if ~isempty(globsol.stseaspos.results(i).Ar)
            res(:,1) = globsol.stseaspos.results(i).Ar;
            sigma_res(:,1) = globsol.stseaspos.results(i).sigma_Ar;
            res_ph(:,1) = globsol.stseaspos.results(i).phaser;
            sigma_res_ph(:,1) = globsol.stseaspos.results(i).sigma_phaser;
        else
            res(1:nwav,1) = 0;
            sigma_res(1:nwav,1) = 0;
            res_ph(1:nwav,1) = 0;
            sigma_res_ph(1:nwav,1) = 0;

        end
        if ~isempty(globsol.stseaspos.results(i).Ae)
            res(:,2) = globsol.stseaspos.results(i).Ae;
            sigma_res(:,2) = globsol.stseaspos.results(i).sigma_Ae;
            res_ph(:,2) = globsol.stseaspos.results(i).phasee;
            sigma_res_ph(:,2) = globsol.stseaspos.results(i).sigma_phasee;

        else
            res(1:nwav,2) = 0;
            sigma_res(1:nwav,2) = 0;
            res_ph(1:nwav,2) = 0;
            sigma_res_ph(1:nwav,2) = 0;

        end
        if ~isempty(globsol.stseaspos.results(i).An)
            res(:,3) = globsol.stseaspos.results(i).An;
            sigma_res(:,3) = globsol.stseaspos.results(i).sigma_An;
            res_ph(:,3) = globsol.stseaspos.results(i).phasen;
            sigma_res_ph(:,3) = globsol.stseaspos.results(i).sigma_phasen;

        else
            res(1:nwav,3) = 0;
            sigma_res(1:nwav,3) = 0;
            res_ph(1:nwav,3) = 0;
            sigma_res_ph(1:nwav,3) = 0;

        end
        
        Pwavres=[Pwav res sigma_res res_ph sigma_res_ph];
        
        fprintf(fid,' %15.4f     %7.3f  %7.3f  %7.3f     %7.3f  %7.3f  %7.3f     %7.1f  %7.1f  %7.1f     %7.1f  %7.1f  %7.1f\n', Pwavres');
        
        clear Pwavres  res sigma_res
        
    end
end


fprintf(fid,'\n\n\n');
fprintf(fid,'\n\n %% Pole tide Love number,  correction in [-], formal error in [-]\n');
fprintf(fid,' %7.4f    %7.4f \n',[globsol.hpole.val globsol.hpole.sigma]');
fprintf(fid,'\n\n %% Pole tide Shida number,  correction in [-], formal error in [-]\n');
fprintf(fid,' %7.4f    %7.4f \n',[globsol.lpole.val globsol.lpole.sigma]');


fprintf(fid,'\n\n\n');
fprintf(fid,'\n\n %% APL regression coefficients\n');
fprintf(fid,' %% station    correction   formal error   total value   reference pressure \n');
fprintf(fid,' %% station    [cm/hPa]     [cm/hPa]       [cm/hPa]      [hPa]\n');
if globsol.antenna.idrg==1
    [rgname,ind]=sortrows(globsol.antenna.rgname);
    a = char(rgname);
    rg     = globsol.antenna.rg(ind,:);
    sigma_rg= globsol.antenna.sigma_rg(ind,:);
    apr_rg = globsol.antenna.apriori_rg(ind,:);
    total = rg+apr_rg(:,2);
    refpres = apr_rg(:,1);
    for i=1:length(globsol.antenna.rg)
        fprintf(fid,'%c%c%c%c%c%c%c%c       %6.3f      %6.3f        %6.3f         %6.3f\n' ,a(i,:), rg(i), sigma_rg(i), total(i), refpres(i));
    end
    clear a
end



%% INTERNAL
if globsol.dlove.id==1
    fprintf(fid,'\n\n____________________________________________________________________\n');
    fprintf(fid,'\n\n %% Constant Love number: h2 = 0.6078,  correction in [-], formal error in [-]\n');
    fprintf(fid,' %2.0f    %7.4f    %7.4f \n',[globsol.dlove.nr_h2 globsol.dlove.val_h2 globsol.dlove.sigma_h2]');

    fprintf(fid,'\n %% Diurnal real Love numbers: corrections in [-], formal errors in [-]\n');
    fprintf(fid,' %2.0f    %7.4f    %7.4f \n',[globsol.dlove.nr_diurnal globsol.dlove.val_diurnal globsol.dlove.sigma_diurnal]');

    fprintf(fid,'\n %% Diurnal imaginary Love numbers: corrections in [-], formal errors in [-]\n');
    fprintf(fid,' %2.0f    %7.4f    %7.4f \n',[globsol.dlove.nr_diurnal_imag globsol.dlove.val_diurnal_imag globsol.dlove.sigma_diurnal_imag]');


    fprintf(fid,'\n %% Long-period real Love numbers: corrections in [-], formal errors in [-]\n');
    fprintf(fid,' %2.0f    %7.4f    %7.4f \n',[globsol.dlove.nr_long globsol.dlove.val_long globsol.dlove.sigma_long]');

    fprintf(fid,'\n %% Long-period imaginary Love numbers: corrections in [-], formal errors in [-]\n');
    fprintf(fid,' %2.0f    %7.4f    %7.4f \n',[globsol.dlove.nr_long_imag globsol.dlove.val_long_imag globsol.dlove.sigma_long_imag]');
end

if globsol.dshida.id==1
    fprintf(fid,'\n\n %% Constant Shida number: l2 = 0.0847,  correction in [-], formal error in [-]\n');
    fprintf(fid,' %2.0f    %7.4f    %7.4f \n',[globsol.dshida.nr_l2 globsol.dshida.val_l2 globsol.dshida.sigma_l2]');

    fprintf(fid,'\n %% Diurnal real Shida numbers: corrections in [-], formal errors in [-]\n');
    fprintf(fid,' %2.0f    %7.4f    %7.4f\n',[globsol.dshida.nr_diurnal globsol.dshida.val_diurnal globsol.dshida.sigma_diurnal]');

    fprintf(fid,'\n %% Diurnal imaginary Shida numbers: corrections in [-], formal errors in [-]\n');
    fprintf(fid,' %2.0f    %7.4f    %7.4f\n',[globsol.dshida.nr_diurnal_imag globsol.dshida.val_diurnal_imag globsol.dshida.sigma_diurnal_imag]');


    fprintf(fid,'\n %% Long-period real Shida numbers: corrections in [-], formal errors in [-]\n');
    fprintf(fid,' %2.0f    %7.4f    %7.4f \n',[globsol.dshida.nr_long globsol.dshida.val_long globsol.dshida.sigma_long]');

    fprintf(fid,'\n %% Long-period imaginary Shida numbers: corrections in [-], formal errors in [-]\n');
    fprintf(fid,' %2.0f    %7.4f    %7.4f \n',[globsol.dshida.nr_long_imag globsol.dshida.val_long_imag globsol.dshida.sigma_long_imag]');
end

if globsol.dNDFWset.id==1
    fprintf(fid,'\n\n %% Nearly Diurnal Free Wobble frequency (from solid Earth tides), correction in [cycle/sid.day], formal error in [cycle/sid.day]\n');
    fprintf(fid,' %10.9f   %10.9f\n',globsol.dNDFWset.val,globsol.dNDFWset.sigma);

    fprintf(fid,'\n %% corresponding Free Core Nutation period (from solid Earth tides), a priori [sid.day], new in [sid.day], formal error in [sid.day]\n');
    fprintf(fid,' %5.2f    %5.2f    %5.2f\n',globsol.dNDFWset.FCN_apr, globsol.dNDFWset.FCN_new ,globsol.dNDFWset.dFCN_sigma);
end


if globsol.daccSSB.id==1
    fprintf(fid,'\n\n %% Acceleration of the Solar-System Barycenter (x,y,z), correction in [cm/sec^2], formal error in [cm/sec^2]\n');
    if globsol.daccSSB.id==1
        fprintf(fid,' %15.11f   %15.11f\n',globsol.daccSSB.val(1),globsol.daccSSB.sigma(1));
        fprintf(fid,' %15.11f   %15.11f\n',globsol.daccSSB.val(2),globsol.daccSSB.sigma(2));
        fprintf(fid,' %15.11f   %15.11f\n',globsol.daccSSB.val(3),globsol.daccSSB.sigma(3));
    end
end

if globsol.gamma.id==1
    fprintf(fid,'\n\n\n');
    fprintf(fid,'\n\n %% Relativistic parameter gamma ,  correction in [-], formal error in [-]\n');
    fprintf(fid,' %10.6f    %10.6f \n',[globsol.gamma.val globsol.gamma.sigma]');
end

fclose(fid);

