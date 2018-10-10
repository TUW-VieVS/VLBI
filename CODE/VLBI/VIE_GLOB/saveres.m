% ************************************************************************
%   Description:
%   This function saves all estimates into a matlab structure 'globsol'
%
%   Input:										
%
%   Output:                
%      globsol              a structure with estimates and all relavant
%                           information about the global adjustment
%
%   External calls: 	
%      globind              					    											
%       
%   Coded for VieVS: 
%   03 Aug 2010 by Hana Spicakova
%
%   Revision: 
%   05 Oct 2010 by Hana Spicakova: globsol.source.datumdef
%   12 Oct 2010 by Hana Spicakova: globsol.antenna.datumdef
%   25 Jun 2012 by Hana Krásná:
%               added: Love and Shida numbers, FCN per from SET, velocity of sources, acceleration of SSB, seasonal
%               variations in station positions (annual, semi-annual), pole
%               tide Love and Shida numbers
%   25 Sep 2012 by Hana Krásná: relativistic parameter gamma 
%   24 Oct 2012 by Hana Krásná: change units of RA in globsol.source.rade
%   and globsol.source.sigma_rade from mas to ms
%   04 Oct 2013 by Hana Krasna: antenna axis offset added
%   06 Dec 2013 by Hana Krasna: APL regression coefficient added
%   25 Aug 2015 by Sigrid Boehm: Tidal ERP coefficients added
% ************************************************************************


function globsol = saveres(parGS,parGS_hl,globsol,refname,refnamec,datumant,scale,...
                           qrefname,datumsouc,fixedsou,sou_dz,RA_all,De_all,eb,...
                           trffile,crffile,aostname,stseasname,P_stseaspos,rgstname,apriori_aplrg)

[g] = globind(parGS);

globsol.antenna.antennas=[''];
globsol.antenna.apriori_cat=[''];
globsol.antenna.datum=[''];

globsol.antenna.id_dxyz=parGS(g.g_coord(1)).id;
globsol.antenna.refname_coord=[''];
globsol.antenna.dxyz=[];
globsol.antenna.sigma_xyz=[];
globsol.antenna.epoch_break=[];
globsol.antenna.dxyzcol=[];

globsol.antenna.id_dvxyz=parGS(g.g_vel(1)).id;
globsol.antenna.dvxyz=[];
globsol.antenna.sigma_vxyz=[];
globsol.antenna.dvxyzcol=[];
globsol.antenna.datumdef=[''];

globsol.antenna.idao=parGS(g.g_ao).id;
globsol.antenna.aoname = [''];
globsol.antenna.ao = [];
globsol.antenna.sigma_ao = [];
globsol.antenna.col_ao = [];

globsol.antenna.idrg=parGS(g.g_aplrg).id;
globsol.antenna.rgname = [''];
globsol.antenna.rg = []; 
globsol.antenna.sigma_rg = []; 
globsol.antenna.col_rg = []; 
globsol.antenna.apriori_rg=[];

globsol.source.id=parGS(g.g_srade(1)).id;
globsol.source.refname=[''];
globsol.source.apriori_cat=[''];
globsol.source.apriori_rade=[];
globsol.source.datum = [''];
globsol.source.fixed_apr = [];
globsol.source.drade=[];
globsol.source.sigma_rade=[];
globsol.source.dradecol=[];
globsol.source.datumdef=[''];

globsol.source.id_vel = parGS(g.g_svrade(1)).id; 
globsol.source.dvrade=[]; 
globsol.source.sigma_vrade=[];
globsol.source.dvradecol=[]; 

globsol.xpol.id=parGS(g.g_eop(1)).id;
globsol.xpol.val=[];
globsol.xpol.sigma=[];
globsol.xpol.mjd=[];
globsol.xpol.col=[];
globsol.ypol.id=parGS(g.g_eop(2)).id;
globsol.ypol.val=[];
globsol.ypol.sigma=[];
globsol.ypol.mjd=[];
globsol.ypol.col=[];
globsol.dut1.id=parGS(g.g_eop(3)).id;
globsol.dut1.val=[];
globsol.dut1.sigma=[];
globsol.dut1.mjd=[];
globsol.dut1.col=[];
globsol.dX.id=parGS(g.g_eop(4)).id;
globsol.dX.val=[];
globsol.dX.sigma=[];
globsol.dX.mjd=[];
globsol.dX.col=[];
globsol.dY.id=parGS(g.g_eop(5)).id;
globsol.dY.val=[];
globsol.dY.sigma=[];
globsol.dY.mjd=[];
globsol.dY.col=[];


globsol.stseaspos.periods_solardays = [];
if parGS(g.g_stseaspos_Ar(1)).id==1 || parGS(g.g_stseaspos_Ae(1)).id==1 || parGS(g.g_stseaspos_An(1)).id==1
    globsol.stseaspos.id=1;
else
    globsol.stseaspos.id=0;
end

for i=1:size(stseasname)
    globsol.stseaspos.results(i).aname=[''];
    globsol.stseaspos.results(i).Acr=[];
    globsol.stseaspos.results(i).Asr=[];
    globsol.stseaspos.results(i).Ar=[];
    globsol.stseaspos.results(i).Ace=[];
    globsol.stseaspos.results(i).Ase=[];
    globsol.stseaspos.results(i).Ae=[];
    globsol.stseaspos.results(i).Acn=[];
    globsol.stseaspos.results(i).Asn=[];
    globsol.stseaspos.results(i).An=[];
    
    globsol.stseaspos.results(i).sigma_Acr=[];
    globsol.stseaspos.results(i).sigma_Asr=[];
    globsol.stseaspos.results(i).sigma_Ar=[];
    globsol.stseaspos.results(i).sigma_Ace=[];
    globsol.stseaspos.results(i).sigma_Ase=[];
    globsol.stseaspos.results(i).sigma_Ae=[];
    globsol.stseaspos.results(i).sigma_Acn=[];
    globsol.stseaspos.results(i).sigma_Asn=[];
    globsol.stseaspos.results(i).sigma_An=[];
end


globsol.hpole.id=parGS(g.g_hpole).id;
globsol.hpole.val=[];
globsol.hpole.sigma=[];
globsol.hpole.col=[];

globsol.lpole.id=parGS(g.g_lpole).id;
globsol.lpole.val=[];
globsol.lpole.sigma=[];
globsol.lpole.col=[];




globsol.dlove.id=parGS(g.g_love).id;
globsol.dlove.val_h2=[];
globsol.dlove.sigma_h2=[];
globsol.dlove.nr_h2=[];
globsol.dlove.col_h2=[];
globsol.dlove.val_diurnal=[];
globsol.dlove.sigma_diurnal=[];
globsol.dlove.nr_diurnal=[];
globsol.dlove.col_diurnal=[];
globsol.dlove.val_diurnal_imag=[];
globsol.dlove.sigma_diurnal_imag=[];
globsol.dlove.nr_diurnal_imag=[];
globsol.dlove.col_diurnal_imag=[];
globsol.dlove.val_long=[];
globsol.dlove.sigma_long=[];
globsol.dlove.nr_long=[];
globsol.dlove.col_long=[];
globsol.dlove.val_long_imag=[];
globsol.dlove.sigma_long_imag=[];
globsol.dlove.nr_long_imag=[];
globsol.dlove.col_long_imag=[];

globsol.dshida.id=parGS(g.g_shida).id;
globsol.dshida.val_l2=[];
globsol.dshida.sigma_l2=[];
globsol.dshida.nr_l2=[];
globsol.dshida.col_l2=[];
globsol.dshida.val_diurnal=[];
globsol.dshida.sigma_diurnal=[];
globsol.dshida.nr_diurnal=[];
globsol.dshida.col_diurnal=[];
globsol.dshida.val_diurnal_imag=[];
globsol.dshida.sigma_diurnal_imag=[];
globsol.dshida.nr_diurnal_imag=[];
globsol.dshida.col_diurnal_imag=[];
globsol.dshida.val_long=[];
globsol.dshida.sigma_long=[];
globsol.dshida.nr_long=[];
globsol.dshida.col_long=[];
globsol.dshida.val_long_imag=[];
globsol.dshida.sigma_long_imag=[];
globsol.dshida.nr_long_imag=[];
globsol.dshida.col_long_imag=[];

globsol.dNDFWset.id=parGS(g.g_FCNset).id;
globsol.dNDFWset.val=[];
globsol.dNDFWset.sigma=[];
globsol.dNDFWset.col=[];

globsol.dNDFWset.FCN_apr=[];
globsol.dNDFWset.FCN_new=[];
globsol.dNDFWset.dFCN_sigma=[];

globsol.daccSSB.id=parGS(g.g_accSSB).id;
globsol.daccSSB.val=[];
globsol.daccSSB.sigma=[];

globsol.gamma.id=parGS(g.g_gamma).id;
globsol.gamma.val=[];
globsol.gamma.sigma=[];
globsol.gamma.col=[];


%%

x = globsol.x;
varpar = globsol.varpar;

globsol.antenna.antennas=refname;
if parGS(g.g_coord(1)).id==1 || parGS(g.g_vel(1)).id==1
    globsol.antenna.apriori_cat=trffile;
    globsol.antenna.datum = datumant;
    globsol.antenna.epoch_break=eb;
 
    globsol.antenna.refname_coord=refnamec;
    globsol.antenna.dxyz=[x(parGS(g.g_coord(1)).newcol),x(parGS(g.g_coord(2)).newcol),x(parGS(g.g_coord(3)).newcol)];
    globsol.antenna.sigma_xyz=[varpar(parGS(g.g_coord(1)).newcol),varpar(parGS(g.g_coord(2)).newcol),varpar(parGS(g.g_coord(3)).newcol)];
    globsol.antenna.dxyzcol=[parGS(g.g_coord(1)).newcol;parGS(g.g_coord(2)).newcol;parGS(g.g_coord(3)).newcol]';
    
    if parGS(g.g_vel(1)).id==1 % velocities
        
        ivx=parGS(g.g_vel(1)).newcol;
        ivy=parGS(g.g_vel(2)).newcol;
        ivz=parGS(g.g_vel(3)).newcol;
        
        globsol.antenna.dvxyz=[x(ivx),x(ivy),x(ivz)]*0.01; %[cm*century/day]=[cm/100y] --> [cm/y]
        globsol.antenna.sigma_vxyz=[varpar(ivx),varpar(ivy),varpar(ivz)]*0.01; %[cm*century/day]=[cm/100y] --> [cm/y]
        globsol.antenna.dvxyzcol=[ivx;ivy;ivz]';
    end
    
    % which datum definition was applied?
    if scale=='0'
        globsol.antenna.datumdef = 'NNT/NNR';
    elseif scale=='1'
        globsol.antenna.datumdef = 'NNT/NNR/NNS';
    end
end

if parGS(g.g_srade(1)).id==1
    globsol.source.refname = qrefname;
    globsol.source.apriori_cat=crffile;
    globsol.source.apriori_rade = [RA_all', De_all'];
    globsol.source.datum = datumsouc;
    globsol.source.fixed_apr = fixedsou; % names of sources fixed to apriori coordinates
    globsol.source.drade = [x(parGS(g.g_srade(1)).newcol)./15,x(parGS(g.g_srade(2)).newcol)]; % [ms, mas]
    globsol.source.sigma_rade = [varpar(parGS(g.g_srade(1)).newcol)./15,varpar(parGS(g.g_srade(2)).newcol)]; % [ms, mas]
    globsol.source.dradecol = [parGS(g.g_srade(1)).newcol;parGS(g.g_srade(2)).newcol]';
    % which datum definition was applied?
    if (sou_dz=='1' & isempty(fixedsou))
        globsol.source.datumdef = 'NNR+dz';
    elseif (sou_dz=='0'& isempty(fixedsou))
        globsol.source.datumdef = 'NNR';
    elseif (isempty(sou_dz) & ~isempty(fixedsou))
        globsol.source.datumdef = 'some fixed sources';
    elseif (sou_dz=='1' & ~isempty(fixedsou))
        globsol.source.datumdef = 'NNR+dz and some fixed sources';
    elseif (sou_dz=='0' & ~isempty(fixedsou))
        globsol.source.datumdef = 'NNR and some fixed sources';
    end
end

if parGS(g.g_svrade(1)).id==1
    globsol.source.dvrade = [x(parGS(g.g_svrade(1)).newcol)./15,x(parGS(g.g_svrade(2)).newcol)]*0.01; %[mas*century/day]=[mas/100y] --> [ms/y],[mas/y];
    globsol.source.sigma_vrade = [varpar(parGS(g.g_svrade(1)).newcol)./15,varpar(parGS(g.g_svrade(2)).newcol)]*0.01; %[mas*century/day]=[mas/100y] --> [ms/y],[mas/y];
    globsol.source.dvradecol = [parGS(g.g_svrade(1)).newcol;parGS(g.g_svrade(2)).newcol]';
end

if parGS(g.g_eop(1)).id==1
   globsol.xpol.val=x(parGS(g.g_eop(1)).newcol);
   globsol.xpol.sigma=varpar(parGS(g.g_eop(1)).newcol);
   globsol.xpol.mjd=parGS(g.g_eop(1)).mjdstep';
   globsol.xpol.col=parGS(g.g_eop(1)).newcol;
end


if parGS(g.g_eop(2)).id==1
   globsol.ypol.val=x(parGS(g.g_eop(2)).newcol);
   globsol.ypol.sigma=varpar(parGS(g.g_eop(2)).newcol);
   globsol.ypol.mjd=parGS(g.g_eop(2)).mjdstep';
   globsol.ypol.col=parGS(g.g_eop(2)).newcol;
end

if parGS(g.g_eop(3)).id==1
   globsol.dut1.val=x(parGS(g.g_eop(3)).newcol)/15; % mas --> ms
   globsol.dut1.sigma=varpar(parGS(g.g_eop(3)).newcol)/15;    % mas --> ms
   globsol.dut1.mjd=parGS(g.g_eop(3)).mjdstep';
   globsol.dut1.col=parGS(g.g_eop(3)).newcol;
end

if parGS(g.g_eop(4)).id==1
   globsol.dX.val=x(parGS(g.g_eop(4)).newcol);
   globsol.dX.sigma=varpar(parGS(g.g_eop(4)).newcol);
   globsol.dX.mjd=parGS(g.g_eop(4)).mjdstep';
   globsol.dX.col=parGS(g.g_eop(4)).newcol;
end

if parGS(g.g_eop(5)).id==1
   globsol.dY.val=x(parGS(g.g_eop(5)).newcol);
   globsol.dY.sigma=varpar(parGS(g.g_eop(5)).newcol);
   globsol.dY.mjd=parGS(g.g_eop(5)).mjdstep';
   globsol.dY.col=parGS(g.g_eop(5)).newcol;
end

if parGS(g.g_ao).id==1
    globsol.antenna.aoname = aostname; 
    globsol.antenna.ao = x(parGS(g.g_ao).newcol);
    globsol.antenna.sigma_ao = varpar(parGS(g.g_ao).newcol);
    globsol.antenna.col_ao = parGS(g.g_ao).newcol;
end


if parGS(g.g_aplrg).id==1
    globsol.antenna.rgname = rgstname; 
    globsol.antenna.rg = x(parGS(g.g_aplrg).newcol); 
    globsol.antenna.sigma_rg = varpar(parGS(g.g_aplrg).newcol); 
    globsol.antenna.col_rg = parGS(g.g_aplrg).newcol; 
    globsol.antenna.apriori_rg = [apriori_aplrg(:,1) apriori_aplrg(:,2).*100]; % hPa cm/hPa
end
if parGS(g.g_stseaspos_Ar(1)).id==1 || parGS(g.g_stseaspos_Ae(1)).id==1 || parGS(g.g_stseaspos_An(1)).id==1
    id=find(cell2mat(parGS(g.g_stseaspos_Ar(1)).spec)==1);
    globsol.stseaspos.periods_solardays = P_stseaspos(id);
    for i=1:size(stseasname)
        globsol.stseaspos.results(i).aname=stseasname(i,:);
    end
    if parGS(g.g_stseaspos_Ar(1)).id==1 
        for i=1:size(stseasname)
            c=[]; s=[]; A=[]; ph=[]; sigma_ph=[];
            c=x(cell2mat(parGS(g.g_stseaspos_Ar(1)).newcol(i)));
            s=x(cell2mat(parGS(g.g_stseaspos_Ar(2)).newcol(i)));
            globsol.stseaspos.results(i).Acr=c;
            globsol.stseaspos.results(i).Asr=s;
            
            sig_c=varpar(cell2mat(parGS(g.g_stseaspos_Ar(1)).newcol(i)));
            sig_s=varpar(cell2mat(parGS(g.g_stseaspos_Ar(2)).newcol(i)));
            globsol.stseaspos.results(i).sigma_Acr=sig_c;
            globsol.stseaspos.results(i).sigma_Asr=sig_s;
            
            A=sqrt(c.^2+s.^2);
            ph=atan2(s,c);
            sigma_ph=sqrt((c./(c.^2+s.^2)).^2 .* sig_s.^2 + (s./(c.^2+s.^2)).^2 .* sig_c.^2);
            globsol.stseaspos.results(i).Ar=A;
            globsol.stseaspos.results(i).sigma_Ar=sqrt(c.^2.*sig_c.^2 + s.^2.*sig_s.^2)./A;
            globsol.stseaspos.results(i).phaser=ph./pi.*180; %rad --> deg
            globsol.stseaspos.results(i).sigma_phaser=sigma_ph./pi.*180;
        end
    end
    if parGS(g.g_stseaspos_Ae(1)).id==1 
        for i=1:size(stseasname)
            c=[]; s=[]; A=[]; ph=[]; sigma_ph=[];
            c=x(cell2mat(parGS(g.g_stseaspos_Ae(1)).newcol(i)));
            s=x(cell2mat(parGS(g.g_stseaspos_Ae(2)).newcol(i)));
            globsol.stseaspos.results(i).Ace=c;
            globsol.stseaspos.results(i).Ase=s;
            
            sig_c=varpar(cell2mat(parGS(g.g_stseaspos_Ae(1)).newcol(i)));
            sig_s=varpar(cell2mat(parGS(g.g_stseaspos_Ae(2)).newcol(i)));
            globsol.stseaspos.results(i).sigma_Ace=sig_c;
            globsol.stseaspos.results(i).sigma_Ase=sig_s;
            
            A=sqrt(c.^2+s.^2);
            ph=atan2(s,c);
            sigma_ph=sqrt((c./(c.^2+s.^2)).^2 .* sig_s.^2 + (s./(c.^2+s.^2)).^2 .* sig_c.^2);
            globsol.stseaspos.results(i).Ae=A;
            globsol.stseaspos.results(i).sigma_Ae=sqrt(c.^2.*sig_c.^2 + s.^2.*sig_s.^2)./A;
            globsol.stseaspos.results(i).phasee=ph./pi.*180; %rad --> deg
            globsol.stseaspos.results(i).sigma_phasee=sigma_ph./pi.*180;


        end
    end
    if parGS(g.g_stseaspos_An(1)).id==1 
        for i=1:size(stseasname)
            c=[]; s=[]; A=[]; ph=[]; sigma_ph=[];
            c=x(cell2mat(parGS(g.g_stseaspos_An(1)).newcol(i)));
            s=x(cell2mat(parGS(g.g_stseaspos_An(2)).newcol(i)));
            globsol.stseaspos.results(i).Acn=c;
            globsol.stseaspos.results(i).Asn=s;
            
            sig_c=varpar(cell2mat(parGS(g.g_stseaspos_An(1)).newcol(i)));
            sig_s=varpar(cell2mat(parGS(g.g_stseaspos_An(2)).newcol(i)));
            globsol.stseaspos.results(i).sigma_Acn=sig_c;
            globsol.stseaspos.results(i).sigma_Asn=sig_s;
            
            A=sqrt(c.^2+s.^2);
            ph=atan2(s,c);
            sigma_ph=sqrt((c./(c.^2+s.^2)).^2 .* sig_s.^2 + (s./(c.^2+s.^2)).^2 .* sig_c.^2);
            globsol.stseaspos.results(i).An=A;
            globsol.stseaspos.results(i).sigma_An=sqrt(c.^2.*sig_c.^2 + s.^2.*sig_s.^2)./A;
            globsol.stseaspos.results(i).phasen=ph./pi.*180; %rad --> deg
            globsol.stseaspos.results(i).sigma_phasen=sigma_ph./pi.*180;
        end
    end
end

if parGS(g.g_hpole).id==1
    globsol.hpole.val=x(parGS(g.g_hpole).newcol);
    globsol.hpole.sigma=varpar(parGS(g.g_hpole).newcol);
    globsol.hpole.col=parGS(g.g_hpole).newcol;
end

if parGS(g.g_lpole).id==1
    globsol.lpole.val=x(parGS(g.g_lpole).newcol);
    globsol.lpole.sigma=varpar(parGS(g.g_lpole).newcol);
    globsol.lpole.col=parGS(g.g_lpole).newcol;
end


if parGS(g.g_love).id==1
    clear i1 i2
    [i1,i2]=find(parGS_hl.love.nr==1);
    globsol.dlove.val_h2=x(parGS(g.g_love).newcol(i2));
    globsol.dlove.sigma_h2=varpar(parGS(g.g_love).newcol(i2));
    globsol.dlove.nr_h2=parGS_hl.love.nr(i2)';
    globsol.dlove.col_h2=parGS(g.g_love).newcol(i2)';
    if isempty(parGS_hl.love.nr_d)==0
        clear i1 i2
        [i1,i2]=find(5<parGS_hl.love.nr & parGS_hl.love.nr<37);
        globsol.dlove.val_diurnal=x(parGS(g.g_love).newcol(i2));
        globsol.dlove.sigma_diurnal=varpar(parGS(g.g_love).newcol(i2));
        globsol.dlove.nr_diurnal=parGS_hl.love.nr(i2)';
        globsol.dlove.col_diurnal=parGS(g.g_love).newcol(i2)';
    end
    if isempty(parGS_hl.love.nr_di)==0
        clear i1 i2
        [i1,i2]=find(36<parGS_hl.love.nr & parGS_hl.love.nr<68);
        globsol.dlove.val_diurnal_imag=x(parGS(g.g_love).newcol(i2));
        globsol.dlove.sigma_diurnal_imag=varpar(parGS(g.g_love).newcol(i2));
        globsol.dlove.nr_diurnal_imag=parGS_hl.love.nr(i2)';
        globsol.dlove.col_diurnal_imag=parGS(g.g_love).newcol(i2)';
    end
    
    if isempty(parGS_hl.love.nr_long)==0
        clear i1 i2
        [i1,i2]=find(67<parGS_hl.love.nr & parGS_hl.love.nr<73);
        globsol.dlove.val_long=x(parGS(g.g_love).newcol(i2));
        globsol.dlove.sigma_long=varpar(parGS(g.g_love).newcol(i2));
        globsol.dlove.nr_long=parGS_hl.love.nr(i2)';
        globsol.dlove.col_long=parGS(g.g_love).newcol(i2)';
    end
    if isempty(parGS_hl.love.nr_longi)==0
        clear i1 i2
        [i1,i2]=find(72<parGS_hl.love.nr & parGS_hl.love.nr<78);
        globsol.dlove.val_long_imag=x(parGS(g.g_love).newcol(i2));
        globsol.dlove.sigma_long_imag=varpar(parGS(g.g_love).newcol(i2));
        globsol.dlove.nr_long_imag=parGS_hl.love.nr(i2)';
        globsol.dlove.col_long_imag=parGS(g.g_love).newcol(i2)';
    end
    
end

if parGS(g.g_shida).id==1
    clear i1 i2
    [i1,i2]=find(parGS_hl.shida.nr==1);
    globsol.dshida.val_l2=x(parGS(g.g_shida).newcol(i2));
    globsol.dshida.sigma_l2=varpar(parGS(g.g_shida).newcol(i2));
    globsol.dshida.nr_l2=parGS_hl.shida.nr(i2)';
    globsol.dshida.col_l2=parGS(g.g_shida).newcol(i2)';
    if isempty(parGS_hl.shida.nr_d)==0
        [i1,i2]=find(7<parGS_hl.shida.nr & parGS_hl.shida.nr<39);
        globsol.dshida.val_diurnal=x(parGS(g.g_shida).newcol(i2));
        globsol.dshida.sigma_diurnal=varpar(parGS(g.g_shida).newcol(i2));
        globsol.dshida.nr_diurnal=parGS_hl.shida.nr(i2)';
        globsol.dshida.col_diurnal=parGS(g.g_shida).newcol(i2)';
    end
    if isempty(parGS_hl.shida.nr_di)==0
        clear i1 i2
        [i1,i2]=find(38<parGS_hl.shida.nr & parGS_hl.shida.nr<70);
        globsol.dshida.val_diurnal_imag=x(parGS(g.g_shida).newcol(i2));
        globsol.dshida.sigma_diurnal_imag=varpar(parGS(g.g_shida).newcol(i2));
        globsol.dshida.nr_diurnal_imag=parGS_hl.shida.nr(i2)';
        globsol.dshida.col_diurnal_imag=parGS(g.g_shida).newcol(i2)';
    end
    
    if isempty(parGS_hl.shida.nr_long)==0
        clear i1 i2
        [i1,i2]=find(69<parGS_hl.shida.nr & parGS_hl.shida.nr<75);
        globsol.dshida.val_long=x(parGS(g.g_shida).newcol(i2));
        globsol.dshida.sigma_long=varpar(parGS(g.g_shida).newcol(i2));
        globsol.dshida.nr_long=parGS_hl.shida.nr(i2)';
        globsol.dshida.col_long=parGS(g.g_shida).newcol(i2)';
    end
    if isempty(parGS_hl.shida.nr_longi)==0
        clear i1 i2
        [i1,i2]=find(74<parGS_hl.shida.nr & parGS_hl.shida.nr<80);
        globsol.dshida.val_long_imag=x(parGS(g.g_shida).newcol(i2));
        globsol.dshida.sigma_long_imag=varpar(parGS(g.g_shida).newcol(i2));
        globsol.dshida.nr_long_imag=parGS_hl.shida.nr(i2)';
        globsol.dshida.col_long_imag=parGS(g.g_shida).newcol(i2)';
    end
    
end

if parGS(g.g_FCNset).id==1
    globsol.dNDFWset.val = x(parGS(g.g_FCNset).newcol);
    globsol.dNDFWset.sigma = varpar(parGS(g.g_FCNset).newcol);
        
    dNDFW=globsol.dNDFWset.val; sNDFW=globsol.dNDFWset.sigma;
    
    NDFW_virt=1.0023181; % FCN
    
    %dFCN_val = -1/(1/(dNDFW+NDFW_virt)-1) + 1/(1/NDFW_virt-1);  % sid. days in CRF
    FCN_apr=1/(1-(NDFW_virt));
    FCN_new=1/(1-(NDFW_virt+dNDFW));
    dFCN_sigma= abs( 1/(1-(sNDFW+NDFW_virt+dNDFW)) - 1/(1-(NDFW_virt+dNDFW)) );  % sid. days in CRF
    
    
    globsol.dNDFWset.FCN_apr = FCN_apr;
    globsol.dNDFWset.FCN_new=FCN_new;
    globsol.dNDFWset.dFCN_sigma=dFCN_sigma;
    
    globsol.dNDFWset.col = parGS(g.g_FCNset).newcol;
      
end


c = 29979245800; % light velocity in cm/sec
if parGS(g.g_accSSB).id==1
    globsol.daccSSB.val = x(parGS(g.g_accSSB).newcol)./c; %cm^2/sec^3 --> cm/sec^2
    globsol.daccSSB.sigma = varpar(parGS(g.g_accSSB).newcol)./c; %cm^2/sec^3 --> cm/sec^2
end


if parGS(g.g_gamma).id==1
    globsol.gamma.val=x(parGS(g.g_gamma).newcol);
    globsol.gamma.sigma=varpar(parGS(g.g_gamma).newcol);
    globsol.gamma.col=parGS(g.g_gamma).newcol;
end

% Recover sine and cosine terms of tidal ERP amplitudes

if parGS(g.g_tidpm).id==1
    path_level='../';
    load([path_level 'DATA/GLOB/tide'],'tide');
    pro = find(tide.gmstpi<2);
    numpro = length(parGS(g.g_tidpm).spectid);
    numret = length(parGS(g.g_tidpm).spectidret);
    globsol.tidAp.val = x(parGS(g.g_tidpm).newcol(1:numpro))*1e3; % mas --> uas
    globsol.tidBp.val = x(parGS(g.g_tidpm).newcol(numpro+1:numpro*2))*1e3; % mas --> uas
    globsol.tidAm.val = x(parGS(g.g_tidpm).newcol(numpro*2+1:numpro*2+numret))*1e3; % mas --> uas
    globsol.tidBm.val = x(parGS(g.g_tidpm).newcol(numpro*2+numret+1:numpro*2+numret*2))*1e3; % mas --> uas
    globsol.tidAp.sigma = varpar(parGS(g.g_tidpm).newcol(1:numpro))*1e3; % mas --> uas
    globsol.tidBp.sigma = varpar(parGS(g.g_tidpm).newcol(numpro+1:numpro*2))*1e3; % mas --> uas
    globsol.tidAm.sigma = varpar(parGS(g.g_tidpm).newcol(numpro*2+1:numpro*2+numret))*1e3; % mas --> uas
    globsol.tidBm.sigma = varpar(parGS(g.g_tidpm).newcol(numpro*2+numret+1:numpro*2+numret*2))*1e3; % mas --> uas
    globsol.tidAp.col = parGS(g.g_tidpm).newcol(1:numpro);
    globsol.tidBp.col = parGS(g.g_tidpm).newcol(numpro+1:numpro*2);
    globsol.tidAm.col = parGS(g.g_tidpm).newcol(numpro*2+1:numpro*2+numret);
    globsol.tidBm.col = parGS(g.g_tidpm).newcol(numpro*2+numret+1:numpro*2+numret*2);
    globsol.tidnum.pro = parGS(g.g_tidpm).spectid;
    globsol.tidnum.ret = parGS(g.g_tidpm).spectid(parGS(g.g_tidpm).spectid>pro(end));
end

if parGS(g.g_tidut).id==1
    numut = length(parGS(g.g_tidut).spectid);
    globsol.tidutc.val = x(parGS(g.g_tidut).newcol(1:numut))/15*1e3; % mas --> us
    globsol.tiduts.val = x(parGS(g.g_tidut).newcol(numut+1:numut*2))/15*1e3; % mas --> us
    globsol.tidutc.sigma = varpar(parGS(g.g_tidut).newcol(1:numut))/15*1e3; % mas --> us
    globsol.tiduts.sigma = varpar(parGS(g.g_tidut).newcol(numut+1:numut*2))/15*1e3; % mas --> us
    globsol.tidutc.col = parGS(g.g_tidut).newcol(1:numut);
    globsol.tiduts.col = parGS(g.g_tidut).newcol(numut+1:numut*2);
    globsol.tidnum.ut1 = parGS(g.g_tidut).spectid;
end

%%