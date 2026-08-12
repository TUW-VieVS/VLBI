% ************************************************************************
%   Description:
%   function for description of the columns in the N_sinex matrix
%
%   Reference: 
%
%   Input:										
%       col_est      vector with the old order of parameters
%       x_           information about the estimated parameters
%       antenna      order of stations
%
%   Output:
%       col_sinex    structure array     information about the columns in
%                                        N_sinex
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   05 Oct 2010 by Hana Spicakova
% ************************************************************************ 

function col_sinex=snx_newcol(col_est,x_,antenna,outsnx, parameter)

    if parameter.lsmopt.stc
        if parameter.lsmopt.stc_all
            cx = "coorx";
            cy = "coory";
            cz = "coorz";
        elseif parameter.lsmopt.stc_sat
            cx = "coorx_sat";
            cy = "coory_sat";
            cz = "coorz_sat";
        elseif parameter.lsmopt.stc_qu
            cx = "coorx_qu";
            cy = "coory_qu";
            cz = "coorz_qu";
        elseif parameter.lsmopt.stc_qs
            cx = ["coorx_sat", "coorx_qu"];
            cy = ["coory_sat", "coory_qu"];
            cz = ["coorz_sat", "coorz_qu"];
        end
    else
        cx = [];
        cy = [];
        cz = [];
    end

    SRP_Names = ["D0"; "Y0"; "B0"; "DC"; "YC"; "BC"; "DS"; "YS"; "BS"];
    ORB_Names = ["sma"; "ecc"; "inc"; "raan"; "argp"; "argl"];

    % x-coordinate
    newcol_x=[];
    for j=1:length(cx)
        clear old
        old=[x_.(cx(j)).col]; 
        for i=1:length(old)
            [a,newcol_x(j,i)]=find(old(i)==col_est);
        end
    end

    % y-coordinate
    newcol_y=[];
    for j=1:length(cy)
        clear old
        old=[x_.(cy(j)).col]; 
        for i=1:length(old)
            [a,newcol_y(j,i)]=find(old(i)==col_est);
        end
    end
    
    % z-coordinate
    newcol_z=[];
    for j=1:length(cz)
        clear old
        old=[x_.(cz(j)).col]; 
        for i=1:length(old)
            [a,newcol_z(j,i)]=find(old(i)==col_est);
        end
    end

    %%
    if outsnx.eop==1
        % xpole
        clear old
        old=[x_.xpol.col]; newcol_xp=[];
        for i=1:length(old)
            [a,newcol_xp(i)]=find(old(i)==col_est);
        end

        % ypole
        clear old
        old=[x_.ypol.col]; newcol_yp=[];
        for i=1:length(old)
            [a,newcol_yp(i)]=find(old(i)==col_est);
        end

        % dut1
        clear old
        old=[x_.dut1.col]; newcol_dut1=[];
        for i=1:length(old)
            [a,newcol_dut1(i)]=find(old(i)==col_est);
        end

         % dX
        clear old
        old=[x_.nutdx.col]; newcol_dX=[];
        for i=1:length(old)
            [a,newcol_dX(i)]=find(old(i)==col_est);
        end

         % dY
        clear old
        old=[x_.nutdy.col]; newcol_dY=[];
        for i=1:length(old)
            [a,newcol_dY(i)]=find(old(i)==col_est);
        end
    end
    %%
    if outsnx.sou==1
         % RA
        clear old
        old=x_.col_soura; newcol_ra=[];
        for i=1:length(old)
            [a,newcol_ra(i)]=find(old(i)==col_est);
        end

        % De
        clear old
        old=x_.col_soude; newcol_de=[];
        for i=1:length(old)
            [a,newcol_de(i)]=find(old(i)==col_est);
        end
    end
    %%
    % zwd
    if outsnx.zwd==1
        for iant = 1 : length(antenna)
            clear old
            old=[x_.zwd(iant).col]; 
            for i=1:length(old)
            	[a,newcol.zwd(iant).col(i)]=find(old(i)==col_est);
            end
        end
    end
    %%
    % troposphere gradients
    if outsnx.tgr==1
        for iant = 1 : length(antenna)
            clear old
            old=[x_.ngr(iant).col]; newcol.ngr=[];
            for i=1:length(old)
                [a,newcol_ngr(iant).col(i)]=find(old(i)==col_est);
            end
        end

        for iant = 1 : length(antenna)
            clear old
            old=[x_.egr(iant).col]; newcol.egr=[];
            for i=1:length(old)
                [a,newcol_egr(iant).col(i)]=find(old(i)==col_est);
            end
        end
    end

    %%
    %Keplerelemente
    if outsnx.orb==1
        % KepEle1
        clear old
        old=[x_.ORB.sma.col]; newcol_sma=[];
        for i=1:length(old)
            [a,newcol_sma(i)]=find(old(i)==col_est);
        end

        % KepEle2
        clear old
        old=[x_.ORB.ecc.col]; newcol_ecc=[];
        for i=1:length(old)
            [a,newcol_ecc(i)]=find(old(i)==col_est);
        end

        % KepEle3
        clear old
        old=[x_.ORB.inc.col]; newcol_inc=[];
        for i=1:length(old)
            [a,newcol_inc(i)]=find(old(i)==col_est);
        end

        % KepEle4
        clear old
        old=[x_.ORB.raan.col]; newcol_raan=[];
        for i=1:length(old)
            [a,newcol_raan(i)]=find(old(i)==col_est);
        end

        % KepEle5
        clear old
        old=[x_.ORB.argp.col]; newcol_argp=[];
        for i=1:length(old)
            [a,newcol_argp(i)]=find(old(i)==col_est);
        end

        % KepEle6
        clear old
        old=[x_.ORB.argl.col]; newcol_argl=[];
        for i=1:length(old)
            [a,newcol_argl(i)]=find(old(i)==col_est);
        end

        % SRP D0
        clear old
        old=[x_.SRP.D0.col]; newcol_srp_D0=[];
        for i=1:length(old)
            [a,newcol_srp_D0(i)]=find(old(i)==col_est);
        end

        % SRP Y0
        clear old
        old=[x_.SRP.Y0.col]; newcol_srp_Y0=[];
        for i=1:length(old)
            [a,newcol_srp_Y0(i)]=find(old(i)==col_est);
        end

        % SRP B0
        clear old
        old=[x_.SRP.B0.col]; newcol_srp_B0=[];
        for i=1:length(old)
            [a,newcol_srp_B0(i)]=find(old(i)==col_est);
        end

        % SRP DC
        clear old
        old=[x_.SRP.DC.col]; newcol_srp_DC=[];
        for i=1:length(old)
            [a,newcol_srp_DC(i)]=find(old(i)==col_est);
        end

        % SRP YC
        clear old
        old=[x_.SRP.YC.col]; newcol_srp_YC=[];
        for i=1:length(old)
            [a,newcol_srp_YC(i)]=find(old(i)==col_est);
        end

        % SRP BC
        clear old
        old=[x_.SRP.BC.col]; newcol_srp_BC=[];
        for i=1:length(old)
            [a,newcol_srp_BC(i)]=find(old(i)==col_est);
        end

        % SRP DS
        clear old
        old=[x_.SRP.DS.col]; newcol_srp_DS=[];
        for i=1:length(old)
            [a,newcol_srp_DS(i)]=find(old(i)==col_est);
        end

        % SRP YS
        clear old
        old=[x_.SRP.YS.col]; newcol_srp_YS=[];
        for i=1:length(old)
            [a,newcol_srp_YS(i)]=find(old(i)==col_est);
        end

        % SRP BS
        clear old
        old=[x_.SRP.BS.col]; newcol_srp_BS=[];
        for i=1:length(old)
            [a,newcol_srp_BS(i)]=find(old(i)==col_est);
        end

        satnames = {};
        for i = 1:length(ORB_Names)
            v = x_.ORB.(ORB_Names{i});
            if ~isempty(v.val)
                for j=1:length(v.val)
                    satnames = [satnames; {v.name}];  % Vektorisierte Extraktion aller Names auf einmal
                end
            end
        end
        for i = 1:length(SRP_Names)
            v = x_.SRP.(SRP_Names{i});
            if ~isempty(v.val)
                for j=1:length(v.val)
                    satnames = [satnames; {v.name}];  % Vektorisierte Extraktion aller Names auf einmal
                end
            end
        end
        satnames = unique(satnames, 'stable'); 
    end

    
    
    %%
    
    for iant = 1 : length(antenna)
        col_sinex.zwd(iant).col=[];
        col_sinex.zwd(iant).mjd=[];
        col_sinex.ngr(iant).col=[];
        col_sinex.ngr(iant).mjd=[];
        col_sinex.egr(iant).col=[];
        col_sinex.egr(iant).mjd=[];
    end
    
    
    for iant = 1 : length(antenna)
        col_sinex.antnames(iant).name = antenna(iant).name;
       
       %zwd
        if outsnx.zwd==1
            col_sinex.zwd(iant).col=newcol.zwd(iant).col;
            col_sinex.zwd(iant).mjd=x_.zwd(iant).mjd;
        end
        
       %ngr, egr
        if outsnx.tgr==1
            col_sinex.ngr(iant).col=newcol_ngr(iant).col;
            col_sinex.ngr(iant).mjd=x_.ngr(iant).mjd;
            col_sinex.egr(iant).col=newcol_egr(iant).col;
            col_sinex.egr(iant).mjd=x_.egr(iant).mjd;
        end
    end
    
    for j=1:length(cx)
        col_sinex.(cx(j))=newcol_x(j,:);
        col_sinex.(cy(j))=newcol_y(j,:);
        col_sinex.(cz(j))=newcol_z(j,:);
    end
    
    
    col_sinex.xp=[];
    col_sinex.yp=[];
    col_sinex.dut1=[];
    col_sinex.dX=[];
    col_sinex.dY=[];
    col_sinex.mjd_xp=[];
    col_sinex.mjd_yp=[];
    col_sinex.mjd_dut1=[];
    col_sinex.mjd_dX=[];
    col_sinex.mjd_dY=[];

    if outsnx.eop==1
        col_sinex.xp=newcol_xp;
        col_sinex.yp=newcol_yp;
        col_sinex.dut1=newcol_dut1;
        col_sinex.dX=newcol_dX;
        col_sinex.dY=newcol_dY;
        col_sinex.mjd_xp=x_.xpol.mjd;
        col_sinex.mjd_yp=x_.ypol.mjd;
        col_sinex.mjd_dut1=x_.dut1.mjd;
        col_sinex.mjd_dX=x_.nutdx.mjd;
        col_sinex.mjd_dY=x_.nutdy.mjd;
    end
    
    col_sinex.sounames=[];
    col_sinex.ra=[];
    col_sinex.de=[];
    if outsnx.sou==1 % sources cannot be reduced in this version!!!
        for isou = 1 : length(x_.source)
           col_sinex.sounames(isou).name = x_.source(isou).name;
        end
        col_sinex.ra=newcol_ra;
        col_sinex.de=newcol_de;
    end

    col_sinex.orb_sma=[];
    col_sinex.orb_ecc=[];
    col_sinex.orb_inc=[];
    col_sinex.orb_raan=[];
    col_sinex.orb_argp=[];
    col_sinex.orb_argl=[];
    col_sinex.srp_D0=[];
    col_sinex.srp_Y0=[];
    col_sinex.srp_B0=[];
    col_sinex.srp_DC=[];
    col_sinex.srp_YC=[];
    col_sinex.srp_BC=[];
    col_sinex.srp_DS=[];
    col_sinex.srp_YS=[];
    col_sinex.srp_BS=[];

    col_sinex.satnames=[];
    if outsnx.orb
        col_sinex.orb_sma=newcol_sma;
        col_sinex.orb_ecc=newcol_ecc;
        col_sinex.orb_inc=newcol_inc;
        col_sinex.orb_raan=newcol_raan;
        col_sinex.orb_argp=newcol_argp;
        col_sinex.orb_argl=newcol_argl;
        col_sinex.srp_D0=newcol_srp_D0;
        col_sinex.srp_Y0=newcol_srp_Y0;
        col_sinex.srp_B0=newcol_srp_B0;
        col_sinex.srp_DC=newcol_srp_DC;
        col_sinex.srp_YC=newcol_srp_YC;
        col_sinex.srp_BC=newcol_srp_BC;
        col_sinex.srp_DS=newcol_srp_DS;
        col_sinex.srp_YS=newcol_srp_YS;
        col_sinex.srp_BS=newcol_srp_BS;

        col_sinex.mjd_orb_sma=x_.ORB.sma.mjd;
        col_sinex.mjd_orb_ecc=x_.ORB.ecc.mjd;
        col_sinex.mjd_orb_inc=x_.ORB.inc.mjd;
        col_sinex.mjd_orb_raan=x_.ORB.raan.mjd;
        col_sinex.mjd_orb_argp=x_.ORB.argp.mjd;
        col_sinex.mjd_orb_argl=x_.ORB.argl.mjd;

        col_sinex.mjd_srp_D0=x_.SRP.D0.mjd;
        col_sinex.mjd_srp_Y0=x_.SRP.Y0.mjd;
        col_sinex.mjd_srp_B0=x_.SRP.B0.mjd;
        col_sinex.mjd_srp_DC=x_.SRP.DC.mjd;
        col_sinex.mjd_srp_YC=x_.SRP.YC.mjd;
        col_sinex.mjd_srp_BC=x_.SRP.BC.mjd;
        col_sinex.mjd_srp_DS=x_.SRP.DS.mjd;
        col_sinex.mjd_srp_YS=x_.SRP.YS.mjd;
        col_sinex.mjd_srp_BS=x_.SRP.BS.mjd;

        col_sinex.satnames = satnames;
    end