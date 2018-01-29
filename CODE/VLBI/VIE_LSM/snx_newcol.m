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

function col_sinex=snx_newcol(col_est,x_,antenna,outsnx)


    % x-coordinate
    clear old
    old=[x_.coorx.col]; newcol_x=[];
    for i=1:length(old)
        [a,newcol_x(i)]=find(old(i)==col_est);
    end

    % y-coordinate
    clear old
    old=[x_.coory.col]; newcol_y=[];
    for i=1:length(old)
        [a,newcol_y(i)]=find(old(i)==col_est);
    end
    
    % z-coordinate
    clear old
    old=[x_.coorz.col]; newcol_z=[];
    for i=1:length(old)
        [a,newcol_z(i)]=find(old(i)==col_est);
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
    
    col_sinex.coorx=newcol_x;
    col_sinex.coory=newcol_y;
    col_sinex.coorz=newcol_z;
    
    
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
     
    
    
    
    
    