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

function col_sinex=newcol_sinex(col_est,x_,antenna)

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
    
    %%
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
    
    %%
    for iant = 1 : length(antenna)
       col_sinex.antnames(iant).name = antenna(iant).name;
    end
    col_sinex.coorx=newcol_x;
    col_sinex.coory=newcol_y;
    col_sinex.coorz=newcol_z;
    
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
    
    col_sinex.sounames=[];
    for isou = 1 : length(x_.source)
       col_sinex.sounames(isou).name = x_.source(isou).name;
    end
    col_sinex.ra=newcol_ra;
    col_sinex.de=newcol_de;

    
    