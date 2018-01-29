% To calculate the residuals it is necessary to calculate again the design
% matrixes for every scan

%   Revision: 
%   19 Jun 2014 by Hana Krasna: sources can be estimated in vie_lsm with
%   NNR condition




disp(' ')
disp('8. RE-CREATION OF DESIGN MATRIXES TO CALCULATE RESIDUALS');

v_total=[];
obs_num=1;
t=([scan.mjd]-mjd0).*24.*60;     % minutes from midnight to the time of the scan
v_real=[];


if opt.basic_outlier == 1
    AQh=[];
end


for iscan=1:ntim  % loop over all scans   
    
   
    %~~~~~~~~~~~~~~~~~~
    % CREATION OF DESIGN MATRICES
  
    [per_stat,T_]=assign_parameters_scan(scan, na,opt,iscan,antenna,t,mjd0,T_); %assign input parameters for each station

     
     
     oc_observ=sparse(oc_observ_all(obs_num:obs_num+scan(iscan).nobs-1,1));
     obs_num=obs_num+scan(iscan).nobs;
 
    %~~~~~~~~~
    % DESIGN MATRICES OF REAL OBSERVATIONS
   
    % Design matrix for pwl clock offsets, rate and quadratic terms
    [Apwclk, Arqclk,num_inter_qclk] = apwq_clk_scan(per_stat,T_,opt,na,t,iscan);  
    % Design matrix for pwl zwd estimates
    [Azwd] = apw_zwd_scan(per_stat,T_,na,t,iscan);
    % Design matrix for north gradients
    [Apwngr] = apw_ngr_scan(per_stat,T_,na);
    % Design matrix for east gradients
    [Apwegr] = apw_egr_scan(per_stat,T_,na);

    if opt.pw_stc==0 %if NO station coord. is estimated
     % Design matrix of the station dx,dy, dz coordinate offset estimate
    [Ax, Ay, Az] = a_xyz_scan(per_stat,na);
    else if opt.pw_stc==1 % if station coord. is estimated as pwl offsets
     % Design matrix for antenna dx,dy,dz coordinates
     [Apwx, Apwy, Apwz] = apw_xyz_scan(per_stat,T_,na);    
        end
    end
    
   
    %~~~~~~~~~
    % DESIGN MATRICES FOR SOURCES COORDINATES
  
    [per_source,T_source] = sourcewisepar_scan(opt,scan,iscan,sources,t,mjd_scan);
   
    if  (opt.global_solve==1 && opt.est_source==1) || (opt.ascii_snx==1 && opt.est_source==1)|| opt.est_sourceNNR==1; % +hana 18Jun14
         [A_ra_glob,A_de_glob]=a_source_scan(per_source);
         if opt.est_sourceNNR==1
             Ara=A_ra_glob;
             Ade=A_de_glob;
         end

    else 
         if opt.pw_sou == 1 % if pwl offsets are estimated
       [Apw_ra,Apw_de] = apw_source_scan(per_source,T_source);
         end
    end
    



    %~~~~~~~~~
    % DESIGN MATRICES FOR EOP
  
    % Number of the stations observing in this scan
    %obs_stat=find(~cellfun(@isempty, {scan(iscan).stat.temp})); 
    
    % Xpol (piecewise) - ahp_xpol
    [Apwxpol,T_xpol] = ahp_xpol_scan_A(scan,iscan,opt,c,rad2mas,t);

    % Ypol (piecewise) - ahp_ypol
    [Apwypol,T_ypol] = ahp_ypol_scan_A(scan,iscan,c,rad2mas,opt,t);

    % dUT1 (piecewise) - ahp_dut1
    [Apwdut1,T_dut1] = ahp_dut1_scan_A(scan,iscan,c,rad2mas,opt,t);

    % Celestial Pole Offsets dx (DEPS) (piecewise) - ahp_nutdx
    [Apwnutdx,T_nutdx] = ahp_nutdx_scan_A(scan,iscan,c,rad2mas,opt,t);

    % Celestial Pole Offsets dy (DPSI) (piecewise) - ahp_nutdy
    [Apwnutdy,T_nutdy] = ahp_nutdy_scan_A(scan,iscan,c,rad2mas,opt,t);
    %TIEMPO(4)=toc;
    
    %~~~~~~~~~
    % ONE DESIGN MATRIX PER SCAN

    A(1).sm= horzcat(Apwclk.sm);  % clk piecewise linear offset  
    A(2).sm= horzcat(Arqclk.sm);  % clk rate and quadratic term    
    A(3).sm= horzcat(Azwd.sm);    % zwd
    A(4).sm= horzcat(Apwngr.sm);  % ngr
    A(5).sm= horzcat(Apwegr.sm);  % egr
    A(6).sm=Apwxpol;              % Xpol
    A(7).sm=Apwypol;              % Ypol
    A(8).sm=Apwdut1;              % dUT1
    A(9).sm=Apwnutdx;             % nutdx
    A(10).sm=Apwnutdy;            % nutdy
   
    % Station coordinates
    if opt.pw_stc==0
       A(13).sm= horzcat(Ax.sm); 
       A(14).sm= horzcat(Ay.sm); 
       A(15).sm= horzcat(Az.sm); 
    else if opt.pw_stc==1
       A(13).sm= horzcat(Apwx.sm); 
       A(14).sm= horzcat(Apwy.sm); 
       A(15).sm= horzcat(Apwz.sm);    
        end
    end
   
    % Sources
    if opt.pw_sou == 1
       A(11).sm= Apw_ra;
       A(12).sm= Apw_de;
    end
    if opt.est_sourceNNR==1
       A(11).sm= Ara;
       A(12).sm= Ade;
    end
    
    clear nmi_observ so  Qll  Apwclk Apwdut1 Apwegr Apwngr
    clear Apwnutdx Apwnutdy Apwxpol Apwypol Arqclk Ax Ay Az Azwd 
    clear Apw_ra Apw_de
    
       
    %~~~~~~~~~
    % CONCATENATING so that A scanwise has the same number of columns as A
    % sessionwise 
   
    [A_session]=concat_sess_scan(scan,iscan,A,opt,t,t_first,na,ns...
        ,num_inter_xpol,num_inter_ypol,num_inter_dut1,num_inter_nutdx...
        ,num_inter_nutdy,num_inter_egr,num_inter_ngr,num_inter_clk...
        ,num_inter_zwd,per_source,num_inter_sou,t_first_sou);
    col_A2=size(A_session(2).sm,2);
    col_A4=size(A_session(4).sm,2);
    col_A5=size(A_session(5).sm,2);
    
    clear A 


    %~~~~~~~~~ 
    % EXCLUDING REFERENCE STATION, PARAMETERS AND SOURCES
    % Excluding the offset parameters of the specified station's zwd, ngr, egr
    % and coordinates
     
     [A_session] = delparam_scan_A(na,opt,A_session,num_inter_zwd,num_inter_ngr,num_inter_egr);
 
    % Excluding the specified (fixed) sources from the design matrix
    if opt.pw_sou == 1;
       [A_session] = delsource_scan_A(opt,A_session,num_inter_sou,ns);    
    end

    % Deleting reference clock
    [A_session] = delref_scan_A(A_session,num_inter_clk,nistat,opt);
    
    % Excluding models
    [A_session] = delmodel_scan_A(opt,A_session);

   
    %~~~~~~~~~
    % OBTAINING ABLK
    
    Ablk=[]; 
   
    % Concatenating A_session to obtain Ablk
    dj=zeros(1,15);
    for i=1:15 
      Ablk=sparse((horzcat(Ablk,A_session(i).sm))); 
      dj(i)=size(A_session(i).sm,2);
    end

    
    % Obtaining residuals
    
    A=vertcat(Ablk,Hblk);
    oc_observ=sparse(vertcat(oc_observ,och_total));
    
    v=A*x-oc_observ;
    v_scan = v(1:scan(1,iscan).nobs,1);
    v_real=vertcat(v_real,v_scan);
   
    
    
    % For basic outlier test
    if opt.basic_outlier == 1
       AQh_scan=Ablk*sqrtmQxx;
       AQh=vertcat(AQh,AQh_scan);
       clear P N_outlier Qxx_outlier qll_scan AQh_scan 
    end
    
    
    clear A_session A v obs_stat oc_observ  per_stat v_scan  Ablk 

end







    