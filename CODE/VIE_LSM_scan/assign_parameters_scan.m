%************************************************************************** 
%   Description:
%   Function to assign the input parameters to each station for one scan,
%   in vie_lsm it is done in the main program.
%   It also determines the estimation intervals (done in stwisepar in vie_lsm)
%
%   Input:
%       'scan'       structure array     (for info. /DOC/scan.doc)
%       'na'         (1,1)                Number of antennas
%       'opt'        structure array     (for info. /DOC/opt.doc)
%       'iscan'      (1,1)                Number of the scan that is being analysed
%       'antenna'    structure array     (for info. /DOC/antenna.doc)
%       't'          (1,num. of scans)    minutes from midnight to the time of the scans
%       'mjd0'       (1,1)                midnight (beginning of the day of the session)
%       'T_'         structure array      estimation intervals of clk, zwd, ngr, egr, xyz
%
%   Output:
%       'per_stat'   structure array     (for info. /DOC/per_stat.doc)
%       'T_'         structure array      estimation intervals of clk, zwd, ngr, egr, xyz
%
%   External calls: 	
%   
%   Coded for VieVS: 
%   4 July 2012 by Claudia Tierno Ros
%
%   Revision: 
%   10 Jan 2014 by Hana Krasna: AO and APL RgC added
%   21 Oct 2015 by Hana Krasna: seasonal variations in station positions (annual, semi-annual) added
%**************************************************************************


function [per_stat,T_]=assign_parameters_scan(scan, na,opt,iscan,antenna,t,mjd0,T_)



for istat = 1:na % number of stations
    
    i = 0; k = 0;
   
    nobs=scan(iscan).nobs;
    per_stat(istat).mjd(1,1:nobs) = 0;
    per_stat(istat).oc_nob(1,1:nobs) = 0;
    per_stat(istat).zd(1,1:nobs) =0;
    per_stat(istat).first(1,1:nobs) = 0;
    per_stat(istat).other(1,1:nobs) = 0;
        
    
    per_stat(istat).int_clk=opt.stat(istat).int_clk;
    per_stat(istat).int_zwd=opt.stat(istat).int_zwd;
    per_stat(istat).int_ngr=opt.stat(istat).int_ngr;
    per_stat(istat).int_egr=opt.stat(istat).int_egr;
    per_stat(istat).int_xyz=opt.stat(istat).int_xyz;
             
    per_stat(istat).name = antenna(istat).name; % The name of the station

    per_stat(istat).no_int_zwd=1; %always 1 because there is one scan 
    per_stat(istat).no_int_clk=1;
    per_stat(istat).no_int_egr=1;
    per_stat(istat).no_int_ngr=1;
    per_stat(istat).no_int_xyz=1;

    if length(scan(iscan).stat)>=istat && ~isempty(scan(iscan).stat(istat).az)
        per_stat(istat).mf=scan(iscan).stat(istat).mfw;
        per_stat(istat).az=scan(iscan).stat(istat).az;
        per_stat(istat).zd=scan(iscan).stat(istat).zd;
        per_stat(istat).x0=scan(iscan).stat(istat).x(1);
        per_stat(istat).y0=scan(iscan).stat(istat).x(2);
        per_stat(istat).z0=scan(iscan).stat(istat).x(3);
     else
        per_stat(istat).mf=0;
        per_stat(istat).az=0;
        per_stat(istat).zd=0;
        per_stat(istat).x0=0;
        per_stat(istat).y0=0;
        per_stat(istat).z0=0;
    end
        
     % preallocating
       per_stat(istat).nob=[1:nobs]; 
       per_stat(istat).total=nobs;
       
       per_stat(istat).mf(1:nobs) = 0;
       per_stat(istat).zd(1:nobs) = 0;
       per_stat(istat).minute=(scan(iscan).mjd-mjd0)*24*60;
     
       
        per_stat(istat).dAO(1,1:nobs) = 0;

        per_stat(istat).pAcr(1,1:nobs)=0;
        per_stat(istat).pAce(1,1:nobs)=0;
        per_stat(istat).pAcn(1,1:nobs)=0;
        per_stat(istat).pAsr(1,1:nobs)=0;
        per_stat(istat).pAse(1,1:nobs)=0;
        per_stat(istat).pAsn(1,1:nobs)=0;
        
        per_stat(istat).drg(1,1:nobs) = 0; 
       
       per_stat(istat).dx(1:nobs) = 0; 
       per_stat(istat).dy(1:nobs) = 0;
       per_stat(istat).dz(1:nobs) = 0;
       
      for iobs = 1:nobs % number of observations per scan 
         i = i + 1; % i : i. observation in the session
         i1 = scan(iscan).obs(iobs).i1;
         i2 = scan(iscan).obs(iobs).i2;
          
       
         if i1 == istat
            k = k + 1; % k : k. observation of the specific station  
            per_stat(istat).mjd(iobs) = scan(iscan).mjd; % The times of scans [day]
            per_stat(istat).oc_nob(iobs) = i;
            per_stat(istat).zd(iobs) = scan(iscan).stat(istat).zd;  % Boehm 21 Aug 2009, 15 Sep 2010
            per_stat(istat).mf(iobs) = scan(iscan).stat(istat).mfw;
             
            per_stat(istat).first(iobs) = -1;
            per_stat(istat).other(iobs) = i2;  % Boehm 21 Aug 2009
            
            per_stat(istat).dAO(iobs) = -scan(iscan).obs(iobs).pAO_st1; % axis offset           
            
            per_stat(istat).pAcr(iobs)=-scan(iscan).obs(iobs).pAcr_st1;
            per_stat(istat).pAce(iobs)=-scan(iscan).obs(iobs).pAce_st1;
            per_stat(istat).pAcn(iobs)=-scan(iscan).obs(iobs).pAcn_st1;
            per_stat(istat).pAsr(iobs)=-scan(iscan).obs(iobs).pAsr_st1;
            per_stat(istat).pAse(iobs)=-scan(iscan).obs(iobs).pAse_st1;
            per_stat(istat).pAsn(iobs)=-scan(iscan).obs(iobs).pAsn_st1;
            
            per_stat(istat).drg(iobs) = -scan(iscan).obs(iobs).prg_st1; % APL Regress. Coef.
            

         
         elseif i2 == istat
            k = k + 1; % k : k. observation of the specific station  
            per_stat(istat).mjd(iobs) = scan(iscan).mjd; % The times of scans [day]
            per_stat(istat).oc_nob(iobs) = i;
            per_stat(istat).zd(iobs) = scan(iscan).stat(istat).zd;  % Boehm 21 Aug 2009, 15 Sep 2010
            per_stat(istat).mf(iobs) = scan(iscan).stat(istat).mfw;
             
            per_stat(istat).first(iobs) = +1;
            per_stat(istat).other(iobs) = i1;  % Boehm 21 Aug 2009
            
            per_stat(istat).dAO(iobs) = scan(iscan).obs(iobs).pAO_st2; % axis offset  
            
            per_stat(istat).pAcr(iobs)=scan(iscan).obs(iobs).pAcr_st2;
            per_stat(istat).pAce(iobs)=scan(iscan).obs(iobs).pAce_st2;
            per_stat(istat).pAcn(iobs)=scan(iscan).obs(iobs).pAcn_st2;
            per_stat(istat).pAsr(iobs)=scan(iscan).obs(iobs).pAsr_st2;
            per_stat(istat).pAse(iobs)=scan(iscan).obs(iobs).pAse_st2;
            per_stat(istat).pAsn(iobs)=scan(iscan).obs(iobs).pAsn_st2;
            
            per_stat(istat).drg(iobs) = scan(iscan).obs(iobs).prg_st2; % APL Regress. Coef.
            
         end
                
         if per_stat(istat).mjd(iobs)~=0
            per_stat(istat).dx(iobs) = -scan(iscan).obs(iobs).pstat1(1); % The partial derivatives of delay wrt antenna coordinates
            per_stat(istat).dy(iobs) = -scan(iscan).obs(iobs).pstat1(2);
            per_stat(istat).dz(iobs) = -scan(iscan).obs(iobs).pstat1(3);        
         end
     end
     
     
  % Estimation intervals   
    if sum(per_stat(istat).mjd)~=0 %if the station is observing in that scan
         
            % determining the ZWD estimation intervals and puting observation times in to them 
            t1 = floor(t(iscan)/per_stat(istat).int_zwd)*per_stat(istat).int_zwd;
            t2=t1+per_stat(istat).int_zwd; 

            T_(istat).zwd = [t1,t2]; % The estimation intervals for zwd

            % determining the CLOCK estimation intervals and puting observation times in to them 
            t1 = floor(t(iscan)/per_stat(istat).int_clk)*per_stat(istat).int_clk;
            t2=t1+per_stat(istat).int_clk; 
            
            T_(istat).clk = [t1,t2]; % The estimation intervals for clk

            % determining the EAST GRADIENTS estimation intervals and puting observation times in to them 
            t1 = floor(t(iscan)/per_stat(istat).int_egr)*per_stat(istat).int_egr;
            t2=t1+per_stat(istat).int_egr; 

            T_(istat).egr = [t1,t2]; % The estimation intervals for east gradients

            % determining the NORTH GRADIENTS estimation intervals and puting observation times in to them 
            t1 = floor(t(iscan)/per_stat(istat).int_ngr)*per_stat(istat).int_ngr;
            t2=t1+per_stat(istat).int_ngr; 

            T_(istat).ngr = [t1,t2]; % The estimation intervals for north gradients

            % determining the COORDINATES estimation intervals and puting observation times in to them 
            t1 = floor(t(iscan)/per_stat(istat).int_xyz)*per_stat(istat).int_xyz;
            t2=t1+per_stat(istat).int_xyz; 

            T_(istat).xyz = [t1,t2]; % The estimation intervals for coordinates

        else %if the station is not observing in that scan
            T_(istat).zwd=[0,0];
            T_(istat).clk=[0,0];
            T_(istat).egr=[0,0];
            T_(istat).ngr=[0,0];
            T_(istat).xyz=[0,0];
            
    end


     
end
