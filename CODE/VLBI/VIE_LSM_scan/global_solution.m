%
% Script to create the design matrices for the global solution in a
% scanwise way
%
%   10 Jan 2014 by Hana Krasna: AO and APL RgC added
%   21 Oct 2015 by Hana Krasna: seasonal variations in station positions (annual, semi-annual) added




%~~~~~~~~~~~
% Initializing some variables
A_love = [];
A_shida = [];
for i = 1 : 15           
    glob_dj(i) = size(A_session(i).sm,2);
end
A_ra_glob = []; A_de_glob = [];
A_vx=[]; A_vy=[]; A_vz=[];
A_FCN=[];
A_Acr=[]; A_Ace=[]; A_Acn=[]; A_Asr=[]; A_Ase=[]; A_Asn=[]; 
A_vra_glob=[]; A_vde_glob=[];
A_ao=[];
A_rg=[];
%~~~~~~~~~~~


%~~~~~~~~~~~
% DESIGN MATRIX FOR SOURCES

 A_ra_glob = []; A_de_glob = [];
    
 if opt.est_source==1
    
    isou=scan(iscan).iso;
    a(isou).ra = []; a(isou).de = [];
    [Ara,Ade] = a_source_scan(per_source);
    a(isou).ra = Ara; a(isou).de = Ade;
    temp_ra_glob = horzcat(A_ra_glob,a(isou).ra); % design matrix for right ascensions of sources 
    temp_de_glob = horzcat(A_de_glob,a(isou).de); % design matrix for declination of sources
    A_ra_glob=zeros(scan(iscan).nobs,ns);
    A_ra_glob(:,scan(iscan).iso)=temp_ra_glob;   
    glob_dj(length(glob_dj)+1) = size(A_ra_glob,2);
    A_de_glob=zeros(scan(iscan).nobs,ns);
    A_de_glob(:,scan(iscan).iso)=temp_de_glob;
    glob_dj(length(glob_dj)+1) = size(A_de_glob,2);
    clear temp_ra_glob temp_de_glob
   
    if opt.est_source_velo ==1
            refvelsou_mjd= 51544; %2000
            opt.refvelsou_mjd=refvelsou_mjd;
            fctvelsou=(mjd0-refvelsou_mjd)/36525;  % [day/century] 
            A_vra_glob = A_ra_glob*fctvelsou;  
            A_vde_glob = A_de_glob*fctvelsou;     
    end

    
 else
    glob_dj(length(glob_dj)+1) = 0;
    glob_dj(length(glob_dj)+1) = 0;    
 end
 

%~~~~~~~~~~~
% DESIGN MATRIX FOR STATIONS VELOCITIES

    if opt.est_vel == 1
        % +hana 18Nov09
        refvel_mjd=modjuldat(opt.refvel);
        opt.refvel_mjd=refvel_mjd;
        fctvel=(mjd0-refvel_mjd)/36525;  % [day/century] +hana 12Nov2010
        vx = A_session(13).sm; A_vx = vertcat(A_vx,vx*fctvel); glob_dj(length(glob_dj)+1) = size(A_vx,2);
        vy = A_session(14).sm; A_vy = vertcat(A_vy,vy*fctvel); glob_dj(length(glob_dj)+1) = size(A_vy,2);
        vz = A_session(15).sm; A_vz = vertcat(A_vz,vz*fctvel); glob_dj(length(glob_dj)+1) = size(A_vz,2);
        clear vx vy vz
        % -hana
    elseif opt.est_vel == 0
        A_vx = []; A_vy = []; A_vz = [];
        glob_dj(length(glob_dj)+1) = 0;
        glob_dj(length(glob_dj)+1) = 0;
        glob_dj(length(glob_dj)+1) = 0; 
    end
    
    
    % A matrix for Axis Offset
    if opt.est_AO == 1
        [A_ao] = a_AO_scan(per_stat,na);
        glob_dj(length(glob_dj)+1) = size(A_ao,2);
    else
        glob_dj(length(glob_dj)+1) = 0;
    end
    
 %  DESIGN MATRICES FOR AMPLITUDES OF SEASONAL VARIATIONS IN THE STATION
%  POSITIONS


nsewa=0;
if opt.est_stsespos ==1
    for i=1:na
        [AAcr,AAce,AAcn,AAsr,AAse,AAsn] = a_stsespos(per_stat,scan(iscan).nobs);
        % Concatenating - for global adjustment
        A_Acr = horzcat(A_Acr,AAcr);
        A_Ace = horzcat(A_Ace,AAce);
        A_Acn = horzcat(A_Acn,AAcn);
        A_Asr = horzcat(A_Asr,AAsr);
        A_Ase = horzcat(A_Ase,AAse);
        A_Asn = horzcat(A_Asn,AAsn);
        nsewa=size(per_stat(i).pAcr,2);
        clear AAcr AAce AAcn AAse AAsn
    end
end
  
%~~~~~~~~~~~
% amplitudes of seasonal variations in station position
  glob_dj(length(glob_dj)+1) = size(A_Acr,2);
  glob_dj(length(glob_dj)+1) = size(A_Ace,2);
  glob_dj(length(glob_dj)+1) = size(A_Acn,2);
  glob_dj(length(glob_dj)+1) = size(A_Asr,2);
  glob_dj(length(glob_dj)+1) = size(A_Ase,2);
  glob_dj(length(glob_dj)+1) = size(A_Asn,2);    


%~~~~~~~~~~~
% Love number for pole tide
  A_hpole=[];
  if opt.est_hpole==1
        for k = 1 : scan(iscan).nobs
            A_hpole(k,1) = scan(iscan).obs(k).phpole(1); % assigning observationwise partial derivatives
        end
  end
  glob_dj(length(glob_dj)+1) = size(A_hpole,2);


%~~~~~~~~~~~    
% Shida number for pole tide
  A_lpole=[];
  if opt.est_lpole==1
        for k = 1 : scan(iscan).nobs
            A_lpole(k,1) = scan(iscan).obs(k).plpole(1); % assigning observationwise partial derivatives
        end
  end
  glob_dj(length(glob_dj)+1) = size(A_lpole,2);
    
  
 % A matrix for APL Regression Coefficients
    if opt.est_rg == 1
        [A_rg] = a_rg_scan(per_stat,na);
        glob_dj(length(glob_dj)+1) = size(A_rg,2);
    else
        glob_dj(length(glob_dj)+1) = 0;
    end

    

%~~~~~~~~~~~
% DESIGN MATRIX FOR LOVE NUMBERS
if opt.est_love == 1
        [A_love] = a_love([scan(iscan).obs],scan(iscan).nobs);
        glob_dj(length(glob_dj)+1) = size(A_love,2);
else
        glob_dj(length(glob_dj)+1) = 0;
end
    
%~~~~~~~~~~~
% DESIGN MATRIX FOR SHIDA NUMBERS
if opt.est_shida == 1
    	[A_shida] = a_shida([scan(iscan).obs],scan(iscan).nobs);
        glob_dj(length(glob_dj)+1) = size(A_shida,2);
else
        glob_dj(length(glob_dj)+1) = 0;
end


%~~~~~~~~~~~
% DESIGN MATRIX FOR FCN PERIOD
if opt.est_FCNset==1
        for k = 1 :scan(iscan).nobs
            A_FCN(k,1) = scan(iscan).obs(k).pFCN(1); % assigning observationwise partial derivatives
        end
end

glob_dj(length(glob_dj)+1) = size(A_FCN,2);



%~~~~~~~~~ 
% DESIGN MATRIX FOR ACCELERATION OF SSB
A_accSSB=[];
if opt.est_accSSB==1
        for k = 1 : scan(iscan).nobs
            % assigning observationwise partial derivatives
            A_accSSB(k,1:3) = scan(iscan).obs(k).pacc;    % [sec^3/cm]     (don't use: .*(100*c)^3; %sec^3/cm --> cm^2)
        end
end
glob_dj(length(glob_dj)+1) = size(A_accSSB,2);

    
%~~~~~~~~~~~   
% velocities of sources added
  glob_dj(length(glob_dj)+1) = size(A_vra_glob,2);
  glob_dj(length(glob_dj)+1) = size(A_vde_glob,2);
    
%~~~~~~~~~~~
% Gamma parameter
  A_gamma=[];
  if opt.est_gamma==1
        for k = 1 : scan(iscan).nobs
            % assigning observationwise partial derivatives
            A_gamma(k,1) = scan(iscan).obs(k).pGamma *c*100;  % [cm] 
        end
  end
  glob_dj(length(glob_dj)+1) = size(A_gamma,2);


  
  


