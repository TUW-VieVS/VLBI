%*************************************************************************
% vie_lsm_scanwise coded for VieVS by Claudia Tierno Ros
% last Revision:
%   Jan 2014 by Hana, bug corrected for LEVEL2 subdirectory, dirpthL2 added
%   as input parameter
%   10 Jan 2014 by Hana Krasna: AO and APL RgC added
%   19 Jun 2014 by Hana Krasna: sources can be estimated in vie_lsm with
%   NNR condition
%   21 Oct 2015 by Hana Krasna: seasonal variations in station positions (annual, semi-annual) added
%   17 Jul 2016 by A. Hellerschmied: Content of "sources" structure (natural sources, quasars) is now stored in the sub-structure "sources.q"
%   23 Feb 2017 by Hana Krasna: param. addnoise added
%   23 Feb 2017 by A. Hellerschmied: addnoise = opt.add_const_noise_cm
%   01 Mar 2017 by A. Hellerschmied: - Call of "write_opt.m" removed along with removing the "manually find clock break" functions
%                                    - The variable "opt" is not global any more
%                                    - General code revision
%*************************************************************************


function vie_lsm_scanwise(antenna,sources,scan,parameter,dirpth,dirpthL2)
disp('---------------------------------------------------------------')
disp('|                   Welcome to VIE_LSM                        |')
disp('|                    scanwise update                          |')
disp('---------------------------------------------------------------')

% Get the sub-structure with the natural sources!
sources = sources.q;

%..........................................................................
disp(' ')
disp('1. PREPARING DATA');
%..........................................................................


if ~isempty(dirpth) && ~exist(['../DATA/LEVEL3/',dirpth],'dir')
    mkdir(['../DATA/LEVEL3/',dirpth])
end

if ~isempty(dirpthL2) && ~exist(['../DATA/LEVEL2/',dirpthL2],'dir')
    mkdir(['../DATA/LEVEL2/',dirpthL2])
end

%..........................................................................
% SOME NEEDED VALUES
%..........................................................................

c = 299792458; % velocity in m/s
rad2mas = (180/pi)*3600*1000; % radian to milli arc second

% Initializing some variables
N_global=0;
b_global=0;
bsnx=0;
Nsnx=0;
col_est=0;
col_A_add=0;
Nt=0;
bt=0;
Pobserv_all=[];
Hblk=[];
och_total=[];


opt=parameter.lsmopt;
if opt.addSnxSource==1
    opt.est_source=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.residuals=1; % 1 = calculate residuals and make outlier test
% 0 = do not calculate residuals or make outlier test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ntim = length(scan);             % number of scans
na = length(antenna);            % number of antennas
ns = length(sources);            % number of sources
nobserv = sum([scan.nobs]);      % number of observations
mjd1 = min([scan.mjd]);          % time of first scan in mjd
mjd0 = floor(mjd1);              % midnight (beginning of the day of the session)
t=([scan.mjd]-mjd0).*24.*60;     % minutes from midnight to the time of the scan
opt.scans_total = ntim;          % total number of scans
opt.midnight = mjd0+1;           % beginning of the 2nd day of session

fname = antenna(1).session;

fprintf('number of scans:    %5d\n',ntim);
fprintf('number of antennas: %5d\n',na);
fprintf('number of sources:  %5d\n',ns);
fprintf('number of obs.:     %5d\n',sum([scan.nobs]));

mjd_scan=zeros(1,ntim); %preallocating
for i=1: ntim
    mjd_scan(i)=scan(i).mjd; %time of the scans
end


%~~~~~~~~~~~~~~~~~
% ADDING IN "sources" TIME AND NUMBER OF THE SCAN OF THE OBSERVATIONS

if opt.pw_sou==1
    for i=1:ns
        sources(i).iscan=[];
    end
    for i=1: ntim
        % adds in 'sources' the number of the scans were the source is being observed
        sources(scan(i).iso).iscan=horzcat(sources(scan(i).iso).iscan,i);
    end
end



%~~~~~~~~~~~~~~~~~
%ASSIGNING VARIABLES FROM 'scan' TO SOURCE 'obs_per_source'
% obs_per_source=zeros(1,ns); %preallocating
% if  (opt.global_solve==1 && opt.est_source==1) || (opt.ascii_snx==1 && opt.est_source==1) || opt.control_gui_vie_lsm==1 || opt.pw_sou == 1 || opt.est_sourceNNR==1 % +hana 05Oct10
for isou = 1 : ns
    i = 0; k = 0;
    for itim = 1 : ntim % number of scans per session
        for iobs = 1 : scan(itim).nobs % number of observations per scan
            i = i + 1;
            iso = scan(itim).iso;
            if iso == isou
                k = k + 1; % k. observation of the specific source
                obs_per_source(isou).mjd(k) = scan(itim).mjd;
                obs_per_source(isou).iso(k) = scan(itim).iso;
                obs_per_source(isou).i1(k) = scan(itim).obs(iobs).i1;
                obs_per_source(isou).i2(k) = scan(itim).obs(iobs).i2;
            end
        end
    end
end

% end

%~~~~~~~~~~~~~~~~~~
% A PRIORI COORDINATES OF STATIONS (1 value per station for the whole session)
% it takes the coord. form scan of the 1st observation of that station
x0=zeros(1,na); %preallocating
y0=zeros(1,na);
z0=zeros(1,na);
for i=1:na
    x0(i)=0;
    y0(i)=0;
    z0(i)=0;
    
    j=1;
    while x0(i)==0
        temp=(~cellfun(@isempty, {scan(j).stat.temp}));
        if length(temp)>=i && temp(i)==1 %is the station observing in scan j?
            x0(i)=scan(j).stat(i).x(1);
            y0(i)=scan(j).stat(i).x(2);
            z0(i)=scan(j).stat(i).x(3);
            clear temp
            j=1;
            break
        else
            j=j+1;
            clear temp
        end
    end
end

clear iobs iso



%..........................................................................
disp(' ')
disp('2. CREATING DEFAULT OPTIONS');
%..........................................................................
if (opt.global_solve==1 && opt.est_source==1) || (opt.ascii_snx==1 && opt.est_source==1) || opt.control_gui_vie_lsm==1 || opt.pw_sou == 1
    opt = lsmopt_scan_global(antenna,na,obs_per_source,ns,sources, opt);
else
    opt = lsmopt_scan(antenna,na,ns,obs_per_source,sources, opt);
end


for istat = 1 : na
    if opt.stat(istat).ref == 1 %gets number of reference station
        nistat = istat;
        break
    end
end




%..........................................................................
disp(' ')
disp('3. CALCULATING NUMBER OF ESTIMATION INTERVALS');
%..........................................................................
% xpol
t20 = ceil(t(ntim)/opt.xpol.int)*opt.xpol.int; % Upper limit of the last estimation interval
t10= floor(t(1)/opt.xpol.int)*opt.xpol.int;    % Lower limit of the first estimation interval
num_inter_xpol=(t20-t10)/opt.xpol.int;         % Number of estimation intervals

% ypol
t20 = ceil(t(ntim)/opt.ypol.int)*opt.ypol.int; % Upper limit of the last estimation interval
t10= floor(t(1)/opt.ypol.int)*opt.ypol.int;    % Lower limit of the first estimation interval
num_inter_ypol=(t20-t10)/opt.ypol.int;         % Number of estimation intervals

% dUT1
t20 = ceil(t(ntim)/opt.dut1.int)*opt.dut1.int; % Upper limit of the last estimation interval
t10= floor(t(1)/opt.dut1.int)*opt.dut1.int;    % Lower limit of the first estimation interval
num_inter_dut1=(t20-t10)/opt.dut1.int;         % Number of estimation intervals

% nutdx
t20 = ceil(t(ntim)/opt.nutdx.int)*opt.nutdx.int; % Upper limit of the last estimation interval
t10= floor(t(1)/opt.nutdx.int)*opt.nutdx.int;    % Lower limit of the first estimation interval
num_inter_nutdx=(t20-t10)/opt.nutdx.int;         % Number of estimation intervals

% nutdy
t20 = ceil(t(ntim)/opt.nutdy.int)*opt.nutdy.int;   % Upper limit of the last estimation interval
t10= floor(t(1)/opt.nutdy.int)*opt.nutdy.int;      % Lower limit of the first estimation interval
num_inter_nutdy=(t20-t10)/opt.nutdy.int;           % Number of estimation intervals

% zwd, pwl clock offsets, ngr and egr
% preallocating
num_inter_zwd=zeros(1,na);
num_inter_clk=zeros(1,na);
num_inter_egr=zeros(1,na);
num_inter_ngr=zeros(1,na);
t_first=zeros(1,na);
t_last=zeros(1,na);

for i=1:na
    nobs=1;
    nobs2=ntim;
    
    %finds 1st observation of station i
    while sum(i==find(~cellfun(@isempty, {scan(nobs).stat.temp})))==0
        nobs=nobs+1;
    end
    
    %finds last observation of station i
    while sum(i==find(~cellfun(@isempty, {scan(nobs2).stat.temp})))==0
        nobs2=nobs2-1;
    end
    
    % zwd
    t1 = floor(t(nobs)/opt.int_zwd)*opt.int_zwd; % Lower limit of the first estimation interval of station i
    t2 = ceil(t(nobs2)/opt.int_zwd)*opt.int_zwd; % Upper limit of the last estimation interval of station i
    num_inter_zwd(i)=(t2-t1)/opt.int_zwd;        % Number of estimation intervals of station i
    
    % pwl clock offsets
    t1 = floor(t(nobs)/opt.int_clk)*opt.int_clk; % Lower limit of the first estimation interval of station i
    t2 = ceil(t(nobs2)/opt.int_clk)*opt.int_clk; % Upper limit of the last estimation interval of station i
    num_inter_clk(i)=(t2-t1)/opt.int_clk;        % Number of estimation intervals of station i
    T_(i).clk_first_obs = t1;
    
    % egr
    t1 = floor(t(nobs)/opt.int_egr)*opt.int_egr; % Lower limit of the first estimation interval of station i
    t2 = ceil(t(nobs2)/opt.int_egr)*opt.int_egr; % Upper limit of the last estimation interval of station i
    num_inter_egr(i)=(t2-t1)/opt.int_egr;        % Number of estimation intervals of station i
    
    % ngr
    t1 = floor(t(nobs)/opt.int_ngr)*opt.int_ngr; % Lower limit of the first estimation interval of station i
    t2 = ceil(t(nobs2)/opt.int_ngr)*opt.int_ngr; % Upper limit of the last estimation interval of station i
    num_inter_ngr(i)=(t2-t1)/opt.int_ngr;        % Number of estimation intervals of station i
    
    
    t_first(i)=t(nobs); % time of the first observation of station i
    t_last(i)=t(nobs2); % time of the last observation of station i
    
end
clear t1 t2 t20 nobs nobs2


% ra and de source coord.
%preallocating
num_inter_sou=zeros(1,ns);
t_first_sou=zeros(1,ns);
t_last_sou=zeros(1,ns);

if opt.pw_sou == 1
    for i=1:ns
        
        % 1st observation of source i
        nobs=sources(i).iscan(1);
        
        % last observation of source i
        nobs2=sources(i).iscan(length(sources(i).iscan));
        
        % number of estimation intervals
        t1 = floor(t(nobs)/opt.sour_int_rade)*opt.sour_int_rade; % Lower limit of the first estimation interval of source i
        t2 = ceil(t(nobs2)/opt.sour_int_rade)*opt.sour_int_rade; % Upper limit of the last estimation interval of source i
        num_inter_sou(i)=(t2-t1)/opt.sour_int_rade;              % Number of estimation intervals for source i
        
        t_first_sou(i)=t(nobs); % time of the first observation of source i
        t_last_sou(i)=t(nobs2); % time of the last observation of source i
        
        nso(i).sources=num_inter_sou(i)+1;
    end
else
    for i=1:ns
        num_inter_sou(i)=0;
        t_first_sou(i)=0;
        nso(i).sources=num_inter_sou(i)+1;
    end
end
clear t1 t2 t20 nobs nobs2



%..........................................................................
disp(' ')
disp('4. FORMING THE REDUCED OBSERVATION VECTOR "oc_observ"');
%..........................................................................

% Forming the reduced observation vector (oc_observ)
temp = [scan.obs];
oc_observ = sparse(([temp.obs]-[temp.com])'.*c*100); % [cm] observed minus computed vector

addnoise = opt.add_const_noise_cm * 1e-2; % m
fprintf('Constant noise added to input sig.: %4.2f cm\n', opt.add_const_noise_cm);

if opt.first ~= 1
    first_solution.ref_st = [];
elseif opt.first == 1
    [oc_observ_all,first_solution,opt] = reduce_oc_scan(nobserv,na,ntim,scan,oc_observ,opt,c,addnoise);
end

if opt.second == 0
    return
end

clear temp


%..........................................................................
disp(' ')
disp('5. CREATION OF DESIGN MATRIXES');
%..........................................................................

obs_num=1;
for iscan=1:ntim  % loop over all scans
    
    % counter
    if mod(iscan,100)==0
        fprintf(['processing scan ', num2str(iscan),' of ',num2str(ntim),'\n'])
    end
    
    %~~~~~~~~~~~~~~~~~~
    % CREATION OF DESIGN MATRICES
    
    [per_stat,T_]=assign_parameters_scan(scan, na,opt,iscan,antenna,t,mjd0,T_); %assign input parameters for each station
    
    nmi_observ = sqrt((addnoise/c)^2.+([scan(1,iscan).obs.sig]).^2); % [seconds]
    
    Pobserv = diag(sparse(1./((nmi_observ.^2).*c^2*100^2))); % [1/cm^2]
    
    Pobserv_all=sparse(blkdiag(Pobserv_all,Pobserv));
    
    oc_observ=oc_observ_all(obs_num:obs_num+scan(iscan).nobs-1,1);
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
    elseif opt.pw_stc==1 % if station coord. is estimated as pwl offsets
        % Design matrix for antenna dx,dy,dz coordinates
        [Apwx, Apwy, Apwz] = apw_xyz_scan(per_stat,T_,na);
    end
    
    
    
    %~~~~~~~~~
    % DESIGN MATRICES FOR SOURCES COORDINATES
    [per_source,T_source] = sourcewisepar_scan(opt,scan,iscan,sources,t,mjd_scan);
    
    if  (opt.global_solve==1 && opt.est_source==1) || (opt.ascii_snx==1 && opt.est_source==1) || (opt.est_sourceNNR==1) % +hana 18Jun14
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
    
    
    % Xpol (piecewise) - ahp_xpol
    [Apwxpol,T_xpol] = ahp_xpol_scan_A(scan,iscan,opt,c,rad2mas,t);
    
    % Ypol (piecewise) - ahp_ypol
    [Apwypol,T_ypol] = ahp_ypol_scan_A(scan,iscan,c,rad2mas,opt,t);
    
    % dUT1 (piecewise) - ahp_dut1
    [Apwdut1,T_dut1] = ahp_dut1_scan_A(scan,iscan,c,rad2mas,opt,t);
    
    % Celestial Pole Offsets dx (DEPS) (piecewise)
    [Apwnutdx,T_nutdx] = ahp_nutdx_scan_A(scan,iscan,c,rad2mas,opt,t);
    
    % Celestial Pole Offsets dy (DPSI) (piecewise)
    [Apwnutdy,T_nutdy] = ahp_nutdy_scan_A(scan,iscan,c,rad2mas,opt,t);
    
    
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
    elseif opt.pw_stc==1
        A(13).sm= horzcat(Apwx.sm);
        A(14).sm= horzcat(Apwy.sm);
        A(15).sm= horzcat(Apwz.sm);
    end
    
    % Sources
    if opt.pw_sou == 1
        A(11).sm= Apw_ra;
        A(12).sm= Apw_de;
    end
    if opt.est_sourceNNR==1
        A(11).sm= Ara;
        A(12).sm= Ade;
        
        for isou=1:length(opt.source)
            ra(isou)=opt.source(isou).ra;
            de(isou)=opt.source(isou).de;
        end
    end
    
    clear nmi_observ so  Qll  Apwclk Apwdut1 Apwegr Apwngr
    clear Apwnutdx Apwnutdy Apwxpol Apwypol Arqclk Ax Ay Az Azwd
    clear Apw_ra Apw_de Ara Ade
    
    
    
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
    if opt.pw_sou == 1
        [A_session] = delsource_scan_A(opt,A_session,num_inter_sou,ns);
    end
    
    % Deleting reference clock
    [A_session] = delref_scan_A(A_session,num_inter_clk,nistat,opt);
    
    % Excluding models
    [A_session] = delmodel_scan_A(opt,A_session);
    
    
    
    %~~~~~~~~~
    % GLOBAL SOLUTION
    if opt.global_solve==1 || opt.ascii_snx ==1
        clear glob_dj
        global_solution;  % Design matrixes for the global solution
    end
    
    
    %~~~~~~~~~
    % OBTAINING ABLK, N AND b
    
    Ablk=[];
    
    % Concatenating A_session to obtain Ablk
    dj=zeros(1,15); %preallocating
    for i=1:15
        Ablk=sparse(horzcat(Ablk,A_session(i).sm));
        dj(i)=size(A_session(i).sm,2);
    end
    
    % Calculating N and b
    N=Ablk'*Pobserv*Ablk;
    b=Ablk'*Pobserv*oc_observ;
    
    % Stacking matrixes (to obtain them for the whole session)
    Nt=Nt+N;
    bt=bt+b;
    
    col_Ablk=size(Ablk,2);
    
    
    %~~~~~~~~~
    % GLOBAL SOLUTION AND SINEX OUTPUT
    if opt.global_solve==1 || opt.ascii_snx ==1
        A_add = horzcat(A_ra_glob,A_de_glob,A_vx,A_vy,A_vz,A_ao,A_love,A_shida,A_FCN,A_accSSB,A_vra_glob,A_vde_glob,...
            A_Acr,A_Ace,A_Acn,A_Asr,A_Ase,A_Asn,A_hpole,A_lpole,A_gamma,A_rg);
        col_A_add=size(A_add,2);
        A_glob = horzcat(Ablk,A_add);
        N_global_scan = A_glob' * Pobserv * A_glob;
        b_global_scan = A_glob' * Pobserv * oc_observ;
        
        N_global=N_global+N_global_scan;
        b_global=b_global+b_global_scan;
        
        % SINEX output
        if opt.ascii_snx ==1
            [col_est] = snx_split_scan(A_session,opt);
            A_sinex=A_glob;
            
            Nsnx_scan=A_sinex' * Pobserv * A_sinex;
            bsnx_scan=A_sinex' * Pobserv * oc_observ;
            
            Nsnx=Nsnx+Nsnx_scan;
            bsnx= bsnx+ bsnx_scan;
            
        end
        clear   bsnx_scan  Nsnx_scan A_sinex A_glob A_add N_global_scan b_global_scan
        clear  A_FCN  A_de_glob A_ra_glob A_glob A_love A_shida A_v* Ablk
    end
    %~~~~~~~~~
    
    clear H N Ph Pobserv b A_session  obs_stat oc_observ och per_stat
    
    
end % loop over scans
clear t


%..........................................................................
% CREATION OF CONSTRAINTS
%..........................................................................

sum_dj = zeros(1,15);
for i = 1 : 15
    sum_dj(i+1) = sum_dj(i) + dj(i);
end

disp(' ')
disp('6. CREATION OF CONSTRAINTS MATRIXES');

Nh=0;
bh=0;

Ngc=0;
bgc=0;

% H(1), H(2), Ph(1), Ph(2), och(1), och(2) clock
[H1,Ph1,H2,Ph2,och1] = hpoc_clk_scan(num_inter_clk,opt,na,c,col_A2);
[H1,Ph1,och1,H2,Ph2] = delref_scan_H(opt,H1,Ph1,H2,Ph2,och1,num_inter_clk,nistat,na,num_inter_qclk); %excludes reference clock
[H1,Ph1,och1,H2,Ph2] = delmodel_scan_const_1_2(opt,H1,Ph1,och1,H2,Ph2); % excludes model
ncol=0;
[Nh,bh,ncol,Ngc,bgc,Hblk,och_total]=create_const(H1,Ph1,och1,Nh,bh,ncol,col_Ablk,col_A_add,Ngc,bgc,opt,col_est,Hblk,och_total);
ncol=ncol+size(H2,2);
nrow=size(H1,1)+size(H2,1);
num_psob(1)=length(find((any(H1,1)==1)));
num_psob(2)=length(find((any(H2,1)==1)));
clear H1 Ph1 H2 Ph2 och1

% H(3), Ph(3), och(3) zwd
[H3,Ph3,och3] = hpoc_zwd_scan(num_inter_zwd,opt,na,c);
[H3,Ph3,och3]=delparam_scan_const_zwd(opt,na,num_inter_zwd,H3,Ph3,och3); %exclude zwd parameters
[Nh,bh,ncol,Ngc,bgc,Hblk,och_total]=create_const(H3,Ph3,och3,Nh,bh,ncol,col_Ablk,col_A_add,Ngc,bgc,opt,col_est,Hblk,och_total);
nrow=nrow+size(H3,1);
num_psob(3)=length(find((any(H3,1)==1)));
clear Ph3 och3 H3

% H(4), Ph(4), och(4) north gradient
[H4,Ph4,och4] = hpoc_ngr_scan(opt,na,col_A4,num_inter_ngr);
[Nh,bh,ncol,Ngc,bgc,Hblk,och_total]=create_const(H4,Ph4,och4,Nh,bh,ncol,col_Ablk,col_A_add,Ngc,bgc,opt,col_est,Hblk,och_total);
nrow=nrow+size(H4,1);
num_psob(4)=length(find((any(H4,1)==1)));
clear Ph4 och4 H4

% H(5), Ph(5), och(5) east gradient
[H5,Ph5,och5] = hpoc_egr_scan(opt,na,col_A5,num_inter_egr);
[Nh,bh,ncol,Ngc,bgc,Hblk,och_total]=create_const(H5,Ph5,och5,Nh,bh,ncol,col_Ablk,col_A_add,Ngc,bgc,opt,col_est,Hblk,och_total);
nrow=nrow+size(H5,1);
num_psob(5)=length(find((any(H5,1)==1)));
clear Ph5 och5 H5

% H(6), Ph(6), och(6) x pol
[H6,Ph6,och6] = ahp_xpol_scan_H(opt,num_inter_xpol);
[H6,Ph6,och6]=delmodel_scan_const_6(opt,H6,Ph6,och6); % Excludes model
[Nh,bh,ncol,Ngc,bgc,Hblk,och_total]=create_const(H6,Ph6,och6,Nh,bh,ncol,col_Ablk,col_A_add,Ngc,bgc,opt,col_est,Hblk,och_total);
nrow=nrow+size(H6,1);
num_psob(6)=length(find((any(H6,1)==1)));
clear Ph6 och6 H6

% H(7), Ph(7), och(7) y pol
[H7,Ph7,och7] = ahp_ypol_scan_H(opt,num_inter_ypol);
[H7,Ph7,och7] = delmodel_scan_const_7(opt,H7,Ph7,och7); %Eliminates model
[Nh,bh,ncol,Ngc,bgc,Hblk,och_total]=create_const(H7,Ph7,och7,Nh,bh,ncol,col_Ablk,col_A_add,Ngc,bgc,opt,col_est,Hblk,och_total);
nrow=nrow+size(H7,1);
num_psob(7)=length(find((any(H7,1)==1)));
clear Ph7 och7 H7

% H(8), Ph(8), och(8) dUT1
[H8,Ph8,och8] = ahp_dut1_scan_H(opt,num_inter_dut1);
[H8,Ph8,och8]=delmodel_scan_const_8(opt,H8,Ph8,och8); %Eliminates model
[Nh,bh,ncol,Ngc,bgc,Hblk,och_total]=create_const(H8,Ph8,och8,Nh,bh,ncol,col_Ablk,col_A_add,Ngc,bgc,opt,col_est,Hblk,och_total);
nrow=nrow+size(H8,1);
num_psob(8)=length(find((any(H8,1)==1)));
clear Ph8 och8 H8

% H(9), Ph(9), och(9) nutation dx
[H9,Ph9,och9] = ahp_nutdx_scan_H(opt,num_inter_nutdx);
[H9,Ph9,och9]=delmodel_scan_const_9(opt,H9,Ph9,och9); % Eliminates model
[Nh,bh,ncol,Ngc,bgc,Hblk,och_total]=create_const(H9,Ph9,och9,Nh,bh,ncol,col_Ablk,col_A_add,Ngc,bgc,opt,col_est,Hblk,och_total);
nrow=nrow+size(H9,1);
num_psob(9)=length(find((any(H9,1)==1)));
clear Ph9 och9 H9

% H(10), Ph(10), och(10) nutation dy
[H10,Ph10,och10] = ahp_nutdy_scan_H(opt,num_inter_nutdy);
[H10,Ph10,och10]=delmodel_scan_const_10(opt,H10,Ph10,och10); % Eliminates model
[Nh,bh,ncol,Ngc,bgc,Hblk,och_total]=create_const(H10,Ph10,och10,Nh,bh,ncol,col_Ablk,col_A_add,Ngc,bgc,opt,col_est,Hblk,och_total);
nrow=nrow+size(H10,1);
num_psob(10)=length(find((any(H10,1)==1)));
clear Ph10 och10 H10


% H(11), Ph(11), och(11) ra of source  H(12), Ph(12), och(12) de of source
sumso.sources(1) = 0;
for isou = 1 : ns
    sumso.sources(isou+1) = sumso.sources(isou) + nso(isou).sources; % total vector of source coor.                                            % estimates after eleminating non-observed sources
end

if opt.pw_sou == 1  && opt.constr_sou == 1  || opt.est_sourceNNR==1
    [H11,H12,Ph11,Ph12,och11,och12] = hpoc_sources_scan(ns,opt,num_inter_sou); % Constraints for source coordinates
    [H11,Ph11,och11,H12,Ph12,och12 ]=delmodel_scan_const_11(ns,opt,num_inter_sou,sumso,H11,Ph11,och11,H12,Ph12,och12); % Eliminates model
    [Nh,bh,ncol,Ngc,bgc,Hblk,och_total]=create_const(H11,Ph11,och11,Nh,bh,ncol,col_Ablk,col_A_add,Ngc,bgc,opt,col_est,Hblk,och_total);
    nrow=nrow+size(H11,1);
    num_psob(11)=length(find((any(H11,1)==1)));
    clear Ph11 och11 H11
    [Nh,bh,ncol,Ngc,bgc,Hblk,och_total]=create_const(H12,Ph12,och12,Nh,bh,ncol,col_Ablk,col_A_add,Ngc,bgc,opt,col_est,Hblk,och_total);
    nrow=nrow+size(H12,1);
    num_psob(12)=length(find((any(H12,1)==1)));
    clear Ph12 och12 H12
end

% H(13), Ph(13), och(13) x station
[H13,H14,H15,Ph13,Ph14,Ph15,och13,och14,och15] = hpoc_xyz_scan(opt,na);
[H13,H14,H15,Ph13,Ph14,Ph15,och13,och14,och15] =...
    delparam_scan_const_xyz(opt,H13,H14,H15,Ph13,Ph14,Ph15,och13,och14,och15,na); %Eliminates parameters
[Nh,bh,ncol,Ngc,bgc,Hblk,och_total]=create_const(H13,Ph13,och13,Nh,bh,ncol,col_Ablk,col_A_add,Ngc,bgc,opt,col_est,Hblk,och_total);
nrow=nrow+size(H13,1);
num_psob(13)=length(find((any(H13,1)==1)));
clear Ph13 och13 H13

%H(14), Ph(14), och(14) y station
[Nh,bh,ncol,Ngc,bgc,Hblk,och_total]=create_const(H14,Ph14,och14,Nh,bh,ncol,col_Ablk,col_A_add,Ngc,bgc,opt,col_est,Hblk,och_total);
nrow=nrow+size(H14,1);
num_psob(14)=length(find((any(H14,1)==1)));
clear Ph14 och14 H14

% H(15), Ph(15), och(15) z station
[Nh,bh,ncol,Ngc,bgc,Hblk,och_total]=create_const(H15,Ph15,och15,Nh,bh,ncol,col_Ablk,col_A_add,Ngc,bgc,opt,col_est,Hblk,och_total);
nrow=nrow+size(H15,1);
num_psob(15)=length(find((any(H15,1)==1)));
clear  Ph15 och15 H15

clear col_est col_red col_Ablk col_A_add  ncol



%..........................................................................
% LEAST SQUARES ADJUSTMENT
%..........................................................................

ess=opt.est_singleses;
if ess==1 % +hana 10Nov10
    
    disp(' ')
    disp('7. LEAST SQUARES ADJUSTMENT');
    
    %~~~~~~~~~
    % ADJUSTMENT
    N=Nt+Nh;
    
    if opt.stc == 1 % if station coordinates will be estimated
        if opt.nnt_stc == 0
            opt.fixed_station = 'all or some of them';
            disp('!!! some or all station coordinates are not estimated (selected as fixed)!!!');
        elseif opt.nnt_stc == 1
            opt.fixed_station = 'none';
            disp('!!! NNT and NNR conditions are introduced to matrix N for station coordinates!!!');
            [N] = nnt_cond_scan(na,x0,y0,z0,opt,sum_dj,N) ; % Add NNT/NNR
        end
    end
    
    % Introducing NNR for source coordinates
    if opt.est_sourceNNR==1
        disp('!!! NNR condition is introduced to matrix N for source coordinates!!!');
        [N] = nnr_cond(ns,ra,de,opt,sum_dj,N);
    end
    
    Qxx=inv(N);
    Qxx = Qxx(1:sum_dj(length(sum_dj)),1:sum_dj(length(sum_dj)));
    b=bt+bh;
    x=Qxx*b;
    
    
    %~~~~~~~~~
    % ACCURACY CRITERIA
    
    lTPl=oc_observ_all'*Pobserv_all*oc_observ_all;
    vTPv=lTPl-x'*b;
    mo = sqrt(vTPv/(nobserv+nrow-length(x)));  % RMS [2] SINEX V2.01 Appendix(2) Eq.(20)
    opt.mo = mo;
    
    mi = mo*sqrt(diag(Qxx)); % std. dev. of estimated parameters [cm,mas]
    fprintf('chi-squared of main solution: %3.4f\n',mo);
    
    % Calculation of residuals
    if opt.basic_outlier==1
        sqrtmQxx=sqrtm(Qxx);
    end
    
    if opt.residuals==1 % calculate residuals and make outlier test??
        calc_residuals
        
        %~~~~~~~~~
        % OUTLIER TEST (begins here)
        
        % Outlier test begins here
        % DETECTING OUTLIER
        if opt.basic_outlier == 1
            qll = diag(inv(Pobserv_all));
        end
        
        
        if opt.simple_outlier == 1
            qvv = 1;
        end
        
        count_out = 0; out_v = [];
        if opt.basic_outlier == 1 || opt.simple_outlier == 1
            for v_i = 1 : length(v_real)
                if opt.basic_outlier == 1
                    qvv = qll(v_i) - AQh(v_i,:)*AQh(v_i,:)';
                end
                if abs(v_real(v_i)) >= opt.par_outlier*mo*sqrt(qvv)
                    count_out = count_out + 1;
                    out_v(count_out) = v_i; % out_v is the vector that contains all outliers in oc_observ
                end
            end
            disp('----------');
            
            disp('detecting outliers');
            fprintf('num. of outliers : %2d\n ',count_out);
        else
            disp('----------');
            disp('outlier detection test was not applied!')
            disp('----------');
        end
        
        
        % Adressing and bookkeeping the outliers as station1-station2-mjd in ASCII file
        % Added on the 4th of November 2009 by Kamil Teke
        out = []; inum = 0; ico = 0;
        if ~isempty(out_v)
            for isc = 1 : length(scan)
                for iobs = 1 : length(scan(isc).obs)
                    inum = inum + 1;
                    for iout = 1 : size(out_v,2)
                        if inum == out_v(iout)
                            ico = ico + 1;
                            stat1_out(ico).name = antenna(scan(isc).obs(iobs).i1).name;
                            stat2_out(ico).name = antenna(scan(isc).obs(iobs).i2).name;
                            mjd_out(ico) = scan(isc).mjd;
                            out(ico,1) = isc; 
                            out(ico,2) = iobs;
                        end
                    end
                end
            end
            yr=str2double(fname(1:2));
            if ~isempty([parameter.lsmopt.dirout])
                outdir = ['../DATA/OUTLIER/',parameter.lsmopt.dirout,'/'];
            else
                outdir = '../DATA/OUTLIER/';
            end
            if isempty(yr)
                outdir=outdir; %non-standard format of NGS-file-name...
            elseif yr>60
                outdir=[outdir,num2str(1900+yr),'/'];
            else
                outdir=[outdir,num2str(2000+yr),'/'];
            end
            if ~exist(outdir,'dir')
                mkdir(outdir);
            end
            fid1 = fopen([outdir,fname,'.OUT'],'a');
            fprintf('writing outliers to file %s%s.OUT\n',outdir,fname);
            for isc_out = 1 : size(out,1)
                fprintf(fid1,'%s %s %4.12f\n',stat1_out(isc_out).name,stat2_out(isc_out).name,mjd_out(isc_out));
            end
            fclose(fid1);
            fprintf('total %2d outlier observations are found but NOT eleminated\n',count_out);
            disp('----------');
        end

        
        % +++ write residuals/outlier info to new variable +++
        % if the res file already exist, ie if it was written in first solution
        resFilename=['../DATA/LEVEL3/', dirpth, '/res_', fname, '.mat'];
        if exist(resFilename, 'file')
            load(resFilename);
        else
            % else: first solution was not run - 'res' variable does not exist
            res.allStatNames={antenna.name};
            res.allSourceNames={sources.name};
            
            lengthOfScans=cellfun(@length, {scan.obs});
            mjdOfObs=zeros(sum(lengthOfScans), 1);
            res.baselineOfObs=zeros(sum(lengthOfScans),2);
            res.source=zeros(sum(lengthOfScans),1);
            
            runningInd=1;
            for iScan=1:size(scan,2)
                % source index of current scan
                iSo=scan(iScan).iso;
                res.source(runningInd:runningInd+lengthOfScans(iScan)-1)=...
                    repmat(iSo, lengthOfScans(iScan),1);
                
                % get station names for each observation of this scan
                res.baselineOfObs(runningInd:runningInd+lengthOfScans(iScan)-1, :)=...
                    [[scan(iScan).obs.i1]', [scan(iScan).obs.i2]'];
                for iObs=1:lengthOfScans(iScan)
                    mjdOfObs(runningInd,1)=scan(iScan).mjd;
                    runningInd=runningInd+1;
                end
            end
            
            % save mjd/baselines if this was not done in first solution
            res.mjd=          mjdOfObs;
        end
        
        % main solution residuals and outliers should be saved in any case
        res.mainVal=v_real;
        res.outlier=out_v;
        
        % save file (again)
        save(resFilename, 'res');
        
        % --- write residuals/outlier info to new variable ---
    end
    
    %~~~~~~~~~
    % DIVIDING THE x VECTOR IN TO SUB-VECTORS
    
    opt.total_est = sum_dj;
    
    fprintf('total number of estimated parameters: %4d\n',sum_dj(end));
    fprintf('total clock offsets: %4d\n',dj(1));
    fprintf('total rate and quad. terms of clock funct.: %4d\n',dj(2));
    fprintf('total zenith wet delay offsets: %4d\n',dj(3));
    fprintf('total tropo. north gradients: %4d\n',dj(4));
    fprintf('total tropo. east gradients: %4d\n',dj(5));
    fprintf('total pole coor. (x-pol) offsets: %4d\n',dj(6));
    fprintf('total pole coor. (y-pol) offsets: %4d\n',dj(7));
    fprintf('total dUT1 offsets: %4d\n',dj(8));
    fprintf('total celestial pole (nutation dx) offsets: %4d\n',dj(9));
    fprintf('total celestial pole (nutation dy) offsets: %4d\n',dj(10));
    fprintf('total right ascension offsets of sources : %4d\n',dj(11));
    fprintf('total declination offsets of sources : %4d\n',dj(12));
    fprintf('antenna coor. dx offsets: %4d\n',dj(13));
    fprintf('antenna coor. dy offsets: %4d\n',dj(14));
    fprintf('antenna coor. dz offsets: %4d\n',dj(15));
    
    [n_]=nscan2lsm(num_inter_clk,num_inter_egr,num_inter_ngr,num_inter_zwd,na,nistat,num_inter_qclk,opt);
    mjd1=mjd_scan(1);
    
    [t, T, tso]=tTscan2lsm(opt, t_first,t_last,na,nistat,num_inter_xpol,num_inter_ypol,num_inter_dut1,num_inter_nutdx,num_inter_nutdy,num_inter_sou,ns);
    
    [x_] = splitx(x,first_solution,mi,na,sum_dj,n_,mjd0,mjd1,t,T,opt,antenna,ns,nso,tso,ess);
    
    [atpa_.mat] = N;
    [atpl_.vec] = b;
    [opt_] = opt;
    
    
    %~~~~~~~~~
    % SAVE FILES
    
    if ~isempty(dirpth) && ~exist(['../DATA/LEVEL3/',dirpth],'dir')
        mkdir(['../DATA/LEVEL3/',dirpth])
    end
    
    disp('----------');
    fprintf('estimated parameters are saved as ../VieVS/DATA/LEVEL3/%s/x_%s.mat\n',dirpth,fname);
    save(['../DATA/LEVEL3/',dirpth,'/x_',fname,'.mat'],'x_');
    
    fprintf('estimation options are saved as ../VieVS/DATA/LEVEL3/%s/opt_%s.mat\n',dirpth,fname);
    save(['../DATA/LEVEL3/',dirpth,'/opt_',fname,'.mat'],'opt_');
    
    fprintf('normal equation matrix is saved as ../VieVS/DATA/LEVEL3/%s/atpa_%s.mat\n',dirpth,fname);
    save(['../DATA/LEVEL3/',dirpth,'/atpa_',fname,'.mat'],'atpa_');
    
    fprintf('right hand side vector is saved as ../VieVS/DATA/LEVEL3/%s/atpl_%s.mat\n',dirpth,fname);
    save(['../DATA/LEVEL3/',dirpth,'/atpl_',fname,'.mat'],'atpl_');
    
    
end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% -------------------------------------------------------------------------
% +Boehm 2 July 2009
% FOR GLOBAL SOLUTION (CLOCK IS FIXED (or NNT) and WITH CONSTRAINTS BETWEEN OFFSETS)

% +hana 05Oct10
% and for N,b for SINEX output (because of source coordinates)

% +hana 10Nov10
% In case, that the single-session solution wasn't applied,
% in the x_ variable will be written only information about
% columns in the N and b and the time information

clear Hblk N Nh Nt Pobserv_all Qxx T_* bh bt num_inter* oc_observ_all och_total

if opt.global_solve == 1 || opt.ascii_snx ==1 % +hana 05Oct10
    disp('----------');
    disp('Preparing Data for Global Solution');
    
    [x_] = splitx(x,first_solution,mi,na,sum_dj,n_,mjd0,mjd1,t,T,opt,antenna,ns,nso,tso,ess);
    
    x_.source='';
    if opt.est_source == 1
        % +hana 21 Jun 2010
        % Assign only the name of the sources to x_.source.name
        for isou = 1 : length(opt.source)
            x_.source(isou).name = opt.source(isou).name;
            x_.source(isou).IERSname = opt.source(isou).IERSname; % hana 24Mar13
        end
        % -hana
    end
    
    
    
    sum_glob_dj(1) = 0;
    for i = 1 : length(glob_dj)
        sum_glob_dj(i+1) = sum_glob_dj(i) + glob_dj(i);
    end
    N_global=N_global+Ngc;
    b_global=b_global+bgc;
    
    IDglobdj=16;
    
    x_.col_soura = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_soude = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_vx    = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_vy    = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_vz    = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_AO    = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    
    nAcr=sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    nAce=sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    nAcn=sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    nAsr=sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    nAse=sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    nAsn=sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    
    % number of seasonal waves (station position variation) : nsewa
    for i=1:na % number of antennas
        x_.col_ssp_Acr(i).col=nAcr(nsewa*i-nsewa+1:nsewa*i);
        x_.col_ssp_Ace(i).col=nAce(nsewa*i-nsewa+1:nsewa*i);
        x_.col_ssp_Acn(i).col=nAcn(nsewa*i-nsewa+1:nsewa*i);
        x_.col_ssp_Asr(i).col=nAsr(nsewa*i-nsewa+1:nsewa*i);
        x_.col_ssp_Ase(i).col=nAse(nsewa*i-nsewa+1:nsewa*i);
        x_.col_ssp_Asn(i).col=nAsn(nsewa*i-nsewa+1:nsewa*i);
    end
    
    x_.col_hpole   = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_lpole   = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    
    x_.col_rg = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    
    x_.col_love  = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_shida = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    
    x_.col_FCNset= sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_accSSB= sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_souvra= sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    x_.col_souvde= sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    
    x_.col_gamma = sum_glob_dj(IDglobdj) + 1 : 1 : sum_glob_dj(IDglobdj+1);  IDglobdj=IDglobdj+1;
    
    
    
    allant= [antenna.name];
    for i=1:size(allant,2)/8
        ant(i,1:8)= allant((i*8)-7:i*8);
    end
    an_glob.name = ant;
    an_glob.in_trf = [antenna.in_trf];
    an_glob.x=[antenna.x];
    an_glob.y=[antenna.y];
    an_glob.z=[antenna.z];
    an_glob.vx=[antenna.vx];
    an_glob.vy=[antenna.vy];
    an_glob.vz=[antenna.vz];
    an_glob.epoch=[antenna.epoch];
    an_glob.start=[antenna.start];
    an_glob.end=[antenna.end];
    an_glob.firstscan_mjd=mjd1;
    % a priori APL Reg Coeff
    for i=1:length(antenna)
        an_glob.aplrg(i,:) = antenna(i).crg; % [ref.pres, a priori RgC]
    end
    
    fprintf('Data for GLOBAL SOLUTION is saved as ../VieVS/DATA/LEVEL2/%s/%s_an_glob.mat\n',dirpth,fname);
    fprintf('Data for GLOBAL SOLUTION is saved as ../VieVS/DATA/LEVEL2/%s/%s_par_glob.mat\n',dirpth,fname);
    fprintf('Data for GLOBAL SOLUTION is saved as ../VieVS/DATA/LEVEL2/%s/%s_Nb_glob.mat\n',dirpth,fname);
    
    opt_.trf=parameter.vie_init.trf;
    opt_.crf=parameter.vie_init.crf;
    
    opt_.total_obs = nobserv;
    opt_.nconstr = nrow;
    opt_.lTPl = lTPl; %l'Pl for global solution
    
    glob1.an = an_glob;
    glob2.x = x_;
    glob2.opt = opt_;
    glob3.N = N_global;
    glob3.b = b_global;
    
    
    if ~isempty(dirpthL2) && ~exist(['../DATA/LEVEL2/',dirpthL2],'dir')
        mkdir(['../DATA/LEVEL2/',dirpthL2])
    end
    
    
    save(['../DATA/LEVEL2/',dirpthL2,'/',fname,'_an_glob.mat'],'glob1');
    save(['../DATA/LEVEL2/',dirpthL2,'/',fname,'_par_glob.mat'],'glob2');
    save(['../DATA/LEVEL2/',dirpthL2,'/',fname,'_Nb_glob.mat'],'glob3');
    
    
    
end


% -------------------------------------------------------------------------
% +hana 16 May 2011
% N and b FOR SINEX OUTPUT from the global matrix (CLOCK IS FIXED (or NNT) and WITH CONSTRAINTS BETWEEN OFFSETS)
% Reduce: always clock
%         optional zwd, trop. gradients
% Keep in sinex file: station coordinates, EOP
% Source coordinates can be writen into sinex file in the N-matrix, bur
% cannot be estimated from the single-session solution

% Constraints of parameters which will be written into SINEX file are
% removed!

if opt.ascii_snx == 1
    
    [col_red col_est] = snx_split(x_,opt.outsnx);
    
    %     addpath ../VIE_MOD_V21
    
    Nsnx=Nsnx+Ngc;
    
    % REDUCTION
    N11=Nsnx(col_est,col_est);
    N12=Nsnx(col_est,col_red);
    N21=Nsnx(col_red,col_est);
    N22=Nsnx(col_red,col_red);
    
    b1=bsnx(col_est);
    b2=bsnx(col_red);
    
    N_sinex = N11 - N12*inv(N22)*N21;
    b_sinex = b1 - N12*inv(N22)*b2;
    
    
    % Reduction of lTPl:
    %lTPl = oc_observ_real'*pobserv*oc_observ_real;
    lTPlreduc = lTPl - (b2'*inv(N22)*b2);
    
    % number of pseudo-observations which were reduced
    % they are added to the number of real observation and the sum is
    % written into the sinex file
    nconstr_red=0;
    if opt.outsnx.clk==0; nconstr_red = nconstr_red + num_psob(1)+num_psob(2); end
    if opt.outsnx.zwd==0; nconstr_red = nconstr_red + num_psob(3); end
    if opt.outsnx.tgr==0; nconstr_red = nconstr_red + num_psob(4)+num_psob(5); end
    if opt.outsnx.xyz==0; nconstr_red = nconstr_red + num_psob(13)+num_psob(14)+num_psob(15); end
    if opt.outsnx.eop==0; nconstr_red = nconstr_red + sum(num_psob(6:10)); end
    
    
    total_est = sum_dj;
    real_obs = nobserv;
    all_obs = real_obs + nconstr_red;
    
    
    % Save info about columns
    col_sinex=snx_newcol(col_est,x_,antenna,opt.outsnx);
    % Save info about statistic
    col_sinex.lTPlreduc = lTPlreduc;
    col_sinex.nr_unknowns = total_est(end);
    col_sinex.nr_obs = all_obs;
    col_sinex.vTPv = vTPv;
    col_sinex.varfac = mo;
    col_sinex.outsnx=opt.outsnx;
    
    % Change units of the N matrix and b vector !!!
    [N_sinex, b_sinex]=snx_changeunits(N_sinex,b_sinex,col_sinex,opt.outsnx);
    
    
    
    
    if exist(['../DATA/LEVEL3/',dirpth,'/SINEX/'])~=7
        mkdir(['../DATA/LEVEL3/',dirpth,'/SINEX/'])
    end
    
    fprintf('\nReduced N and b for SINEX output are saved in ../VieVS/DATA/LEVEL3/%s/SINEX/ \n',dirpth);
    save(['../DATA/LEVEL3/',dirpth,'/SINEX/N_sinex_',fname,'.mat'],'N_sinex');
    save(['../DATA/LEVEL3/',dirpth,'/SINEX/b_sinex_',fname,'.mat'],'b_sinex');
    save(['../DATA/LEVEL3/',dirpth,'/SINEX/col_sinex_',fname,'.mat'],'col_sinex');
    
    % create an ascii sinex file in DATA/SNX
    fprintf('\nWriting SINEX file ... \n');
    write_sinex_vievs(fname, [dirpth '/']);
    fprintf('\nSINEX file is saved in ../VieVS/DATA/SNX/%s.SNX \n\n',fname);
    
end


% -------------------------------------------------------------------------
disp('----------');
disp('8. vie_lsm IS COMPLETED!');
disp('----------');
% -------------------------------------------------------------------------

