% ************************************************************************
%   Description:
%   In vie_glob the global adjustment of the VLBI session is done. In this
%   version you can estimate station coordinates and velocities, source
%   coordinates and EOP.
%
%   References:
%
%   Input:
%      LEVEL2 data created with vie_lsm
%       - sessionname_an_glob.mat
%       - sessionname_Nb_glob.mat
%       - sessionname_par_glob.mat
%
%   Output:
%      the output data is stored
%       - in the structure array: globsol_*.mat
%       - in the txt file: glob_results_*.txt
%       in VieVS/OUT/GLOB/*/
%
%   External calls:
%
%
%   Coded for VieVS:
%   05 Aug 2010 by Hana Spicakova
%
%   Revision:
%   07 Sep 2010 by Hana Spicakova:  an apriori catalogue also in .txt-format
%                                   can be read in (coded by T.Nilsson)
%   08 Sep 2010 by Sun Jing: modification to make the program compatible
%                                   with Linux system
%   04 Oct 2010 by Hana Spicakova:  datum definition for source coordinates - choice between NNR+dz (3 rotations and
%                                   translation in declination) and only NNR
%   12 Oct 2010 by Hana Spicakova:  datum definition for station
%                                   coordinates - choice between NNT/NNR and NNT/NNR/NNS
%   14 Oct 2010 by Hana Spicakova:  bug in the computation of standard
%                                   deviations was fixed
%   15 Dec 2010 by Hana Spicakova:  empty spaces in antbr fixed (it concerned stations DSS13 and GOLDMARS)
%   15 Dec 2010 by Hana Spicakova:  names of antennas are written into globsol.antenna.an_order,
%                                   even if the station coordinates were not estimated
%   12 Jan 2011 by Hana Spicakova:  a bug corrected which occured if a
%                                   break in external file is before beginning of the used catalogue
%   July 2011   by Hana Spicakova:  major changes
%               Major parts of VIE_GLOB were newly written to allow the session-wise reduction for station,
%               which observation period was not sufficient to estimate reliable coordinates and velocities.
%               Also sources can be session-wise reduced, if they show a significant
%               positional instability in either RA or De.
%               If you are interested in reduced parameters, i.e., these
%               stations and sources coordinates, but also EOP, zwd or troposheric
%               gradients, a function "backward_solution.m" was
%               added to VieVS/OUT/GLOB/. For using it you have to store
%               the session-wise N matrices and b vectors. This can be
%               chosen in the 1st GUI of vie_glob.
%  16 Sep 2011  by Hana Spicakova:  bug fixed for the case, if the single session solution in
%               vie_lsm is not done
%  16 Sep 2011  input for function refname_ant corrected
%  16 Sep 2011  by Hana Spicakova: bug by estimation of EOP with different
%               resolution was fixed
%  10 May 2012  by Hana Kr�sn�: changed to be compatibel with the GUI 2.0
%  20 Jun 2012  by Hana Kr�sn�: 
%               added: Love and Shida numbers, FCN per from SET, velocity of sources, acceleration of SSB, seasonal
%               variations in station positions (annual, semi-annual), pole
%               tide Love and Shida numbers
%  07 Aug 2012  by Hana Kr�sn�: compatible with superstation file
%  25 Sep 2012  by Hana Kr�sn�: relativistic parameter gamma 
%  23 Oct 2012  by Hana Krasna: TRF and CRF catalogue are stored in an
%               official format
%  23 Oct 2012  by Hana Krasna: The N and b, and covariance matrix of the global solution can be
%               saved in SINEX files (currently it is working only for
%               the simulateous estimation of TRF (coor+vel) and CRF (coor))
%  24 Oct 2012 by Hana Krasna: the units of the RA estimates change from
%               mas to ms. Change only in globsol.source and in the
%               .txt output
%  20 Dec 2012 by Hana Krasna: correction of station epochs if velocity
%               fixed
%  24 Mar 2013 by Hana Krasna: changed to IERS name of sources and changes
%              according to supersouce file
%  21 Jun 2013 by Hana Krasna: it is possible to distinquish between
%              intervals for velocity constraints at one station with breaks
%  19 Sep 2013 by Hana Krasna: small bug corrected which caused the program
%              to crash if TRF was fixed
%  04 Oct 2013 by Hana Krasna: estimation of antenna axis offset added
%  06 Dec 2013 by Hana Krasna: APL regression coefficient added
%  24 Aug 2015 by Sigrid Boehm: tidal ERP coefficients added
%  02 Oct 2015 by A. Hellerschmied: Error exception added for plotting the station network (It crashes, if certain MATLAB "mapping" toolboxes are not installed)
%  13 Oct 2015 by D. Mayer: Bug fix and minor changes
%  21 Oct 2015 by Hana Krasna: seasonal variations in station positions (annual, semi-annual), pole
%               tide Love and Shida numbers, APL regression coefficient
%               added
%  14 Dec 2015 by S. B�hm: Bug-fix: If sessions with too high RMS were found, the program stopped
%  21 Dec 2015 by Hana Krasna: Error at call of function "par_newcol" fixed (wrong input arguments)
%  17 Jan 2017 by Hana Krasna: Seasonal harmonic signal and pole tide Lova and Shida numbers added to the GUI
%  19 Feb 2017 by Hana Krasna: APL regression coeff. added to the GUI
%  21 Feb 2017 by David Mayer: Added possibility to remove RMS check
%  21 Feb 2017 by David Mayer: Added possibility to fix EOP for specified sessions
%  21 Feb 2017 by David Mayer: Added possibility to save intermediate results (workspace)
%  24 Mar 2017 by Hana Krasna: minor bug - if bred2 is empty it is set to zero: bred =[] --> bred=0;
%  10 Apr 2017 by David Mayer: Added IVS source names 
%  21 Jun 2017 by David Mayer: Adjusted saving of globsol for really big global solutions
%  21 Jun 2017 by David Mayer: The fixing EOP option now ignores the version number
%  08 Sep 2017 by David Mayer: Source name bug fix
%  22 Jan 2018 by Hana Krasna: Backward solution is called by vie_glob
% ************************************************************************

function vie_glob
% clear all
% close all
% clc
save_intermediate_results_flag = 0;
path_level='../';

% % Read which parametres are to be estimated ('paramGS.m')
% guiglob
% uiwait(guiglob)

load([path_level 'DATA/GLOB/parGS'],'parGS');
load([path_level 'DATA/GLOB/pathGS'],'pathGS');

[sizeparGS] = size(parGS);
IDglob = sizeparGS(2)+1;

% add parameters which are not in GUI yet
parGS(IDglob).name = 'love'; parGS(IDglob).id=0; IDglob=IDglob+1;
parGS(IDglob).name = 'shida'; parGS(IDglob).id=0; IDglob=IDglob+1;
parGS(IDglob).name = 'FCNset'; parGS(IDglob).id=0; IDglob=IDglob+1;
parGS(IDglob).name = 'accSSB'; parGS(IDglob).id=0; IDglob=IDglob+1;

parGS(IDglob).name = 'svra'; parGS(IDglob).id=0; IDglob=IDglob+1;
parGS(IDglob).name = 'svde'; parGS(IDglob).id=parGS(IDglob-1).id; IDglob=IDglob+1;
parGS(IDglob).name = 'gamma'; parGS(IDglob).id=0; IDglob=IDglob+1;

parGS(IDglob).name = 'bdco'; parGS(IDglob).id=2; IDglob=IDglob+1; % reduce (2) or delete (0) bas-dep clock offset

[g] = globind(parGS);


if parGS(g.g_coord(1)).id==0 && parGS(g.g_vel(1)).id==1
    parGS(g.g_vel(1)).id=0;
    fprintf('\n Station velocities cannot be estimated without station coordinates. Station velocities will be fixed.\n')
end

if parGS(g.g_coord(1)).id==0 && parGS(g.g_stseaspos_Ar(1)).id==1 ||...
        parGS(g.g_coord(1)).id==0 && parGS(g.g_stseaspos_Ae(1)).id==1 ||...
        parGS(g.g_coord(1)).id==0 && parGS(g.g_stseaspos_An(1)).id==1
    parGS(g.g_stseaspos_Ar(1)).id=0;
    parGS(g.g_stseaspos_Ae(1)).id=0;
    parGS(g.g_stseaspos_An(1)).id=0;
    fprintf('\n Seasonal variations of station position cannot be estimated if the stations are fixed.\n')
end


if parGS(g.g_srade(1)).id==0 && parGS(g.g_svrade(1)).id==1
    parGS(g.g_svrade(1)).id=0;
    fprintf('\n Source velocities cannot be estimated without source coordinates. Source velocities will be fixed.\n')
end



save([path_level 'DATA/GLOB/parGS'],'parGS');


if pathGS.bckwrdsol==1
    if exist([path_level 'OUT/GLOB/BACKWARD_DATA/' pathGS.out ])~=7
        mkdir([path_level 'OUT/GLOB/BACKWARD_DATA/' pathGS.out ])
    end
end

% Seasonal variations in station positions
% sin and cos amplitudes are always estimated together
% for all stations
% R E N
stseaspos_spec=[ 1 1 1    %Ssa
    1 1 1];  %365.25 solar days

freq_cpsd= [
    0.005461              %Ssa - semi-annual
    0.002730375267416];   %365.25 solar days - annual
% !!!!! has to be consistent with VIE_MOD/stseasonal.m !!!!!!!
sol2sid=366.2422/365.2422; %1.002737909
Psd = 1./freq_cpsd; % sid day
P_stseaspos = Psd./sol2sid; % solar day


parGS(g.g_stseaspos_Ar(1)).spec = {stseaspos_spec(:,1)};
parGS(g.g_stseaspos_Ae(1)).spec = {stseaspos_spec(:,2)};
parGS(g.g_stseaspos_An(1)).spec = {stseaspos_spec(:,3)};

% Tidal ERP coefficients
% initialize special tide parameters 'spectid,-ret' for tidal ERP terms
if parGS(g.g_tidpm).id == 1
    parGS = spec_tiderp(parGS,g,path_level);
end



% Love and Shida numbers (solid Earth tides)
% GUI for Love and Shida numbers
if (parGS(g.g_love).id==1 || parGS(g.g_shida).id==1)
    guiglob_hl
    uiwait(guiglob_hl)
end

parGS_hl=[];
if parGS(g.g_love).id==1 || parGS(g.g_shida).id==1
    load([path_level 'DATA/GLOB/parGS_hl'],'parGS_hl');
    if parGS(g.g_love).id==1 && isempty(parGS_hl.love.nr)
        parGS(g.g_love).id=0;
        fprintf('\n Love numbers won''t be estimated! \n You didn''t specify any! \n')
    end
    if parGS(g.g_shida).id==1 && isempty(parGS_hl.shida.nr)
        parGS(g.g_shida).id=0;
        fprintf('\n Shida numbers won''t be estimated! \n You didn''t specify any! \n')
    end
end



pathGS.path_out = [path_level 'OUT/GLOB/'];

dir_in=pathGS.L2;

% antenna breaks from an external file (provided by IVS analysis coordinator (A.Nothnagel, Bonn))
% 'VLBI-DISCONT.txt' converted to a mat.file: antbr
load ([path_level 'DATA/GLOB/TRF/DISCONT/' pathGS.discont])    % antbr
% in case that station is missing in the lower part of the file
% 'VLBI-DISCONT.txt', empty space is changed to 0.
for i=1:length(antbr)
    if isempty(antbr(i).break)
        antbr(i).break=0;
    end
end


% stations, for which constant velocity will be estimated, although there are
% discontinuities in the position
velconst=[''];
if parGS(g.g_vel(1)).id==1
    file_velconst = [path_level 'DATA/GLOB/TRF/VELOC/' pathGS.velconst];
    fid=fopen(file_velconst);
    %     if fid~=-1
    %         velconst=textread(file_velconst, '%8c', 'whitespace','','commentstyle','matlab');
    %         fclose(fid);
    %     end
    if fid~=-1
        velconst_lines=textread(file_velconst, '%s', 'delimiter', '\n', 'whitespace','','commentstyle','matlab');
        fclose(fid);
    else
        velconst_lines=[];
        velconst_interv=[];
    end
    for i=1:length(velconst_lines)
        velconst(i,:)=velconst_lines{i}(1:8);
        if isempty (str2num(velconst_lines{i}(9:end)))
            velconst_interv{i}='0';
        else
            velconst_interv{i}=velconst_lines{i}(9:end);
        end
    end
    
    
    % velocity ties between different stations
    file_velties=[path_level 'DATA/GLOB/TRF/VELOC/TIES/' pathGS.velties];
end
if parGS(g.g_coord(1)).id==1
    % path to the file with station names for NNT/NNR condition
    file_datumsta=[path_level 'DATA/GLOB/TRF/DATUM/' pathGS.datumsta];
end

if parGS(g.g_srade(1)).id==1
    % path to the file with source names for NNR condition
    file_datumsou=[path_level 'DATA/GLOB/CRF/DATUM/' pathGS.datumsou];
    
    % path to the file with source names which will be fixed to apriori
    % coordinates
    file_fixedsou=[path_level 'DATA/GLOB/CRF/FIXED_SOURCES/' pathGS.fixedsou];
    
    % path to the file with source names which will be session-wise reduced
    file_soured  =[path_level 'DATA/GLOB/CRF/REDUCE/' pathGS.soured];
end

if parGS(g.g_ao).id==1
    file_ao = [path_level 'DATA/GLOB/TRF/AO/' pathGS.ao];
end

if parGS(g.g_stseaspos_Ar(1)).id==1 || parGS(g.g_stseaspos_Ae(1)).id==1 || parGS(g.g_stseaspos_An(1)).id==1
    file_stseason = [path_level 'DATA/GLOB/TRF/STSEASON/' pathGS.stseason];
end



if parGS(g.g_aplrg).id==1
    file_aplrg = [path_level 'DATA/GLOB/TRF/APLRG/' pathGS.aplrg];
end

fprintf('---------------------------------------------------------------\n')
fprintf('|                   Welcome to VIE_GLOB!                      |\n')
fprintf('---------------------------------------------------------------\n\n')

%  path = [pathGS.path_in dir_in '\'];
%  dir_all = dir([path '*_par_glob.mat']);
%  lse = length(dir_all);


% changed by Sun J. to make it compatible with LINUX
% +sun
path = [pathGS.path_in dir_in];
tmp = dir(path);
lse = 0;
for it = 1 : length(tmp)
    if ~isempty(strfind(tmp(it).name,'_par_glob.mat'))    
        lse = lse + 1;
        dir_all(lse)=tmp(it);
    end
end
path = [path '/'];
% -sun



if lse==0
    fprintf('\n There are no files in the %s directory!\n The program will stop!\n\n', path)
end

aposteriori=[];
for ise = 1:lse
    load ([path dir_all(ise).name]);
    sesname=dir_all(ise).name;
    ses(ise)={[sesname(1:strfind(sesname,'_par_glob.mat')-1)]};
    if isfield(glob2.opt, 'mo') %failes if only N matrices are calculated in previous run
        if ~isempty(glob2.opt.mo)
            aposteriori(ise) = glob2.opt.mo;
        end
    end
end

% Delete sessions with big a posteriori sigma
rms_check = 1;
maxRMS=pathGS.maxRMS;
badses='';
badses_mo=[];
if ~isempty(aposteriori)
    if rms_check
		[a,bad]=find(aposteriori > maxRMS | isnan(aposteriori));
    else
        bad = [];
    end
	badses=ses(bad);
    for i = 1:length(bad)
        fprintf('\n Session %s has a high RMS and will be removed \n', badses{i});
    end
    badses_mo=aposteriori(bad);
    ses(bad)=[];
end

lse=size(ses,2);

for ise=1:lse
    load ([path ses{ise} '_par_glob.mat']);
    lTPl(ise)=glob2.opt.lTPl;
    nobserv(ise)=glob2.opt.total_obs;
    nconstr(ise)=glob2.opt.nconstr;
end


%%

fprintf(' Computing ... \n\n')

for i=1:length(parGS)
    parGS(i).oldcol=[];
    parGS(i).mjd=[];
    parGS(i).mjdi=[];
    parGS(i).newcol=[];
    parGS(i).mjdstep=[];
end

%% Station positions and velocities
if save_intermediate_results_flag
	save('workspace1.mat');
	%load('workspace1.mat');
end
fprintf(' 1 ... Checking the station/source names \n\n')

% Get all names and coordinates of the stations
[refname_all,refantbr_all,antactiv_all,mjd_all] = refname_ant(path, ses, antbr, maxRMS, dir_in);

% Read in stations which will be reduced
anames_fo=[''];  % fo = "few observations"
file_antred=([path_level 'DATA/GLOB/TRF/REDUCE/' pathGS.antred]);

% fid=fopen(file_antred);
% if fid~=-1
%     [antred_lines] =textread(file_antred, '%s', 'delimiter', '\n', 'whitespace','','commentstyle','matlab');
%     fclose(fid);
% end
%
% for i=1:length(antred_lines)
%     anames_fo(i,:)=antred_lines{i}(1:8);
%     if isempty (str2num(antred_lines{i}(9:end)))
%         anames_fo_interv{i}='0';
%     else
%         anames_fo_interv{i}=antred_lines{i}(9:end); % intervals which should be reduced
%     end
% end

fid=fopen(file_antred);
if fid~=-1
    anames_fo=textread(file_antred, '%8c', 'whitespace','','commentstyle','matlab');
    fclose(fid);
end



% Leave only stations going to the global adjustment
% Delete names of stations which will be reduced, i.e., which won't be in the
% global adjustment from: refname, refantbr, antactiv
if ~isempty(anames_fo)
    [refname, refantbr, antactiv]=refname_wored(anames_fo,refname_all,refantbr_all,antactiv_all);
else
    refname=refname_all;
    refantbr=refantbr_all;
    antactiv=antactiv_all;
end


% names of the stations where the seasonal position variation will be
% estimated
stseasname=[''];
if parGS(g.g_stseaspos_Ar(1)).id==1 || parGS(g.g_stseaspos_Ae(1)).id==1 || parGS(g.g_stseaspos_An(1)).id==1
    %stseasname_orig=refname;
    % stseasname_orig=['WETTZELL']
    
    % as datum stations
    fid=fopen(file_datumsta);
    if fid~=-1
        stseasname_orig=textread(file_stseason, '%8c', 'whitespace','','commentstyle','matlab');
        fclose(fid);
    end
    
    % check if these stations are really in the global adjustment
    k=0;
    for i=1:size(stseasname_orig,1)
        f=find(strcmp(cellstr(stseasname_orig(i,:)),cellstr(refname))==1);
        if ~isempty(f)
            k=k+1;
            stseasname(k,:)=stseasname_orig(i,:);
        end
    end
end
lnstseasvar=size(stseasname,1);



% stations where axis offset will be estimated
aostname=[''];
aostname_orig=[''];
if parGS(g.g_ao).id==1
    
    fid=fopen(file_ao);
    if fid~=-1
        aostname_orig=textread(file_ao, '%8c', 'whitespace','','commentstyle','matlab');
        fclose(fid);
    end
    
    % check if these stations are really in the global adjustment
    k=0;
    for i=1:size(aostname_orig,1)
        f=find(strcmp(cellstr(aostname_orig(i,:)),cellstr(refname))==1);
        if ~isempty(f)
            k=k+1;
            aostname(k,:)=aostname_orig(i,:);
        end
    end
end
lnao=size(aostname,1);


% stations where APL RgC will be estimated
rgstname=[''];
rgstname_orig=[''];
if parGS(g.g_aplrg).id==1
    
    fid=fopen(file_aplrg);
    if fid~=-1
        rgstname_orig=textread(file_aplrg, '%8c', 'whitespace','','commentstyle','matlab');
        fclose(fid);
    end
    
    % check if these stations are really in the global adjustment
    k=0;
    for i=1:size(rgstname_orig,1)
        f=find(strcmp(cellstr(rgstname_orig(i,:)),cellstr(refname))==1);
        if ~isempty(f)
            k=k+1;
            rgstname(k,:)=rgstname_orig(i,:);
        end
    end
end
lrg=size(rgstname,1);


ln=size(antactiv,1);
antactiv(ln+1,:)=mjd_all;

%  plot the activity of the antennas
[antactiv1,idses]=sortrows(antactiv',ln+1);

% save for SINEX
antactiv_refnameOrder = antactiv1';


antactiv1(:,ln+1)=[];
clear antactiv
antactiv=antactiv1';
[antactiv_plot,id_ant]=sortrows(antactiv,-[1:size(antactiv,2)]);
antname_plot=refname(id_ant,:);

plot_antactiv(antname_plot,antactiv_plot)


clear antactiv1 r i j k

% session order according to time "ses_time" ( according to alphabet "ses")
ses_time=ses(idses);

%%

k=0; % refnamec: names of antenna, also multiple if there are breaks in the position, that will be considered
for i = 1:length(refantbr)
    [xx,inb]=find(refantbr(i).interv > 0);
    for j=1:length(inb)
        refnamec(k+j,:)=refantbr(i).name;
    end
    k=size(refnamec,1);
end

%%

clear inant fs

lnc=0; % lnc: "how many columns for i.e. x-coordinates" will be in the N-matrix (consideration of breaks)

% number of columns for station coordinates in the new N-matrix
if parGS(g.g_coord(1)).id==1
    lnc=size(refnamec,1);
    parGS(g.g_coord(1)).newcol=1:lnc; % columns in the N-matrix for antenna coordinates x
    parGS(g.g_coord(2)).newcol=lnc+1:2*lnc;     % y
    parGS(g.g_coord(3)).newcol=2*lnc+1:3*lnc;   % z
    ncoord=3*lnc; % number of coordinates
end
if parGS(g.g_coord(1)).id==0  % antenna coordinates
    ncoord=0;
end

% number of columns for station velocities in the new N-matrix
nveloc=0;
lnv=0;
if parGS(g.g_vel(1)).id==1
    lnv=lnc;
    parGS(g.g_vel(1)).newcol=ncoord+1:ncoord+lnv;           % vx
    parGS(g.g_vel(2)).newcol=ncoord+lnv+1:ncoord+2*lnv;     % vy
    parGS(g.g_vel(3)).newcol=ncoord+2*lnv+1:ncoord+3*lnv;   % vz
    nveloc=3*lnv;
end


%% source coordinates

nsou=0;
nvsou=0;
lns=0;
qrefname.IERS=[''];
fixedsou=[''];
reducsou=[''];

if parGS(g.g_srade(1)).id==1
    % fixed sources
    fid=fopen(file_fixedsou);
    if fid~=-1
        fixedsou=textread(file_fixedsou, '%8c', 'whitespace','','commentstyle','matlab');
        fclose(fid);
    end
    % reduced sources
    fid=fopen(file_soured);
    if fid~=-1
        reducsou=textread(file_soured, '%8c', 'whitespace','','commentstyle','matlab');
        fclose(fid);
    end
    
    
    % create qrefname: names of the sources in the adjustment
    [qrefname, souactiv] = refname_sou(path, ses, fixedsou, reducsou);
    
    lns = size(qrefname.IERS,1); % final number of sources in the global adjustment
    
    % for SINEX
    souactiv_qrefnameOrder1 = souactiv;
    souactiv_qrefnameOrder1(lns+1,:)=mjd_all;
    souactiv_qrefnameOrder1=sortrows(souactiv_qrefnameOrder1',lns+1);
    souactiv_qrefnameOrder=souactiv_qrefnameOrder1';
    
    % Figure 3: plot activity of the sources
    [souname_plot,souactiv_plot] = plot_souactiv(qrefname.IERS,souactiv,mjd_all);
    
    
    % number of columns for source coordinates in the new N-matrix
    parGS(g.g_srade(1)).newcol = ncoord+nveloc+1 : ncoord+nveloc+lns;           % RA
    parGS(g.g_srade(2)).newcol = ncoord+nveloc+lns+1 : ncoord+nveloc+2*lns;     % De
    nsou = 2*lns;
    
    if parGS(g.g_svrade(1)).id==1
        parGS(g.g_svrade(1)).newcol = ncoord+nveloc+nsou+1 : ncoord+nveloc+nsou+lns;     % vel RA
        parGS(g.g_svrade(2)).newcol = ncoord+nveloc+nsou+lns+1 : ncoord+nveloc+nsou+2*lns;     % vel De
        nvsou = 2*lns;
    end
    
    
end

clear nrq



% Next parameters: EOP
[lmjd_eop,llove,lshida,lFCNset,laccSSB,lstseaspos1,lhpole,llpole,lgamma,...
    ltidpm,ltidut,parGS] = numpar_next(parGS,parGS_hl,lse,ses,path);
lstseaspos=lstseaspos1*lnstseasvar; % for each station is a set of amplitudes


%%  Loop to get the new position of the estimated parametres
if save_intermediate_results_flag
	save('workspace2.mat');
	%load('workspace2.mat');
end
fprintf(' 2 ... Allocation of the parameters in the new N matrix \n\n')

%------------------------------special EOP-------------------
special_EOP = 0;
if special_EOP
    special_EOP_file = 'fix_EOP_for_single_baseline_sessions.txt';
    fid_special_EOP = fopen(['../DATA/GLOB/EOP/' special_EOP_file]);
    special_EOP_sessions = textscan(fid_special_EOP,'%s');
    special_EOP_sessions = char(special_EOP_sessions{1});
    special_EOP_sessions = cellstr(special_EOP_sessions(:,1:end-5));
    fclose(fid_special_EOP);
    reduced_flag = 0;
end
%------------------------------------------------------------
for ise=1:lse
	fprintf(1, ['load session ' ses{ise} '\n'])
    load ([path ses{ise} '_par_glob.mat']);
    load ([path ses{ise} '_an_glob.mat']);
    
    %------------------------------special EOP-------------------
    if special_EOP
        index_xpol = find(strcmp({parGS.name},'xpol'));
        if reduced_flag
            parGS(index_xpol).id = 2; %xpol
            parGS(index_xpol+1).id = 2; %ypol
            parGS(index_xpol+2).id = 2; %dut1
            parGS(index_xpol+3).id = 2; %nutdx
            parGS(index_xpol+4).id = 2; %nutdy
            reduced_flag = 0;
        end
        
        if length(ses{ise})==9
            ases = ses{ise};
        else
            ases = ses{ise}(1:end-5);
        end
            
        
%         if (parGS(index_xpol).id == 2) && (sum(strcmp(ses(ise,1:end-5),special_EOP_sessions))>0)
        if (parGS(index_xpol).id == 2) && (sum(strcmp(ases,special_EOP_sessions))>0)
            parGS(index_xpol).id = 0; %xpol
            parGS(index_xpol+1).id = 0; %ypol
            parGS(index_xpol+2).id = 0; %dut1
            parGS(index_xpol+3).id = 0; %nutdx
            parGS(index_xpol+4).id = 0; %nutdy
            reduced_flag = 1;
        end
        
    end
    %------------------------------------------------------------
    
    antenna = glob1.an;
    x_ = glob2.x;
    
    % columns and mjd of parameters in this session in the old N matrix
    [parGS] = par_oldcol(x_,parGS);

    
    % stop vie_glob if parameters which should be estimated are not in the old
    % N-matrices
    for i=1:length(parGS)
        if parGS(i).id==1  & isempty(parGS(i).oldcol)
            fprintf('\n Parameter %s is not included in the N matrix in the %1.0f. session! \n',parGS(i).name,ise)
            fprintf(' VIE_GLOB will stop!\n')
            return
        end
        if parGS(i).id==2 & isempty(parGS(i).oldcol)
            fprintf('\n WARNING! You wanted to reduce %s, but it is not included in the N matrix in the %1.0f. session! \n',parGS(i).name,ise)
        end
    end
    
    
    % columns of antennas, which position will be
    % session-wise reduced
    antfo_old=[];  % indices of stations to be reduced (wrt refname)
    nfo=size(anames_fo,1);
    for i=1:nfo
        [inant,xx]=find(strcmp(cellstr(anames_fo(i,1:8)),cellstr(antenna.name)) == 1);
        antfo_old=[antfo_old; inant];
    end
    antfo_old=unique(antfo_old);
    
    % order of the sources in this session
    qnames=[''];
    if parGS(g.g_srade(1)).id==1
        for iso=1:length(glob2.x.source)
            qnames(iso,:) = glob2.x.source(iso).IERSname;
        end
    end
    
    % session-wise reduced sources
    reducsou_old=[];  % indices of sources to be reduced (wrt qnames)
    nfo=size(reducsou,1);
    for i=1:nfo
        [insou,xx]=find(strcmp(cellstr(reducsou(i,1:8)),cellstr(qnames)) == 1);
        reducsou_old=[reducsou_old; insou];
    end
    reducsou_old=unique(reducsou_old);
    
    % NEW COLUMNS for the parameters
    [ar,actp,parGS] = par_newcol(parGS,parGS_hl,antenna,refantbr,refnamec,qnames,qrefname.IERS,aostname,rgstname,ise,lnc,lnv,lns,lnao,...
        llove,lshida,lFCNset,stseasname,lse,nvsou,laccSSB,lhpole,llpole,lgamma,lrg,...
        ltidpm,ltidut);

    parsplit(ise).ar=ar; % actual-reference
    
    
    % parameters, which will be SESSION-WISE REDUCED
    [redpos_all redcoor redsou] = redpar(parGS,antfo_old,reducsou_old);
    parsplit(ise).redpos=redpos_all; % all columns in the old N matrix with parameters, which will be reduced
    parsplit(ise).redcoor=redcoor; %columns in the old N matrix for stations - needed for backwards solution
    parsplit(ise).redsou=redsou;  %columns in the old N matrix for sources - needed for backwards solution
    
    if pathGS.bckwrdsol==1
        save([path_level 'OUT/GLOB/BACKWARD_DATA/' pathGS.out '/parGS_' pathGS.L2 '_' num2str(ise)],'parGS')
    end
    
    
    % store reference pressure and a priori APL RgC for stations in the
    % global adjustment
    aname=antenna.name;
    for i=1:size(aname,1)
        f=find(strcmp(cellstr(aname(i,:)),cellstr(rgstname))==1);
        if ~isempty(f)
            id=f;
            apriori_aplrg(id,:)=antenna.aplrg(i,:);
        end
    end
    
    clear glob2 qnames
end

if pathGS.bckwrdsol==1
    save([path_level 'OUT/GLOB/BACKWARD_DATA/' pathGS.out '/parsplit_' dir_in],'parsplit')
    save([path_level 'OUT/GLOB/BACKWARD_DATA/' pathGS.out '/paths_' dir_in],'pathGS')
end

for i = 1:length(parGS)
    if isnumeric(parGS(i).newcol) % the amplitudes of seasonal station positions variations are stored as cells (they are not necessary for this loop)
        parGS(i).newcol=unique(parGS(i).newcol);
        parGS(i).newcol=nonzeros(parGS(i).newcol)';
    end
end

% aesp --> all estimated parameters: size of the new N-matrix
aesp = ncoord+nveloc+nsou+nvsou+lmjd_eop+lnao+llove+lshida+lFCNset+laccSSB+lstseaspos+lhpole+llpole+lgamma+lrg+...
    ltidpm+ltidut;

%% Sorting and splitting the N-matrix and b-vector
if save_intermediate_results_flag
	save('workspace3.mat');
	%load('workspace3.mat');
end
fprintf(' 3 ... Sorting and splitting of the old N matrices \n\n')


redparam=0;
for ise=1:lse
    load ([path ses{ise} '_Nb_glob.mat']);
    fprintf('       session %4.0f out of %4.0f -- %s\n',ise,lse,ses{ise})
    
    Norig = sparse(glob3.N);
    borig = sparse(glob3.b);
    clear glob3
    
    ar = parsplit(ise).ar; % estimated parameters
    redpos = parsplit(ise).redpos;    % position of columns of the parameters, which are to be reduced
    
    %% Splitting and division of the N-matrix for reduction
    
    N11=zeros(aesp);
    N12=zeros(aesp,length(redpos));
    
    N11(ar(:,2),ar(:,2))=Norig(ar(:,1),ar(:,1));
    N12(ar(:,2),:)=Norig(ar(:,1),redpos);
    N21=N12';
    N22=Norig(redpos,redpos);
    
    b1=zeros(aesp,1);
    b1(ar(:,2))=borig(ar(:,1));
    b2=borig(redpos);
    
    
    %% REDUCTION
    
    redcoor = parsplit(ise).redcoor;
    redsou = parsplit(ise).redsou;
    
    %% constraints for antenna coordinates
    
    if ~isempty(redcoor)
        
        [xx idco]=intersect(redpos,redcoor);
        clear xx
        
        nco=length(redcoor);
        
        H=zeros(nco,length(redpos));
        for i=1:nco
            H(i,idco(i))=1;
        end
        
        cons(1:nco) = 1/4; % weight
        Ph=diag(cons);
        HTPH =  H'*Ph*H;
        N22 = N22 + HTPH;
        
        w(1:nco)=0; %0 mas - ideal value
        b2(idco) = b2(idco) + Ph*w';
        
        clear Ph w H HTPH cons
    end
    
    %% constraints for source coordinates
    if ~isempty(redsou)
        
        [xx idsr]=intersect(redpos,redsou);
        clear xx
        
        nsr=length(redsou);
        
        H=zeros(nsr,length(redpos));
        for i=1:nsr
            H(i,idsr(i))=1;
        end
        
        cons(1:nsr) = 1/4; % weight
        Ph=diag(cons);
        HTPH =  H'*Ph*H;
        N22 = N22 + HTPH;
        
        w(1:nsr)=0; %0 mas - ideal walue
        b2(idsr) = b2(idsr) + Ph*w';
        
        clear Ph w H HTPH cons
        
    end
    
    %%
    
    invN22 = inv(N22);
    Nred = N11 - N12*invN22*N21;
    
    bred2 = N12*invN22*b2;
    if isempty(bred2); bred2=0; end

    bred = b1 -  bred2;
    
    % %    size_N22=size(N22); rank_N22=rank(N22); if size_N22>rank_N22; ise; end
    
    
    %% STACKING
 
    
    b2invN22b2 = b2'*invN22*b2;
    if isempty(b2invN22b2); b2invN22b2=0; end
    
    
    if ise==1
        N_stc = Nred;
        b_stc = bred;
        lTPl_reduc_stc = lTPl(1) - b2invN22b2;
        
    else
        N_stc = N_stc + Nred;
        b_stc = b_stc + bred;
        lTPl_reduc_stc = lTPl_reduc_stc + (lTPl(ise) - b2invN22b2);
        
    end
    
    
    if pathGS.bckwrdsol==1
        save([path_level 'OUT/GLOB/BACKWARD_DATA/' pathGS.out '/b2_' pathGS.L2 '_' num2str(ise)],'b2')
        save([path_level 'OUT/GLOB/BACKWARD_DATA/' pathGS.out '/N22_' pathGS.L2 '_' num2str(ise)],'N22')
        save([path_level 'OUT/GLOB/BACKWARD_DATA/' pathGS.out '/N21_' pathGS.L2 '_' num2str(ise)],'N21')
    end
    
    
    redparam = redparam+length(b2);
    
    clear b1 b2 b N N11 N12 N21 N22 invN22
    clear Norig borig ar redpos
    
    
end


Nfree = N_stc;
bfree = b_stc;

clear N_stc b_stc

%%
if save_intermediate_results_flag
	save('workspace4.mat');
	%load('workspace4.mat');
end
if parGS(g.g_coord(1)).id==1 || parGS(g.g_srade(1)).id==1
    fprintf('\n 4 ... Preparing B matrix for datum definition \n\n')
end
%% Stations for NNT/NNR condition
datumantc=[''];

if parGS(g.g_coord(1)).id==1
    [inidantc,datumantc]=datum(file_datumsta,refnamec);
end

if parGS(g.g_vel(1)).id==1
    inidantv=inidantc;
    datumantv=datumantc;
end

%%
% The apriori catalogue coordinates for NNT/NNR condition
if parGS(g.g_coord(1)).id==1
    
    % load the apriori trf catalogue (and rename it to 'trf.' structure)
    load ([path ses{1} '_par_glob.mat']);
    
    trffile=glob2.opt.trf;
    
    % COPIED FROM VIE_INIT.MAT (March 2013)
    % (1) TRF: load superstation file
    % if a superstation file is given (==~manual textfile is given) - load that
    if strcmp(trffile{1}(end-3:end), '.mat')
        trf=load(trffile{1});
        nam=fieldnames(trf);
        trf=eval(['trf.' nam{1}]);
    else % a manual trf file is given -> BUT load default superstat file in any case!!
        trf=load([path_level 'TRF/superstation.mat']);
        nam=fieldnames(trf);
        trf=eval(['trf.' nam{1}]);
    end
    
    % (2) TRF: load manual TRF (if chosen)
    if strcmpi(trffile{1}(end-3:end), '.txt')
        % put in superstation file (variable trf)
        trf=manualTrfToSuperstatTrf(trf,trffile{1}); % function from VIE_INIT
        
        % write "trfname" to trffile{2}
        trffile{2}='manualTrf';
    end
    
    
    % intervals for station coordinates + a priori coordinates in these
    % intervals, needed for NNT/NNR condition
    [X0,Y0,Z0,brstart,brend,refantbr] = refantbr_ext(refantbr,refname,trf,trffile,path,ses);
end

slc=0; Bc=[];
scale=pathGS.datumst_scale;
if parGS(g.g_coord(1)).id==1
    Re=6371000;
    X0r=X0/Re; Y0r=Y0/Re; Z0r=Z0/Re;
    excidantc=setdiff(1:size(refnamec,1),inidantc); %indices of excluded stations from the NNT/NNR condition (wrt refnamec)
    Bc=nntnnr(X0r,Y0r,Z0r,lnc,excidantc,scale);
    slc=3;
end

slv=0; Bv=[];
if parGS(g.g_vel(1)).id==1
    Bv=Bc;
    slv=3; %  vx,vy,vz
end

% plot station in the solution (included+excluded in the nnr/nnt condition)
if parGS(g.g_coord(1)).id==1
    try
        plot_ant(X0,Y0,Z0,excidantc,refnamec)
    catch plot_network_error
        fprintf(2, getReport(plot_network_error, 'extended', 'hyperlinks', 'off')); % display error message
    end
end


%% velocity ties
B_vties=[];
if parGS(g.g_vel(1)).id==1
    velties=read_velocity_ties(file_velties);
    B_vties=create_B_vties(refnamec,velties,velconst,velconst_interv);
end

%% NNR condition for sources
datumsouc=[''];
excidsouc=[];
RA_all=[]; De_all=[];
sls=0;
crffile=[''];
sou_dz=pathGS.datumsou_dz; % '1' NNR+dz; '0' NNR
Bq=[];

if parGS(g.g_srade(1)).id==1
    
    load ([path ses{1} '_par_glob.mat']);
    crffile = glob2.opt.crf;
    
    % COPIED FROM VIE_INIT.MAT (March 2013)
    % CRF:
    % (1) CRF: load supersource file
    if strcmp(crffile{1}(end-3:end), '.mat')
        crf=load(crffile{1});
        nam=fieldnames(crf);
        crf=eval(['crf.', nam{1}]);
    else
        fprintf('not working yet (manual CRF!)\n')
    end
    
    
    
    for i=1:size(qrefname.IERS,1)
        qname=qrefname.IERS(i,:);
        
        for k = 1:length(crf)
            r(k)= strcmp(cellstr(qname),cellstr(crf(k).IERSname));
        end
        SouId=find(r==1);
        
        
        % save source apriori coordinates from the crf catalogue
        if isempty(SouId)==0
            if ~isempty(crf(SouId).(crffile{2})) % if it is in chosen catalogue
                RA_all(i) = crf(SouId).(crffile{2}).ra;
                De_all(i) = crf(SouId).(crffile{2}).de;
            elseif ~isempty(crf(SouId).vievsCrf) % backup vievsTRF catalogue
                RA_all(i) = crf(SouId).vievsCrf.ra;
                De_all(i) = crf(SouId).vievsCrf.de;
            end
        else
            fprintf('\n The source %s was not found in the supersource file!\n', qname)
        end
    end
    
    if isempty(pathGS.datumsou)==0
        [inidsouc,datumsouc]=datum(file_datumsou,qrefname.IERS);
        excidsouc=setdiff(1:size(qrefname.IERS,1),inidsouc); %indices of excluded sources from the NNR condition
        
        Bq=nnr_source(RA_all,De_all,lns,excidsouc,sou_dz);
        sls=3;
        
        % plot sources in the solution (included+excluded in the NNR condition)
        plot_sou(RA_all,De_all,excidsouc);
        
    end
end


slvs=0; Bvq=[];
if parGS(g.g_svrade(1)).id==1
    Bvq=Bq;
    slvs=3;
end


%% Apply B matrices (datum definition)

clear B
sBc1=size(Bc,1);
sBc2=size(Bc,2);
sBv1=size(Bv,1);
sBv2=size(Bv,2);
sBq1=size(Bq,1);
sBq2=size(Bq,2);

sBvq1=size(Bvq,1);
sBvq2=size(Bvq,2);
sBvt1=size(B_vties,1);
sBvt2=size(B_vties,2);









Bcvqv = [Bc                zeros(sBc1,sBv2)  zeros(sBc1,sBq2)   zeros(sBc1,sBvq2)
    zeros(sBv1,sBc2)  Bv                zeros(sBv1,sBq2)   zeros(sBv1,sBvq2)
    zeros(sBq1,sBc2)  zeros(sBq1,sBv2)  Bq                 zeros(sBq1,sBvq2)
    zeros(sBvq1,sBc2)  zeros(sBvq1,sBv2)  zeros(sBvq1,sBq2)  Bvq
    zeros(sBvt1,sBc2) B_vties           zeros(sBvt1,sBq2)  zeros(sBvt1,sBvq2)
    ];

B = [Bcvqv zeros(sBc1+sBv1+sBq1+sBvq1+sBvt1,lmjd_eop+lnao+llove+lshida+...
    lFCNset+laccSSB+lstseaspos+lhpole+llpole+lgamma+lrg+...
    ltidpm+ltidut)];

b=[bfree; zeros(sBc1+sBv1+sBq1+sBvq1+sBvt1,1)];


N = [Nfree B'
    B     zeros(size(B,1))];



%   size_N=size(N);
%   Nxx=full(N);
%   rank_N=rank(Nxx);

%% Final solution

invN=inv(N);
x = invN*b;

%% ACCURACY CRITERIA


vTPv = lTPl_reduc_stc - x'*b;

nobserv_all = sum(nobserv);
uparam = length(b)+redparam;
nconstr_all = sum(nconstr);


% aposteriori variance of unit weight
sigma0_2 = vTPv/(nobserv_all - uparam + nconstr_all + size(B,1));
sigma0 = sqrt(sigma0_2);


% Covariance matrix of the estimated parameters (squared values)
Q=sigma0_2*invN;

% aposterioti variance of estimated parameters
diagQ = diag(Q);
varpar=sqrt(diagQ(1:aesp));
varpar=full(varpar);


%%
eb=[];
if parGS(g.g_coord(1)).id==1
    if parGS(g.g_vel(1)).id==1
        epoch(1,1:size(refnamec,1)) = glob2.opt.refvel_mjd; % epoch is chosen in vie_lsm
    else
        epoch = fepochs(refantbr, refnamec, brstart, brend);
    end
    eb=[epoch; brstart; brend]';
end

clear globsol


x=full(x(1:aesp));

% SAVING into structure "globsol"
globsol.Nfree = Nfree;
globsol.bfree = bfree;
globsol.N = N;
globsol.b = b;
globsol.x = x;
globsol.sigma0 = sigma0;
globsol.Q = Q;
globsol.varpar = varpar;
globsol.nr_of_conditions = size(B,1);
globsol.sessions = ses;
globsol.maxRMS = maxRMS;
globsol.numObs = nobserv_all + nconstr_all; % SINEX
globsol.numEst = uparam; % SINEX
globsol.vTPv = vTPv; % SINEX
globsol.varfac = sigma0_2; % SINEX

if parGS(g.g_coord(1)).id==0
    trffile='';
end
if parGS(g.g_ao).id==0
    aostname='';
end
if parGS(g.g_aplrg).id==0
    rgstname='';
    apriori_aplrg=[];
end


globsol = saveres(parGS,parGS_hl,globsol,refname,refnamec,datumantc,scale,...
    qrefname,datumsouc,fixedsou,sou_dz,RA_all,De_all,eb,...
    trffile,crffile,aostname,stseasname,P_stseaspos,rgstname,apriori_aplrg);


if parGS(g.g_coord(1)).id==1
    globsol.antenna.apr_pos = refantbr;
    globsol.antenna.antactiv=antactiv_refnameOrder; %for SINEX
end
if parGS(g.g_srade(1)).id==1
    globsol.source.souactiv=souactiv_qrefnameOrder; %for SINEX
end

if save_intermediate_results_flag
	save('workspace5.mat');
	%load('workspace5.mat');
end

% CREATES TXT file with estimates and saves it in [../../OUT/GLOB/_ESTIMATES/pathGS.out/]
createTXT(globsol,ses_time,badses,badses_mo,pathGS)

% Writes text file with tidal ERP amplitudes and formal errors in
% [../../OUT/GLOB/_ESTIMATES/pathGS.out/]
if parGS(g.g_tidpm).id == 1
    tiderpTXT(globsol,ses_time,pathGS)
end

%% Estimation of a new trf/crf

if (parGS(g.g_coord(1)).id==1 || parGS(g.g_vel(1)).id==1)
    % a file with the new coordinates will be created in '/TRF_new/'
    newtrf(refantbr,globsol,pathGS,trf);
end

if (parGS(g.g_srade(1)).id==1)
    newcrf(globsol,pathGS,crf);
%    newcrf_binary(globsol,pathGS); % binary format of the crf catalogue
end

%% Save TRF and CRF in the sinex format
% currently it is working only for the TRF (pos+vel) and CRF (pos)
% Uncomment if you want to create the sinex files:

% write_sinex_vie_glob_CRFTRF(globsol,pathGS,trf,crf);



%%
save([pathGS.path_out '_ESTIMATES/' pathGS.out '/globsol_' dir_in],'globsol','-v7.3');


%% Save plots into ../../OUT/GLOB/_PLOTS

if exist([pathGS.path_out '_PLOTS/' pathGS.out])~=7
    mkdir([pathGS.path_out '_PLOTS/' pathGS.out])
end
%orient landscape
print('-f1', '-depsc' ,'-r500',[pathGS.path_out '_PLOTS/' pathGS.out '/ant_activity_' dir_in]);
if parGS(g.g_coord(1)).id==1
    print('-f2', '-depsc' ,'-r500',[pathGS.path_out '_PLOTS/' pathGS.out '/ant_map_' dir_in]);
end
if parGS(g.g_srade(1)).id==1
    print('-f3', '-depsc' ,'-r500',[pathGS.path_out '_PLOTS/' pathGS.out '/sou_activity_' dir_in]);
    if ~isempty(sou_dz) % plot only if NNR (+dz) condition on source coordinates was applied
        print('-f4', '-depsc' ,'-r500',[pathGS.path_out '_PLOTS/' pathGS.out '/sou_map_' dir_in]);
    end
end

%% Run the backward solution
if pathGS.bckwrdsol==1
    fprintf('\n 5 ... Running the backward solution \n\n')
    backward_solution(dir_in, pathGS.out)
end

%%

fprintf('\n\n Done! \n');

fprintf('\n Estimates in TXT format are stored in VieVS/OUT/GLOB/_ESTIMATES/%s/glob_results_%s.txt \n',pathGS.out ,dir_in);
fprintf('\n Figures in EPS format are stored in VieVS/OUT/GLOB/_PLOTS/%s/... \n',pathGS.out);






