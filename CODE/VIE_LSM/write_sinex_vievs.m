% This function writes SINEX output from mat files created by Vienna VLBI
% software (LEVEL3 = final results).
%
%
% created   Matthias Madzak, 28.6.2010
%
% modified
%   10 Feb 2011 by Hana Spicakova: Revised.
%   06 May 2011 by Matthias Madzak: Added revision parts from Hanas version
%       of this function.
%   10 May 2011 by Matthias Madzak: Epoch for stat coords
%       (vievs) is midnight. 
%   19 May 2011 by Hana Spicakova: troposphere gradients and zwd added
%   19 May 2011 by Hana Spicakova: function works also if station
%       coordinates are fixed
%   19 May 2011 by Hana Spicakova: input for this function is name of the
%       session (not process list), because write_sinex_vievs.m is used as internal
%       function of vie_lsm.m
%   03 June 2011 by Hana Spicakova: a problem with matlab function textscan
%      in Matlab version 7.10. solved
%   17 June 2011 by Hana Spicakova: problems with blank space in the name of station solved (NRAO85 3 --> NRAO85_3)
%   17 June 2011 by Hana Spicakova: if a station is missing in ns-codes.txt, it will be refered to it as:
%      siteCode=9999, domes=99999S000, description=antenna(i).name;
%   17 June 2011 by Hana Spicakova: fclose for ns-codes.txt was missing 
%   05 July 2011 by Hana Spicakova: FILE/COMMENT block according to
%      template from Axel Nothnagel - Version 2.0
%   10 July 2012 by Hana Krasna: compatibility of non-tidal APL with
%   version2.0
%   10 July 2012 by Hana Krasna: hydrology loading added
%   30 Oct 2012 by Hana Krasna: ns_codes.txt changed to superstation file
%   30 Oct 2012 by Hana Krasna: the total tropospheric zenith delay added
%   18 Dec 2012 by Hana Krasna: changes in header according to the new
%       representation of constraints in the GUI 2.1
%   30 Apr 2013 by Hana Krasna: TROWET changed to TROTOT
%   18 Dec 2013 by Hana Krasna: ICRF designation for sources added
%   18 Dec 2013 by Hana Krasna: change in the header
%   20 Dec 2013 by Hana Krasna: change to IERS name
%   10 Feb 2014 by Hana Krasna: information about tidal UT corrections
%       added
%   19 Jun 2014 by Hana Krasna: sources can be estimated in vie_lsm with
%       NNR condition
%   17 Jul 2014 by Hana Krasna: description added for the case when sources
%       are estimated with NNR
%   13 Jan 2015 by Matthias Madzak: Comment responsible person can be 
%       changed using input
%   03 Feb 2015 by Andreas Hellerschmied: SINEX files are now saved in the
%       defined LEVEL3 sub-directory (in /DATA/SNX/<LEVEL3_SubDir>)  
%   09 Dec 2015 by A. Hellerschmied: Bug-fix: Now "find(strcmp(...))" is used instead of "strfind(...)"
%   13 June 2016 by Daniel Landskron: Bug fix with Station comments
%   15 July 2016 by David Mayer: IVS source names are now written to the
%   SNX file instead of IERS names
%   26 Sep 2016 by Hana Krasna: update related changes in supersource file,
%   change of the IVS name format string to char
%   25 Jan 2017 by Daniel Landskron: adapted for splitted mf and gradients
%   06 Jul 2018 by Daniel Landskron: VMF3 added to the troposphere models 
%   23 Aug 2018 by Daniel Landskron: adapted for vgosDB
%   29 Nov 2018 by Daniel Landskron: structure of observation restrictions standardized
%   06 Mar 2019 by Daniel Landskron: suffix checkbox added to the sinex files
%   27 Jul 2019 by Daniel Landskron: zwet parameter added to scan structure
%
% call this function:
%
%  - output files are in LEVEL3 folder:
%       write_sinex_vievs('C:/process_list.mat')
%
%  - there is a subfolder between:
%       write_sinex_vievs('C:/process_list.mat', subfolder) or
%       write_sinex_vievs('C:/process_list.mat', 'CONT08')


function write_sinex_vievs(fname, varargin)

ispc = false; % quick fix --> should be properly implemented for early matlab versions
% if subfolder is given as an input


if nargin>1     % Subfolder defined 
    subfolder=varargin{1};
else            % No subfolder defined
    subfolder='/';
end

if nargin>2     % Email Adress, name and suffix defined
    firstname=varargin{2};
    if isempty(firstname)
        firstname='TU';
    end
    
    lastname=varargin{3};
    if isempty(lastname)
        lastname='Wien';
    end
    
    email=varargin{4};
    if isempty(email)
        email='vlbi@geo.tuwien.ac.at';
    end
    
    suffix=varargin{5};
    if ~isempty(suffix)
       suffix = ['_' suffix];
    end
end


% you can load here a real process list, if you want to use this function as
% external one (not as a part of vie_lsm)
% load('c:/vievs/WORK/process_list.mat');
process_list=['xxxx/' fname];

% get size of process list
sizePl=size(process_list);

% if the last row in process_list is empty (' '), delete this last row
if strcmp(process_list(sizePl(1)), ' ')
    process_list(sizePl(1),:)=[];
end



%% calling this function as a script - just for testing...
% clc
% clear all
% 
%  global a_tidefree f_tidefree
%  a_tidefree = 6378136.6; %m      Equatorial radius of the Earth
%  f_tidefree = 1/298.25642;     % Flattening factor of the Earth
%  
%  load('process_list.mat');
%  subfolder='test'
% 
 
%%

% for all sessions in the process list
for pl=1:size(process_list,1)
    files.scan=['../DATA/LEVEL3/', subfolder, '/', process_list(pl,6:end), '_scan.mat'];
    % x_ --> later, is only loaded when parameters were estimated!
    % opt_ --> later, is only loaded when parameters were estimated!
    if nargin>=7
        N_sinex = varargin{6};
        b_sinex = varargin{7};
        col_sinex = varargin{8};
    else
        files.N_sinex=['../DATA/LEVEL3/', subfolder, '/SINEX/N_sinex_', process_list(pl,6:end), '.mat'];
        files.b_sinex=['../DATA/LEVEL3/', subfolder, '/SINEX/b_sinex_', process_list(pl,6:end), '.mat'];
        files.col_sinex=['../DATA/LEVEL3/', subfolder, '/SINEX/col_sinex_', process_list(pl,6:end), '.mat'];
    end
    files.parameter=['../DATA/LEVEL3/', subfolder, '/', process_list(pl,6:end), '_parameter.mat'];
    files.antenna=['../DATA/LEVEL3/', subfolder, '/', process_list(pl,6:end), '_antenna.mat'];
    files.sources=['../DATA/LEVEL3/', subfolder, '/', process_list(pl,6:end), '_sources.mat'];
    
    snxFile=['../DATA/SNX/', subfolder, process_list(pl,6:end), suffix, '.snx'];
    files.superstation=('../TRF/superstation.mat');
    
    
    % check if SNX Folder exist
    snxfolder=['../DATA/SNX/', subfolder];
    if ~exist(snxfolder, 'dir')
        mkdir(snxfolder); % if not, create it
    end
    clear snxfolder

    % check if all files exist
    fieldn=fieldnames(files);
    for k=1:length(fieldn)
        curFileVar=['files.', fieldn{k}];
        if ~exist(eval(curFileVar), 'file')
            fprintf('\nFile ''%s'' not found\n', eval(curFileVar))
            fileNotFound=1;
            break   % if not, stop
        else
            fileNotFound=0;
        end
    end
    if fileNotFound==1
        break
    end

    % load all files in struct 'files'
    for k=1:length(fieldn)
        curFileVar=['files.', fieldn{k}];
        load(eval(curFileVar));
    end
    clear files curFileVar fieldn

    
    outsnx=col_sinex.outsnx;
    
       
    
    % load x_ file (since it's existance depends on estimation)
    % if no estimation was done: x_ = []; 

    if parameter.lsmopt.est_singleses==0
        opt_=[];
        x_=[];
    else
        opt_=['../DATA/LEVEL3/', subfolder, '/opt_', process_list(pl,6:end), '.mat'];
        x_=['../DATA/LEVEL3/', subfolder, '/x_', process_list(pl,6:end), '.mat'];
        load(x_);
        load(opt_);
    end

    % define space between two blocks
    blockSpace={'*', '* -----------------------------------------------------------------------------', '*'};

    % open file ('a': writing, append, 'w': discharge content if any)
    fid=fopen(snxFile, 'w');
    
   
    % add site code and itrfXYZ to x_.antenna struct
    for k=1:size(antenna, 2)
        % get index of current station in itrf struct
        aname=antenna(k).name;
        
        idblank=find(aname(1:max(find(aname(:)~=' ')))==' '); %hana Jun11
        aname(idblank)='_';
        
%         a=find(strcmp({superstations.name}, aname));
        antennaIdx = find(strcmp({superstations.name}, aname));
%         antennaIdx=~cellfun('isempty',a);
        if ~isempty(antennaIdx)
            if length(antennaIdx) > 1
                fprintf('More than one antenna entries of superstation file match antenna(k).name')
            end
%             antennaIdx=find(antennaIdx);
            
            %SiteCode
            curSiteCode=str2double(superstations(antennaIdx).CDP);
            if isnan(curSiteCode) % if siteCode in File = '----'
                antenna(k).siteCode=0;
            else
                antenna(k).siteCode=curSiteCode;
            end
            % Domes
            antenna(k).domes=superstations(antennaIdx).domes;
            % PointCode
            antenna(k).pointCode='A';
            
            %Description (if too long: only 22 characters)
            if length(superstations(antennaIdx).comments) <= 22
                antenna(k).description = superstations(antennaIdx).comments;
            else
                antenna(k).description = superstations(antennaIdx).comments(1:22);
            end
        else
            antenna(k).siteCode=9999;     %hana Jun11
            antenna(k).pointCode='A';
            antenna(k).domes=cellstr('99999S000');
            antenna(k).description=antenna(k).name;
        end
    end
   
%%
    
    
    % add 0 as site codes for stations which do not occur in itrf
    %if sum(cellfun('isempty',{x_.antenna.siteCode})) > 0
    %    x_.antenna(cellfun('isempty',{x_.antenna.siteCode})).siteCode=0;
    %end

    %[pointCode{1:size(x_.antenna, 2), 1}]=deal('A');
    %[obsCode{1:size(x_.antenna, 2), 1}]=deal('R');          % R=VLBI
    %[statDescription{1:size(antenna, 2), 1}]=deal('');
    
    % get eops  
    % get mjds of min intervals of eops (max number of these mjds)
    maxNumMjdIdx=find([size(col_sinex.mjd_xp,2), size(col_sinex.mjd_yp,2), size(col_sinex.mjd_dut1,2), size(col_sinex.mjd_dX,2), size(col_sinex.mjd_dY,2)]==max([size(col_sinex.mjd_xp,2), size(col_sinex.mjd_yp,2), size(col_sinex.mjd_dut1,2), size(col_sinex.mjd_dX,2), size(col_sinex.mjd_dY,2)]));
    switch maxNumMjdIdx(1)
        case 1
            maxMjdEop=col_sinex.mjd_xp;
        case 2
            maxMjdEop=col_sinex.mjd_yp;
        case 3
            maxMjdEop=col_sinex.mjd_dut1;
        case 4
            maxMjdEop=col_sinex.mjd_dX;
        case 5
            maxMjdEop=col_sinex.mjd_dY;
    end
    if outsnx.eop==1 
        allNeededMjds=sort(unique([col_sinex.mjd_xp, col_sinex.mjd_yp, col_sinex.mjd_dut1, col_sinex.mjd_dX, col_sinex.mjd_dY]));
        eop_matrix=eop_out2matrix(x_, opt_, parameter, allNeededMjds);
        eop_matrix=eop_matrix';
        % get dates of eop epochs (=first column in eop matrix)
        eopDates=mjd2yydoysecod(eop_matrix(:,1));
        eopDatesYrStr=num2str(eopDates(:,1));
    end
    numXpol=size(col_sinex.xp,2);
    numYpol=size(col_sinex.yp,2);
    numDut1=size(col_sinex.dut1,2);
    numXnut=size(col_sinex.dX,2);
    numYnut=size(col_sinex.dY,2);
    %numEstEOPs=numXpol+numYpol+numDut1+numXnut+numYnut;
    %numEOPs=size(eop_matrix,1);
    numStat=size(antenna, 2);
    numSou=length(col_sinex.sounames);
    
    %% write headerline
    %if existsnx==0
    snxFormat=2.10;         % format of sinex
    creatingAgency='VIE';   % agency creating this file
    curDate=clock;          % current date and time
    doy=datenum(curDate-[curDate(1),0,0,curDate(4:end)]);  % day of year of current time
    secod=curDate(4)*60*60+curDate(5)*60+curDate(6);       % seconds of day of current time
    yrStr=num2str(curDate(1));
    dataAgency='VIE';       % agency providing the data
    dataStartEnd=mjd2yydoysecod([min([antenna.firstObsMjd]), max([antenna.lastObsMjd])]);
    sessionTimeConsumption=max([antenna.lastObsMjd])-min([antenna.firstObsMjd]);    % needed later when writing epochs where parameters are valid 
    dataStartEndYrStr=num2str(dataStartEnd(:,1));
    obsCode='R';            % observation code (R=VLBI)
    numEstPar=length(b_sinex);
    constrCode=num2str(2);  % 0=tight, 1=significant, 2=unconstrained
    
    solcontents=[''];
    if outsnx.xyz==1; solcontents = ['S ']; end
    if outsnx.eop==1; solcontents =[solcontents 'E ']; end
    if outsnx.sou==1; solcontents =[solcontents 'C ']; end
    if outsnx.zwd==1 || outsnx.tgr==1; solcontents =[solcontents 'T']; end
    
    %                    2.10 VIE   08:  025 :09862  TUW   00: 000  :00000    00:000   :00000               
    fprintf(fid, '%%=SNX %4.2f %3s %2s:%03.0f:%05.0f %3s %02s:%03.0f:%05.0f %02s:%03.0f:%05.0f %1s %05.0f, %1s, %s\n', snxFormat, creatingAgency, yrStr(3:4), doy, secod, dataAgency, dataStartEndYrStr(5:2:end), dataStartEnd(1,2), dataStartEnd(1,3), dataStartEndYrStr(6:2:end), dataStartEnd(2,2), dataStartEnd(2,3), obsCode, numEstPar, constrCode, solcontents);

    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});

    %% write 4 FILE/REFERENCE block if tickbox is ticked
    blockName='FILE/REFERENCE';

    infoType={'DESCRIPTION', 'OUTPUT', 'CONTACT', 'SOFTWARE'};
    info={'TU Wien', 'LEVEL3 output of Vienna VLBI Software', [firstname, ' ', lastname, ' <',email, '>'], 'VieVS - Vienna VLBI Software'};
    frBlock=[infoType; info];

    %write data
    fprintf(fid, '+%s\n', blockName);
    fprintf(fid, ' %-18s%-60s\n', frBlock{:});
    fprintf(fid, '-%s\n', blockName);

    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});
    
    %% write 5 FILE/COMMENT block
    blockName='FILE/COMMENT';
    %comment={'Constraints for EOP are always relative.'};
    
    if strcmp(parameter.vie_init.zhd,'in situ'); zhd=[' b) use ' parameter.data_type ' file content'];
    elseif strcmp(parameter.vie_init.zhd,'no'); zhd=' c) other (no)';
    elseif strcmp(parameter.vie_init.zhd,'vmf3'); zhd=' c) other (VMF3)';
    elseif strcmp(parameter.vie_init.zhd,'vmf1'); zhd=' c) other (VMF1)';
    elseif strcmp(parameter.vie_init.zhd,'gpt3'); zhd=' c) other (GPT3 function)';
    else zhd=' c) other'; 
    end
    
    if strcmp(parameter.vie_init.tp,'in situ'); tp=[' b) use ' parameter.data_type ' file content'];
    elseif strcmp(parameter.vie_init.tp,'vmf3'); tp=' c) other (VMF3)';
    elseif strcmp(parameter.vie_init.tp,'vmf1'); tp=' c) other (VMF1)';
    elseif strcmp(parameter.vie_init.tp,'gpt3'); tp=' c) other (GPT3 function)';
    else tp=' c) other';
    end 
    
    
    if parameter.vie_mod.cts==1; cts='IERS2010';
    else cts='none';
    end
    if parameter.vie_mod.cto==1; cto=parameter.vie_mod.ocm;
    else cto=' none';
    end
    if parameter.vie_mod.cta==1
        if strcmp(parameter.vie_mod.ctam,'s12_cm_noib_leonid.mat')
            cta =' Ponte and Ray (2002) model; Leonid Petrov';% '(http://gemini.gsfc.nasa.gov/aplo/aplo_s1_s2_noib_1.0x1.0deg.nc)';
        else cta=parameter.vie_mod.ctam;
        end
    else cta=' none'; 
    end
    if parameter.vie_mod.cnta==1; cnta= parameter.vie_mod.cntam;
    else cnta=' none';
    end
    if parameter.vie_mod.chl==1; chl= parameter.vie_mod.chlm;
    else chl=' none';
    end
    if parameter.vie_mod.ctp==1
        if strcmp(parameter.vie_mod.ctpm,'linear'); ctp=' linear IERS2003';
        else ctp=' cubic IERS2010'; % strcmp(parameter.vie_mod.ctpm,'cubic'
        end
    else ctp=' none';
    end
    if parameter.vie_mod.linear==1; eopIntrpl=' linear';
    else eopIntrpl=' lagrange';
    end
    
    if parameter.vie_mod.dXdY==0
        eop_dXdY = ' - but nutation offsets dX, dY were set to zero!';
    else eop_dXdY = '';
    end
    if parameter.vie_mod.eophf==1
        hfERP = [' ocean tides: ' parameter.vie_mod.eopoc ','];
        if parameter.vie_mod.lib_pm ==1
            hfERP=[hfERP ' libration in polar motion,'];
        end
        if parameter.vie_mod.lib_ut ==1
            hfERP=[hfERP ' libration in UT1'];
        end  
    else hfERP='';
    end
    
    if parameter.vie_mod.tidalUT==1
        tidalUT = ' UT1S all constituents';
        if parameter.vie_mod.tidalUT35==1
            tidalUT =[''];
            tidalUT = ' UT1R < 35d';
        end
    else
        tidalUT = ' not applied';
    end
    
    if parameter.vie_mod.therm == 1
        therm = ['IVS antenna thermal expansion model of Nothnagel (2008) using scan-wise temperatures from ' parameter.data_type ' file'];
    else
        therm = 'none';
    end
    
    
    if parameter.lsmopt.pw_clk==3
        clk=' b) other (piece-wise linear offsets + one rate + one squared rate per station';
    elseif parameter.lsmopt.pw_clk==2
        clk=' b) other (piece-wise linear offsets + one rate per station';
    elseif parameter.lsmopt.pw_clk==1
        clk=' b) other (piece-wise linear offsets per station';
    else %parameter.lsmopt.pw_clk==0
        clk=' b) other (clocks were not estimated)'; 
    end
    clk_int=[''];
    if parameter.lsmopt.pw_clk~=0
        clk_int=[' with ' num2str(parameter.lsmopt.int_clk) ' min interval and constraints ' num2str(parameter.lsmopt.coef_clk) ' cm after ' num2str(parameter.lsmopt.int_clk) ' min' ];
    end
    if isfield(parameter.lsmopt, 'bdco_nrlist')
        bdco=num2str(size(parameter.lsmopt.bdco_nrlist,1));
    else
        bdco = '0';
    end
    if parameter.lsmopt.pw_zwd==1
        zwd=[' zenith wet delay: ' num2str(parameter.lsmopt.int_zwd) ' min offsets with rel. constr. ' num2str(parameter.lsmopt.coef_zwd) ' cm after '  num2str(parameter.lsmopt.int_zwd) ' min' ];
    else
        zwd=' was not estimated';
    end
    
    % north and east gradients have to have the same parameterization for
    % the sinex output in this version
    if parameter.lsmopt.pw_ngr==1 && parameter.lsmopt.pw_egr==1
        if parameter.lsmopt.constr_rel_ngr==1 && parameter.lsmopt.constr_abs_ngr==0
            tgr=[' ' num2str(parameter.lsmopt.int_ngr) ' min offsets with rel. constr. ' num2str(parameter.lsmopt.coef_rel_ngr) ' cm after ' num2str(parameter.lsmopt.int_ngr) ' min'];
        elseif parameter.lsmopt.constr_rel_ngr==1 && parameter.lsmopt.constr_abs_ngr==1
            tgr=[' ' num2str(parameter.lsmopt.int_ngr) ' min offsets with rel. constr. ' num2str(parameter.lsmopt.coef_rel_ngr) ' cm after ' num2str(parameter.lsmopt.int_ngr) ' min and abs. constr. ' num2str(parameter.lsmopt.coef_abs_ngr) ' cm'];
        elseif parameter.lsmopt.constr_rel_ngr==0 && parameter.lsmopt.constr_abs_ngr==1
            tgr=[' ' num2str(parameter.lsmopt.int_ngr) ' min offsets with abs. constr. ' num2str(parameter.lsmopt.coef_abs_ngr) ' cm'];
        elseif parameter.lsmopt.constr_rel_ngr==0
            tgr=[' ' num2str(parameter.lsmopt.int_ngr) ' min offsets without any constraints'];
        end
    else
        tgr=' were not estimated';
    end
    
    outlr='';
    if parameter.vie_init.rm_outlier==1
        outlr=' + with residuals larger than 5*aposteriori variance of unit weight)'; % usual number which is used for outlier detection (in previous run of VieVS by creating the outlier file)
    end
    
    
    if parameter.lsmopt.est_sourceNNR ==1
        estsou = ' c) other (NNR condition on sources in the a priori catalogue and constraints of 10 mas on all sources)';
    else
        estsou = ' a) fix all positions to their a prioris';
    end
    cutoff = parameter.obs_restrictions.cut_off_elev/pi*180; % rad --> deg
        
    comment={'Analysis description - template by Axel Nothnagel';
        'Version 2.0 (2011-05-26)'; 
        '(http://vlbi.geod.uni-bonn.de/IVS-AC/Docs/Analysis_description.txt)';
        '1. Origin of input data';
        [' c) ' parameter.data_type ' file from IVS Data Center'];
        '3. Origin of meteorological data';
        [' zenith hydrostatic delay: ' zhd];
        [' temperature: ' tp];
        '4. Use numerical weather model';
        ' a) none';
        '5. Setting of mapping functions';
        ' b) other (created by VieVS)';
        '6. Origin of cable cal data';
        [' b) Use ' parameter.data_type ' file content'];
        '7. Application of cable cal flags';
        [' a) all on by default (data from ' parameter.data_type ' file is used as it is)'];
        '8. Selection/settings of geophysical models';
        '8.1 A priori Earth orientation';
        ['- A priori precession/nutation model: ' parameter.vie_mod.nutmod];
        ['- A priori EOP: ' parameter.vie_mod.EOPfile eop_dXdY];
        ['- A priori short-period tidal variations in polar motion and UT1:'];
        hfERP;
        [' EOP interpolation: ' eopIntrpl];
        ['- Tidal UT variations (rg_zont2): ' tidalUT];
        '8.2 A priori geophysical models';
        ['- Solid Earth tides: ' cts];
        ['- Tidal ocean loading: ' cto];
        ['- Tidal atmosphere loading: ' cta];
        ['- Non-tidal atmosphere loading: ' cnta];
        ['- Pole tide: ' ctp];
        ['- Hydrology loading: ' chl];
        ['- Antenna thermal deformation: ' therm];
        ['- Axis offsets: taken from the IVS Analysis Coordinator''s file antenna-info.txt'];
        '8.3 Others'
        ['- A priori TRF: ' parameter.vie_init.trf{2}];
        ['- A priori CRF: ' parameter.vie_init.crf{2}];
        ['- A priori hydrostatic troposphere gradients: ' parameter.vie_mod.apgm_h];
        ['- A priori wet troposphere gradients: ' parameter.vie_mod.apgm_w];
        ['- Hydrostatic mapping functions: ' parameter.vie_mod.mfh];
        ['- Wet mapping functions: ' parameter.vie_mod.mfw];
        '9. Parameterization';
        'Informations about constraints refer only to the SOLUTION/ESTIMATES block and';
        'pre-reduced parameters. The SOLUTION/NORMAL_EQUATION_MATRIX block is without';
        'constraints.';
        '9.1 Clocks';
        clk;
        clk_int;
        [' baseline dependent clock offsets: ' bdco ' baselines'];
        '9.2 Atmospheres';
        ' b) other: piece-wise linear function:';
        zwd;
        ' troposphere horizontal gradients: ';
        tgr;
        '10. Outlier elimination';
        [' c) other (we do not use outliers with quality flag higher than ' num2str(parameter.obs_restrictions.Qlim)];
        outlr;
        '11. Handling of source positions';
        estsou;
        '12. Generation of pre-reduced datum-free normal equations';
        ' b) Program is based on authentic classical least squares adjustment and matrix';
        ' elements of interest are squeezed out by inversion steps.';
        ' vT*P*v = lT*P*l_red - xT*b';
        ' lT*P*l_red = lT*P*l - b2T*inv(N22)*b2';
        '';
        ['Cut-off elevation angle: ' num2str(cutoff) ' degree'];
        };
    
    
    %write data
    fprintf(fid, '+%s\n', blockName);
    fprintf(fid, ' %-79s\n', comment{:});
    fprintf(fid, '-%s\n', blockName);

    
    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});
    
    %% write 9. NUTATION/DATA block if tickbox is ticked
    blockName='NUTATION/DATA';

    nutationModel=parameter.vie_mod.nutmod;
    nutationDescription=' ';

    % write data to file
    fprintf(fid, '+%s\n', blockName);
    fprintf(fid, ' %-8s %-70s\n', nutationModel, nutationDescription);
    fprintf(fid, '-%s\n', blockName);

    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});

    %% write PRECESSION/DATA block if tickbox is ticked
    blockName='PRECESSION/DATA';

    precessionModel=parameter.vie_mod.nutmod;
    precessionDescription=' ';

    % write data to file
    fprintf(fid, '+%s\n', blockName);
    fprintf(fid, ' %-8s %-70s\n', precessionModel, precessionDescription);
    fprintf(fid, '-%s\n', blockName);

    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});

    %% write 11. SOURCE/ID block if tickbox is ticked
    blockName='SOURCE/ID';
    
    % write blockstart line to file
    fprintf(fid, '+%s\n', blockName);
    
    useIVSnames = false;
    % write description line
    if useIVSnames
        fprintf(fid, '%s\n', '*Code IVS_des ICRF_designator  Comments');
    else
        fprintf(fid, '%s\n', '*Code IERS_des ICRF_designator  Comments');
    end
    
    % for all sources
    for k=1:numSou
        if useIVSnames
            if isempty(sources.q(k).IVSname)
                fprintf(fid, ' %04s %-8s %-16s %-68s\n', num2str(k), sources.q(k).IERSname, sources.q(k).ICRFdes, num2str(sources.q(k).numobs));
            else
                fprintf(fid, ' %04s %-8s %-16s %-68s\n', num2str(k), sources.q(k).IVSname, sources.q(k).ICRFdes, num2str(sources.q(k).numobs));
            end
        else
            fprintf(fid, ' %04s %-8s %-16s %-68s\n', num2str(k), sources.q(k).IERSname, sources.q(k).ICRFdes, num2str(sources.q(k).numobs));
        end
        sources.q(k).siteCode=num2str(k);
    end
    
    % write blockend line to file
    fprintf(fid, '-%s\n', blockName);
    
    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});

    %% write 12. SITE/ID block if tickbox is ticked
    blockName='SITE/ID';
    writeFormat=' %4s %2s %9s %1s %-22s %3.0f %02.0f %04.1f %3.0f %02.0f %04.1f %7.1f\n';

    % get lat and lon from extimated x,y,z
    [approxLat,approxLon,approxH]=xyz2ell([[antenna.x]', [antenna.y]', [antenna.z]']);

    % rad -> �
    approxLat=approxLat*180/pi;
    approxLon=approxLon*180/pi;

    % [-180, 180] -> [0, 360]
    approxLon(approxLon<0)=approxLon(approxLon<0)+360;

    % � -> dms
    %dmsLat=degrees2dms(approxLat); not in basic Matlab
    degrees = fix(approxLat);
    minutes = abs(fix((approxLat-degrees).*60));
    seconds = (abs(approxLat-degrees).*60 - minutes).*60;
    dmsLat = [degrees,minutes,seconds];

    %dmsLon=degrees2dms(approxLon); not in basic Matlab
    degrees = fix(approxLon);
    minutes = abs(fix((approxLon-degrees).*60));
    seconds = (abs(approxLon-degrees).*60 - minutes).*60;
    dmsLon = [degrees,minutes,seconds];
    clearvars degrees minutes seconds;

    % write blockstart line to file
    fprintf(fid, '+%s\n', blockName);

    % write description line
    fprintf(fid, '%s\n', '*CODE PT DOMES____ T STATION_DESCRIPTION___ APPROX_LON_ APPROX_LAT_ APP_H__');
    % write data
    for k=1:size(antenna, 2)  
        %                                  -4s                     %-2s                 %9s             %1s            %-22s               %3.0f         %02.0f      %04.1f        %3.0f        %02.0f         %04.1f    7.1f\n';
        fprintf(fid, writeFormat, num2str(antenna(k).siteCode), antenna(k).pointCode, antenna(k).domes(:), obsCode, antenna(k).description, dmsLon(k,1), dmsLon(k,2), dmsLon(k,3), dmsLat(k,1), dmsLat(k,2), dmsLat(k,3), approxH(k));
        %fprintf(fid, '\n');
    end

    % write blockend line to file
    fprintf(fid, '-%s\n', blockName);

    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});

    %% write 17. SITE/ECCENTRICITY block if tickbox is ticked
    blockName='SITE/ECCENTRICITY';
    writeFormat=' %4.0f %2s %4.0f %1s %2s:%03.0f:%05.0f %2s:%03.0f:%05.0f %3s %8.4f %8.4f %8.4f\n';
    % write blockstart line to file
    fprintf(fid, '+%s\n', blockName);

    % write info line
    fprintf(fid, '*Code PT SBIN T Data_Start__ Data_End____ typ Up/X____ North/Y_ East/Z__\n');
    
    % get number of scans
    %numScan=size(scan,2);
    
    % for all stations
    for k=1:numStat
        %curPtCode='A';
        soln=1; % solution number
        
        time=mjd2yydoysecod([antenna(k).firstObsMjd, antenna(k).lastObsMjd, (antenna(k).firstObsMjd+antenna(k).lastObsMjd)/2]);
        timeYrStr=num2str(time(:,1));
        
        % write data
        fprintf(fid, writeFormat, antenna(k).siteCode, antenna(k).pointCode, soln, obsCode, timeYrStr(7:3:end), time(1,2), time(1,3), timeYrStr(8:3:end), time(2,2), time(2,3), antenna(k).ecctype, antenna(k).c_ecc(1), antenna(k).c_ecc(2), antenna(k).c_ecc(3));
    end

    % write blockend line to file
    fprintf(fid, '-%s\n', blockName);

    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});

    %% write 20. SOLUTION/EPOCHS block if tickbox is ticked
    blockName='SOLUTION/EPOCHS';
    writeFormat=' %4.0f %2s %4.0f %1s %2s:%03.0f:%05.0f %2s:%03.0f:%05.0f %2s:%03.0f:%05.0f\n';

    % write blockstart line to file
    fprintf(fid, '+%s\n', blockName);

    % write info line
    fprintf(fid, '*Code PT SOLN T Data_start__ Data_end____ Mean_epoch__\n');

    % write data
    % for each station: one line (since we have only one combination of point
    % code: A / solution ID: 1 / ObservationCode: R=VLBI)
    
    for k=1:numStat
        % get time of first, last and mean mjd of observation
        time=mjd2yydoysecod([antenna(k).firstObsMjd, antenna(k).lastObsMjd, (antenna(k).firstObsMjd+antenna(k).lastObsMjd)/2]);

        % year: float 2008 -> string 08 (again for first, last and median scan)
        timeYrStr=num2str(time(:,1));
        solID=1;
        %obsCode='R'; % already defined before

        fprintf(fid, writeFormat, antenna(k).siteCode, antenna(k).pointCode, solID, obsCode, timeYrStr(7:3:end), time(1,2), time(1,3), timeYrStr(8:3:end), time(2,2), time(2,3), timeYrStr(9:3:end), time(3,2), time(3,3));
        clear time
    end

    % write blockend line to file
    fprintf(fid, '-%s\n', blockName);

    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});
    
    
    %% write 22 SOLUTION/STATISTICS block
    % define name of block and output format
    blockName='SOLUTION/STATISTICS';
    writeFormat=' %-30s %22.15f\n';
    
   
    %if estimation was carried out
    if parameter.lsmopt.est_singleses==1
        infoType={'NUMBER OF OBSERVATIONS', 'NUMBER OF UNKNOWNS', 'SQUARE SUM OF RESIDUALS (VTPV)', 'VARIANCE FACTOR' , 'WEIGHTED SQUARE SUM OF O-C'}; % hana
        info={col_sinex.nr_obs, col_sinex.nr_unknowns, col_sinex.vTPv, col_sinex.varfac, col_sinex.lTPlreduc}; % hana

    else
        infoType={'NUMBER OF OBSERVATIONS', 'NUMBER OF UNKNOWNS', 'WEIGHTED SQUARE SUM OF O-C'};
        info={col_sinex.nr_obs, col_sinex.nr_unknowns, col_sinex.lTPlreduc}; % hana
    end
    frBlock=[infoType; info];

    %write data
    fprintf(fid, '+%s\n', blockName);
    fprintf(fid, writeFormat, frBlock{:});
    fprintf(fid, '-%s\n', blockName);

    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});
    
    
    
 %% Compute (interpolate) apriori hydrostatic zenith delay   
 if outsnx.zwd==1
    for k=1:numStat
        clear mjdZdry
        clear mjdZwet
        % Interpolate (linear) the apriori (hydrostatic) delay from scan. to
        % the estimated mjds
        iscan1=1;
        for iscan =1:length(scan)
            if length(scan(iscan).stat)>=k
                if ~isempty(scan(iscan).stat(k).zdry)
                    mjdZdry(iscan1,:) = [scan(iscan).mjd scan(iscan).stat(k).zdry];
                    mjdZwet(iscan1,:) = [scan(iscan).mjd scan(iscan).stat(k).zwet];
                    iscan1=iscan1+1;
                end
            end
        end

        % check if there are not two scans with the same date
        for inu=1:size(mjdZdry,1)-1
            nul=mjdZdry(inu,1)- mjdZdry(inu+1,1);
            if nul==0
                 mjdZdry(inu+1,1) = mjdZdry(inu+1,1) + 1e-10;
                 mjdZwet(inu+1,1) = mjdZwet(inu+1,1) + 1e-10;
            end
        end

        
        for imjd = 1:length(col_sinex.zwd(k).mjd)
            antenna(k).zdry(imjd) = interp1(mjdZdry(:,1),mjdZdry(:,2),col_sinex.zwd(k).mjd(imjd),'linear','extrap');
            antenna(k).zwet(imjd) = interp1(mjdZwet(:,1),mjdZwet(:,2),col_sinex.zwd(k).mjd(imjd),'linear','extrap');
        end
    end
 end
    
   
        
    %% write 23 SOLUTION/ESTIMATES block if tickbox is ticked
    
    % if session lasts one (0.9) day
    if sessionTimeConsumption>0.9
        midmjd=floor(max([antenna.lastObsMjd])-0.5)+0.5;        %hana
        aprDate=mjd2yydoysecod(midmjd);
        aprDateYrStr=num2str(aprDate(1));
    else %session lasts less than one (0.9) day, eg. one hour
        midmjd=floor(max([antenna.lastObsMjd]))+0.5;          % hana jun11
        aprDate=mjd2yydoysecod(midmjd);
        aprDateYrStr=num2str(aprDate(1));
    end
    
    
    % if option was chosen in gui
    if parameter.lsmopt.est_singleses==1
        
        % define name of block and output format
        blockName='SOLUTION/ESTIMATES';
        writeFormat=' %5.0f %-6s %04s %2s %4s %2s:%03.0f:%05.0f %-4s %1.0f %21s %11s\n';
        %writeFormat=' %5.0f %-6s %4s %2s %4s %2s:%03.0f:%05.0f %-4s %1.0f %21.13d %11.4d\n';
        if ispc, formatEstVal='%22.14e'; formatStDev='%12.5e'; else
                 formatEstVal='%21.14e'; formatStDev='%11.5e';
        end
        constrEstCoord=1;
        curIndex=1;

        % write blockstart line to file
        fprintf(fid, '+%s\n', blockName);

        % write info line
        fprintf(fid, '*Index Type__ Code Pt Soln Ref_Epoch___ Unit S Estimated_Value______ Std_Dev____\n');

        % write data  

        if outsnx.xyz==1
        % write x,y,z coordinates
            cpsd_all = cPostSeismDeform(midmjd,antenna);
            for k=1:numStat
                % Solution ID   
                soln=num2str(1);


                % calculate apriori values
                %                  coordx   +     vx      * (   mean scan mjd - itrf epoch mjd )/ diff(mjd)->years                                
                antenna(k).aprX=antenna(k).x+antenna(k).vx*(midmjd-antenna(k).epoch)/365.25 + cpsd_all(1,1,k);
                antenna(k).aprY=antenna(k).y+antenna(k).vy*(midmjd-antenna(k).epoch)/365.25 + cpsd_all(2,1,k);
                antenna(k).aprZ=antenna(k).z+antenna(k).vz*(midmjd-antenna(k).epoch)/365.25 + cpsd_all(3,1,k);


                % calculate total estimated values
                antenna(k).totX=antenna(k).aprX+x_.coorx(k).val/100;
                antenna(k).totY=antenna(k).aprY+x_.coory(k).val/100;
                antenna(k).totZ=antenna(k).aprZ+x_.coorz(k).val/100;


                % preparation of values for format e+02 instead of e+002
                totXstring=sprintf(formatEstVal, antenna(k).totX);
                if ispc, totXstring = strrep(totXstring, 'e+0', 'e+'); totXstring = strrep(totXstring, 'e-0', 'e-');   end
                totYstring=sprintf(formatEstVal, antenna(k).totY);
                if ispc, totYstring = strrep(totYstring, 'e+0', 'e+'); totYstring = strrep(totYstring, 'e-0', 'e-'); end
                totZstring=sprintf(formatEstVal, antenna(k).totZ);
                if ispc, totZstring = strrep(totZstring, 'e+0', 'e+'); totZstring = strrep(totZstring, 'e-0', 'e-'); end

                totXstdString=sprintf(formatStDev, x_.coorx(k).mx/100);
                if ispc, totXstdString = strrep(totXstdString, 'e+0', 'e+'); totXstdString = strrep(totXstdString, 'e-0', 'e-');   end
                totYstdString=sprintf(formatStDev, x_.coory(k).mx/100);
                if ispc, totYstdString = strrep(totYstdString, 'e+0', 'e+'); totYstdString = strrep(totYstdString, 'e-0', 'e-');  end
                totZstdString=sprintf(formatStDev, x_.coorz(k).mx/100);
                if ispc, totZstdString = strrep(totZstdString, 'e+0', 'e+'); totZstdString = strrep(totZstdString, 'e-0', 'e-');  end

                %                   Code  PT  1   yr:time(2) time3  m     1   estimate  stdev      
    %             fprintf(fid, format, curIndex, 'STAX',   num2str(antenna(k).siteCode), antenna(k).pointCode, soln, statTimeYrStr(7:3:end), statTime(1,2), statTime(1,3), 'm', constrEstCoord, totXstring, totXstdString);
    %             fprintf(fid, format, curIndex+1, 'STAY', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, statTimeYrStr(8:3:end), statTime(2,2), statTime(2,3), 'm', constrEstCoord, totYstring, totYstdString);
    %             fprintf(fid, format, curIndex+2, 'STAZ', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, statTimeYrStr(9:3:end), statTime(3,2), statTime(3,3), 'm', constrEstCoord, totZstring, totZstdString);

                fprintf(fid, writeFormat, curIndex,  'STAX', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm', constrEstCoord, totXstring, totXstdString);
                fprintf(fid, writeFormat, curIndex+1,'STAY', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm', constrEstCoord, totYstring, totYstdString);
                fprintf(fid, writeFormat, curIndex+2,'STAZ', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm', constrEstCoord, totZstring, totZstdString);

             curIndex=curIndex+3;                                   
            end
        end
        
        if outsnx.eop==1
            % get eops --> moved to above because we need number of total estimated
            % parameters in header line!

            % write EOPs 
            soln=num2str(1);

            % write xpol

            % constraint
            if parameter.lsmopt.xpol.constrain==1 && parameter.lsmopt.xpol.coef<1
                constrEstXpol=0;
            elseif parameter.lsmopt.xpol.constrain==0
                constrEstXpol=2;
            else %parameter.lsmopt.xpol.constrain==1 && parameter.lsmopt.xpol.coef>1
                constrEstXpol=1;
            end

            for k=1:numXpol
                curLine=eop_matrix(:,1)==col_sinex.mjd_xp(k); % index of current xpol in eop_matrix

                % make correct format
                totXpo=sprintf(formatEstVal, eop_matrix(curLine,2));
                if ispc, totXpo = strrep(totXpo, 'e+0', 'e+'); totXpo = strrep(totXpo, 'e-0', 'e-');  end
                totXpoStdDev=sprintf(formatStDev, eop_matrix(curLine,17));
                if ispc, totXpoStdDev = strrep(totXpoStdDev, 'e+0', 'e+'); totXpoStdDev = strrep(totXpoStdDev, 'e-0', 'e-'); end

                fprintf(fid, writeFormat, curIndex, 'XPO', '----', '--', soln, eopDatesYrStr(curLine,3:4), eopDates(curLine,2), eopDates(curLine,3), 'mas', constrEstXpol, totXpo, totXpoStdDev);  
                curIndex=curIndex+1;
            end

            % write ypol

            % constraint
            if parameter.lsmopt.ypol.constrain==1 && parameter.lsmopt.ypol.coef<1
                constrEstYpol=0;
            elseif parameter.lsmopt.ypol.constrain==0
                constrEstYpol=2;
            else %parameter.lsmopt.ypol.constrain==1 && parameter.lsmopt.ypol.coef>1
                constrEstYpol=1;
            end

            for k=1:numYpol
                curLine=eop_matrix(:,1)==col_sinex.mjd_yp(k); % index of current xpol in eop_matrix
                % make correct format
                totYpo=sprintf(formatEstVal, eop_matrix(curLine,3));
                if ispc, totYpo = strrep(totYpo, 'e+0', 'e+'); totYpo = strrep(totYpo, 'e-0', 'e-');   end
                totYpoStdDev=sprintf(formatStDev, eop_matrix(curLine,18));
                if ispc, totYpoStdDev = strrep(totYpoStdDev, 'e+0', 'e+'); totYpoStdDev = strrep(totYpoStdDev, 'e-0', 'e-');  end

                fprintf(fid, writeFormat, curIndex, 'YPO', '----', '--', soln, eopDatesYrStr(curLine,3:4), eopDates(curLine,2), eopDates(curLine,3), 'mas', constrEstYpol, totYpo, totYpoStdDev);  
                curIndex=curIndex+1;
            end

            % write dUT1

            % constraint
            if parameter.lsmopt.dut1.constrain==1 && parameter.lsmopt.dut1.coef<1
                constrEstDut1=0;
            elseif parameter.lsmopt.dut1.constrain==0
                constrEstDut1=2;
            else %parameter.lsmopt.dut1.constrain==1 && parameter.lsmopt.dut1.coef>1
                constrEstDut1=1;
            end

            for k=1:numDut1
                curLine=eop_matrix(:,1)==col_sinex.mjd_dut1(k); % index of current xpol in eop_matrix
                % make correct format
                totDut1=sprintf(formatEstVal, eop_matrix(curLine,4));
                if ispc, totDut1 = strrep(totDut1, 'e+0', 'e+'); totDut1 = strrep(totDut1, 'e-0', 'e-');   end
                totDut1StdDev=sprintf(formatStDev, eop_matrix(curLine,19));
                if ispc, totDut1StdDev = strrep(totDut1StdDev, 'e+0', 'e+'); totDut1StdDev = strrep(totDut1StdDev, 'e-0', 'e-');   end

                fprintf(fid, writeFormat, curIndex, 'UT', '----', '--', soln, eopDatesYrStr(curLine,3:4), eopDates(curLine,2), eopDates(curLine,3), 'ms', constrEstDut1, totDut1, totDut1StdDev);  
                curIndex=curIndex+1;
            end

            % write NUTx

            % constraint
            if parameter.lsmopt.nutdx.constrain==1 && parameter.lsmopt.nutdx.coef<1
                constrEstNutx=0;
            elseif parameter.lsmopt.nutdx.constrain==0
                constrEstNutx=2;
            else %parameter.lsmopt.nutdx.constrain==1 && parameter.lsmopt.nutdx.coef>1
                constrEstNutx=1;
            end

            for k=1:numXnut
                curLine=eop_matrix(:,1)==col_sinex.mjd_dX(k); % index of current xpol in eop_matrix
                % make correct format
                totdX=sprintf(formatEstVal, eop_matrix(curLine,5));
                if ispc, totdX = strrep(totdX, 'e+0', 'e+'); totdX = strrep(totdX, 'e-0', 'e-');   end
                totdXStdDev=sprintf(formatStDev, eop_matrix(curLine,20));
                if ispc, totdXStdDev = strrep(totdXStdDev, 'e+0', 'e+'); totdXStdDev = strrep(totdXStdDev, 'e-0', 'e-');   end

                fprintf(fid, writeFormat, curIndex, 'NUT_X', '----', '--', soln, eopDatesYrStr(curLine,3:4), eopDates(curLine,2), eopDates(curLine,3), 'mas', constrEstNutx, totdX, totdXStdDev);  
                curIndex=curIndex+1;
            end

            %question: nut: 3x at midnight, xpol(and rest): between 17:00  and 17:00 every hour (or so)
            %what are the estimates from eop_matrix in the column of nutation (do they correspond to the three midnight mjds or the hourly epochs)?
            % probable solution: make unique (all mjds) to get eop estimates
            % (in eop_matrix) so that all existing eops are interpolated at all
            % availoable epochs (the ones that will not be used are nnecessary
            % anyhow...
            % write NUTy


            % write NUTy

            % constraint
            if parameter.lsmopt.nutdy.constrain==1 && parameter.lsmopt.nutdy.coef<1
                constrEstNuty=0;
            elseif parameter.lsmopt.nutdy.constrain==0
                constrEstNuty=2;
            else %parameter.lsmopt.nutdy.constrain==1 && parameter.lsmopt.nutdy.coef>1
                constrEstNuty=1;
            end

            for k=1:numYnut
                curLine=eop_matrix(:,1)==col_sinex.mjd_dY(k); % index of current xpol in eop_matrix
                % make correct format
                totdY=sprintf(formatEstVal, eop_matrix(curLine,6));
                if ispc, totdY = strrep(totdY, 'e+0', 'e+'); totdY = strrep(totdY, 'e-0', 'e-');   end
                totdYStdDev=sprintf(formatStDev, eop_matrix(curLine,21));
                if ispc, totdYStdDev = strrep(totdYStdDev, 'e+0', 'e+'); totdYStdDev = strrep(totdYStdDev, 'e-0', 'e-');  end

                fprintf(fid, writeFormat, curIndex, 'NUT_Y', '----', '--', soln, eopDatesYrStr(curLine,3:4), eopDates(curLine,2), eopDates(curLine,3), 'mas', constrEstNuty, totdY, totdYStdDev);  
                curIndex=curIndex+1;
            end
        end
        
        
        if outsnx.sou==1
            for k=1:numSou
                % Solution ID   
                soln=num2str(1);
                
                % calculate total estimated values
                if ~isempty(x_.soura(k).val)
                    sources.q(k).totRa = sources.q(k).ra2000 +x_.soura(k).val*15/1000/3600/180*pi;
                    sources.q(k).totDe = sources.q(k).de2000 +x_.soude(k).val/1000/3600/180*pi;

                    sources.q(k).totRaStd = x_.soura(k).mx*15/1000/3600/180*pi;
                    sources.q(k).totDeStd = x_.soude(k).mx/1000/3600/180*pi;
                else             
                    sources.q(k).totRa = 0;
                    sources.q(k).totDe = 0;
                    sources.q(k).totRaStd = 0;
                    sources.q(k).totDeStd = 0;
                end
                
                if x_.soura(k).inNNR == 1
                    sources.q(k).constrEstSou=1;    % 1 - NNR, 2 - not in the NNR
                else
                    sources.q(k).constrEstSou=2;    % 1 - NNR, 2 - not in the NNR
                end
                
                % preparation of values for format e+02 instead of e+002
                totRastring=sprintf(formatEstVal, sources.q(k).totRa);
                if ispc, totRastring = strrep(totRastring, 'e+0', 'e+'); totRastring = strrep(totRastring, 'e-0', 'e-');   end
                totDestring=sprintf(formatEstVal, sources.q(k).totDe);
                if ispc, totDestring = strrep(totDestring, 'e+0', 'e+'); totDestring = strrep(totDestring, 'e-0', 'e-'); end

                totRastdString=sprintf(formatStDev, sources.q(k).totRaStd);
                if ispc, totRastdString = strrep(totRastdString, 'e+0', 'e+'); totRastdString = strrep(totRastdString, 'e-0', 'e-');   end
                totDestdString=sprintf(formatStDev,sources.q(k).totDeStd);
                if ispc, totDestdString = strrep(totDestdString, 'e+0', 'e+'); totDestdString = strrep(totDestdString, 'e-0', 'e-');  end

                fprintf(fid, writeFormat, curIndex,  'RS_RA', num2str(sources.q(k).siteCode), '--', soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'rad', sources.q(k).constrEstSou, totRastring, totRastdString);
                fprintf(fid, writeFormat, curIndex+1,'RS_DE', num2str(sources.q(k).siteCode), '--', soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'rad', sources.q(k).constrEstSou, totDestring, totDestdString);

                curIndex=curIndex+2;
            end
        end
        
        
        % write troposphere zenith total delay
        if outsnx.zwd==1
            if parameter.lsmopt.constr_zwd==1 && parameter.lsmopt.coef_zwd<0.1
                constrEstZwd=0;
            elseif parameter.lsmopt.constr_zwd==0
                constrEstZwd=2;
            else %parameter.lsmopt.constr_zwd==1 && parameter.lsmopt.coef_zwd>0.1
                constrEstZwd=1;
            end

            for k=1:numStat
                clear zwdTime zwdTimeYrStr totZtd1 totZtdStdDev1    
                zwdTime=mjd2yydoysecod(col_sinex.zwd(k).mjd);
                zwdTimeYrStr=num2str(zwdTime(:,1));
                numZwd=length(col_sinex.zwd(k).col);
                
                aprZtd = antenna(k).zdry + antenna(k).zwet;
                
                totZtd1 = aprZtd +  x_.zwd(k).val./100; %m
                totZtdStdDev1 = x_.zwd(k).mx./100; %m
                for j=1:numZwd 
                    totZtd=sprintf(formatEstVal, totZtd1(j));
                    if ispc, totZtd = strrep(totZtd, 'e+0', 'e+'); totZtd = strrep(totZtd, 'e-0', 'e-');   end
                    totZtdStdDev=sprintf(formatStDev, totZtdStdDev1(j));
                    if ispc, totZtdStdDev = strrep(totZtdStdDev, 'e+0', 'e+'); totZtdStdDev = strrep(totZtdStdDev, 'e-0', 'e-');  end

                    fprintf(fid, writeFormat, curIndex, 'TROTOT', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, zwdTimeYrStr(j,3:end), zwdTime(j,2), zwdTime(j,3), 'm', constrEstZwd, totZtd, totZtdStdDev);
                    curIndex=curIndex+1;
                end
            end
        end
        
        % write troposphere north and east gradients
        if outsnx.tgr==1
            if parameter.lsmopt.constr_abs_ngr==1 && parameter.lsmopt.coef_abs_ngr<0.6
                constrEstNgr=0;
            elseif parameter.lsmopt.constr_abs_ngr==0
                constrEstNgr=2;
            else % loose
                constrEstNgr=1;
            end

            if parameter.lsmopt.constr_abs_egr==1 && parameter.lsmopt.coef_abs_egr<0.6
                constrEstEgr=0;
            elseif parameter.lsmopt.constr_abs_egr==0
                constrEstEgr=2;
            else % loose
                constrEstEgr=1;
            end


            for k=1:numStat
                clear ngrTime ngrTimeYrStr
                ngrTime=mjd2yydoysecod(col_sinex.ngr(k).mjd);
                ngrTimeYrStr=num2str(ngrTime(:,1));
                numNgr=length(col_sinex.ngr(k).col);

                aprNgr=antenna(k).apriori_ngr/1000;
                totNgr1=aprNgr + x_.ngr(k).val./100;
                totNgrStdDev1= x_.ngr(k).mx./100;

                aprEgr=antenna(k).apriori_egr/1000;
                totEgr1=aprEgr + x_.egr(k).val./100;
                totEgrStdDev1= x_.egr(k).mx./100;

                for j=1:numNgr 
                    totNgr=sprintf(formatEstVal, totNgr1(j));
                    if ispc, totNgr = strrep(totNgr, 'e+0', 'e+'); totNgr = strrep(totNgr, 'e-0', 'e-');   end
                    totNgrStdDev=sprintf(formatStDev, totNgrStdDev1(j));
                    if ispc, totNgrStdDev = strrep(totNgrStdDev, 'e+0', 'e+'); totNgrStdDev = strrep(totNgrStdDev, 'e-0', 'e-');  end

                    totEgr=sprintf(formatEstVal, totEgr1(j));
                    if ispc, totEgr = strrep(totEgr, 'e+0', 'e+'); totEgr = strrep(totEgr, 'e-0', 'e-');   end
                    totEgrStdDev=sprintf(formatStDev, totEgrStdDev1(j));
                    if ispc, totEgrStdDev = strrep(totEgrStdDev, 'e+0', 'e+'); totEgrStdDev = strrep(totEgrStdDev, 'e-0', 'e-');  end

                    fprintf(fid, writeFormat, curIndex, 'TGNTOT', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, ngrTimeYrStr(j,3:end), ngrTime(j,2), ngrTime(j,3), 'm', constrEstNgr, totNgr, totNgrStdDev);
                    fprintf(fid, writeFormat, curIndex+1, 'TGETOT', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, ngrTimeYrStr(j,3:end), ngrTime(j,2), ngrTime(j,3), 'm', constrEstEgr, totEgr, totEgrStdDev);
                    curIndex=curIndex+2;
                end
            end
        end
        
        
        % write blockend line to file
        fprintf(fid, '-%s\n', blockName);

        % write headerspace
        fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});
    
    else %if writing estimates was not chosen in gui
        % only calculate a priori values for a priori block
        cpsd_all = cPostSeismDeform(midmjd,antenna);
        for k=1:numStat % hana
            %                  coordx   +     vx      * (  mean scan mjd - itrf epoch mjd )/ diff(mjd)->years                                
                antenna(k).aprX=antenna(k).x+antenna(k).vx*(midmjd-antenna(k).epoch)/365.25 + cpsd_all(1,1,k);
                antenna(k).aprY=antenna(k).y+antenna(k).vy*(midmjd-antenna(k).epoch)/365.25 + cpsd_all(2,1,k);
                antenna(k).aprZ=antenna(k).z+antenna(k).vz*(midmjd-antenna(k).epoch)/365.25 + cpsd_all(3,1,k);
        end

    end

    %% write 24. SOLUTION/APRIORI block if tickbox is ticked
    blockName='SOLUTION/APRIORI';
    writeFormat=' %5.0f %-6s %04s %2s %4s %2s:%03.0f:%05.0f %-4s %1.0f %21s %11s\n';
    if ispc, formatAprVal='%22.14e'; formatStDev='%12.5e'; else
             formatAprVal='%21.14e'; formatStDev='%11.5e';
    end

    curIndex=1;
    constrApr=2;
    
    % write blockstart line to file
    fprintf(fid, '+%s\n', blockName);

    % write info line
    fprintf(fid, '*Index Type__ CODE PT Soln Ref_epoch___ Unit S Apriori_value________ Std_Dev____\n');

    % write data
    if outsnx.xyz==1
        % x,y,z coordinates
        for k=1:numStat
            % solution ID
            soln=num2str(1);

            antenna(k).aprX_sigma=0;
            antenna(k).aprY_sigma=0;
            antenna(k).aprZ_sigma=0;

            % get format e+02
            aprX=sprintf(formatAprVal, antenna(k).aprX);
            if ispc, aprX = strrep(aprX, 'e+0', 'e+'); aprX = strrep(aprX, 'e-0', 'e-');   end
            aprY=sprintf(formatAprVal, antenna(k).aprY);
            if ispc, aprY = strrep(aprY, 'e+0', 'e+'); strrep(aprY, 'e-0', 'e-');  end
            aprZ=sprintf(formatAprVal, antenna(k).aprZ);
            if ispc, aprZ = strrep(aprZ, 'e+0', 'e+'); strrep(aprZ, 'e-0', 'e-');  end

            aprX_sigma=sprintf(formatStDev, antenna(k).aprX_sigma);
            if ispc, aprX_sigma = strrep(aprX_sigma, 'e+0', 'e+'); aprX_sigma = strrep(aprX_sigma, 'e-0', 'e-');  end
            aprY_sigma=sprintf(formatStDev, antenna(k).aprY_sigma);
            if ispc, aprY_sigma = strrep(aprY_sigma, 'e+0', 'e+'); aprY_sigma = strrep(aprY_sigma, 'e-0', 'e-');   end
            aprZ_sigma=sprintf(formatStDev, antenna(k).aprZ_sigma);
            if ispc, aprZ_sigma = strrep(aprZ_sigma, 'e+0', 'e+'); aprZ_sigma = strrep(aprZ_sigma, 'e-0', 'e-');  end

            % write one line for each x,y,z
            fprintf(fid, writeFormat, curIndex,  'STAX', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm', constrApr, aprX, aprX_sigma);
            fprintf(fid, writeFormat, curIndex+1,'STAY', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm', constrApr, aprY, aprY_sigma);
            fprintf(fid, writeFormat, curIndex+2,'STAZ', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm', constrApr, aprZ, aprZ_sigma);
            curIndex=curIndex+3;
        end
    end
    
    if outsnx.eop==1
        % get apriori dates
        aprDates=mjd2yydoysecod(eop_matrix(:,1));
        aprDatesYrStr=num2str(aprDates(:,1));

        % write apriori EOPs 
        numEOPs=size(eop_matrix,1);

        % write xpol
        for k=1:numXpol
            curLine=eop_matrix(:,1)==col_sinex.mjd_xp(k); % index of current xpol in eop_matrix

            % make correct format
            aprXpo=sprintf(formatAprVal, eop_matrix(curLine,7));
            if ispc, aprXpo = strrep(aprXpo, 'e+0', 'e+'); aprXpo = strrep(aprXpo, 'e-0', 'e-');  end
            aprXpoStdDev=sprintf(formatStDev, 0);
            if ispc, aprXpoStdDev = strrep(aprXpoStdDev, 'e+0', 'e+'); aprXpoStdDev = strrep(aprXpoStdDev, 'e-0', 'e-');   end

            fprintf(fid, writeFormat, curIndex, 'XPO', '----', '--', soln, aprDatesYrStr(curLine,3:4), aprDates(curLine,2), aprDates(curLine,3), 'mas', constrApr, aprXpo, aprXpoStdDev);  
            curIndex=curIndex+1;        
        end

        % write ypol
        for k=1:numYpol
            curLine=eop_matrix(:,1)==col_sinex.mjd_yp(k); % index of current xpol in eop_matrix
            % make correct format
            aprYpo=sprintf(formatAprVal, eop_matrix(curLine,8));
            if ispc, aprYpo = strrep(aprYpo, 'e+0', 'e+'); aprYpo = strrep(aprYpo, 'e-0', 'e-');  end
            aprYpoStdDev=sprintf(formatStDev, 0);
            if ispc, aprYpoStdDev = strrep(aprYpoStdDev, 'e+0', 'e+'); aprYpoStdDev = strrep(aprYpoStdDev, 'e-0', 'e-');  end

            fprintf(fid, writeFormat, curIndex, 'YPO', '----', '--', soln, aprDatesYrStr(curLine,3:4), aprDates(curLine,2), aprDates(curLine,3), 'mas', constrApr, aprYpo, aprYpoStdDev);  
            curIndex=curIndex+1;
        end

        % write dUT1
        for k=1:numDut1
            curLine=eop_matrix(:,1)==col_sinex.mjd_dut1(k); % index of current xpol in eop_matrix
            % make correct format
            aprDut1=sprintf(formatAprVal, eop_matrix(curLine,9));
            if ispc, aprDut1 = strrep(aprDut1, 'e+0', 'e+'); aprDut1 = strrep(aprDut1, 'e-0', 'e-');   end
            aprDut1StdDev=sprintf(formatStDev, 0);
            if ispc, aprDut1StdDev = strrep(aprDut1StdDev, 'e+0', 'e+'); aprDut1StdDev = strrep(aprDut1StdDev, 'e-0', 'e-');  end

            fprintf(fid, writeFormat, curIndex, 'UT', '----', '--', soln, aprDatesYrStr(curLine,3:4), aprDates(curLine,2), aprDates(curLine,3), 'ms', constrApr, aprDut1, aprDut1StdDev);  
            curIndex=curIndex+1;
        end

        % write NUTx
        for k=1:numXnut
            curLine=eop_matrix(:,1)==col_sinex.mjd_dX(k); % index of current xpol in eop_matrix
            % make correct format
            aprdX=sprintf(formatAprVal, eop_matrix(curLine,10));
            if ispc, aprdX = strrep(aprdX, 'e+0', 'e+'); aprdX = strrep(aprdX, 'e-0', 'e-');   end
            aprdXStdDev=sprintf(formatStDev, 0);
            if ispc, aprdXStdDev = strrep(aprdXStdDev, 'e+0', 'e+'); aprdXStdDev = strrep(aprdXStdDev, 'e-0', 'e-');  end

            fprintf(fid, writeFormat, curIndex, 'NUT_X', '----', '--', soln, aprDatesYrStr(curLine,3:4), aprDates(curLine,2), aprDates(curLine,3), 'mas', constrApr, aprdX, aprdXStdDev);  
            curIndex=curIndex+1;
        end

        % write NUTy
        for k=1:numYnut
            curLine=eop_matrix(:,1)==col_sinex.mjd_dY(k); % index of current xpol in eop_matrix
            % make correct format
            aprdY=sprintf(formatAprVal, eop_matrix(curLine,11));
            if ispc, aprdY = strrep(aprdY, 'e+0', 'e+'); aprdY = strrep(aprdY, 'e-0', 'e-');  end
            aprdYStdDev=sprintf(formatStDev, 0);
            if ispc, aprdYStdDev = strrep(aprdYStdDev, 'e+0', 'e+'); aprdYStdDev = strrep(aprdYStdDev, 'e-0', 'e-');  end

            fprintf(fid, writeFormat, curIndex, 'NUT_Y', '----', '--', soln, aprDatesYrStr(curLine,3:4), aprDates(curLine,2), aprDates(curLine,3), 'mas', constrApr, aprdY, aprdYStdDev);  
            curIndex=curIndex+1;
        end
    end
    
    
    % ra, de sources
    if outsnx.sou==1    
        for k=1:numSou
            % solution ID
            soln=num2str(1);


            sources.q(k).aprra_sigma=0;
            sources.q(k).aprde_sigma=0;

            % get format e+02
            aprra=sprintf(formatAprVal, sources.q(k).ra2000);
            if ispc, aprra = strrep(aprra, 'e+0', 'e+'); aprra = strrep(aprra, 'e-0', 'e-');   end
            aprde=sprintf(formatAprVal, sources.q(k).de2000);
            if ispc, aprde = strrep(aprde, 'e+0', 'e+'); aprde = strrep(aprde, 'e-0', 'e-');  end

            aprra_sigma=sprintf(formatStDev, sources.q(k).aprra_sigma);
            if ispc, aprra_sigma = strrep(aprra_sigma, 'e+0', 'e+'); aprra_sigma = strrep(aprra_sigma, 'e-0', 'e-');  end
            aprde_sigma=sprintf(formatStDev, sources.q(k).aprde_sigma);
            if ispc, aprde_sigma = strrep(aprde_sigma, 'e+0', 'e+'); aprde_sigma = strrep(aprde_sigma, 'e-0', 'e-');   end

            % write one line for each RA De
            fprintf(fid, writeFormat, curIndex,  'RS_RA',num2str(sources.q(k).siteCode), '--', soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'rad', constrApr, aprra, aprra_sigma);
            fprintf(fid, writeFormat, curIndex+1,'RS_DE',num2str(sources.q(k).siteCode), '--', soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'rad', constrApr, aprde, aprde_sigma);
            curIndex=curIndex+2;
        end
    end
    
    % write troposphere zenith wet delay
    if outsnx.zwd==1
        for k=1:numStat
            clear zwdTime zwdTimeYrStr
            zwdTime=mjd2yydoysecod(col_sinex.zwd(k).mjd);
            zwdTimeYrStr=num2str(zwdTime(:,1));
            numZwd=length(col_sinex.zwd(k).col);
                   
            for j=1:numZwd 
                aprZwd=sprintf(formatAprVal, antenna(k).zdry(j));
                if ispc, aprZwd = strrep(aprZwd, 'e+0', 'e+'); aprZwd = strrep(aprZwd, 'e-0', 'e-');   end
                aprZwdStdDev=sprintf(formatStDev, 0);
                if ispc, aprZwdStdDev = strrep(aprZwdStdDev, 'e+0', 'e+'); aprZwdStdDev = strrep(aprZwdStdDev, 'e-0', 'e-');  end

                fprintf(fid, writeFormat, curIndex, 'TROTOT', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, zwdTimeYrStr(j,3:end), zwdTime(j,2), zwdTime(j,3), 'm', constrApr, aprZwd, aprZwdStdDev);
                curIndex=curIndex+1;
            end
        end
    end
    
    
    % write tropospheric north and east gradients
    if outsnx.tgr==1
        for k=1:numStat
            clear ngrTime ngrTimeYrStr
            ngrTime=mjd2yydoysecod(col_sinex.ngr(k).mjd);
            ngrTimeYrStr=num2str(ngrTime(:,1));
            numNgr=length(col_sinex.ngr(k).col);
            for j=1:numNgr 
                aprNgr=sprintf(formatAprVal, antenna(k).apriori_ngr/1000);
                if ispc, aprNgr = strrep(aprNgr, 'e+0', 'e+'); aprNgr = strrep(aprNgr, 'e-0', 'e-');   end
                aprEgr=sprintf(formatAprVal, antenna(k).apriori_egr/1000);
                if ispc, aprEgr = strrep(aprEgr, 'e+0', 'e+'); aprEgr = strrep(aprEgr, 'e-0', 'e-');   end

                aprNgrStdDev=sprintf(formatStDev, 0);
                if ispc, aprNgrStdDev = strrep(aprNgrStdDev, 'e+0', 'e+'); aprNgrStdDev = strrep(aprNgrStdDev, 'e-0', 'e-');  end
                aprEgrStdDev=sprintf(formatStDev, 0);
                if ispc, aprEgrStdDev = strrep(aprEgrStdDev, 'e+0', 'e+'); aprEgrStdDev = strrep(aprEgrStdDev, 'e-0', 'e-');  end

                fprintf(fid, writeFormat, curIndex, 'TGNTOT', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, ngrTimeYrStr(j,3:end), ngrTime(j,2), ngrTime(j,3), 'm', constrApr, aprNgr, aprNgrStdDev);
                fprintf(fid, writeFormat, curIndex+1, 'TGETOT', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, ngrTimeYrStr(j,3:end), ngrTime(j,2), ngrTime(j,3), 'm', constrApr, aprEgr, aprEgrStdDev);
                curIndex=curIndex+2;
            end
        end
    end

    % write blockend line to file
    fprintf(fid, '-%s\n', blockName);

    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});

    %% write 27. SOLUTION/NORMAL_EQUATION_VECTOR block if tickbox is ticked
    blockName='SOLUTION/NORMAL_EQUATION_VECTOR';
    writeFormat=' %5.0f %-6s %04s %2s %4s %2s:%03.0f:%05.0f %4s %1.0f %21s\n';
    if ispc, formatVectorValue='%22.14e'; else
             formatVectorValue='%21.14e'; 
    end
    constrNevCoord=2;
    % define varible index (one for each estimate)
    curIndex=1;

    % write blockstart line to file
    fprintf(fid, '+%s\n', blockName);

    % write info line
    fprintf(fid, '*Index Type__ Code Pt Soln Ref_Epoch___ Unit S RightHandSideVector_b\n');

    % write b vector to file
    % write b vector for station coordinates
    if outsnx.xyz==1
        for k=1:numStat
            % make proper format (e+05 instead of e+005)
            stax_b=sprintf(formatVectorValue, b_sinex(col_sinex.coorx(k)));
            if ispc, stax_b = strrep(stax_b, 'e+0', 'e+'); stax_b = strrep(stax_b, 'e-0', 'e-');  end
            stay_b=sprintf(formatVectorValue, b_sinex(col_sinex.coory(k)));
            if ispc, stay_b = strrep(stay_b, 'e+0', 'e+'); stay_b = strrep(stay_b, 'e-0', 'e-');   end
            staz_b=sprintf(formatVectorValue, b_sinex(col_sinex.coorz(k)));
            if ispc, staz_b = strrep(staz_b, 'e+0', 'e+'); staz_b = strrep(staz_b, 'e-0', 'e-');   end

            fprintf(fid, writeFormat, curIndex, 'STAX',   num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm', constrNevCoord, stax_b);
            fprintf(fid, writeFormat, curIndex+1, 'STAY', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm', constrNevCoord, stay_b);
            fprintf(fid, writeFormat, curIndex+2, 'STAZ', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'm', constrNevCoord, staz_b);
            curIndex=curIndex+3;
        end
    end
    
    
    % write b vector for eops 
    soln=num2str(1);
    
    if outsnx.eop==1
        % constraint
        constrNevXpol=2;
        % write xpol
        for k=1:numXpol
            curLine=eop_matrix(:,1)==col_sinex.mjd_xp(k); % index of current xpol in eop_matrix

            % make correct format
            xpo_b=sprintf(formatVectorValue, b_sinex(col_sinex.xp(k)));
            if ispc, xpo_b = strrep(xpo_b, 'e+0', 'e+'); xpo_b = strrep(xpo_b, 'e-0', 'e-');    end

            fprintf(fid, writeFormat, curIndex, 'XPO', '----', '--', soln, eopDatesYrStr(curLine,3:4), eopDates(curLine,2), eopDates(curLine,3), 'mas', constrNevXpol, xpo_b);  
            curIndex=curIndex+1;   
        end

        % write ypol
        % constraint
        constrNevYpol=2;
        for k=1:numYpol
            curLine=eop_matrix(:,1)==col_sinex.mjd_yp(k); % index of current xpol in eop_matrix

            % make correct format
            ypo_b=sprintf(formatVectorValue, b_sinex(col_sinex.yp(k)));
            if ispc, ypo_b = strrep(ypo_b, 'e+0', 'e+'); ypo_b = strrep(ypo_b, 'e-0', 'e-');  end

            fprintf(fid, writeFormat, curIndex, 'YPO', '----', '--', soln, eopDatesYrStr(curLine,3:4), eopDates(curLine,2), eopDates(curLine,3), 'mas', constrNevYpol, ypo_b);  
            curIndex=curIndex+1;
        end

        % write dUT1
        % constraint
        constrNevDut1=2;
        for k=1:numDut1
            curLine=eop_matrix(:,1)==col_sinex.mjd_dut1(k); % index of current xpol in eop_matrix

            % make correct format
            dut1_b=sprintf(formatVectorValue, b_sinex(col_sinex.dut1(k)));
            if ispc, dut1_b = strrep(dut1_b, 'e+0', 'e+'); dut1_b = strrep(dut1_b, 'e-0', 'e-');  end

            fprintf(fid, writeFormat, curIndex, 'UT', '----', '--', soln, eopDatesYrStr(curLine,3:4), eopDates(curLine,2), eopDates(curLine,3), 'ms', constrNevDut1, dut1_b);  
            curIndex=curIndex+1;
        end

        % write NUTx
        % constraint
        constrNevNutx=2;
        for k=1:numXnut
            curLine=eop_matrix(:,1)==col_sinex.mjd_dX(k); % index of current xpol in eop_matrix

            % make correct format
            dX_b=sprintf(formatVectorValue, b_sinex(col_sinex.dX(k)));
            if ispc, dX_b = strrep(dX_b, 'e+0', 'e+'); dX_b = strrep(dX_b, 'e-0', 'e-');   end

            fprintf(fid, writeFormat, curIndex, 'NUT_X', '----', '--', soln, eopDatesYrStr(curLine,3:4), eopDates(curLine,2), eopDates(curLine,3), 'mas', constrNevNutx, dX_b);  
            curIndex=curIndex+1;
        end

        % write NUTy
        % constraint
        constrNevNuty=2;
        for k=1:numYnut
            curLine=eop_matrix(:,1)==col_sinex.mjd_dY(k); % index of current xpol in eop_matrix
            % make correct format
            dY_b=sprintf(formatVectorValue, b_sinex(col_sinex.dY(k)));
            if ispc, dY_b = strrep(dY_b, 'e+0', 'e+'); dY_b = strrep(dY_b, 'e-0', 'e-');   end

            fprintf(fid, writeFormat, curIndex, 'NUT_Y', '----', '--', soln, eopDatesYrStr(curLine,3:4), eopDates(curLine,2), eopDates(curLine,3), 'mas', constrNevNuty, dY_b);  
            curIndex=curIndex+1;
        end
    end
    
    % write b vector for sources
    if outsnx.sou==1
        constrNevSou=2;
        for k=1:numSou
            % make proper format (e+05 instead of e+005)
            ra_b=sprintf(formatVectorValue, b_sinex(col_sinex.ra(k)));
            if ispc, ra_b = strrep(ra_b, 'e+0', 'e+'); ra_b = strrep(ra_b, 'e-0', 'e-');  end
            de_b=sprintf(formatVectorValue, b_sinex(col_sinex.de(k)));
            if ispc, de_b = strrep(de_b, 'e+0', 'e+'); de_b = strrep(de_b, 'e-0', 'e-');   end

            fprintf(fid, writeFormat, curIndex,  'RS_RA', num2str(sources.q(k).siteCode), '--', soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'rad', constrNevSou, ra_b);
            fprintf(fid, writeFormat, curIndex+1,'RS_DE', num2str(sources.q(k).siteCode), '--', soln, aprDateYrStr(3:end), aprDate(1,2), aprDate(1,3), 'rad', constrNevSou, de_b);
            curIndex=curIndex+2;
        end
    end
    
    
    % write troposphere zenith wet delay
    if outsnx.zwd==1
        constrZwd=2;
        for k=1:numStat
            clear zwdTime zwdTimeYrStr
            zwdTime=mjd2yydoysecod(col_sinex.zwd(k).mjd);
            zwdTimeYrStr=num2str(zwdTime(:,1));
            numZwd=length(col_sinex.zwd(k).col);
            for j=1:numZwd 
                zwd_b=sprintf(formatVectorValue, b_sinex(col_sinex.zwd(k).col(j)));
                if ispc, zwd_b = strrep(zwd_b, 'e+0', 'e+'); zwd_b = strrep(zwd_b, 'e-0', 'e-');   end
                fprintf(fid, writeFormat, curIndex, 'TROTOT', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, zwdTimeYrStr(j,3:end), zwdTime(j,2), zwdTime(j,3), 'm', constrZwd, zwd_b);
                curIndex=curIndex+1;
            end
        end
    end
    
    % write troposphere north and east gradients
    if outsnx.tgr==1
        constrNgr=2;
        constrEgr=2;
        for k=1:numStat
            clear ngrTime ngrTimeYrStr
            ngrTime=mjd2yydoysecod(col_sinex.ngr(k).mjd);
            ngrTimeYrStr=num2str(ngrTime(:,1));
            numNgr=length(col_sinex.ngr(k).col);
            for j=1:numNgr 
                ngr_b=sprintf(formatVectorValue, b_sinex(col_sinex.ngr(k).col(j)));
                if ispc, ngr_b = strrep(ngr_b, 'e+0', 'e+'); ngr_b = strrep(ngr_b, 'e-0', 'e-');   end
                egr_b=sprintf(formatVectorValue, b_sinex(col_sinex.egr(k).col(j)));
                if ispc, egr_b = strrep(egr_b, 'e+0', 'e+'); egr_b = strrep(egr_b, 'e-0', 'e-');   end
                fprintf(fid, writeFormat, curIndex, 'TGNTOT', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, ngrTimeYrStr(j,3:end), ngrTime(j,2), ngrTime(j,3), 'm', constrNgr, ngr_b);
                fprintf(fid, writeFormat, curIndex+1, 'TGETOT', num2str(antenna(k).siteCode), antenna(k).pointCode, soln, ngrTimeYrStr(j,3:end), ngrTime(j,2), ngrTime(j,3), 'm', constrEgr, egr_b);
                curIndex=curIndex+2;
            end
        end
    end

    
    % write blockend line to file
    fprintf(fid, '-%s\n', blockName);

    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});

    
    
    
    
    %% write 28. SOLUTION/NORMAL_EQUATION_MATRIX block if tickbox is ticked
    blockName='SOLUTION/NORMAL_EQUATION_MATRIX L';
    %writeFormat=' %5.0f %5.0f %s\n';
    formatNormEquInd=' %5.0f %5.0f';
    if ispc, formatNormEquVal='  %21.13e';  else
             formatNormEquVal='  %20.13e';
    end

    % write blockstart line to file
    fprintf(fid, '+%s\n', blockName);

    % write info line
    fprintf(fid, '*Row__ Col__ Norm_Equ_Matrix_Value Norm_Equ_Matrix_Valu2 Norm_Equ_Matrix_Valu3');
    
    % get max number of estimated parameter of eops (e.g.: 25 values for
    % xpo/ypo/ut1, 2 values for NutationsX/Y -->25
    maxNumEops=max([numXpol, numYpol, numDut1, numXnut, numYnut]);
    tmpMatEop=zeros(maxNumEops, 5);
    
    % get all indices for N Matrix in one array
    tmpMatCoor=[col_sinex.coorx; col_sinex.coory; col_sinex.coorz];
    tmpMatEop(1:numXpol, 1)=col_sinex.xp;
    tmpMatEop(1:numYpol, 2)=col_sinex.yp;
    tmpMatEop(1:numDut1, 3)=col_sinex.dut1;
    tmpMatEop(1:numXnut, 4)=col_sinex.dX;
    tmpMatEop(1:numYnut, 5)=col_sinex.dY;
    tmpMatSou=[col_sinex.ra; col_sinex.de];
    tmpMatZwd=[col_sinex.zwd.col];
    tmpMatTgr=[col_sinex.ngr.col; col_sinex.egr.col];
    
    % make one vector out of it...
    tmpMat=[tmpMatCoor(:); tmpMatEop(:); tmpMatSou(:); tmpMatZwd(:); tmpMatTgr(:)];
    
    % ... and delete zeros
    tmpMat(tmpMat==0)=[];
    
    % reorder the normal equation matrix
    % preallocate
    N=zeros(size(N_sinex, 1), size(N_sinex, 2));

    
    numStatEst=numStat;
    if outsnx.xyz==0; numStatEst=0; end % if station coordinates are fixed
    
    for col=1:size(N,1)
        for stat=1:numStatEst
            N((stat-1)*3+1,col)=N_sinex(col_sinex.coorx(stat), tmpMat(col));
            N((stat-1)*3+2,col)=N_sinex(col_sinex.coory(stat), tmpMat(col));
            N((stat-1)*3+3,col)=N_sinex(col_sinex.coorz(stat), tmpMat(col));
        end
        for xpo=1:numXpol
            N(numStatEst*3+ xpo, col)=N_sinex(col_sinex.xp(xpo), tmpMat(col));
        end
        for ypo=1:numYpol
            N(numStatEst*3+numXpol +ypo, col)=N_sinex(col_sinex.yp(ypo), tmpMat(col));
        end
        for dut1=1:numDut1
            N(numStatEst*3+numXpol+numYpol +dut1, col)=N_sinex(col_sinex.dut1(dut1), tmpMat(col));
        end
        for dX=1:numXnut
            N(numStatEst*3+numXpol+numYpol+numDut1 +dX, col)=N_sinex(col_sinex.dX(dX), tmpMat(col));
        end
        for dY=1:numYnut
            N(numStatEst*3+numXpol+numYpol+numDut1+numXnut +dY, col)=N_sinex(col_sinex.dY(dY), tmpMat(col));
        end
        for sou=1:numSou
            N(numStatEst*3+numXpol+numYpol+numDut1+numXnut+numYnut +(sou-1)*2+1, col)=N_sinex(col_sinex.ra(sou), tmpMat(col));
            N(numStatEst*3+numXpol+numYpol+numDut1+numXnut+numYnut +(sou-1)*2+2, col)=N_sinex(col_sinex.de(sou), tmpMat(col));
        end
        numZwdAll=length(tmpMatZwd);
        for zwd=1:numZwdAll
            N(numStatEst*3+numXpol+numYpol+numDut1+numXnut+numYnut+numSou*2 +zwd, col)=N_sinex(tmpMatZwd(zwd), tmpMat(col));
        end
        tmpMatNgr=[col_sinex.ngr.col]; tmpMatEgr=[col_sinex.egr.col];
        numNgrAll=size(tmpMatTgr,2);
        for tgr=1:numNgrAll
            N(numStatEst*3+numXpol+numYpol+numDut1+numXnut+numYnut+numSou*2+numZwdAll +(tgr-1)*2+1, col)=N_sinex(tmpMatNgr(tgr), tmpMat(col)); %north
            N(numStatEst*3+numXpol+numYpol+numDut1+numXnut+numYnut+numSou*2+numZwdAll +(tgr-1)*2+2, col)=N_sinex(tmpMatEgr(tgr), tmpMat(col)); %east
        end
        
%         for eop=1:numEOPs
%             N(numStatEst*3+ (eop-1)*5+1,col)=N_sinex(col_sinex.xp(eop), tmpMat(col));
%             N(numStatEst*3+ (eop-1)*5+2,col)=N_sinex(col_sinex.yp(eop), tmpMat(col));
%             N(numStatEst*3+ (eop-1)*5+3,col)=N_sinex(col_sinex.dut1(eop), tmpMat(col));
%             N(numStatEst*3+ (eop-1)*5+4,col)=N_sinex(col_sinex.dX(eop), tmpMat(col));
%             N(numStatEst*3+ (eop-1)*5+5,col)=N_sinex(col_sinex.dY(eop), tmpMat(col));
%         end
    end

    % write normal equation matrix
    for row=1:size(N,1)
        for col=1:row
            if mod(col,3)==1
                fprintf(fid, '\n');
                fprintf(fid, formatNormEquInd, row, col);
            end
            % make proper format
            Nstring=sprintf(formatNormEquVal, N(row,col));
            if ispc, Nstring = strrep(Nstring, 'e+0', 'e+'); Nstring = strrep(Nstring, 'e-0', 'e-');   end
               fprintf(fid, '%21s', Nstring);
        end
    end
            
%     writeFormat=' %5.0f %5.0f %s\n';
%     formatNormEquVal='%22.14e';
%     for row=1:size(N,1)
%         for col=1:row
%             % make proper format
%             Nstring=sprintf(formatNormEquVal, N(row, col));
%             if ispc, Nstring = strrep(Nstring, 'e+0', 'e+'); Nstring = strrep(Nstring, 'e-0', 'e-');  end
%             
%             fprintf(fid, writeFormat, row, col, Nstring);
%         end
%     end

    % write blockend line to file
    fprintf(fid, '\n-%s\n', blockName);

    % write headerspace
    fprintf(fid, '%s\n%s\n%s\n', blockSpace{:});


    %% write endsinex line
    fprintf(fid, '%%ENDSNX');
    
    %% delete variables
    clear obsCode
    fclose(fid);
    
    
    fprintf('%s%s\n','SINEX file is saved as ',snxFile);
    
end %for - process_list


