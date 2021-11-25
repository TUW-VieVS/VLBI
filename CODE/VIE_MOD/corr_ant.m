% ************************************************************************
%   Description:
%   Loads data for station corrections and writes it to the antenna 
%   structure. 
%  
%   Input:										
%      mjd1                Modified Julian Date of first observation for 
%                          the scan; gives information of the year of 
%                          observation[d]
%      tim      (7,1)      year,month,day,hour,min,sec,doy
%      antenna             VieVS antenna structure array
%      parameter           VieVS parameter structure array
%      flagmess            VieVS error flag structure
% 
%   Output:
%      antenna             VieVS anntenna structure array with station
%                          correction data
%      flagmess            VieVS error flag structure
% 
%   Data directories:
%      ../OTIDE/
%      ../ATIDE/
%      ../ATM/
%      ../THERMDEF/TDEF.mat
%
%   External calls: 	
%		xyz2ell.m
%		ren2xyz.m
%		call_ren2xyz.m
%      %grepstat.m, read_atm.m                 			
%       greapload.m
%		datachecking.m
%       
%   Coded for VieVS: 
%   23 Nov 2009 by Hana Spicakova
%
%   Revision: 
%   03 Feb 2010 by Lucia Plank: thermdef
%   05 Feb 2010 by Tobias Nilsson: Now reads atm/vm1 data from next year for
%       sessions ending close to the end of the year.
%   12 Mar 2010 by Lucia Plank: eccentricities added
%   30 Mar 2010 by Lucia Plank: tim instead of jul2dat.m
%   26 May 2010 by Lucia Plank: flagmessage
%   31 Aug 2010 by Sigrid Böhm and Tobias Nilsson: Blanks in station names
%   10 Nov 2011 by Hana Spicakova: change between more APL series
%   10 Nov 2011 by Hana Spicakova: ocean pole tide loading added (IERS2010)
%   27 Mar 2012 by Hana Spicakova: hydrology loading
%   08 Nov 2012 by Hana Kr�sn�: changed to superstation file (ocean tidal loading,
%                  atmosphere tidal loading, ocean pole tide loading, thermal deformation coeeficients)
%   05 Dec 2013 by Hana Krasna: APL regression coefficients added 
%   23 Jan 2014 by Hana Krasna: GIA uplift rates as apriori added
%   12 Jun 2014 by Hana Krasna: works now also without ocean loading
%                  corrections
%   01 Jul 2014 by Hana Krasna: works now also without oc. pole tide loading
%                  (should not be needed as the corrections are interpolated)
%   04 Nov 2015 by A. Girdiuk: thermal antenna deformation was included in antenna structure in VIE_INIT
%							   there does not include anymore
%   26 Jan 2016 by D. Landskron: slight changes made to make it faster
%   27 Apr 2016 by D. Landskron: non-tidal ATM now uses only ATM data of 1 year, as this is enough
%   29 Apr 2016 by M. Madzak: Added post-seismic deformation to be loaded
%                             from superstat file (eg ITRF2014a)
%   17 May 2016 by A.Girdiuk: loops are removed. Necessity and usage of next year data are checked.
%							  External grepstat.m and grepstat70.m functions
%                             are replaced by greploaddat.m
%							  datachecking.m is added to help when next year 
%							  includes some equaled data to previous year.
%  24 May 2016 by A.Girdiuk : bug-fix in case of no loads or vmf1 is available
%  23 Jun 2016 by A.Girdiuk : bug-fix in loadings
%  08 Aug 2016 by A.Girdiuk : new hydrology loading provider is added
%  24 Jan 2017 by D. Landskron: zhd and zwd also written to parameter.vm1
%  22 Feb 2017 by A. Hellerschmied: Manual TRF file support established
%  23 Feb 2017 by A. Hellerschmied: flagmess is added as in- and output argument (not global any more!)
%  13 Dec 2017 by D. Landskron: VMF1 input changed from .mat to .vmf1_r
%  11 Jan 2018 by D. Landskron: VM1 folder moved to TRP and renamed to VMF1
%  07 Feb 2018 by D. Landskron: bug-fix with VM1 close to New Year
%  06 Jul 2018 by D. Landskron: vm1 renamed to vmf1 and VMF3 added to the troposphere models 
%  06 Jun 2019 by D. Landskron: reading of non-tidal APL coefficients adapted to ascii files
% *************************************************************************
%
function [antenna, flagmess] = corr_ant(time_lim, tim, antenna, parameter, flagmess)

mjd1=time_lim(1);
mjd2=time_lim(2);

%  READ EXTERNAL FILES

trffile = parameter.vie_init.trf;
if strcmp(trffile{1}(end-3:end), '.mat')
    trf=load(trffile{1});
    nam=fieldnames(trf);
    trf=eval(['trf.' nam{1}]);
else % a manual trf file is given -> BUT load default superstat file in any case!!
    trf=load('../TRF/superstation.mat');
    nam=fieldnames(trf);
    trf=eval(['trf.' nam{1}]);
end


iye  = tim(1);
idoy = tim(7);

if (idoy>360) || (idoy<6) % There is a chance that we need ATM/VMF1 data from next year.
  numyrs=2;
else
  numyrs=1;
end
nant=length(antenna);

% Loop over antennas to read corrections from external files
for ist=1:nant
    
    aname = antenna(ist).name;
    ant = [antenna(ist).x, antenna(ist).y, antenna(ist).z];    
    
    IDsuper = antenna(ist).IDsuper;
    if length(trf)<IDsuper % case of a new station in manual trf file
        IDsuper = [];
    end
    
    %%% Eccentricities
    if strcmp(antenna(ist).ecctype,'NEU')
        [phi,lam]=xyz2ell(ant);
        antenna(ist).c_ecc = ren2xyz([antenna(ist).ecc(3),...
                                      antenna(ist).ecc(2),...
                                      antenna(ist).ecc(1)],phi,lam);
    else
        antenna(ist).c_ecc = antenna(ist).ecc;
    end
    
    %%% Thermal antenna deformation %%%
    % if no coefficients --> send error message
    if isempty (antenna(ist).thermal)
        flagmess.thermal(ist)=1;
    end

        
    %%%% underline/blank in station name %%%%%%%%  
    idblank = aname(1:find(aname~=' ', 1, 'last' ))==' ';
    aname(idblank)='_';


    %%% Tidal Ocean Loading, cto %%%
    if parameter.vie_mod.cto == 1 
        % Get apriori tidal ocean loading amplitudes and phases for stations in the session      
        
        if ~isempty(trf(IDsuper).ocean_loading)
             antenna(ist).cto = trf(IDsuper).ocean_loading.(parameter.vie_mod.ocm);
        else
            antenna(ist).cto=[];
        end
        % If no ocean loading data for the station, send a message
        if isempty (antenna(ist).cto)
            flagmess.cto(ist)=1;
        end
    end 
    
    
    %%% Tidal Atmosphere Loading, cta %%%
    if parameter.vie_mod.cta == 1 
        % Get apriori tidal ocean loading amplitudes and phases for stations in the session
        model = parameter.vie_mod.ctam;
        antenna(ist).cta = trf(IDsuper).atmosphere_tidal_loading.(model);
        
        % If no atmosphere loading data for the station, send a message
        if isempty (antenna(ist).cta)
            flagmess.cta(ist)=1;
        end
    end
    
    %%% Non tidal atmosphere loading  %%% cnta
    if parameter.vie_mod.cnta == 1
        
        antenna(ist).cnta_dx=[];
        for idyr=0:numyrs-1
            fil=num2str(sprintf(['../ATM/' parameter.vie_mod.cntam '/y%4d.apl_r'],iye+idyr));
            if exist(fil,'file')
                fid = fopen(fil);
                apl_data = textscan(fid,'%s%f%f%f%f','CommentStyle','!');
                fclose(fid);
                if ~isempty(find(strcmpi(apl_data{1},strtrim(aname))))
                   
                    % reduce VMF3 data to lines that contain the respective station
                    wantedLines = ismember(apl_data{1,1}, strtrim(aname));
                    for k=1:length(apl_data)
                        apl_data{k} = apl_data{k}(wantedLines);
                    end
                    
                    % make VMF3 data smaller by extracting only data on the respective day + 1 day before and after
                    wantedLines = apl_data{2}>mjd1-1.25 & apl_data{2}<mjd2+1;
                    for k=1:length(apl_data)
                        apl_data{k} = apl_data{k}(wantedLines);
                    end
                    
                    if ~isempty(apl_data)
                        % Save into antenna.cnta_dx
                        matm2 = [apl_data{2} apl_data{3} apl_data{4} apl_data{5}];
                        [nta_dxyz] = call_ren2xyz(matm2,ant);
                        antenna(ist).cnta_dx = [antenna(ist).cnta_dx ; nta_dxyz];  % = [station,tmjd,ah,aw,zhd,zwd]
                    end
                end
            end
        end
        
        % Reduce data, if it contains 2 years
        if isempty(antenna(ist).cnta_dx)
            flagmess.cnta(ist)=1;
        else
            [antenna(ist).cnta_dx,flag]=datachecking(antenna(ist).cnta_dx,mjd1,mjd2);
            if ~flag
                flagmess.cnta(ist)=1;
                antenna(ist).cnta_dx=[];
                fprintf('%s : not enough non-tidal atm at session time: %f to %f\n',aname,mjd1,mjd2)
            end
        end

    end
    
    
    %%% APL with regression coefficients, crg %%%
    if parameter.vie_mod.crg == 1 
        model = parameter.vie_mod.crgm;
        antenna(ist).crg = trf(IDsuper).aplrg.(model); %[p0 RG[m/hPa]]
        
        % If no APL RG at the station, send a message
        if isempty (antenna(ist).crg)
            flagmess.crg(ist)=1;
        end
    else
        refpres=antenna(ist).gpt3.p; % reference pressure taken from GPT3
        antenna(ist).crg = [refpres 0];
    end
    
    
    % GIA vertical velocity
    if parameter.vie_mod.gia == 1
        model = parameter.vie_mod.giam;
        
        gia_vren = [0 0 0 0]; %mjdref r e n
        
        gia_vren(1) = trf(IDsuper).gia.refEpoch;
        gia_vren(2) = trf(IDsuper).gia.(model); % up velocity [m/y]
        gia_dvx = call_ren2xyz(gia_vren,ant);
        
        antenna(ist).gia_dvx = gia_dvx; 
        
        % If no GIA data, send a message
        if length(antenna(ist).gia_dvx)<4
            flagmess.gia(ist)=1;
        end

    end 
    

    %%% Hydrology loading %%%
     if parameter.vie_mod.chl == 1
		antenna(ist).chl_dx=[];
		if strcmp(parameter.vie_mod.chlm,'GSFC')
			numyrs_hyd = numyrs+1;
	        iye_hyd = iye-1;
    	    for idyr=0:numyrs_hyd
    	        fil=num2str(sprintf(['../HYDLO/' parameter.vie_mod.chlm '/%4d_CMTE_HYDLO.mat'],iye_hyd+idyr));
    	        if exist(fil,'file')
    	        	load (fil);
    	            if ~isempty(find(strcmpi({hydlo.ivsname},aname)))
    	            	mhyd = greploaddata(hydlo,aname,mjd1,mjd2); 
    	                if ~isempty(mhyd)
    	                    [hyd_dxyz] = call_ren2xyz(mhyd,ant);
    	                    antenna(ist).chl_dx = [antenna(ist).chl_dx; hyd_dxyz];
    	                end
    	            end
    	        end
    	    end
        
		else
			fil=sprintf(['../HYDLO/' parameter.vie_mod.chlm '/erahyd.mat']);
			load (fil);
            if ~isempty(find(strcmpi({erahyd.ivsname},aname)))
                mload = greploaddata(erahyd,aname,mjd1,mjd2);
                mload(:,2:4) = mload(:,2:4)/1000; % North,East,Up: [mm] => [m]
                mload(:,2:4) = [mload(:,4) mload(:,3) mload(:,2)];
                if ~isempty(mload)
                    [hyd_dxyz] = call_ren2xyz(mload,ant);
                    antenna(ist).chl_dx = [antenna(ist).chl_dx; hyd_dxyz];
                end
            end
		end        
        % Save into antenna.chl_dx
		if isempty (antenna(ist).chl_dx)
    		flagmess.chl(ist)=1;
		else
    		[antenna(ist).chl_dx,flag]=datachecking(antenna(ist).chl_dx,mjd1,mjd2);
    	    if ~flag
    	    	flagmess.chl(ist)=1;
    	        antenna(ist).chl_dx=[];
    	        fprintf('%s : not enough hydlo at session time: %f to %f\n',aname,mjd1,mjd2)
			end
		end
     end
        
    %%% Ocean pole tide loading %%%
    if parameter.vie_mod.ctop == 1 
        if ~isempty(trf(IDsuper).oceanPoleTideLoading)
           antenna(ist).opl = trf(IDsuper).oceanPoleTideLoading.desai;  %u_r^R, u_r^I, u_n^R, u_n^I, u_e^R, u_e^I
        %   antenna(ist).opl = trf(IDsuper).oceanPoleTideLoading.(parameter.vie_mod.ctopm);  %u_r^R, u_r^I, u_n^R, u_n^I, u_e^R, u_e^I
        else
            antenna(ist).opl =[];
        end
      
        % If no OPL coefficients for the station, send a message
        if isempty (antenna(ist).opl)
            flagmess.ctop(ist)=1;
        end
    end
    
    
    %%% Vienna Mapping Function 1, VMF3 %%%
    if strcmp(parameter.vie_init.zhd,'vmf3')   ||   strcmp(parameter.vie_init.zwd,'vmf3')   ||   strcmp(parameter.vie_mod.mfh,'vmf3')   ||   strcmp(parameter.vie_mod.mfw,'vmf3')
        
        antenna(ist).vmf3=[];
        for idyr=0:numyrs-1
            fil=num2str(sprintf('../TRP/VMF3/y%4d.vmf3_r',iye+idyr));
            if exist(fil,'file')
                fid = fopen(fil);
                vmf3_data = textscan(fid,'%s%f%f%f%f%f%f%f%f','CommentStyle','#');
                fclose(fid);
                if ~isempty(find(strcmpi(vmf3_data{1},strtrim(aname))))
                   
                    % reduce VMF3 data to lines that contain the respective station
                    wantedLines = ismember(vmf3_data{1,1}, strtrim(aname));
                    for k=1:length(vmf3_data)
                        vmf3_data{k} = vmf3_data{k}(wantedLines);
                    end
                    
                    % make VMF3 data smaller by extracting only data on the respective day + 1 day before and after
                    wantedLines = vmf3_data{2}>mjd1-1.25 & vmf3_data{2}<mjd2+1;
                    for k=1:length(vmf3_data)
                        vmf3_data{k} = vmf3_data{k}(wantedLines);
                    end
                    
                    if ~isempty(vmf3_data)
                        % Save into antenna.vmf3
                        antenna(ist).vmf3 = [antenna(ist).vmf3 ; vmf3_data{2},vmf3_data{3},vmf3_data{4},vmf3_data{5},vmf3_data{6}];  % = [station,tmjd,ah,aw,zhd,zwd]
                    end
                end
            end
        end
        
        % Reduce data, if it contains 2 years
        if isempty(antenna(ist).vmf3)
            flagmess.vmf3(ist)=1;
        else
            [antenna(ist).vmf3,flag]=datachecking(antenna(ist).vmf3,mjd1,mjd2);
            if ~flag
                flagmess.vmf3(ist)=1;
                antenna(ist).vmf3=[];
                fprintf('%s : not enough vmf3 at session time: %f to %f\n',aname,mjd1,mjd2)
            end
        end
        
    end
    
    
    %%% Vienna Mapping Function 1, VMF1 %%%
    if strcmp(parameter.vie_init.zhd,'vmf1')   ||   strcmp(parameter.vie_init.zwd,'vmf1')   ||   strcmp(parameter.vie_mod.mfh,'vmf1')   ||   strcmp(parameter.vie_mod.mfw,'vmf1')
        
        antenna(ist).vmf1=[];
        for idyr=0:numyrs-1
            fil=num2str(sprintf('../TRP/VMF1/y%4d.vmf1_r',iye+idyr));
            if exist(fil,'file')
                fid = fopen(fil);
                vmf1_data = textscan(fid,'%s%f%f%f%f%f%f%f%f%f%f','CommentStyle','#');
                fclose(fid);
                if ~isempty(find(strcmpi(vmf1_data{1},strtrim(aname))))
                   
                    % reduce VMF1 data to lines that contain the respective station
                    wantedLines = ismember(vmf1_data{1,1}, strtrim(aname));
                    for k=1:length(vmf1_data)
                        vmf1_data{k} = vmf1_data{k}(wantedLines);
                    end
                    
                    % make VMF1 data smaller by extracting only data on the respective day + 1 day before and after
                    wantedLines = vmf1_data{2}>mjd1-1.25 & vmf1_data{2}<mjd2+1;
                    for k=1:length(vmf1_data)
                        vmf1_data{k} = vmf1_data{k}(wantedLines);
                    end
                    
                    if ~isempty(vmf1_data)
                        % Save into antenna.vmf1
                        antenna(ist).vmf1 = [antenna(ist).vmf1;vmf1_data{2},vmf1_data{3},vmf1_data{4},vmf1_data{5},vmf1_data{6}];  % = [station,tmjd,ah,aw,zhd,zwd]
                    end
                end
            end
        end
        
        % Reduce data, if it contains 2 years
        if isempty(antenna(ist).vmf1)
            flagmess.vmf1(ist)=1;
        else
            [antenna(ist).vmf1,flag]=datachecking(antenna(ist).vmf1,mjd1,mjd2);
            if ~flag
                flagmess.vmf1(ist)=1;
                antenna(ist).vmf1=[];
                fprintf('%s : not enough vmf1 at session time: %f to %f\n',aname,mjd1,mjd2)
            end
        end
        
    end
    
    
    %%% Post-seismic deformation
    if ~strcmp(trffile{2}, 'manualTrf')
        if isfield(trf(IDsuper).(trffile{2}), 'psd')
            antenna(ist).psd=trf(IDsuper).(trffile{2}).psd;
        elseif strcmp(trffile{2}, 'vievsTrf')
            if isfield(trf(IDsuper).itrf2014, 'psd')
                antenna(ist).psd=trf(IDsuper).itrf2014.psd; % apply psd from ITRF2014 to vievsTRF
                fprintf('\n ITRF2014 psd applied at station %8s !!!\n',antenna(ist).name);
            end
        end
    end
    
    
end % 1:nant
