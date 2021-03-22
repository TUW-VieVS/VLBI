% This function creates a superstation (.mat) file where all available
% information for all (IVS) VLBI stations is written to.
% 7.10.2011, Matthias Madzak
% 14 Jan 2014, APL Regression Coefficients added by Hana
% 23 Jan 2014, GIA Uplift Rates added by Hana
% 03 Jun 2014, skip the second entry for Efelsberg in ns-codes.txt by Hana
% 05 Jun 2014, VieTRF10a changed to VieTRF13 by Hana
% 05 Jun 2014, check for _ --> '' added for more files by Hana
% 18 Jun 2014, Vienna tidal atmosphere loading added by Hana
% 12 Nov 2015, blank lines in cs_codes are left out by Matthias
% 10 May 2016, M. Madzak, A. Girdiuk: ITRF2014 and VTRF2014 added
% 11 May 2016, A. Girdiuk: bug-fix
% 07 Jun 2017, H. Krasna: bug fixed, read the station names from the user TRF with 8 chars
% 20 Oct 2017, H. Krasna: DTRF2014 added
 

function [varargout]=mk_superstatFile(inFiles, outFile)

% ==============
% 1. get inFiles
% ==============
idiF=1;
nsCodesFile=inFiles(idiF).name; idiF=idiF+1;
antennaInfoFile=inFiles(idiF).name; idiF=idiF+1;
eccdatFile=inFiles(idiF).name; idiF=idiF+1;
gravdefFile=inFiles(idiF).name; idiF=idiF+1;
itrf2014File=inFiles(idiF).name; idiF=idiF+1;
dtrf2014File=inFiles(idiF).name; idiF=idiF+1;
vtrf2014File=inFiles(idiF).name; idiF=idiF+1;
ivstrf2014bFile=inFiles(idiF).name; idiF=idiF+1;
vieTrf13File=inFiles(idiF).name; idiF=idiF+1;
s123viennaFile=inFiles(idiF).name; idiF=idiF+1;
s12vandamFile=inFiles(idiF).name; idiF=idiF+1;
s12gsfcFile=inFiles(idiF).name; idiF=idiF+1;
s12UserOwnFile=inFiles(idiF).name; 
s12UserOwnFieldname=inFiles(idiF).fieldName; idiF=idiF+1;
oceanLoadingTPXO72File=inFiles(idiF).name;
oceanLoadingTPXO72Fieldname=inFiles(idiF).fieldName; idiF=idiF+1;
oceanLoadingGOT00File=inFiles(idiF).name;
oceanLoadingGOT00Fieldname=inFiles(idiF).fieldName; idiF=idiF+1;
oceanLoadingEOT08aFile=inFiles(idiF).name;
oceanLoadingEOT08aFieldname=inFiles(idiF).fieldName; idiF=idiF+1;
oceanLoadingFES2004File=inFiles(idiF).name;
oceanLoadingFES2004Fieldname=inFiles(idiF).fieldName; idiF=idiF+1;
oceanLoadingAG06File=inFiles(idiF).name;
oceanLoadingAG06Fieldname=inFiles(idiF).fieldName; idiF=idiF+1;
oceanLoadingUserOwnFilename=inFiles(idiF).name;
oceanLoadingUserOwnFieldname=inFiles(idiF).fieldName; idiF=idiF+1;
opoleloadcoefcmcorFile=inFiles(idiF).name; idiF=idiF+1;
opoleloadUserOwnFilename=inFiles(idiF).name; 
opoleloadUserOwnFieldname=inFiles(idiF).fieldName; idiF=idiF+1;
vievsTrfFile=inFiles(idiF).name; idiF=idiF+1;
userOwnTrfFile=inFiles(idiF).name;

% Flags, for the availability of TRF files
flag_trf_itrf2014 = true;
flag_trf_dtrf2014 = true;
flag_trf_vtrf2014 = true;
flag_trf_vieTrf13 = true;
flag_trf_ivstrf2014b = true;

aplrgFilename1 = '../TRF/data/APL_RC_external1.txt'; % added by Hana 01/2014
aplrgFieldname1 = 'p0_RCexternal';
aplrgFilename2 = '../TRF/data/APL_RC_vievs1.txt';
aplrgFieldname2 = 'p0_RCvievs';


giaFilename1 = '../TRF/data/gia_uplift_rates_ICE_5G_VM2_2012.txt';
giaFieldname1 = 'ICE_5G_VM2_2012';
giaMJDFilename1 = '../TRF/data/gia_mean_mjd_8413.txt';
giaMJDFieldname1 = 'refEpoch';

% atmo tidal loading grid file from vandam
atmoTidalLoadVandamGridFile='../TRF/data/s1_s2_def_cm.dat';
atmoTidalLoadViennaGridFile='../TRF/data/s1_s2_s3_cm_noib_grid.dat';

break0=struct('start', [],  'epoch', [],  'end', [], 'x', [], 'y', [],... 
    'z', [], 'x_sigma', [], 'y_sigma', [], 'z_sigma', [], 'vx', [], ...
    'vy', [], 'vz', [], 'vx_sigma', [], 'vy_sigma', [], 'vz_sigma', []);

% preallocate ns_codes
atm0=struct('vienna', [], 'gsfc', [], 'vandam', []);

ns_codes=struct('code', [], 'name', [],'domes', [], 'CDP', [], 'comments', [], ...
    'antenna_info', [], 'ecc', [], ...
    'ocean_loading', [],  ...
    'itrf2014', [], 'dtrf2014', [], ...
    'vtrf2014', [], 'ivsTrf2014b', [], 'VieTRF13', [], 'vievsTrf', [], ...
    'atmosphere_tidal_loading', [], 'oceanPoleTideLoading', [], 'aplrg', [], ...
    'gia', [], 'gravdef', []);

%% ============
% 2. read files
% =============

% First check, if vievsTRF is available!
% - It is required in any case, because it is used as backup TRF file!
if isempty(vievsTrfFile)
    varargout{1} = 'The loction of the vievsTRF.txt (mandatory!) is not specified! Define the file location and run the program again.';
    return;
end

% -----------------------
% 2.1 Antenna information
% -----------------------
fprintf('2.1 Antenna information\n')

% ----------------
% 2.1.1 ns-codes.txt
% ----------------

fprintf('2.1.1 Loading ns_codes.txt\n\n');

if ~isempty(nsCodesFile)
	% open file
	fid=fopen(nsCodesFile);
    if (fid < 0)
        varargout{1} = ['ERROR: Cannot open the file: ', nsCodesFile];
        return;
    end
	% preallocate ns_codes struct

	i=1;

	% read
	while ~feof(fid)
		line=fgetl(fid);
		% if we have a blank line, leave that line
		if isempty(line)
			continue
		end
		
		line=[line,'                                                                                                                                                         '];
		if line(1) ~= '*' & strcmp(line(2:3),'Ef')==0       % skip the second entry for Efelsberg
			ns_codes(i).code = line(1:3);
			ns_codes(i).name = line(5:12);
			ns_codes(i).domes = line(14:22);
			ns_codes(i).CDP = line(24:27);
			ns_codes(i).comments = line(29:length(line));
			i=i+1; % increase index
			
		end
	end
	
	% close file
	fclose(fid);

	% get number of sites
    number_of_site=length(ns_codes);	
else
    fprintf('ns_codes.txt is not available\n\n');
    varargout{1} = 'ns_codes.txt is not specified!';
    return;
end





% --------------------
% 2.1.2 antenna-info.txt
% --------------------

fprintf('\n2.1.2 Loading antenna-info.txt\n\n');


if ~isempty(antennaInfoFile)
	% open file
	fid=fopen(antennaInfoFile);
    if (fid < 0)
        varargout{1} = ['ERROR: Cannot open the file: ', antennaInfoFile];
        return;
    end

	j=0;
	while ~feof(fid)
		
		line = [fgetl(fid),'                                                                                                                                                              '];
		if strcmp(line(1:12), 'ANTENNA_INFO')
			
			% find ivsname in ns_codes (there is no DOMES number)
			
			iStat=find(strcmpi({ns_codes.name}, line(15:22)));
			
			if ~isempty(iStat)
				iStat=iStat(1);
				ns_codes(iStat).antenna_info.focus            = line(25:31);
				ns_codes(iStat).antenna_info.mount            = line(33:39);
				ns_codes(iStat).antenna_info.flag             = line(41:47);
				ns_codes(iStat).antenna_info.measuretype      = line(49:55);
				ns_codes(iStat).antenna_info.reftemp          = str2double(line(58:61));
				ns_codes(iStat).antenna_info.reftemp_sin      = str2double(line(63:66));
				ns_codes(iStat).antenna_info.reftemp_cos      = str2double(line(68:71));
				ns_codes(iStat).antenna_info.refpres          = str2double(line(73:78));
				ns_codes(iStat).antenna_info.diam             = str2double(line(81:85));
				ns_codes(iStat).antenna_info.found_h          = str2double(line(87:93));
				ns_codes(iStat).antenna_info.found_d          = str2double(line(95:100));
				ns_codes(iStat).antenna_info.found_texp       = str2double(line(103:109));
				ns_codes(iStat).antenna_info.fixedaxis        = str2double(line(112:118));
				ns_codes(iStat).antenna_info.fixedaxis_texp   = str2double(line(120:126));
				ns_codes(iStat).antenna_info.axis_offset      = str2double(line(129:135));
				ns_codes(iStat).antenna_info.axis_offset_texp = str2double(line(137:143));
				ns_codes(iStat).antenna_info.avertex          = str2double(line(146:152));
				ns_codes(iStat).antenna_info.avertextexp      = str2double(line(154:160));
				ns_codes(iStat).antenna_info.subref_h         = str2double(line(163:169));
				ns_codes(iStat).antenna_info.subref_h_texp    = str2double(line(171:177));
			else
				fprintf('%s (antenna-info.txt) not found in ns_codes\n', line(15:22));
				
			end
		end
	end
	fclose(fid);
else
    fprintf('antenna-info is not available\n\n');
    varargout{1} = 'antenna-info.txt is not specified!';
    return;
end


% ---------------
% 2.1.3 ECCDAT.dat
% ---------------

fprintf('\n2.1.3 Loading ECCDAT.ecc\n\n');

if ~isempty(eccdatFile)
	fid=fopen(eccdatFile);
    if (fid < 0)
        varargout{1} = ['ERROR: Cannot open the file: ', eccdatFile];
        return;
    end

	while ~feof(fid)
		line = [fgetl(fid),'                                                                                                          '];
		
		% if we have data line
		if strcmp(line(1), ' ')
		
			% find station in ns_codes
			iStat=find(strcmpi({ns_codes.name}, line(3:10)));
			
			% if it is empty, try to find stationname with '_' instead of ' '
			if isempty(iStat)
% 				iStat=find(strcmpi({ns_codes.name}, strrep(line(3:10), ' ', '_'))); % does not work for 'JPL MV1 '
                iStat=find(strcmpi({ns_codes.name}, [strrep(line(3:9), ' ', '_') line(10)]));
			end
			
			

			if ~isempty(iStat)
				if ~isfield(ns_codes(iStat), 'ecc')
					curBreak=1;
				else
					% get number of break (if there was already the stat found)
					if isfield(ns_codes(iStat(1)).ecc, 'break')
						curBreak=length(ns_codes(iStat).ecc.break)+1;
					else
						curBreak=1;
					end
				end
				iStat=iStat(1);
				ns_codes(iStat).ecc.break(curBreak).starting       = line(18:33);
				ns_codes(iStat).ecc.break(curBreak).ending         = line(36:51);
				ns_codes(iStat).ecc.break(curBreak).FCE            = str2double(line(54:63));
				ns_codes(iStat).ecc.break(curBreak).SCE            = str2double(line(65:74));
				ns_codes(iStat).ecc.break(curBreak).TCE            = str2double(line(76:85));
				ns_codes(iStat).ecc.break(curBreak).type_e         = line(88:90);
			else
				fprintf('%s (station in ECCDAT.ecc) not found in ns_codes ('' '' -> ''_'' checked)\n', line(3:10));
			end
		end
	end
	fclose(fid);
else
    fprintf('ECCDAT.ecc is not available\n\n');
    varargout{1} = 'ECCDAT.ecc is not specified!';
    return;
end

% ---------------
% 2.1.4 Gravitational deformation
% ---------------

fprintf('\n2.1.4 Gravitational deformation\n\n');

if ~isempty(gravdefFile)
	fid=fopen(gravdefFile);
    if (fid < 0)
        fprintf(['ERROR: Cannot open the file ', gravdefFile]);
        varargout{1} = ['ERROR: Cannot open the file ', gravdefFile];
        return;
    end

    [ns_codes] = gravdef_parser(ns_codes, gravdefFile);
end
% Markus end


% -----------------
% 2.2 Ocean loading
% -----------------
fprintf('\n2.2 Ocean loading\n')


%fprintf('\n2.2.1 Loading ocean loading files\n\n');
lines2print={'Loading ocean_loading_TPXO72.txt', 'Loading ocean_loading_GOT00.txt',...
    'Loading ocean_loading_EOT08a.txt', 'Loading ocean_loading_FES2004.txt',...
    'Loading ocean_loading_AG06.txt'};
nFiles=size(lines2print,2);

files2load={oceanLoadingTPXO72File, oceanLoadingGOT00File,...
    oceanLoadingEOT08aFile, oceanLoadingFES2004File, ...
    oceanLoadingAG06File};
ol_fieldnames={oceanLoadingTPXO72Fieldname, oceanLoadingGOT00Fieldname,...
    oceanLoadingEOT08aFieldname, oceanLoadingFES2004Fieldname, ...
    oceanLoadingAG06Fieldname};

% if the user own file is given - append to cell strings
if ~isempty(oceanLoadingUserOwnFilename)
    lines2print(nFiles+1)={'Loading ocean_loading user own'};
    files2load(nFiles+1)={oceanLoadingUserOwnFilename};
    ol_fieldnames(nFiles+1)={oceanLoadingUserOwnFieldname};
    nFiles=nFiles+1;
end

for k=1:nFiles % 6 ocean loading files (last is optional - user own
    fprintf(' %s\n', lines2print{k});
    
    % if no file was given -> write to command window
    if isempty(files2load{k})
        fprintf('  No input file!\n')
    else
		% if textfile
		if strcmp(files2load{k}(end-3:end),'.TXT') || strcmp(files2load{k}(end-3:end),'.txt')
			fid=fopen(files2load{k});
            if (fid < 0)
               varargout{1} = ['ERROR: Cannot open the file: ', files2load{k}];
               return;
            end


			while ~feof(fid)
				line=[fgetl(fid), '        '];
				if strcmpi(line,  '        ')
					continue % empty lines are skipped!
				end
				
				% get data line
				if strcmp(line(1:2), '  ') && isnan(str2double(line(3:8)))
					% try to find station of line in ns_codes
					iStat=find(strcmpi({ns_codes.name}, line(3:10)));

					% if there was nothing found-> try to find '_' instad of ' '
					if isempty(iStat)
						iStat=find(strcmpi({ns_codes.name}, strrep(line(3:10), ' ', '_')));
					end

					% only if we have now found something
					if ~isempty(iStat)
						iStat=iStat(1);
						% preallocate
						ns_codes(iStat).ocean_loading.(ol_fieldnames{k})=zeros(6,11);
						
						% get following line (which is usually a $$ - comment - line)
						line=[fgetl(fid), '  '];

						% remove $$ lines (and empty lines!)
						while strcmp(line(1:2), '$$') || strcmpi(line,'  ')
							line=[fgetl(fid), '  '];
						end

						% get data
						ns_codes(iStat).ocean_loading.(ol_fieldnames{k})(1,:)=[str2num(line(2:8)) str2num(line(9:15)) str2num(line(16:22)) str2num(line(23:29)) str2num(line(30:36)) str2num(line(37:43)) str2num(line(44:50)) str2num(line(51:57))  str2num(line(58:64)) str2num(line(65:71)) str2num(line(72:78))]; % Richard Ray format
						line = [fgetl(fid),'                                                                                                          '];                   
						ns_codes(iStat).ocean_loading.(ol_fieldnames{k})(2,:)=str2num(line);
						line = [fgetl(fid),'                                                                                                          '];
						ns_codes(iStat).ocean_loading.(ol_fieldnames{k})(3,:)=str2num(line);
						line = [fgetl(fid),'                                                                                                          '];
						ns_codes(iStat).ocean_loading.(ol_fieldnames{k})(4,:)=str2num(line);
							line = [fgetl(fid),'                                                                                                          '];
						ns_codes(iStat).ocean_loading.(ol_fieldnames{k})(5,:)=str2num(line);

						line = [fgetl(fid),'                                                                                                          '];
						ns_codes(iStat).ocean_loading.(ol_fieldnames{k})(6,:)=str2num(line);
					else
						fprintf('%s (ocean_loading_%s.TXT) not found in ns_codes ('' '' -> ''_'' checked)\n', line(3:10), ol_fieldnames{k});
					end   
				end % we are at data line
			end % while
			fclose(fid);
		elseif strcmp(files2load{k}(end-3:end),'.mat')
			tempp=load(files2load{k});
            fprintf('ocean_loading_%s is introduced\n',ol_fieldnames{k})
            fn=fieldnames(tempp); 
            oceanLOAD=tempp.(fn{1});
            for i=1:length(oceanLOAD)
                if ~isempty(find(strcmpi({ns_codes.name}, oceanLOAD(i).ivsname)))
            		iStat=find(strcmpi({ns_codes.name}, oceanLOAD(i).ivsname));
                    ns_codes(iStat).ocean_loading.(ol_fieldnames{1,k})=oceanLOAD(i).load;
                end
            end
        end
    end
    
end


% --------------------------------
% 2.3 Terrestrial reference frames
% --------------------------------
fprintf('\n2.3 Terrestrial reference frames\n')


% --------------------------
% 2.3.1 ITRF2014-IVS-TRF.snx
% --------------------------
fprintf('\n2.3.1 ITRF2014-IVS-TRF.snx\n\n');
if exist(itrf2014File, 'file')
	itrf2014File_name = 'itrf2014';
	[ns_codes] = trf_by_snx_reader(ns_codes,itrf2014File,itrf2014File_name,break0);

    setStartEndEpoch('itrf2014');

    % add post-seismic deformation
    itrf2014psdfile='../TRF/data/ITRF2014-psd-vlbi.dat';
    if exist(itrf2014psdfile, 'file')
    % load file
        fid=fopen(itrf2014psdfile);
        if (fid < 0)
           varargout{1} = ['ERROR: Cannot open the file: ', itrf2014psdfile];
           return;
        end
    
        curLinePart=0;
        while ~feof(fid)
            curl=fgetl(fid); curLinePart=curLinePart+1;
        
            if curLinePart==4
                curLinePart=1;
            end
        
            if curLinePart==1
                code=str2double(curl(2:5));
                pt=curl(7:8);
                domes=curl(10:18);
                ep=date2mjd([str2double(curl(20:21)) 1 1])+...
                 str2double(curl(23:25))-1+str2double(curl(27:31))/60/60/24;
                e=sscanf(curl(35:min([72,length(curl)])), '%f')';
            elseif curLinePart==2
                n=sscanf(curl(35:min([72,length(curl)])), '%f')';
            else % last (UP) line
                u=sscanf(curl(35:min([72,length(curl)])), '%f')';
            
            % find proper station
                foundStationLog=strcmpi(domes,{ns_codes.domes});
                if sum(foundStationLog)==0
                % try to find code
                    if ~isnan(code) % just to be sure that not '----' or so is compared (and eventually found!!)
                        foundStationLog=strcmpi(num2str(code),{ns_codes.CDP});
                    end
                end
            
                % if (at least now) found
                if sum(foundStationLog)>0
                    foundStationInd=find(foundStationLog);
                    for iFoundStat=1:length(foundStationInd)
                        curFoundStat=foundStationInd(iFoundStat);
                        psdBreak=1;
                        if isfield(ns_codes(curFoundStat).itrf2014, 'psd')
                            psdBreak=length(ns_codes(curFoundStat).itrf2014.psd)+1;
                        end
                        ns_codes(curFoundStat).itrf2014.psd(psdBreak).epoch=ep;
                        ns_codes(curFoundStat).itrf2014.psd(psdBreak).e=e;
                        ns_codes(curFoundStat).itrf2014.psd(psdBreak).n=n;
                        ns_codes(curFoundStat).itrf2014.psd(psdBreak).u=u;
                    
                    end
                else
                    fprintf('Station %s (in ITRF2014 psd file) not found in ns_codes!\n',...
                        domes);
                end
            end
        end
        

        fclose(fid);    
    else
        fprintf('Warning: psd file of ITRF2014 not found\n(%s)\npause 2sek\n',...
            itrf2014psdfile);
        pause(2);
    end
else
    fprintf('ITRF2014-IVS-TRF.snx is not available\n\n');
    flag_trf_itrf2014 = false;
end

% --------------------------
% 2.3.2 DTRF2014
% --------------------------

fprintf('\n2.3.2 DTRF2014_VLBI.snx\n\n');
if exist(dtrf2014File, 'file')
	dtrf2014File_name = 'dtrf2014';
	[ns_codes] = trf_by_snx_reader(ns_codes,dtrf2014File,dtrf2014File_name,break0);

    setStartEndEpoch('dtrf2014');
else
    fprintf('DTRF2014_VLBI.snx is not available\n\n');
    flag_trf_dtrf2014 = false;
end


% ------------------------
% 2.3.3 VTRF2014_final.snx
% ------------------------
fprintf('\n2.3.3 Loading VTRF2014_final.snx\n\n');
if exist(vtrf2014File, 'file')

	vtrf2014File_name = 'vtrf2014';
	[ns_codes] = trf_by_snx_reader(ns_codes,vtrf2014File,vtrf2014File_name,break0);


    setStartEndEpoch('vtrf2014');
else
    fprintf('VTRF2014_final.snx is not available\n\n');
    flag_trf_vtrf2014 = false;
end
        
% --------------------------
% 2.3.4 IVS_TRF2014b.SSC.txt
% --------------------------
fprintf('\n2.3.4 Loading IVS_TRF2014b.SSC.txt\n\n');

if ~isempty(ivstrf2014bFile)
	ep=2005;
    [ns_codes] = ITRF_VLBI_SSC_reader(ivstrf2014bFile,ep,ns_codes,'ivsTrf2014b',break0);
else
    fprintf('IVS_TRF2014b.SSC.txt is not available\n\n');
    flag_trf_ivstrf2014b = false;
end

% --------------------
% 2.3.5 vieTrf13.txt
% --------------------

fprintf('\n2.3.5 Loading VieTRF13.txt\n\n');

if ~isempty(vieTrf13File)
	% open file -> read date -> close file
	fid=fopen(vieTrf13File);
    if (fid < 0)
       varargout{1} = ['ERROR: Cannot open the file: ', vieTrf13File];
       return;
    end
	vieTrf13Data=textscan(fid, '%8s      %13.4f   %13.4f   %13.4f        %7.4f     %7.4f     %7.4f      %5.0f   %5.0f   %5.0f', 'commentstyle', '%', 'delimiter', '||');
	fclose(fid);

	% for all lines of the vieTrf file
	for iLine=1:size(vieTrf13Data{1},1)
		
		% find index of current station (try also _ instead of ' '  and
		% vice-versa)
		iStat=strcmp({ns_codes.name}, vieTrf13Data{1}{iLine})|...
			strcmp({ns_codes.name}, strrep(vieTrf13Data{1}{iLine}, ' ', '_'))|...
			strcmp({ns_codes.name}, strrep(vieTrf13Data{1}{iLine}, '_', ' '));
		% if not found -> append to end of ns_codes
		if sum(iStat)==0
			indInNsCodes=size(ns_codes,2)+1;
		else
			indInNsCodes=find(iStat);
		end
		
		% for all indices where the name was found (always? only one!)
		for iStatFound=1:length(indInNsCodes)
			curInd=indInNsCodes(iStatFound);
			
			% +++ get break number +++
			% if we append station to ns_codes, it's also a new break
			if curInd>size(ns_codes,2)
				curBreak=1;
			% if there exists ns_codes(x).vieTrf (means there is already a
			% break
			elseif isfield(ns_codes(curInd), 'VieTRF13')
				if ~isempty(ns_codes(curInd).VieTRF13)
					curBreak=size(ns_codes(curInd).VieTRF13.break,2)+1;
				else
					curBreak=1;
				end
			else
				curBreak=1;
			end
			% --- get break number ---
			
			% preallocate break (first break only)
			if curBreak==1
				ns_codes(curInd).VieTRF13.break(1)=break0;
			end
			
			% write coordinates and velocities to ns_codes
			ns_codes(curInd).VieTRF13.break(curBreak).start=vieTrf13Data{9}(iLine);
			ns_codes(curInd).VieTRF13.break(curBreak).epoch=vieTrf13Data{8}(iLine);
			ns_codes(curInd).VieTRF13.break(curBreak).end=vieTrf13Data{10}(iLine);
			ns_codes(curInd).VieTRF13.break(curBreak).x=vieTrf13Data{2}(iLine);
			ns_codes(curInd).VieTRF13.break(curBreak).y=vieTrf13Data{3}(iLine);
			ns_codes(curInd).VieTRF13.break(curBreak).z=vieTrf13Data{4}(iLine);
			ns_codes(curInd).VieTRF13.break(curBreak).vx=vieTrf13Data{5}(iLine);
			ns_codes(curInd).VieTRF13.break(curBreak).vy=vieTrf13Data{6}(iLine);
			ns_codes(curInd).VieTRF13.break(curBreak).vz=vieTrf13Data{7}(iLine);
			
		end    
	end
else
    fprintf('VieTRF13.txt is not available\n\n');
    flag_trf_vieTrf13 = false;
end     


% -------------------
% 2.3.6 VieVS TRF file
% -------------------
% This TRF is the backup-TRF for VieVS and it is a mixture of the vieTrf
% and (if no coordinates are available there) the VTRF2008. Also in this
% backup-TRF the estimated coordinates after an earthquake are written to!

fprintf('\n2.3.6 Loading VieVS TRF file\n\n');

% read data using textscan
fid=fopen(vievsTrfFile);
if (fid < 0)
   varargout{1} = ['ERROR: Cannot open the file: ', vievsTrf];
   return;
end
vievsTrf=textscan(fid, '%8s %20f %16f %16f %15f %12f %12f %11f %8f %8f %f %s', 'commentstyle', '%', 'delimiter', '|');
fclose(fid);
nStat=size(vievsTrf{1},1); % get number of stations

% see if last column has same number of entries (otherwise add empty string)
if size(vievsTrf{12},1)<nStat
    vievsTrf{12}{nStat}='';
end

% add all stations to ns_codes struct
for k=1:nStat
    % 1. find station (of current line in vievstrf file) in ns_codes ...
    if sum(strcmpi({ns_codes.name}, vievsTrf{1}{k}))>0
        indStat=find(strcmpi({ns_codes.name}, vievsTrf{1}{k}));
    else
        % try to find the station with '_' instead of ' '
        if sum(strcmpi({ns_codes.name}, strrep(vievsTrf{1}{k}, ' ', '_')))>0
            indStat=find(strcmpi({ns_codes.name}, strrep(vievsTrf{1}{k}, ' ', '_')));
        else
            % if not found, add new entry to ns_codes
            fprintf('station %s (in vievsTRF) was not found in ns_codes -> writing station to new entry!\n', vievsTrf{1}{k});
            indStat=size(ns_codes,2)+1;
        end
    end
    
    % 2. if station has already been inserted in ns_codes -> add new break; else simplty add another break entry
    if size(ns_codes,2)>=indStat % ohterwise i need to create new entry (new station!)
        if isfield(ns_codes(indStat), 'vievsTrf') % needed for k==1: no vievsTrf field exists yet
            indStat=indStat(1);
            if isfield(ns_codes(indStat).vievsTrf, 'break')
                newBreakEntry=size(ns_codes(indStat).vievsTrf.break,2)+1; % 'break' field exists; get index of new break values
            else
                % nothing was written to current station
            %    ns_codes(indStat).vievsTrf.name=vievsTrf{1}{k}; % write station name only once
                newBreakEntry=1;
            end
        else
            % nothing was written to current station
        %    ns_codes(indStat).vievsTrf.name=vievsTrf{1}{k}; % write station name only once
            newBreakEntry=1;
        end
    else
        newBreakEntry=1;
    end

    % 3. write x,y,z, vx,vy,vz of current break to ns_codes
    ns_codes(indStat).vievsTrf.break(newBreakEntry).x=vievsTrf{2}(k);
    ns_codes(indStat).vievsTrf.break(newBreakEntry).y=vievsTrf{3}(k);
    ns_codes(indStat).vievsTrf.break(newBreakEntry).z=vievsTrf{4}(k);
    
    ns_codes(indStat).vievsTrf.break(newBreakEntry).vx=vievsTrf{5}(k);
    ns_codes(indStat).vievsTrf.break(newBreakEntry).vy=vievsTrf{6}(k);
    ns_codes(indStat).vievsTrf.break(newBreakEntry).vz=vievsTrf{7}(k);
    
    ns_codes(indStat).vievsTrf.break(newBreakEntry).epoch=vievsTrf{8}(k);
    ns_codes(indStat).vievsTrf.break(newBreakEntry).start=vievsTrf{9}(k);
    ns_codes(indStat).vievsTrf.break(newBreakEntry).end=vievsTrf{10}(k);
    
    % also write ".indatum" to struct (0 for not in datum, eg earthquake), otherwise 1
    vievsTrf{11}(isnan(vievsTrf{11}))=1;
    ns_codes(indStat).vievsTrf.break(newBreakEntry).indatum=vievsTrf{11}(k);
    ns_codes(indStat).vievsTrf.break(newBreakEntry).comment=vievsTrf{12}{k};

    % if there is no name in ns_codes for current station (usually when it did not exist before), add name
    if isempty(ns_codes(indStat).name)
        ns_codes(indStat).name=vievsTrf{1}{k};
    end
    
end


% since we want all ns-codes stations to have vievsTrf coordinates: Write
% message if not.
noVievsTrfCoords=cellfun(@isempty, {ns_codes.vievsTrf});

if sum(noVievsTrfCoords)>0 % if there are ns_codes stations with no vievsTrf coordinates
    
    % get indices of stations where we don't have vievsTrf coordinates.
    ind=find(noVievsTrfCoords);
    
    % write all those stations as one line
    for k=1:sum(noVievsTrfCoords)
        if ~isempty(ns_codes(ind(k)).VieTRF13)
            % try to get VieTRF13 coordinates
            ns_codes(ind(k)).vievsTrf.break=ns_codes(ind(k)).VieTRF13.break;
        end
        
    end
end

% -----------------------
% 2.3.9 User own TRF file
% -----------------------
if ~isempty(userOwnTrfFile)
    fprintf('\n2.3.9 User own TRF file (%s)\n\n', userOwnTrfFile);

    if strcmp(userOwnTrfFile(end-3:end),'.TXT') || strcmp(userOwnTrfFile(end-3:end),'.txt')
        % open file -> read date -> close file
        fid=fopen(userOwnTrfFile);
        if (fid < 0)
            varargout{1} = ['ERROR: Cannot open the file: ', userOwnTrfFile];
            return;
        end
        userOwnTrfData=textscan(fid, '%8c %f %f %f %f %f %f %f %f %f', 'commentstyle', '%');
        fclose(fid);
    elseif strcmp(userOwnTrfFile(end-3:end),'.mat')
        tempt=load(userOwnTrfFile);
		fnt=fieldnames(tempt); 
        position=tempt.(fnt{1});
        % format / stationname x y z vx vy vz epoch start end
        % station names
        userOwnTrfData{1}={position.ivsname}.';
        % coordinates
        userOwnTrfData{2}={position.x}.';
        userOwnTrfData{3}={position.y}.';
        userOwnTrfData{4}={position.z}.';
        % epoch
        userOwnTrfData{8}={position.mjd}.';
        % velocity
		if ~isfield(position,'vx')
			userOwnTrfData{5}={};
		else
			userOwnTrfData{5}={position.vx};
		end
		if ~isfield(position,'vy')
			userOwnTrfData{6}={};
		else
			userOwnTrfData{5}={position.vy};
		end
		if ~isfield(position,'vz')
			userOwnTrfData{7}={};
		else
			userOwnTrfData{5}={position.vz};
		end
        % star and end date
		if ~isfield(position,'start')
        	userOwnTrfData{9}={};
		else
			userOwnTrfData{5}={position.start};
		end
		if ~isfield(position,'end')
			userOwnTrfData{10}={};
		else
			userOwnTrfData{5}={position.end};
		end
    end
    [ns_codes] = ownTRFread(ns_codes,userOwnTrfData);
end


% ----------------------
% 2.4 Atmosphere loading
% ----------------------
fprintf('\n2.4 Atmosphere loading\n')

% ----------------------------
% 2.4.1 s1_s2_s3_cm_noib_vlbi.dat (VIENNA)
% ----------------------------
fprintf('\n2.4.1 Loading s1_s2_s3_cm_noib_vlbi.dat (VIENNA)\n\n');
% for i=1:length(ns_codes)
%     ns_codes(i).atmosphere_tidal_loading.vienna=[];
% end
if isempty(s123viennaFile)
    varargout{1} = 'ERROR: Location of s1_s2_s3_cm_noib_vlbi.dat not specified!';
    return;
end
fid = fopen(s123viennaFile);
if (fid < 0)
    varargout{1} = ['ERROR: Cannot open the file: ', s123viennaFile];
    return;
end


while ~feof(fid)
    line = [fgetl(fid),'                                                                                                          '];
    linum=str2num(line(9:end));
    % find station index
    iStat=find(strcmpi({ns_codes.name}, line(1:8)));
    
    for iCorrectStat=1:size(iStat,2)
        % get current index of station (could be more than one I want to
        % write atmo loading values to
        curIStat=iStat(iCorrectStat);
    
        % write values to struct
        ns_codes(curIStat).atmosphere_tidal_loading.vienna=zeros(1,12);
        ns_codes(curIStat).atmosphere_tidal_loading.vienna=linum([4 5 6 7 16 17 18 19 10 11 12 13]); %RNE 
    end
end
fclose(fid);






% ---------------------------
% 2.4.2 s12_cm_noib_gsfc.dat
% ---------------------------

fprintf('\n2.4.2 Loading s12_cm_noib_gsfc.dat\n\n');

if isempty(s123viennaFile)
    varargout{1} = 'ERROR: Location of s12_cm_noib_gsfc.dat not specified!';
    return;
end
fid = fopen(s12gsfcFile);
if (fid < 0)
    varargout{1} = ['ERROR: Cannot open the file: ', s12gsfcFile];
    return;
end

while ~feof(fid)
    
    % get line of textfile
    line = [fgetl(fid),'                                                                                                          '];
    
    % find station index
    iStat=find(strcmpi({ns_codes.name}, line(1:8)));
    
    for iCorrectStat=1:size(iStat,2)
        % get current index of station (could be more than one I want to
        % write atmo loading values to
        curIStat=iStat(iCorrectStat);
    
        % write values to struct
        ns_codes(curIStat).atmosphere_tidal_loading.gsfc=zeros(1,12);
        ns_codes(curIStat).atmosphere_tidal_loading.gsfc=str2num(line(34:end));
            
    end
end
fclose(fid);

% ----------------------------
% 2.4.3 s12_cm_noib_vandam.dat
% ----------------------------

fprintf('\n2.4.3 Loading s12_cm_noib_vandam.dat\n\n');

if isempty(s123viennaFile)
    varargout{1} = 'ERROR: Location of s12_cm_noib_vandam.dat not specified!';
    return;
end
fid = fopen(s12vandamFile);
if (fid < 0)
    varargout{1} = ['ERROR: Cannot open the file: ', s12vandamFile];
    return;
end

while ~feof(fid)
    line = [fgetl(fid),'                                                                                                          '];
    
    % find station index
    iStat=find(strcmpi({ns_codes.name}, line(1:8)));
    
    for iCorrectStat=1:size(iStat,2)
        % get current index of station (could be more than one I want to
        % write atmo loading values to
        curIStat=iStat(iCorrectStat);
    
        % write values to struct
        ns_codes(curIStat).atmosphere_tidal_loading.vandam=zeros(1,12);
        ns_codes(curIStat).atmosphere_tidal_loading.vandam=str2num(line(34:end));
    end
end
fclose(fid);



% ---------------------------------
% 2.4.4 User own atmosphere loading
% ---------------------------------

if ~isempty(s12UserOwnFile)
    fprintf('\n2.4.4 User own atmosphere loading data (%s)\n\n', s12UserOwnFile);
    fid = fopen(s12UserOwnFile);
    if (fid < 0)
        varargout{1} = ['ERROR: Cannot open the file: ', s12UserOwnFile];
        return;
    end

    while ~feof(fid)
        line = [fgetl(fid),'                                                                                                          '];

        % find station index
        iStat=find(~cellfun(@isempty, strfind({ns_codes.name}, line(1:8))));

        for iCorrectStat=1:size(iStat,2)
            % get current index of station (could be more than one I want to
            % write atmo loading values to
            curIStat=iStat(iCorrectStat);

            % write values to struct
            ns_codes(curIStat).atmosphere_tidal_loading.(s12UserOwnFieldname)=zeros(1,12);
            ns_codes(curIStat).atmosphere_tidal_loading.(s12UserOwnFieldname)=str2num(line(34:end));
        end
    end
    fclose(fid);
end


% ---------------------------
% 2.5 Ocean Pole Tide Loading
% ---------------------------
fprintf('2.5 Ocean pole tide loading\n')


% ------------------------------------------
% 2.5.1 Ocean Pole Tide Loading, model desai
% ------------------------------------------

fprintf('\n2.5.1 Ocean Pole Tide Loading (model by desai)\n\n');
if ~isempty(opoleloadUserOwnFilename)
    fprintf('\n2.5.2 Ocean Pole Tide Loading (user own - %s)\n\n', opoleloadUserOwnFilename);
    userOwnOptl=1;
else
    userOwnOptl=0;
end

% open file
if isempty(opoleloadcoefcmcorFile)
    varargout{1} = 'ERROR: Ocean pole tide loading file is not specified!';
    return;
end
fid=fopen(opoleloadcoefcmcorFile);
if (fid < 0)
    varargout{1} = ['ERROR: Cannot open file: ', opoleloadcoefcmcorFile];
    return;
end

% get data
oceanPoleTideCoefs=textscan(fid, '%f %f %f %f %f %f %f %f', 'headerlines', 14);

% close file
fclose(fid);


% see if this works faster without a loop:
latsForGrid=unique(oceanPoleTideCoefs{2});
lonsForGrid=unique(oceanPoleTideCoefs{1});

% get grid out of these vectors
[lon_meshgrid, lat_meshgrid]=meshgrid(lonsForGrid, latsForGrid);

% preallocate
empty_grid=zeros(size(lon_meshgrid));
urR=empty_grid'; urI=empty_grid';
unR=empty_grid'; unI=empty_grid';
ueR=empty_grid'; ueI=empty_grid';

urR(:)=oceanPoleTideCoefs{3};
urI(:)=oceanPoleTideCoefs{4};
unR(:)=oceanPoleTideCoefs{5};
unI(:)=oceanPoleTideCoefs{6};
ueR(:)=oceanPoleTideCoefs{7};
ueI(:)=oceanPoleTideCoefs{8};
urR=urR';
urI=urI';
unR=unR';
unI=unI';
ueR=ueR';
ueI=ueI';

if userOwnOptl
    fid2=fopen(oceanLoadingUserOwnFilename);
    oceanPoleTideCoefs_userOwn=textscan(fid2, '%f %f %f %f %f %f %f %f', 'headerlines', 14);
    fclose(fid2);
    
    latsForGrid_userOwn=unique(oceanPoleTideCoefs_userOwn{2});
    lonsForGrid_userOwn=unique(oceanPoleTideCoefs_userOwn{1}); 
    
    [lon_meshgrid_userOwn, lat_meshgrid_userOwn]=meshgrid(lonsForGrid_userOwn, latsForGrid_userOwn);
    
    % preallocate user own
    empty_grid_userOwn=zeros(size(lon_meshgrid_userOwn));
    urR_userOwn=empty_grid_userOwn'; urI_userOwn=empty_grid_userOwn';
    unR_userOwn=empty_grid_userOwn'; unI_userOwn=empty_grid_userOwn';
    ueR_userOwn=empty_grid_userOwn'; ueI_userOwn=empty_grid_userOwn';

    urR_userOwn(:)=oceanPoleTideCoefs_userOwn{3};
    urI_userOwn(:)=oceanPoleTideCoefs_userOwn{4};
    unR_userOwn(:)=oceanPoleTideCoefs_userOwn{5};
    unI_userOwn(:)=oceanPoleTideCoefs_userOwn{6};
    ueR_userOwn(:)=oceanPoleTideCoefs_userOwn{7};
    ueI_userOwn(:)=oceanPoleTideCoefs_userOwn{8};
    urR_userOwn=urR_userOwn';
    urI_userOwn=urI_userOwn';
    unR_userOwn=unR_userOwn';
    unI_userOwn=unI_userOwn';
    ueR_userOwn=ueR_userOwn';
    ueI_userOwn=ueI_userOwn';

end
% Parameter des GRS80
a = 6378137;
b = 6356752.3141;

% for all stations -> interpolate values and write to struct
for k=1:size(ns_codes,2)
    curLon=[];
    curLat=[];
    coordsFound=0;
    % get lon and lat of current station
    % see if we have ITRF2014 coordinates
    if ~isempty(ns_codes(k).itrf2014) 
        coordsFound=1;
        curX=ns_codes(k).itrf2014.break(size(ns_codes(k).itrf2014.break,2)).x;
        curY=ns_codes(k).itrf2014.break(size(ns_codes(k).itrf2014.break,2)).y;
        curZ=ns_codes(k).itrf2014.break(size(ns_codes(k).itrf2014.break,2)).z;
    elseif ~isempty(ns_codes(k).vtrf2014) % if we have VTRF 2014 coords
        coordsFound=1;
        curX=ns_codes(k).vtrf2014.break(size(ns_codes(k).vtrf2014.break,2)).x;
        curY=ns_codes(k).vtrf2014.break(size(ns_codes(k).vtrf2014.break,2)).y;
        curZ=ns_codes(k).vtrf2014.break(size(ns_codes(k).vtrf2014.break,2)).z;
    elseif ~isempty(ns_codes(k).dtrf2014) % if we have DTRF 2014 coords
        coordsFound=1;
        curX=ns_codes(k).dtrf2014.break(size(ns_codes(k).dtrf2014.break,2)).x;
        curY=ns_codes(k).dtrf2014.break(size(ns_codes(k).dtrf2014.break,2)).y;
        curZ=ns_codes(k).dtrf2014.break(size(ns_codes(k).dtrf2014.break,2)).z;
    elseif ~isempty(ns_codes(k).vievsTrf)
        coordsFound=1;
        curX=ns_codes(k).vievsTrf.break(size(ns_codes(k).vievsTrf.break,2)).x;
        curY=ns_codes(k).vievsTrf.break(size(ns_codes(k).vievsTrf.break,2)).y;
        curZ=ns_codes(k).vievsTrf.break(size(ns_codes(k).vievsTrf.break,2)).z;
    end
    
    % if we have found coords -> interpolate and save to struct
    if coordsFound==1
        % if we have xyz found, calc lat/lono from it
        if isempty(curLon)
            [curLat, curLon, curHeight] = kart2ell([curX, curY, curZ],a,b);
            curLat=curLat*180/pi;
            curLon=curLon*180/pi;
            if curLon<0
                curLon=curLon+360;
            end
        end
        
        % interpolate ocean pole tide loading coeffs to station lon/lat
        curUrR=interp2(lon_meshgrid, lat_meshgrid, urR, curLon, curLat);
        curUrI=interp2(lon_meshgrid, lat_meshgrid, urI, curLon, curLat);
        curUnR=interp2(lon_meshgrid, lat_meshgrid, unR, curLon, curLat);
        curUnI=interp2(lon_meshgrid, lat_meshgrid, unI, curLon, curLat);
        curUeR=interp2(lon_meshgrid, lat_meshgrid, ueR, curLon, curLat);
        curUeI=interp2(lon_meshgrid, lat_meshgrid, ueI, curLon, curLat);
 
        allVarsOneVector=[curUrR, curUrI, curUnR, curUnI, curUeR, curUeI];
       
        if userOwnOptl
            curUrR_userOwn=interp2(lon_meshgrid_userOwn, lat_meshgrid_userOwn, urR_userOwn, curLon, curLat);
            curUrI_userOwn=interp2(lon_meshgrid_userOwn, lat_meshgrid_userOwn, urI_userOwn, curLon, curLat);
            curUnR_userOwn=interp2(lon_meshgrid_userOwn, lat_meshgrid_userOwn, unR_userOwn, curLon, curLat);
            curUnI_userOwn=interp2(lon_meshgrid_userOwn, lat_meshgrid_userOwn, unI_userOwn, curLon, curLat);
            curUeR_userOwn=interp2(lon_meshgrid_userOwn, lat_meshgrid_userOwn, ueR_userOwn, curLon, curLat);
            curUeI_userOwn=interp2(lon_meshgrid_userOwn, lat_meshgrid_userOwn, ueI_userOwn, curLon, curLat);
            
             allVarsOneVector_userOwn=[curUrR_userOwn, curUrI_userOwn, ...
                curUnR_userOwn, curUnI_userOwn, curUeR_userOwn, curUeI_userOwn];
        end
        
        if isnan(sum(allVarsOneVector))
            fprintf('warning: one value of ocean pole tide loading is NaN\n')
%             keyboard;
        end
        
        % write data to ns_codes struct
        ns_codes(k).oceanPoleTideLoading.desai=allVarsOneVector;
        if userOwnOptl
            ns_codes(k).oceanPoleTideLoading.(opoleloadUserOwnFieldname)=allVarsOneVector_userOwn;
        end
    else
        curN='--------';
        curCode='--';
        if ~isempty(ns_codes(k).name); curN=ns_codes(k).name; end
        if ~isempty(ns_codes(k).code); curCode=strtrim(ns_codes(k).code); end
        fprintf('%s (%s): No coords found -> No Ocean Pole Tide Loading calculated!\n', curN, curCode);
    end
end


% ---------------------------
% 2.6 APL Regression Coefficients
% ---------------------------
% added by Hana 01/2014
fprintf('2.6 APL Regression Coefficients\n')


% file 1
[Tlscp, RefPres, RgC] = textread(aplrgFilename1, '%8c %f %f', 'commentstyle', 'matlab');

for ian = 1:size(Tlscp,1)
    % find telescope name in ns_codes
    iStat=[];
    if ~strcmp(Tlscp(ian,:),'--------')
        iStat=find(strcmpi({ns_codes.name}, Tlscp(ian,:)));
    end 
    if ~isempty(iStat)
        ns_codes(iStat(1)).aplrg.(aplrgFieldname1) = [RefPres(ian), RgC(ian)];
    end
end

% file 2
Tlscp=[]; RefPres=[]; RgC=[];
[Tlscp, RefPres, RgC] = textread(aplrgFilename2, '%8c %f %f', 'commentstyle', 'matlab');

for ian = 1:size(Tlscp,1)
    % find telescope name in ns_codes
    iStat=[];
    if ~strcmp(Tlscp(ian,:),'--------')
        iStat=find(strcmpi({ns_codes.name}, Tlscp(ian,:)));
    end 
    if ~isempty(iStat)
        ns_codes(iStat(1)).aplrg.(aplrgFieldname2) = [RefPres(ian), RgC(ian)];
    end
end



% ---------------------------
% 2.7 GIA Uplift Rates
% ---------------------------
% added by Hana 01/2014
fprintf('2.7 GIA Uplift Rates\n')


for k=1:length(ns_codes)
    ns_codes(k).gia.(giaMJDFieldname1) = [];
    ns_codes(k).gia.(giaFieldname1) = [];
end



% file with the refence epoch
Tlscp=[];
[Tlscp, GiaMJD0] = textread(giaMJDFilename1, '%8c %f', 'commentstyle', 'matlab');

for ian = 1:size(Tlscp,1)
    % find telescope name in ns_codes
    iStat=[];
    if ~strcmp(Tlscp(ian,:),'--------')
        iStat= find(strcmp({ns_codes.name}, Tlscp(ian,:))|...
            strcmp({ns_codes.name}, strrep(Tlscp(ian,:), ' ', '_'))|...
            strcmp({ns_codes.name}, strrep(Tlscp(ian,:), '_', ' ')));
    end 
    if ~isempty(iStat)
        ns_codes(iStat(1)).gia.(giaMJDFieldname1) = GiaMJD0(ian);
    else
        ns_codes(iStat(1)).gia.(giaMJDFieldname1) = [];
    end
end


% file with GIA uplift rates
Tlscp=[];
[Tlscp, xx1, xx2, GiaUp] = textread(giaFilename1, '%8c %f %f %f', 'commentstyle', 'matlab');
xx1=[];
xx2=[];

for ian = 1:size(Tlscp,1)
    % find telescope name in ns_codes
    iStat=[];
    if ~strcmp(Tlscp(ian,:),'--------')
%         iStat=find(~cellfun(@isempty, strfind({ns_codes.name}, Tlscp(ian,:))))
        
        iStat= find(strcmp({ns_codes.name}, Tlscp(ian,:))|...
            strcmp({ns_codes.name}, strrep(Tlscp(ian,:), ' ', '_'))|...
            strcmp({ns_codes.name}, strrep(Tlscp(ian,:), '_', ' ')));
        
        
    end 
    if ~isempty(iStat)
        ns_codes(iStat(1)).gia.(giaFieldname1) = GiaUp(ian)/1000; %m/y
        if isempty(ns_codes(iStat(1)).gia.(giaMJDFieldname1))
            ns_codes(iStat(1)).gia.(giaMJDFieldname1) = 51544; %mjd; % write the reference epoch 1.1.2000 if there is no other information
        end
    else
        ns_codes(iStat(1)).gia.(giaFieldname1) = [];
    end
end



%%


%% ========================================================================
% 3. get atmosphere ocean loading from grids for stations that were not in
% the default files
% =========================================================================

fprintf('\n3. Interpolating (linear) atmospheric loading grids (vandam,\ngsfc) for stations not in list\n\n');

% get logicals for ns_codes-entries that have no vandam/gsfc/dudy field
a={ns_codes(:).atmosphere_tidal_loading};
noVienna=~cellfun(@(x) isfield(x, 'vienna'), a);
noVandam=~cellfun(@(x) isfield(x, 'vandam'), a);
noGsfc=~cellfun(@(x) isfield(x, 'gsfc'), a);
noAtmoTidLoad=cellfun(@isempty, a);

% if there are stations with no vandam/gsfc atmo loading data
if sum(noVandam|noGsfc|noVienna|noAtmoTidLoad)>0
    % get indices where we have to do something
    indDoAtmoLoading=find(noVienna|noVandam|noGsfc|noAtmoTidLoad);
    
    % only do once: load grid data to
    needToGetViennaGridToProperFormat=1;
    needToGetVandamGridToProperFormat=1;
    needToGetGsfcGridToProperFormat=1;
    
    % for all stations where at least one loading is missing
    for k=1:size(indDoAtmoLoading,2)
        % get index in ns_codes for station I want to calc atmo ocean load
        indInNs_codes=indDoAtmoLoading(k);
        
%         if indInNs_codes==235
%             keyboard;
%         end
        
        % try to get coords
        
        % define boolean if we have found coordinates
        coordsYes=0;

        % preallocate coordinates
        curX=0;
        curY=0;
        curZ=0;

        % if we have an entry for ITRF2014 - use these coordinates
        if ~isempty(ns_codes(indInNs_codes).itrf2014)
            curX=ns_codes(indInNs_codes).itrf2014.break(size(ns_codes(indInNs_codes).itrf2014.break,2)).x;
            curY=ns_codes(indInNs_codes).itrf2014.break(size(ns_codes(indInNs_codes).itrf2014.break,2)).y;
            curZ=ns_codes(indInNs_codes).itrf2014.break(size(ns_codes(indInNs_codes).itrf2014.break,2)).z;
            coordsYes=1;
        elseif ~isempty(ns_codes(indInNs_codes).vtrf2014) % next try: maybe there are coordinates available for vtrf2014:
            curX=ns_codes(indInNs_codes).vtrf2014.break(size(ns_codes(indInNs_codes).vtrf2014.break,2)).x;
            curY=ns_codes(indInNs_codes).vtrf2014.break(size(ns_codes(indInNs_codes).vtrf2014.break,2)).y;
            curZ=ns_codes(indInNs_codes).vtrf2014.break(size(ns_codes(indInNs_codes).vtrf2014.break,2)).z;
            coordsYes=1;
        elseif ~isempty(ns_codes(indInNs_codes).vievsTrf)
            curX=ns_codes(indInNs_codes).vievsTrf.break(size(ns_codes(indInNs_codes).vievsTrf.break,2)).x;
            curY=ns_codes(indInNs_codes).vievsTrf.break(size(ns_codes(indInNs_codes).vievsTrf.break,2)).y;
            curZ=ns_codes(indInNs_codes).vievsTrf.break(size(ns_codes(indInNs_codes).vievsTrf.break,2)).z;
            coordsYes=1;
        end

        % if we have coordinates -> get data
        if coordsYes==1
            % xyz -> lat lon h
            [curLat, curLon, curH]=xyz2ell([curX, curY, curZ]);
            
            % [-180,180] -> [0,360]
            curLon(curLon<0)=curLon(curLon<0)+2*pi;
            
            curLat=curLat*180/pi;
            curLon=curLon*180/pi;
            
            % if we need vienna data
            % ======================
            if noVienna(indInNs_codes)==1
                if ~exist(atmoTidalLoadViennaGridFile, 'file')
                    fprintf('Vienna grid file (%s) does not exist\n-> get it!\n', atmoTidalLoadViennaGridFile);
                else
                    % vienna grid file exists:
                    if needToGetViennaGridToProperFormat==1
                        % get data
                        fid=fopen(atmoTidalLoadViennaGridFile);
                        if (fid < 0)
                            varargout{1} = ['ERROR: Cannot open the file: ', atmoTidalLoadViennaGridFile];
                            return;
                        end
                        viennaGrid=cell2mat(textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f'));
                        fclose(fid);
                        lonV=0.5:359.5;
                        latV=89.5:-1:-89.5;
                        [vieLonMesh,vieLatMesh]=meshgrid(lonV,latV);
                        needToGetViennaGridToProperFormat=0; % do not load it again!
                        viennaData=zeros(length(latV),length(lonV),18); % 18 parameters
                        for aplViePar=1:18
                            viennaData(:,:,aplViePar)=...
                                reshape(viennaGrid(:,2+aplViePar),length(lonV),length(latV))'; % same order as textfile (columns) 
                           
                        end
                        aplVieInd=[4 5 6 7 16 17 18 19 10 11 12 13]-3;
                    end
                    % interpolate values
                    temp_vec=zeros(1,12);
                    for vieAplInd=1:12
                        temp_vec(1,vieAplInd)=...
                            interp2(vieLonMesh,vieLatMesh,viennaData(:,:,aplVieInd(vieAplInd)), curLon, curLat, 'spline'); % spline is needed because curLon/curLat might lie outside of mesh -> linear would return NaN!
                    end
                    ns_codes(indInNs_codes).atmosphere_tidal_loading.vienna=temp_vec; %RNE 
                end
            end
            
            % if we need vandam data
            % ======================
            if noVandam(indInNs_codes)==1
                % do once: loading data
                if ~exist(atmoTidalLoadVandamGridFile, 'file')
                    % url where file can be downloaded
                    url='http://geophy.uni.lu/applications/atm1/download/s1_s2_def_cm.dat.Z';

                    % write message to command window
                    fprintf('Downloading vandam grid file\n');

                    % download .Z file
                    urlwrite(url, '../TRF/data/s1_s2_def_cm.dat.Z');

                    % write message to user to extract file
                    fprintf('Grid file (s1_s2_def_cm.dat.Z) for atmospheric tidal loading not found!\nIt was downloaded to ../TRF/data.\nPlease extract manually and re-run this program (or type ''return'').\n');
                    keyboard;
                end
                % if we need to get the vandam grid into proper format
                % (only do once)
                if needToGetVandamGridToProperFormat==1
                    % get data
                    fid=fopen('../TRF/data/s1_s2_def_cm.dat');
                    if (fid < 0)
                        varargout{1} = ['ERROR: Cannot open the file: ../TRF/data/s1_s2_def_cm.dat'];
                        return;
                    end
                    vandamGrid=textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f');
                    fclose(fid);

                    % create grid for interpolation
                    [meshX,meshY] = meshgrid(0:360,90:-1:-90);

                    % create data grid
                    data=zeros(size(meshX));

                    % done (only do once): load grid vandam data
                    needToGetVandamGridToProperFormat=0;

                end
                % interpolate using matlab's 3d spline interpolation
                % for all 12 columns in grid file (get t1, t2, t3, ... t12)
                ns_codes(indInNs_codes).atmosphere_tidal_loading.vandam=zeros(1,12);
                for iCol=3:14

                    data(:)=vandamGrid{iCol};
                    ns_codes(indInNs_codes).atmosphere_tidal_loading.vandam(iCol-2)=interp2(meshX, meshY, data, curLon, curLat, 'linear');
                end
                
            end % if we need vandam data
            
            % if we need Gsfc data
            % ======================
            if noGsfc(indInNs_codes)==1
                % do once: loading data
                if ~exist('../TRF/data/aplo_s1_s2_noib_1.0x1.0deg.nc', 'file')
                    % download
                    url='http://gemini.gsfc.nasa.gov/aplo/aplo_s1_s2_noib_1.0x1.0deg.nc';
                    
                    % write message to command window
                    fprintf('Downloading gsfc grid file (aplo_s1_s2_noib_1.0x1.0deg.nc)\n');

                    % download .Z file
                    urlwrite(url, '../TRF/data/aplo_s1_s2_noib_1.0x1.0deg.nc');
                end
                
                % load data to proper format ONLY ONCE
                if needToGetGsfcGridToProperFormat ==1
                    [mesh_lon, mesh_lat, a{1}, a{2}, a{3}, ...
                        a{4}, a{5}, a{6}, a{7}, a{8}, a{9}, ...
                        a{10}, a{11}, a{12}]=getGsfcsAploGrid('../TRF/data/aplo_s1_s2_noib_1.0x1.0deg.nc');
                    needToGetGsfcGridToProperFormat=0;
                end
                
                % interpolate nad write to struct
                ns_codes(indInNs_codes).atmosphere_tidal_loading.gsfc=zeros(1,12);
                for ind_k=1:12
                    ns_codes(indInNs_codes).atmosphere_tidal_loading.gsfc(ind_k)=interp2(mesh_lon, mesh_lat, a{ind_k}, curLon, curLat, 'linear');
                end
                
            end % if we need gsfc loading correction
        end % =if we have coordinates found
        
    end % for all stations we need loading corrections
end % if we need any loading for any station (which was not in file) at all
   
     
%% ===========================================
% 4. write to command window for ocean loading
% ============================================
fprintf('\n4.1 Writing stations for ocean tide loading (oso) which HAVE NO OCEAN TIDE LOADING\nfor full list (all stations) see below\n\n');
osoFormat='%-24s %16f %15f %15f\n';

tmp_superstatCell=struct2cell(ns_codes);
% get field-index (actually: logicals) of ocean loading in superstation struct
oceanFieldIndex=strcmp(fieldnames(ns_codes), 'ocean_loading');
% get stations which have no ocean loading corrections
noOceanLoadCorr=cellfun('isempty', tmp_superstatCell(oceanFieldIndex,:,:));
noOceanLoadCorr=noOceanLoadCorr(:);
noOceanLoadCorrInd=find(noOceanLoadCorr);

% get number of blocks
nBlocks=ceil(sum(noOceanLoadCorr)/100);

for iBlock=1:nBlocks
    % get first and last station of current block
    firstStat=1+(iBlock-1)*100;
    lastStat=iBlock*100;
    if lastStat>sum(noOceanLoadCorr)
        lastStat=sum(noOceanLoadCorr);
    end
    
    % write header lines
    fprintf('//Stations %1.0f - %1.0f (if coords available) which don''t have any ocean tide loading corrections\n// -> http://holt.oso.chalmers.se/loading//)\n', firstStat, lastStat);
    fprintf('//Name of station_______| Longitude (deg)| Latitude (deg)| Height (m)  	  \n//Name of station_______|     X (m)      |     Y (m)     |     Z (m)      \n');
    for k=firstStat:lastStat
        curInd=noOceanLoadCorrInd(k);
        curX=0;
        curY=0;
        curZ=0;
        coordsYes=0;
        
        trfsToGetCoords={'itrf2014', 'vtrf2014', 'dtrf2014', 'vievsTrf', 'VieTRF13'}; % get xyz from one of these
        for kTrf=1:length(trfsToGetCoords)
            if ~isempty(ns_codes(curInd).(trfsToGetCoords{kTrf}))
                curX=ns_codes(curInd).(trfsToGetCoords{kTrf}).break(size(ns_codes(curInd).(trfsToGetCoords{kTrf}).break,2)).x;
                curY=ns_codes(curInd).(trfsToGetCoords{kTrf}).break(size(ns_codes(curInd).(trfsToGetCoords{kTrf}).break,2)).y;
                curZ=ns_codes(curInd).(trfsToGetCoords{kTrf}).break(size(ns_codes(curInd).(trfsToGetCoords{kTrf}).break,2)).z;
                coordsYes=1;
                break
            end
        end
        
        if coordsYes==1
            fprintf(osoFormat, ns_codes(curInd).name, curX, curY, curZ);
        end
    end
end
        
% 4.2 Stations which have ocean loading corrections for some (but not all)
% models
fprintf('\n4.2 Writing stations which have ocean loading corrections \nfor some (but not all) models\n\n');
for iStat=1:length(ns_codes)
    if ~isempty(ns_codes(iStat).ocean_loading)
        modelsEmpty=''; % string containing all empty models
        emptyModelFound=0;
        for iModel=1:length(ol_fieldnames)
            % just to be sure: add field if not exists...
            if ~isfield(ns_codes(iStat).ocean_loading, ol_fieldnames{iModel})
            	ns_codes(iStat).ocean_loading.(ol_fieldnames{iModel})=[];
            end
            if isempty(ns_codes(iStat).ocean_loading.(ol_fieldnames{iModel}))
                if emptyModelFound==1
                    modelsEmpty=[modelsEmpty,', '];
                end
                modelsEmpty=[modelsEmpty, ol_fieldnames{iModel}];
                emptyModelFound=1;
            end
   
        end
        if emptyModelFound==1
            fprintf('%-8s %s\n', ns_codes(iStat).name, modelsEmpty);
        end
    end
end
        


fprintf('\n4.3 Writing ALL stations for ocean tide loading (oso)\n\n');

% write text to command window as would be needed for ocean loading http://froste.oso.chalmers.se/loading//
% plot stations in blocks of 100
% preallocate cell array for stations where no coordinates were found
noCoordsStats=cell(size(ns_codes,2),4);

% get number of blocks (eg. 211 stations: 3 blocks)
nBlocks=ceil(size(ns_codes,2)/100);
for iBlock=1:nBlocks
    % get first and last station of block (1st block: 1 - 100, 101 - 200,
    % ...)
    firstStat=1+(iBlock-1)*100;
    lastStat=iBlock*100;
    % if number of station is larger -> use this number instead
    if lastStat>size(ns_codes,2)
        lastStat=size(ns_codes,2);
    end
    % write headerlines to command window
    fprintf('//Stations %1.0f - %1.0f (if coords available) for tidal ocean loading\n// -> http://froste.oso.chalmers.se/loading//)\n', firstStat, lastStat);
    fprintf('//Name of station_______| Longitude (deg)| Latitude (deg)| Height (m)  	  \n//Name of station_______|     X (m)      |     Y (m)     |     Z (m)      \n');
    
    for k=firstStat:lastStat
        % get coordinates (depending on what is available -> if new station: maybe we just have a our own coordinates (and not TRFs are available)
        curX=0;
        curY=0;
        curZ=0;
        coordsYes=0;
        
        trfsToGetCoords={'itrf2014', 'vtrf2014', 'dtrf2014', 'vievsTrf', 'VieTRF13'}; % get xyz from one of these
        for kTrf=1:length(trfsToGetCoords)
            if ~isempty(ns_codes(k).(trfsToGetCoords{kTrf}))
                curX=ns_codes(k).(trfsToGetCoords{kTrf}).break(size(ns_codes(k).(trfsToGetCoords{kTrf}).break,2)).x;
                curY=ns_codes(k).(trfsToGetCoords{kTrf}).break(size(ns_codes(k).(trfsToGetCoords{kTrf}).break,2)).y;
                curZ=ns_codes(k).(trfsToGetCoords{kTrf}).break(size(ns_codes(k).(trfsToGetCoords{kTrf}).break,2)).z;
                coordsYes=1;
                break
            end
        end

        if coordsYes==0
            % write all stations which have no coordinates to this array for error message
            noCoordsStats{k,1}=ns_codes(k).name;
            noCoordsStats{k,2}=ns_codes(k).domes;
            noCoordsStats{k,3}=ns_codes(k).comments;
            
        end
        noCoordsStats{k,4}=coordsYes;    
        if coordsYes==1
            fprintf(osoFormat, ns_codes(k).name, curX, curY, curZ);
        end
    end
end
   
% delete empty entries in noCoordsStats
% noCoordsStatsTest=noCoordsStats;
% noCoordsStatsTest(cell2mat(noCoordsStatsTest(:,4))==1,:)=[];

noCoordsStats(cellfun(@isempty, noCoordsStats(:,1)),:)=[];

% write stations with no coords to command window (for information)
fprintf('\nFollowing stations have no coordinates (neither ITRF2014, vievsTrf(backup),\n VieTRF13, VTRF2014, nor DTRF2014) and therefore no ocean loading\ncorrections can be calculated for them:\n\n  NAME     DOMES      COMMENT\n');
for k=1:size(noCoordsStats,1)
    fprintf('%4.0f %-8s  %9s  %30s\n', k,noCoordsStats{k,1}, noCoordsStats{k,2}, noCoordsStats{k,3});
end


% =============================================
% 5. save to disk and define output of function
% =============================================

fprintf('\n5. Saving superstations struct\n\n');


% Delete fields for unavailable TRFs:



if ~flag_trf_itrf2014
    ns_codes = rmfield(ns_codes, 'itrf2014');
end
if ~flag_trf_dtrf2014
    ns_codes = rmfield(ns_codes, 'dtrf2014');
end
if ~flag_trf_vtrf2014
    ns_codes = rmfield(ns_codes, 'vtrf2014');
end
if ~flag_trf_vieTrf13
    ns_codes = rmfield(ns_codes, 'VieTRF13');
end
if ~flag_trf_ivstrf2014b
    ns_codes = rmfield(ns_codes, 'ivsTrf2014b');
end


% save the main variable under name superstations
superstations=ns_codes;

% save to output mat file
save(outFile, 'superstations');

% make temporary message
message=sprintf('Superstations file was successfully written to %s', outFile);

% if output is wanted = error/successful message
if nargout>0
    varargout{1}=message;
else 
    % printf message to screen
    fprintf(message);
end

% --------------------------------------------------
% START OF NESTED FUNCTIONS
% --------------------------------------------------

    % This function sets all first start epochs to 0 and all end epochs to
    % 99999.
    % This is needed for eh ITRF2014 and VTRF2014 since they have start and
    % end epoch equal to first/last observation of that station (and
    % therefore no later VLBI session could be analyzed with that values)
    % INPUT
    %  trfname  string
    function setStartEndEpoch(trfname)
        for k=1:length(ns_codes)
            if ~isempty(ns_codes(k).(trfname))
                ns_codes(k).(trfname).break(1).start=0;
                ns_codes(k).(trfname).break(length(ns_codes(k).(trfname).break)).end=99999;
            end
        end
        
    end

end %  main function mk_superstatFile
