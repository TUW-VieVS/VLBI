% ------------------------------------------------------------------------
% [coerpr,trpr,hrpr,nvar,nrpr,nqrpr,tbnd,nsat,orbpar,sclpar,locq,satnum,narc,anltyp,source] = readrpr (filename);
%
% Purpose: Read Formatted STD File
%
% Input:   filename          - GNSS precise orbit file name
% Output:  coerpr(i,j,k,m,n) - Polynomial coefficients
%                              i=1,...,nint - integration interval
%                              j=1,...,nsat - satellite nr
%                              k=1,...,nvar - parameter
%                              m=1,...,nq+1 - degree
%                              n=1,2,3,     - component
%          trpr(i)           - epoch expansion point (mjd)
%                              i=1,...,nrpr - integration interval
%          hrpr(i)           - integraiton interval (sec)
%                              i=1,...,nrpr - integration interval
%          nvar              - number of parameters (var eqns)
%          nrpr              - number of integration intervals
%          nqrpr             - polynomial degree
%          tbnd(i)           - integration interval boundaries
%                              i=1,...,nint+1
%          nsat              - number oft satellites
%          orbpar(i,j)       - orbital parameters
%                              i=1,...,nsat - satellite nr
%                              j=1: semimajor axis (m)
%                                2: eccentricity
%                                3: inclination (rad)
%                                4: RA of ascending node (rad)
%                                5: argument of perigee (rad)
%                                6: argument of latitude (rad)
%                                7: direct rpr D0
%                                8: direct rpr Y0
%                                9: direct rpr B0
%                               10: direct rpr DC
%                               11: direct rpr YC
%                               12: direct rpr BC
%                               13: direct rpr DS
%                               14: direct rpr YS
%                               15: direct rpr BS
%          sclpar(i,j)       - parameter scaling
%                              i=1,...,nsat - satellite nr
%                              j=1,...,nvar - parameters
%          locq(i,j,k)       - parameter description                   
%                              i=1,...,nsat - satellite nr
%                              j=1,...,mvar - parameters
%                              k=1,...,6 - description id
%          satnum(m)      - satellite number
%                           m=1,...,nsat - satellite nr
%                            nn: GPS
%                           1nn: GLONASS
%                           2nn: Galileo
%                           3nn: SBAS
%                           4nn: BeiDou
%                           5nn: QZSS
%          narc           - number of acrs (only 1 arc currently supported)
%          anltyp         - model description
%          source         - source of elements
%
% Author:  Urs Hugentobler, FESG
% Date  :  27-11-2022
% Changes: 
% -------------------------------------------------------------------------

function [coerpr,trpr,hrpr,nvar,nrpr,nqrpr,tbnd,nsat,orbpar,sclpar,locq,satnum,narc,anltyp,source] = readrpr (filename)

	format compact;

	fid=fopen(filename,'r');

	% number of arcs
	tline = fgetl(fid);
	if (strcmp(tline,'#P'))
		error('*** Format with stoch pulses not supported') 
	end
	narc  = str2double(tline);
	% read arcs
	if (narc > 1)
		error('*** narc > 1 not implemented') 
	end

	for iarc = 1:narc
		tline = fgetl(fid);
		nsat  = str2double(tline(1:7));
		nrpr  = str2double(tline(8:14));
		nqrpr = str2double(tline(15:21));
		% satnr
		nlin  = ceil(nsat/24);
		num   = nsat;
		satnum = [];
		for i=1:nlin
			tline = fgetl(fid);
			for j = 1:min(num,24)
				sat = str2double(tline((j-1)*3+1:j*3));
				satnum = [satnum,sat];
			end
			num = num-24;
		end
		source = fgetl(fid);

		% boundary epochs
		tline  = fgetl(fid);
		tarpr  = str2double(tline(1:22));
		tbrpr  = str2double(tline(23:45));
		zeron  = str2double(tline(46:68));
		if (zeron ~= 12) 
			error('*** Old format not supported')
		end

		% number of equations
		tline  = fgetl(fid);
		nvar   = str2double(tline(1:22));
		anltyp = tline(24:31);
		
		% integration intervals
		tbnd   = zeros(1,nrpr+1);
		for i=1:nrpr+1
			tline   = fgetl(fid);
			tbnd(i) = str2double(tline(1:22));
		end
		
		% parameters
		orbpar = zeros(nsat,nvar);
		sclpar = zeros(nsat,nvar);
		locq   = zeros(nsat,nvar,6);
		for i = 1:nsat
			for j = 1:nvar
				tline  = fgetl(fid);
				orbpar(i,j) = str2double(tline(1:22));
				sclpar(i,j) = str2double(tline(24:45));
				for k = 1:6
					locq(i,j,k) = str2double(tline(46+(k-1)*6:50+(k-1)*6));
				end
			end
		end
		
		% polynomial coefficients
		trpr   = zeros(1,nrpr);
		hrpr   = zeros(1,nrpr);
		coerpr = zeros(nrpr,nqrpr+1,3,nsat);
		for i = 1:nrpr
			tline = fgetl(fid);
			trpr(i) = str2double(tline(1:22));
			hrpr(i) = str2double(tline(24:45));
			for j=1:nsat
				for k=1:nvar
					for m=1:nqrpr+1
						tline = fgetl(fid);
						coerpr(i,j,k,m,1) = str2double(tline(1:22));
						coerpr(i,j,k,m,2) = str2double(tline(24:45));
						coerpr(i,j,k,m,3) = str2double(tline(47:68));
					end
				end
			end
		end
	end
	fclose(fid);
end