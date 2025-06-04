% ------------------------------------------------------------------------
% [coeff,t0,hstep,nint,nq,nsat,tbound,tosc,oscele,satnum,narc,descr,source] = readstd (filename);
%
% Purpose: Read Formatted STD File
%
% Input:   filename       - GNSS precise orbit file name
% Output:  coeff(i,j,k,m) - Polynomial coefficients
%                           i=1,...,nint - integration interval
%                           j=1,...,nq+1 - degree
%                           k=1,2,3,     - component
%                           m=1,...,nsat - satellite nr
%          t0(i)          - epoch expansion point (mjd)
%                           i=1,...,nint - integration interval
%          hstep(i)       - integraiton interval (sec)
%                           i=1,...,nint - integration interval
%          nint           - number of integration intervals
%          nq             - polynomial degree
%          nsat           - number of satellites
%          tbound(i)      - integration interval boundaries
%                           i=1,...,nint+1
%          tosc           - osculation epoch (mjd)
%          oscele(m,n)    - esculating elements
%                           m=1,...,nsat - satellite nr
%                           n=1: semimajor axis (m)
%                             2: eccentricity
%                             3: inclination (rad)
%                             4: RA of ascending node (rad)
%                             5: argument of perigee (rad)
%                             6: argument of latitude (rad)
%                             7: perigee passage time (sec after tosc)
%          satnum(m)      - satellite number
%                           m=1,...,nsat - satellite nr
%                            nn: GPS
%                           1nn: GLONASS
%                           2nn: Galileo
%                           3nn: SBAS
%                           4nn: BeiDou
%                           5nn: QZSS
%          narc           - number of acrs (only 1 arc currently supported)
%          descr(i,j)     - model description
%          source(i)      - source of elements
%
% Author:  Urs Hugentobler, FESG
% Date  :  27-11-2022
% Changes: 
% -------------------------------------------------------------------------
function [coeff,t0,hstep,nint,nq,nsat,tbound,tosc,oscele,satnum,narc,descr,source] = readstd (filename)

	format compact;

	fid=fopen(filename,'r');

	% number of arcs
	tline = fgetl(fid);
	narc  = str2num(tline);
	% read arcs
	if (narc > 1); error('narc > 1 not implemented'); end;

	% read new format
	if (narc < 0) 
		tline = fgetl(fid);
		iftm  = str2num(tline(1:4));
		narc  = str2num(tline(5:8));
	% read model descriptor
		tline = fgetl(fid);
		nlin  = str2num(tline);
		descr = [];
		for i=1:nlin
			tline = fgetl(fid);
			descr = [descr;tline];
		end
	end

	for iarc = 1:narc
		tline = fgetl(fid);
		nsat  = str2num(tline(1:7));
		nint  = str2num(tline(8:14));
		nq    = str2num(tline(15:21));
		% satnr
		nlin  = ceil(nsat/24);
		num   = nsat;
		satnum = [];
		for i=1:nlin
			tline = fgetl(fid);
			for j = 1:min(num,24);
				sat = str2num(tline((j-1)*3+1:j*3));
				satnum = [satnum,sat];
			end
			num = num-24;
		end
		source = fgetl(fid);

		% epochs of intervals
		tline = fgetl(fid);
		tosc  = str2num(tline(1:26));
		tbound = zeros(1,nint+1);
		for i = 1:nint+1
			tline = fgetl(fid);
			tbound(i) = str2num(tline);
		end
		
		% osculating elements
		oscele = zeros(nsat,7);
		for i=1:nsat
			for j=1:7
				tline = fgetl(fid);
				oscele(i,j) = str2num(tline);
			end
		end
		
		% polynomial coefficients
		t0    = zeros(1,nint);
		hstep = zeros(1,nint);
		coeff = zeros(nint,nq+1,3,nsat);
		for i = 1:nint
			tline = fgetl(fid);
			t0(i)    = str2num(tline(1:26));
			hstep(i) = str2num(tline(27:52));
			for j=1:nq+1
				for k=1:nsat
					tline = fgetl(fid);
					coeff(i,j,1,k) = str2num(tline(1:26));
					coeff(i,j,2,k) = str2num(tline(27:52));
					coeff(i,j,3,k) = str2num(tline(53:78));
				end
			end
		end
	end

	fclose(fid);
end