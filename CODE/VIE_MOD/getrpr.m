% ------------------------------------------------------------------------
% [drdpar,dvdpar] = getrpr(prn,t,coerpr,trpr,hrpr,tbnd,sclpar,satnum);
%
% Purpose: Read Formatted STD File
%
% Input:   filname           - GNSS precise orbit file name
%          prn               - satellite number
%          t                 - epoch of pos, vel (mjd)
%          coerpr(i,j,k,m,n) - Polynomial coefficients
%                              i=1,...,nint - integration interval
%                              j=1,...,nsat - satellite nr
%                              k=1,...,nvar - parameter
%                              m=1,...,nq+1 - degree
%                              n=1,2,3,     - component
%          trpr(i)           - epoch expansion point (mjd)
%                              i=1,...,nrpr - integration interval
%          hrpr(i)           - integraiton interval (sec)
%                              i=1,...,nrpr - integration interval
%          tbnd(i)           - integration interval boundaries
%                              i=1,...,nint+1
%          sclpar(i,j)       - parameter scaling
%                              i=1,...,nsat - satellite nr
%                              j=1,...,nvar - parameters
%          satnum(m)      - satellite number
%                           m=1,...,nsat - satellite nr
%                            nn: GPS
%                           1nn: GLONASS
%                           2nn: Galileo
%                           3nn: SBAS
%                           4nn: BeiDou
%                           5nn: QZSS
%
% Author:  Urs Hugentobler, FESG
% Date  :  27-11-2022
% Changes: 
% -------------------------------------------------------------------------

function [drdpar,dvdpar] = getrpr(prn,t,coerpr,trpr,hrpr,tbnd,sclpar,satnum)

	format compact;

	nint = size(coerpr,1);
	nvar = size(coerpr,3);
	nq   = size(coerpr,4)-1; 

	% find satellite
	isat = find(satnum==prn);
	if (isempty(isat))
		error(sprintf('*** Satellite %d not found',prn))
	end

	% find integration interval
	eps = 2/1440;
	if (t < tbnd(1)-eps)
		error(sprintf('*** Epoch earlier than first epoch %f',tbd(1)))
	end
	if (t > tbnd(end)+eps)
		error(sprintf('*** Epoch later than last epoch %f',tbd(end)))
	end
	ix   = find(tbnd-eps < t);
	int  = min(ix(end),nint);

	% polynomial argument
	ti = (t-trpr(int))*86400/hrpr(int);

	% position
	drdpar = squeeze(coerpr(int,isat,:,end,:));
	for iq=nq:-1:1
		drdpar = drdpar*ti+ squeeze(coerpr(int,isat,:,iq,:));
	end

	% velocity
	dvdpar = squeeze(coerpr(int,isat,:,end,:))*nq;
	for iq=nq:-1:2
		dvdpar = dvdpar*ti + squeeze(coerpr(int,isat,:,iq,:));
	end
	dvdpar = dvdpar/hrpr(int);

	% scaling
	for k = 1:3
		drdpar(:,k) = drdpar(:,k).*sclpar(isat,:)';
		dvdpar(:,k) = dvdpar(:,k).*sclpar(isat,:)';
	end
end