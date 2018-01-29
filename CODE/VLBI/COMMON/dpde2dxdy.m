function [dX,dY]=dpde2dxdy(mjd,dpsi,deps)
%-------------------------------------------------------------------------C
%									  C
%       Conversion of offsets dpsi, deps to dX, dY                        C
%									  C
%-------------------------------------------------------------------------C
%	 								  C
%	References :                                                      C
%		Lieske et al. 1977, A&A 58, pp. 1-16                      C
%               IERS Conventions 2000, chap. 5                            C
%	 	Herring et al. 2002, JGR 107, B4			  C
%	 								  C
%-------------------------------------------------------------------------C
%	 								  C
%	Input :                                                           C
%		mjd : time in mjd                                         C
%		dpsi, deps : offsets in arcsec ref. to IAU2000A           C
%	Output :                                                          C
%               dX, dY : offsets in arcsec ref. to IAU2000A	          C
%	 								  C
%-------------------------------------------------------------------------C
%	 								  C
%	Subroutine interface : none					  C
%	 								  C
%-------------------------------------------------------------------------C
%	 								  C
%	Author : Sebastien Lambert (Sebastien.Lambert@obspm.fr)		  C
%	Last modified : 25/11/02	 				  C
%	 								  C
%-------------------------------------------------------------------------C

	
a2r=4.84813681109535993e-6;
r2a=2.06264806247096355e+5;
	
%	Julian centuries from J2000.0

	jc=(mjd-51544.5)/36525;

%	Lieske et al. (1977) expressions in arcsec and converted in radians
%	and improved numerical value of psiA from Herring et al. (2002)

	eps0=84381.448;

        psiA=5038.47875*jc-1.07259*jc.^2-0.001147*jc.^3;

        chiA=10.5526*jc-2.38064*jc.^2-0.001125*jc.^3;

        epsA=eps0-46.8402*jc-0.00059*jc.^2+0.001813*jc.^3;

        eps0=eps0*a2r;
        psiA=psiA*a2r;
        chiA=chiA*a2r;
        epsA=epsA*a2r;
	
%	Chapter 5, equation (23)

	dX=dpsi*sin(epsA)+deps*(psiA*cos(epsA)-chiA);
	dY=deps-dpsi*sin(epsA)*(psiA*cos(epsA)-chiA);
     

