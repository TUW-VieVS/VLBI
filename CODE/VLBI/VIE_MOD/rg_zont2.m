% ************************************************************************
%   Description:
%   Long term Tidal variations in the Earth rotation. To obtain variations
%   free from tidal effects, the correction should be substracted from the
%   observed value. 
%   Translated subroutine RG_ZONT.F
% 
%   Reference: 
%   IERS Conv. 2010, chapter 8 
%
%   Input:										
%      mjd     (n,1)           time in mjd  TT [d]
%      opt35                   2 options
%                    1 ......  UT1R <35d
%                    2 ......  UT1S all periods 5d-18.6yr
%
%   Output:
%      corr    (n,1)           correction on UT1 in [s]
% 
%   External calls: 	
%      fund_arg.m  
%
%   Coded for VieVS: 
%   02 Nov 2012 by Lucia Plank
%
%   Revision: 
%
% *************************************************************************
%
% *  - - - - - - - - - - -
% *   R G _ Z O N T 2
% *  - - - - - - - - - - -
% *
% *  This routine is part of the International Earth Rotation and
% *  Reference Systems Service (IERS) Conventions software collection.
% *
% *  This subroutine evaluates the effects of zonal Earth tides on the
% *  rotation of the Earth.  The model used is a combination of Yoder
% *  et al. (1981) elastic body tide, Wahr and Bergen (1986) inelastic
% *  body tide, and Kantha et al. (1998) ocean tide models 
% *  as recommended by the IERS Conventions (2010).  Refer to
% *  Chapter 8 pp. 105 - 106.  The latest version of the model is located
% *  at http://tai.bipm.org/iers/convupdt/convupdt_c8.html.
% *
% *  In general, Class 1, 2, and 3 models represent physical effects that
% *  act on geodetic parameters while canonical models provide lower-level
% *  representations or basic computations that are used by Class 1, 2, or
% *  3 models.
% * 
% *  Status:  Class 3 model
% *
% *     Class 1 models are those recommended to be used a priori in the
% *     reduction of raw space geodetic data in order to determine
% *     geodetic parameter estimates.
% *     Class 2 models are those that eliminate an observational
% *     singularity and are purely conventional in nature.
% *     Class 3 models are those that are not required as either Class
% *     1 or 2.
% *     Canonical models are accepted as is and cannot be classified as
% *     a Class 1, 2, or 3 model.
% *
% *  Given:
% *     T           d      TT, Julian centuries since J2000 (Note 1)
% *
% *  Returned:
% *     DUT         d      Effect on UT1 (Note 2)
% *     DLOD        d      Effect on excess length of day (LOD) (Note 3)
% *     DOMEGA      d      Effect on rotational speed (Note 4)
% *
% *  Notes:
% *
% *  1) Though T is strictly TDB, it is usually more convenient to use
% *     TT, which makes no significant difference.  Julian centuries since
% *     J2000 is (JD - 2451545.0)/36525.
% *
% *  2) The expression used is as adopted in IERS Conventions (2010).
% *     DUT is expressed in seconds and is double precision.
% *
% *  3) The expression used is as adopted in IERS Conventions (2010).
% *     DLOD is the excess in LOD and is expressed in seconds per day
% *     and is double precision.  The phrase 'per day' is generally
% *     understood, so it has been omitted commonly in speech and
% *     literature.  
% *     See: Stephenson, F. R., Morrison, L. V., Whitrow, G. J., 1984,
% *     "Long-Term Changes in the Rotation of the Earth: 700 B. C. to
% *     A. D. 1980 [and Discussion]", Phil. Trans. Roy. Soc. of London.
% *     Series A, 313, pp. 47 - 70.
% * 
% *  4) The expression used is as adopted in IERS Conventions (2010).
% *     Rotational speed is expressed in radians per second and is
% *     double precision.
% *  
% *  Called:
% *     FUNDARG      Computation of the fundamental lunisolar arguments
% *
% *  Test case:
% *     given input: T = .07995893223819302 Julian centuries since J2000
% *                  (MJD = 54465)
% *     expected output: DUT    =  7.983287678576557467E-002 seconds
% *                      DLOD   =  5.035331113978199288E-005 seconds / day
% *                      DOMEGA = -4.249711616463017E-014 radians / second
% *
% *  References:
% *
% *     Yoder, C. F., Williams, J. G., and Parke, M. E., (1981),
% *     "Tidal Variations of Earth Rotation," J. Geophys. Res., 86,
% *     pp. 881 - 891.
% *
% *     Wahr, J. and Bergen, Z., (1986), "The effects of mantle 
% *     anelasticity on nutations, Earth tides, and tidal variations
% *     in rotation rate," Geophys. J. Roy. astr. Soc., 87, pp. 633 - 668.
% *
% *     Kantha, L. H., Stewart, J. S., and Desai, S. D., (1998), "Long-
% *     period lunar fortnightly and monthly ocean tides," J. Geophys.
% *     Res., 103, pp. 12639 - 12647.
% *
% *     Gross, R. S., (2009), "Ocean tidal effects on Earth rotation,"
% *     J. Geodyn., 48(3-5), pp. 219 - 225.
% * 
% *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
% *     IERS Technical Note No. 36, BKG (2010)
% *
% *  Revisions:  
% *  2008 January 18 B.E. Stetzler  Initial changes to header
% *               and used 2PI instead of PI as parameter
% *  2008 January 25 B.E. Stetzler Additional changes to header
% *  2008 February 21 B.E. Stetzler Definition of (excess) LOD clarified
% *  2008 March   12 B.E. Stetzler Applied changes to wording of notes.
% *  2008 March   14 B.E. Stetzler Further changes applied to code.
% *  2008 April   03 B.E. Stetzler Provided example test case
% *  2009 February 11 B.E. Stetzler Updated test case due to changes made
% *                                 to FUNDARG.F subroutine
% *  2009 April   10 B.E. Stetzler DLOD corrected to say it is expressed
% *                                in seconds per day
% *  2009 May     04 B.E. Stetzler Code formatting changes based on 
% *                                client recommendations
% *  2009 May     07 B.E. Stetzler Updated test case due to above changes
% *  2010 February 19 B.E. Stetzler Replaced Conventions 2003 recommended
% *                                 model with Conventions 2010 model
% *  2010 February 22 B.E. Stetzler Provided example test case
% *  2010 February 23 B.E. Stetzler Updated values to two decimal places
% *  2010 February 23 B.E. Stetzler Split fundamental arguments and
% *                                 coefficients for four decimal place
% *                                 precision
% *  2010 February 25 B.E. Stetzler Recalculation of fundamental arguments
% *  2010 March    01 B.E. Stetzler Updated table values to four decimal
% *                                 places and double precision
% *  2010 March    12 B.E. Stetzler Applied changes to wording of notes.
% *  2010 March    22 B.E. Stetzler Corrected DOMEGA output for test case
% *  2011 December 20 B.E. Stetzler Corrected typo in the sine coefficient
% *                                 of the 365.26 period for delta LOD.
% *                                 Noted by J. Lundberg. Test case output
% *                                 for DLOD updated.
% *-----------------------------------------------------------------------

function [corr]=rg_zont2(mjd,opt35)
	 	
% 	Mean angular velocity of the Earth in rad/s

%	Om = 7.292115D-5;
	
%	Conversion radians <--> arcseconds

%        a2r=4.84813681109535993D-6;
%        r2a=2.06264806247096355D+5;
	
%	Julian centuries from J2000.0

	jc = (mjd-51544.5D0)./36525.D0;
	
%	Delaunay arguments from IERS Conventions (1996)
% + usage of fund_arg.m --> small differences in the dalaunay arguments +
% + due to change of addition-sequences for numerical reasons           +

    ARG = fund_arg(jc,1)';
    
% * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% *  --------------------------------------------------
% *  Tables of multiples of arguments and coefficients
% *  --------------------------------------------------

TAB1=[...

          1,  0,  2,  2,  2,    -0.0235D0,0.0000D0
          2,  0,  2,  0,  1,    -0.0404D0,0.0000D0
          2,  0,  2,  0,  2,    -0.0987D0,0.0000D0
          0,  0,  2,  2,  1,    -0.0508D0,0.0000D0
          0,  0,  2,  2,  2,    -0.1231D0,0.0000D0
          1,  0,  2,  0,  0,    -0.0385D0,0.0000D0
          1,  0,  2,  0,  1,    -0.4108D0,0.0000D0
          1,  0,  2,  0,  2,    -0.9926D0,0.0000D0
          3,  0,  0,  0,  0,    -0.0179D0,0.0000D0
         -1,  0,  2,  2,  1,    -0.0818D0,0.0000D0
         -1,  0,  2,  2,  2,    -0.1974D0,0.0000D0
          1,  0,  0,  2,  0,    -0.0761D0,0.0000D0
          2,  0,  2, -2,  2,     0.0216D0,0.0000D0
          0,  1,  2,  0,  2,     0.0254D0,0.0000D0
          0,  0,  2,  0,  0,    -0.2989D0,0.0000D0
          0,  0,  2,  0,  1,    -3.1873D0,0.2010D0
          0,  0,  2,  0,  2,    -7.8468D0,0.5320D0
          2,  0,  0,  0, -1,     0.0216D0,0.0000D0
          2,  0,  0,  0,  0,    -0.3384D0,0.0000D0
          2,  0,  0,  0,  1      0.0179D0,0.0000D0
          0, -1,  2,  0,  2,    -0.0244D0,0.0000D0
          0,  0,  0,  2, -1,     0.0470D0,0.0000D0
          0,  0,  0,  2,  0,    -0.7341D0,0.0000D0
          0,  0,  0,  2,  1,    -0.0526D0,0.0000D0
          0, -1,  0,  2,  0,    -0.0508D0,0.0000D0
          1,  0,  2, -2,  1,     0.0498D0,0.0000D0
          1,  0,  2, -2,  2,     0.1006D0,0.0000D0
          1,  1,  0,  0,  0,     0.0395D0,0.0000D0
         -1,  0,  2,  0,  0,     0.0470D0,0.0000D0
         -1,  0,  2,  0,  1,     0.1767D0,0.0000D0
         -1,  0,  2,  0,  2,     0.4352D0,0.0000D0
          1,  0,  0,  0, -1,     0.5339D0,0.0000D0
          1,  0,  0,  0,  0,    -8.4046D0,0.2500D0
          1,  0,  0,  0,  1,     0.5443D0,0.0000D0
          0,  0,  0,  1,  0,     0.0470D0,0.0000D0
          1, -1,  0,  0,  0,    -0.0555D0,0.0000D0
         -1,  0,  0,  2, -1,     0.1175D0,0.0000D0
         -1,  0,  0,  2,  0,    -1.8236D0,0.0000D0
         -1,  0,  0,  2,  1,     0.1316D0,0.0000D0
          1,  0, -2,  2, -1      0.0179D0,0.0000D0
         -1, -1,  0,  2,  0,    -0.0855D0,0.0000D0
         % stop for periods < 35d
          0,  2,  2, -2,  2,    -0.0573D0,0.0000D0
          0,  1,  2, -2,  1,     0.0329D0,0.0000D0
          0,  1,  2, -2,  2,    -1.8847D0,0.0000D0
          0,  0,  2, -2,  0,     0.2510D0,0.0000D0
          0,  0,  2, -2,  1,     1.1703D0,0.0000D0
          0,  0,  2, -2,  2,   -49.7174D0,0.4330D0
          0,  2,  0,  0,  0,    -0.1936D0,0.0000D0
          2,  0,  0, -2, -1,     0.0489D0,0.0000D0
          2,  0,  0, -2,  0,    -0.5471D0,0.0000D0
          2,  0,  0, -2,  1,     0.0367D0,0.0000D0
          0, -1,  2, -2,  1,    -0.0451D0,0.0000D0
          0,  1,  0,  0, -1,     0.0921D0,0.0000D0
          0, -1,  2, -2,  2,     0.8281D0,0.0000D0
          0,  1,  0,  0,  0,   -15.8887D0,0.1530D0
          0,  1,  0,  0,  1,    -0.1382D0,0.0000D0
          1,  0,  0, -1,  0,     0.0348D0,0.0000D0
          2,  0, -2,  0,  0,    -0.1372D0,0.0000D0
         -2,  0,  2,  0,  1,     0.4211D0,0.0000D0
         -1,  1,  0,  1,  0     -0.0404D0,0.0000D0
          0,  0,  0,  0,  2,     7.8998D0,0.0000D0
          0,  0,  0,  0,  1  -1617.2681D0,0.0000D0 ];
	

%	Implementation for time jc in second
	
%	select case (choice)
if opt35==1
    TAB1(42:end,:)=[];
    disp('UT1R')
else
    disp('UT1S') 
end
%	Correction to remove from UT1 in s

	phi=TAB1(:,1:5)*ARG;
    corr=(TAB1(:,6))'*sin(phi)+(TAB1(:,7))'*cos(phi);
	corr=corr'*1.D-4; %[s]

% *  Finished.
% 
% *+----------------------------------------------------------------------
% *
% *  Copyright (C) 2008
% *  IERS Conventions Center
% *
% *  ==================================
% *  IERS Conventions Software License
% *  ==================================
% *
% *  NOTICE TO USER:
% *
% *  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
% *  WHICH APPLY TO ITS USE.
% *
% *  1. The Software is provided by the IERS Conventions Center ("the
% *     Center").
% *
% *  2. Permission is granted to anyone to use the Software for any
% *     purpose, including commercial applications, free of charge,
% *     subject to the conditions and restrictions listed below.
% *
% *  3. You (the user) may adapt the Software and its algorithms for your
% *     own purposes and you may distribute the resulting "derived work"
% *     to others, provided that the derived work complies with the
% *     following requirements:
% *
% *     a) Your work shall be clearly identified so that it cannot be
% *        mistaken for IERS Conventions software and that it has been
% *        neither distributed by nor endorsed by the Center.
% *
% *     b) Your work (including source code) must contain descriptions of
% *        how the derived work is based upon and/or differs from the
% *        original Software.
% *
% *     c) The name(s) of all modified routine(s) that you distribute
% *        shall be changed.
% * 
% *     d) The origin of the IERS Conventions components of your derived
% *        work must not be misrepresented; you must not claim that you
% *        wrote the original Software.
% *
% *     e) The source code must be included for all routine(s) that you
% *        distribute.  This notice must be reproduced intact in any
% *        source distribution. 
% *
% *  4. In any published work produced by the user and which includes
% *     results achieved by using the Software, you shall acknowledge
% *     that the Software was used in obtaining those results.
% *
% *  5. The Software is provided to the user "as is" and the Center makes
% *     no warranty as to its use or performance.   The Center does not
% *     and cannot warrant the performance or results which the user may
% *     obtain by using the Software.  The Center makes no warranties,
% *     express or implied, as to non-infringement of third party rights,
% *     merchantability, or fitness for any particular purpose.  In no
% *     event will the Center be liable to the user for any consequential,
% *     incidental, or special damages, including any lost profits or lost
% *     savings, even if a Center representative has been advised of such
% *     damages, or for any claim by any third party.
% *
% *  Correspondence concerning IERS Conventions software should be
% *  addressed as follows:
% *
% *                     Gerard Petit
% *     Internet email: gpetit[at]bipm.org
% *     Postal address: IERS Conventions Center
% *                     Time, frequency and gravimetry section, BIPM
% *                     Pavillon de Breteuil
% *                     92312 Sevres  FRANCE
% *
% *     or
% *
% *                     Brian Luzum
% *     Internet email: brian.luzum[at]usno.navy.mil
% *     Postal address: IERS Conventions Center
% *                     Earth Orientation Department
% *                     3450 Massachusetts Ave, NW
% *                     Washington, DC 20392
% *
% *
% *-----------------------------------------------------------------------
