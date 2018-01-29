%       SUBROUTINE FCNNUT ( MJD,X,Y,DX,DY )
% *+
% *  - - - - - - - - - - -
% *   F C N N U T
% *  - - - - - - - - - - -
% *
% *  This routine is part of the International Earth Rotation and
% *  Reference Systems Service (IERS) Conventions software collection.
% *
% *  This subroutine computes the effects of the free core nutation.  
% *  Please note that the table is updated each year (see Note 4).
% *  The parameter N needs to be incremented for each additional
% *  year in the table.
% *
% *  In general, Class 1, 2, and 3 models represent physical effects that
% *  act on geodetic parameters while canonical models provide lower-level
% *  representations or basic computations that are used by Class 1, 2, or
% *  3 models.
% * 
% *  Status: Class 3 model
% *
% *     Class 1 models are those recommended to be used a priori in the
% *     reduction of raw space geodetic data in order to determine
% *     geodetic parameter estimates.
% *     Class 2 models are those that eliminate an observational
% *     singularity and are purely conventional in nature.
% *     Class 3 models are those that are not required as either Class
% *     1 or 2.
% *     Canonical models are accepted as is and cannot be classified as a
% *     Class 1, 2, or 3 model.
% *
% *  Given:
% *     mjd           d      Modified Julian Date, TDB (Note 1)
% *
% *  Returned:
% *     X             d      CIP offset x component, in microas (Note 2)
% *     Y             d      CIP offset y component, in microas (Note 2)
% *     dX            d      Uncertainty of x component, in microas (Note 3)
% *     dY            d      Uncertainty of y component, in microas (Note 3) 
% *
% *  Notes:
% *
% *  1) Though the Modified Julian Date (MJD) is strictly TDB, it is
% *     usually more convenient to use TT, which makes no significant
% *     difference.
% *
% *  1) CIP is the Celestial Intermediate Pole.  The expression
% *     used is given in microarcseconds.
% *  
% *  2) The expression used is given in microarcseconds.
% *
% *  3) The updated table is maintained at the website
% *     http://syrte.obspm.fr/~lambert/fcn/.
% *
% *  Test case:
% *     given input: MJD = 54790D0   Modified Julian Date, TDB
% *                  
% *     expected output:  X = -176.0443531306245006D0 microarcseconds
% *                       Y = -93.62265658430092685D0 microarcseconds
% *                       dX = 3.723278688524589874D0 microarcseconds
% *                       dY = 3.723278688524589874D0 microarcseconds
% *
% *  References:
% *
% *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
% *     IERS Technical Note No. 36, BKG (2010)
% *
% *  Revisions:
% *  2007 August   27 S. Lambert    Original code
% *  2008 November 19 B.E.Stetzler  Added header and copyright
% *  2008 November 20 B.E.Stetzler  provided a test case 
% *  2009 July     10 B.E.Stetzler  Updated parameter N to 26, 
% *                                 updated table to 2009, and corrected
% *                                 test case results
% *  2009 August   18 B.E.Stetzler  Capitalized all variables for FORTRAN
% *                                 77 compatibility
% *  2010 July      2 B.E.Stetzler  Updated parameter N to 27, 
% *                                 updated table to 2010, and corrected
% *                                 test case results
% *  2011 July     13 B.E.Stetzler  Updated parameter N to 28, 
% *                                 updated table to 2011, and corrected
% *                                 test case results
% *-----------------------------------------------------------------------

% Coded for VieVS: Hana Spicakova 04.11.2011

function [Qfcn, dQfcn, dQfcnAc, dQfcnAs] = fcnnut(tt)

% FCN parameters
sd=1.002737909; % sid./sol. days

fr=1.0023181; % cpsd FCN
PER=1/(1-fr)/sd;              %! period in solar days
% PER = -430.21;  % orig                      %! period in solar days
%%
sig=(2*pi/PER); %rad/day
t=(tt - 51544.5); %day

                  
% *       Block data of amplitudes for X (microas)
        
% TABLE =...
% [ 45700.D0,     4.55D0,   -36.58D0,    19.72D0,% ! 1984.0
% 46066.D0,  -141.82D0,  -105.35D0,    11.12D0, %! 1985.0
% 46431.D0,  -246.56D0,  -170.21D0,     9.47D0, %! 1986.0
% 46796.D0,  -281.89D0,  -159.24D0,     8.65D0, %! 1987.0
% 47161.D0,  -255.05D0,   -43.58D0,     8.11D0, %! 1988.0
% 47527.D0,  -210.46D0,   -88.56D0,     7.31D0, %! 1989.0
% 47892.D0,  -187.79D0,   -57.35D0,     6.41D0, %! 1990.0
% 48257.D0,  -163.01D0,    26.26D0,     5.52D0, %! 1991.0
% 48622.D0,  -145.55D0,    44.64D0,     4.80D0, %! 1992.0
% 48988.D0,  -146.73D0,    51.48D0,     6.02D0, %! 1993.0
% 49353.D0,  -113.67D0,    13.06D0,     9.99D0, %! 1994.0
% 49718.D0,   -87.14D0,     4.38D0,     8.68D0, %! 1995.0
% 50083.D0,   -88.62D0,     3.15D0,     8.07D0, %! 1996.0
% 50449.D0,   -95.52D0,    33.20D0,     4.44D0, %! 1997.0
% 50814.D0,   -69.26D0,    26.92D0,     3.41D0, %! 1998.0
% 51179.D0,   -43.92D0,   -14.58D0,     3.45D0, %! 1999.0
% 51544.D0,     6.19D0,   -81.41D0,     3.22D0, %! 2000.0
% 51910.D0,    69.27D0,  -133.07D0,     2.80D0, %! 2001.0
% 52275.D0,    86.87D0,  -128.07D0,     2.69D0, %! 2002.0
% 52640.D0,   111.42D0,   -43.33D0,     2.54D0, %! 2003.0
% 53005.D0,   114.60D0,     0.21D0,     2.49D0, %! 2004.0
% 53371.D0,   131.44D0,    -4.14D0,     2.68D0, %! 2005.0
% 53736.D0,   155.41D0,    28.97D0,     2.16D0, %! 2006.0
% 54101.D0,   158.77D0,    59.07D0,     1.84D0, %! 2007.0
% 54466.D0,   155.49D0,   100.93D0,     1.72D0, %! 2008.0
% 54832.D0,   142.30D0,   142.93D0,     1.88D0, %! 2009.0
% 55197.D0,    35.98D0,   184.00D0,     1.96D0, %! 2010.0
% 55562.D0,    23.45D0,   221.84D0,     2.08D0];% ! 2011.0
  

% from http://syrte.obspm.fr/~lambert/fcn/table.txt , 2016-06-11
TABLE =...
[45700.D0,  -219.93D0,  -155.91D0,     8.32D0,% ! 1984.0 
46066.D0,  -245.97D0,  -149.26D0,     6.76D0, %! 1985.0 
46431.D0,  -247.42D0,   -97.46D0,     5.75D0, %! 1986.0 
46796.D0,  -232.93D0,  -118.49D0,     4.95D0, %! 1987.0 
47161.D0,  -218.13D0,   -76.40D0,     4.27D0, %! 1988.0 
47527.D0,  -200.95D0,   -43.62D0,     3.71D0, %! 1989.0 
47892.D0,  -184.49D0,   -14.01D0,     3.32D0, %! 1990.0 
48257.D0,  -171.29D0,    -0.10D0,     3.34D0, %! 1991.0 
48622.D0,  -158.71D0,    13.29D0,     3.37D0, %! 1992.0 
48988.D0,  -142.91D0,     7.30D0,     3.36D0, %! 1993.0 
49353.D0,  -135.67D0,    27.63D0,     3.39D0, %! 1994.0 
49718.D0,  -113.62D0,    32.95D0,     2.94D0, %! 1995.0 
50083.D0,   -89.39D0,    29.53D0,     2.68D0, %! 1996.0 
50449.D0,   -67.61D0,     4.93D0,     2.53D0, %! 1997.0 
50814.D0,   -39.00D0,   -21.69D0,     2.21D0, %! 1998.0 
51179.D0,    -2.66D0,   -56.95D0,     1.93D0, %! 1999.0 
51544.D0,    16.85D0,   -66.32D0,     1.76D0, %! 2000.0 
51910.D0,    44.47D0,   -54.65D0,     1.58D0, %! 2001.0 
52275.D0,    64.48D0,   -56.76D0,     1.52D0, %! 2002.0 
52640.D0,    88.76D0,   -57.87D0,     1.49D0, %! 2003.0 
53005.D0,   113.38D0,   -30.49D0,     1.34D0, %! 2004.0 
53371.D0,   134.47D0,    -0.56D0,     1.24D0, %! 2005.0 
53736.D0,   145.18D0,    49.72D0,     1.11D0, %! 2006.0 
54101.D0,   144.64D0,    73.19D0,     1.10D0, %! 2007.0 
54466.D0,   120.61D0,   106.58D0,     1.04D0, %! 2008.0 
54832.D0,    72.15D0,   182.43D0,     0.84D0, %! 2009.0 
55197.D0,    76.37D0,   208.46D0,     0.78D0, %! 2010.0 
55562.D0,    79.31D0,   217.89D0,     0.78D0, %! 2011.0 
55927.D0,    80.66D0,   220.96D0,     0.76D0, %! 2012.0 
56293.D0,    63.49D0,   227.45D0,     0.80D0, %! 2013.0 
56658.D0,    61.63D0,   232.20D0,     0.82D0, %! 2014.0 
57023.D0,    69.78D0,   235.49D0,     0.88D0];% ! 2015.0 


% *       Amplitudes extracted from the table
N=size(TABLE,1);
DATE=TABLE(:,1);
XC=TABLE(:,2)./1000; %mas
XS=TABLE(:,3)./1000; %mas

%*       Prediction of the amplitude at the input date

if (tt<=DATE(1)) 
   AXC=XC(1);
   AXS=XS(1);
elseif (tt>=DATE(N))
   AXC=XC(N);
   AXS=XS(N);
else
   for I=1:N-1
      if (tt>=DATE(I) && tt<DATE(I+1))
         T=tt-DATE(I);
         DT=DATE(I+1)-DATE(I);
         DAXC=XC(I+1)-XC(I);
         DAXS=XS(I+1)-XS(I);
         AXC=XC(I)+(DAXC/DT)*T;
         AXS=XS(I)+(DAXS/DT)*T;
      end
   end
end
        
%*     Computation of X and Y

% AXC=0; %mas
% AXS=0;

AXC=AXC/1000/3600/180*pi; %rad
AXS=AXS/1000/3600/180*pi; %rad

% Xfcn=AXC*cos(sig*t)-AXS*sin(sig*t); %rad
% Yfcn=AXS*cos(sig*t)+AXC*sin(sig*t); %rad

% Qfcn= [1  0      Xfcn
%        0  1      Yfcn
%        -Xfcn -Yfcn 1];

% don't apply FCN a priori
Qfcn= [1  0 0
       0  1 0
       0  0 1];


%% partials

xix=-AXC*sin(sig*t)-AXS*cos(sig*t);
xiy=-AXS*sin(sig*t)+AXC*cos(sig*t);
dPHI=-2*pi*sd*t;

% period (fix amplitudes to Lambert's model)
dQfcn = [ 0                0   dPHI*xix
         0                0   dPHI*xiy
         -dPHI*xix  -dPHI*xiy 0]; %[day]

% amplitude
dQfcnAc = [ 0            0    cos(sig*t)
            0            0    sin(sig*t)
            -cos(sig*t) -sin(sig*t) 0 ]; %[-]

        
dQfcnAs = [ 0            0    -sin(sig*t)
            0            0    cos(sig*t)
            sin(sig*t) -cos(sig*t) 0 ]; %[-]


