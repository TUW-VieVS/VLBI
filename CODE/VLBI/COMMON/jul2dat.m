% C *******************************************************************
% C     SUBRUTINA JULDAT
% C *******************************************************************
% C
% C USAGE : call juldat (iy,id,ih,im,s,idj,utj,ipaso)
% C
% C PROGRAMMER: N. ZARRAOA
% C
% C ARGUMENTS: IY    - Year
% C            ID    - Day of year
% C            IH    - Hour
% C            IM    - Minutes
% C            S     - Seconds
% C            IDJ   - Modified Julian Day (integer)
% C            UTJ   - Fraction of Julian Day (rad)
% C            IPASO - If ipaso.eq.0  then   ydhms --------> Julian date
% C                    if ipaso.ne.0  then  julian date ---> ydhms
% C
% C REQUIRED ROUTINES : None
% C
% C LASTEST REVISION : September 20, 1989
% C LAST REVISION: 12, January, 2001 (Oleg Titov)  Problem 2001 has been fixed
% C***************************************************************************
%%
% 2009/01 Spicakova Hana

% translated OCCAM/STATION/arc.f/SUBROUTINE juldat

% ipaso=1
% PASS FROM JULIAN DATE TO YEAR + DOY

% external functions:  ----
   
 % input:   mjd        Modified Julian Day

 % output:  yr         year
 %          doy        day of year
 %          h          hour
 %          m          minutes
 %          s          seconds

%%

function [yr,doy,h,m,s] = jul2dat(mjd)

j20=51544;  % mjd for year=2000

mjd_int = fix(mjd);

ic0=mjd_int-j20;
ic=ic0/1461;
ic=fix(ic);

if (ic0<0)
    ic=ic-1;
end

ir=ic0-ic*1461;
yr=2000+ic*4;
yr=fix(yr);

if (ir<366)
    doy = ir+1;
end

if (ir >= 366)
    ir=ir-366;
    iy1=ir/365;
    iy1=fix(iy1);  
    yr=yr+1+iy1;
    doy=ir-iy1*365+1;
end

frd = (mjd-mjd_int);

horas=(frd)*24;
h=fix(horas+0.0001);
pmin=(horas-h)*60;
m=fix(pmin+0.0001);
s=(pmin-m)*60;

