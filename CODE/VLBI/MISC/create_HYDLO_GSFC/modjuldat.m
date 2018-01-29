function tmjd = modjuldat(y,m,d,h,mi,s)

% function to convert date to modified julian date
% input parameters:
% y   years, e.g. 2007
% m   months
% d   days
% h   hours
% mi  minutes
% s   seconds
% 
% Johannes Boehm, 2001-03-27
% mod. 2009-02-03 Johannes Boehm 
% mod. 2009-09-30 Lucia Plank

if (nargin < 6)
    s = 0;
 if (nargin < 5)
     mi = 0;
  if (nargin < 4)
      h = 0;
   if (nargin < 3)
       d = 1;
    if (nargin < 2)
       m = 1;
    end
   end
  end
 end
end

if (m<=2)
   m = m+12;
   y = y-1;
end

% if date is before Oct 4, 1582
if (y<=1582)&&(m<=10)&&(d<=4)
   b = -2;
% if date is after Oct 4, 1582   
else
   b = floor(y/400)-floor(y/100);
end

jd = floor(365.25*y)-2400000.5;

tmjd=jd+floor(30.6001*(m+1))+b+1720996.5+d+h/24+mi/1440+s/86400;
