function tmjd = modjuldat(y,m,d,h,mi,s)

% function to convert date to modified julian date
% input parameters:
% y   years, e.g. 2007 (n x 1 vector)
% m   months (n x 1 vector)
% d   days (n x 1 vector)
% h   hours (n x 1 vector)
% mi  minutes (n x 1 vector)
% s   seconds (n x 1 vector)
%
% Created:
% Johannes Boehm, 2001-03-27
%
% Changes:
% - mod. 2009-02-03 Johannes Boehm
% - mod. 2009-09-30 Lucia Plank
% - 2016-11-28, A. Hellerschmied: Now vecotr input is supported. All input arguments have to be scalars or (n x 1) vectors
% - 2017-03-24, A. Hellerschmied: Ther was a bug in case month and day is not specified as input argument (were set to "0"!)
%                                  - Now months and days are to to "1" (January 1st) if not specified

n_mjd = size(y, 1);

% Check input:
% if size(m,1)~=n_mjd || size(d,1)~=n_mjd || size(h,1)~=n_mjd || size(mi,1)~=n_mjd || size(s,1)~=n_mjd
%     error('Invalid input: Input arguments are not equally sized.');
% end

if (nargin < 6)
    s = zeros(n_mjd, 1);
    if (nargin < 5)
        mi = zeros(n_mjd, 1);
        if (nargin < 4)
            h = zeros(n_mjd, 1);
            if (nargin < 3)
                d = ones(n_mjd, 1);
                if (nargin < 2)
                    m = ones(n_mjd, 1);
                end
            end
        end
    end
end

ind = m <= 2;
m(ind) = m(ind) + 12;
y(ind) = y(ind) - 1;

% if date is before Oct 4, 1582
b = zeros(n_mjd, 1);
ind = (y <= 1582) & (m <= 10) & (d <= 4);
b(ind)  = -2;
% if date is after Oct 4, 1582
b(~ind) = floor(y(~ind)/400)-floor(y(~ind)/100);

jd      = floor(365.25.*y)-2400000.5;
tmjd    = jd + floor(30.6001.*(m+1)) + b + 1720996.5 + d + h./24 + mi./1440 + s./86400;
