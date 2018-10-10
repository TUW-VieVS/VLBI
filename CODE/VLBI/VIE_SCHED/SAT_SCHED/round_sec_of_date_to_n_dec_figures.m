% -------------------------------------------------------------------------
%
%                              round_sec_of_date_to_n_dec_figures
%
%   Rounds the seconds of a given date to a arbitrary number of decimal
%   places. Hour, minute, day of year and year are adjusted, if necessary.
%
%   Author: 
%       2013-11-14: Andreas Hellerschmied (heller182@gmx.at)
%   
%   changes       :
%           
%
%   inputs        :
%   - year          : Year
%   - days_of_year  : Number of days in given year
%   - hr            : Hours
%   - min           : Minutes
%   - sec           : Seconds (float)
%   - dec_places    : Number of decimals
%     
%
%   outputs       :
%   - year          : Year
%   - days_of_year  : Number of days in given year
%   - hr            : Hours
%   - min           : Minutes
%   - sec           : Seconds, rounded to "dec_places" number of decimal 
%                     places (float).
%    
%
%   locals        :
% 
%
%   coupling      :
%   - leapyear      : Finds leap years.
%   
%
%   references    :
%
%-------------------------------------------------------------------------

function [year, days_of_year, hr, min, sec] = round_sec_of_date_to_n_dec_figures(year, days_of_year, hr, min, sec, dec_places)

    temp = (round(sec * 10^dec_places) / (10^dec_places));

    if (temp == 60)
       sec = 0.0;
       min = min + 1;
       if (min == 60)
           min = 0;
           hr = hr + 1;
           if (hr == 24)
              hr = 0;
              days_of_year = days_of_year + 1;
              
              if (leapyear(year)) % year has 366 days
                  if (days_of_year == 367)
                      year = year + 1;
                      days_of_year = 1;
                  end
              else % year has 365 days
                  if (days_of_year == 366)
                      year = year + 1;
                      days_of_year = 1;
                  end
              end
           end
       end
    else
        sec = temp;
    end

return;

