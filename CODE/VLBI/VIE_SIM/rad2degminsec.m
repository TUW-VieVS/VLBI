function [dd,mm,ss] = rad2degminsec(ang)

ang_deg = abs(ang*180/pi);
dd = floor(ang_deg);
ang_min = (ang_deg - dd)*60;
mm = floor(ang_min);
ss = (ang_min - mm)*60;
if (ang < 0.0)
    dd = dd*(-1);
end