function [hh,mm,ss] = rad2houminsec(ang)

ang_deg = ang*180/(pi*15);
hh = floor(ang_deg);
ang_min = (ang_deg - hh)*60;
mm = floor(ang_min);
ss = (ang_min - mm)*60;