function TDB = tt2tdb(TT)

% conversion TT (given in mjd[days] in TDB (mjd[days]) time
% vectorized

% reference: Seidelmann and Fukushima 1992, 'Why new time scales', table 1
% L. Plank, Aug. 2009


% elapsed time from J2000.0 in Julian centuries
T   = (TT - 51544.5)./36525;


P = 0.0016568 * sin(deg2rad(35999.37 * T + 357.5)) + ...
   (0.0000224 * sin(deg2rad(32964.5  * T + 246))   + ...
   (0.0000138 * sin(deg2rad(71998.7  * T + 355))   + ...
   (0.0000048 * sin(deg2rad( 3034.9  * T +  25))   + ...
   (0.0000047 * sin(deg2rad(34777.3  * T + 230)))))); %[sec]

TDB = TT + P./86400; 

% modest

g = deg2rad(357.528 + 35999.05 * T);
TDB= TT + (0.001658*sin(g+0.0167*sin(g)))/86400;