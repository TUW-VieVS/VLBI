function [tau,uu,vv]=modDelay(mod, tele1, tele2, freq, rfreq, ra, de, mjd)
%Usage: [tau,uu,vv]=modDelay(mod, tele1, tele2, freq, rfreq, RA, DEC, jd);
%
% Generates a multiband delay from a model (of Gaussian components) for a
% given pair of telescopes, observing frequencies, source and observation
% date. 
%
% mod is a matrix of parameters for source component models. Each row is a
% separate component of the form: 
% [Amp, FWHMmajor, FWHMminor, majAxisAngle, RAoffset, Decoffset]
% where Amp is amplitude in arbitrary units, FWHMmajor and FWHMminor are
% the FWHM of the major and minor axes on the sky (in mas) and majAxisAngle
% is the position angle of the major axis (in degrees, with zero indicating
% alignment in RA). RAoffset and Decoffset are offsets of the centroid from
% the phase centre (in mas). 
% tele1 and tele2 are vectors of the xyz positions of the two telescopes
% involved in the observation (in MHz)
% freq is a vector of the observing frequency bands to be used in the
% multiband delay estimation (in MHz).
% rfreq is the reference frequency for the observation which is used in
% estimation of the uu,vv tracks. 
% RA and DEC are the j2000 Right Ascension and Declination of the source in
% radians
% mjd is the Modified Julian Day of the observation 
%
% tau is the estimated delay due to structure (in ps)
% uu and vv are the baseline vector of the two telescopes during the
% observation (in Megalambda).
% CHANGES:
% 22. Jan 14 by Lucia Plank: using the code of Vie_Sched for the
%               calculation of u, v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

constants
global c 

lambda=c/rfreq/1E6; % Wavelength of reference frequency.
bl=tele1-tele2;
blx=bl(1);
bly=bl(2);
blz=bl(3);

twopi = 2 * pi;

gmst = tgmst(mjd);
ha = gmst - ra;
while (ha > pi)
    ha = ha - twopi;
end
while (ha < -pi)
    ha = ha + twopi;
end
u = (blx * sin(ha) + bly * cos(ha))/lambda;
%v = blz + sin(de) * (-blx * cos(ha) + bly * sin(ha))/lambda;%
v = (blz*cos(de) + sin(de) * (-blx * cos(ha) + bly * sin(ha)))/lambda;

uu=u/1E6;
vv=v/1E6;

for k=1:length(freq)
    Vis(k,:)=modVis(mod, uu*freq(k)/rfreq, vv*freq(k)/rfreq);
end
for k=1:length(uu) 
    P=polyfit((freq(:)-mean(freq))*1E6, angle(Vis(:,k)), 1);
    tau(k)=P(1)/ 2 / pi * 1E12 ; 
end


