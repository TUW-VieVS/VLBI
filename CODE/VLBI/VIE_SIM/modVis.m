function Vis=modVis(mod, uu, vv);
%Usage: Vis=modVis(mod, uu, vv);
%
% Generates complex visibilities from a model (Gaussian) on a given u,v track 
% mod is a vector of parameters for the model
% mod=[Amp, FWHMmajor, FWHMminor, majAxisAngle, RAoffset, Decoffset]
% where Amp is amplitude in arbitratu units, FWHMmajor and FWHMminor are
% the FWHM of the major and minor axes on the sky (in mas) and majAxisAngle
% is the position angle of the major axis (in degrees, with zero indicating
% alignment in RA). RAoffset and Decoffset are offsets of the centroid from
% the phase centre (in mas). uu and vv are in units of Megalambda. 
Vis=0;
for j=1:length(mod(:,1))
    %offsetscale is the conversion factor between the offset in mas and the
    %fringe spacing in Megalambda
    offsetscale=1E-6 * 180/pi * 3600E3;
    %FWHMscale is ~109 which is the conversion factor between the
    %FWHM in the image plane (in mas) and the FWHM in the visibility plane
    %(in Megalambda) (182, taken from Miriad manual and in Briggs, Schwab and Swanek 1999),
    %divided by sqrt(log(2))*2 (the conversion between FWHM and c in a Gaussian profile. 
    fwhmscale = 182 / (sqrt(log(2))*2);
    Vis=Vis + (mod(j,1) * exp(-( ((cos(pi/180*(90-mod(j,4))) * uu - sin(pi/180*(90-mod(j,4))) * vv )/fwhmscale*mod(j,3) ).^2 + ( (sin(pi/180*(90-mod(j,4))) * uu + cos(pi/180*(90-mod(j,4))) * vv )/fwhmscale*mod(j,2) ).^2) ) .* exp( i * 2* pi / offsetscale * (mod(j,5)*uu + mod(j,6)*vv)));
end
%Vis=mod(1) * exp(-( ((cos(pi/180*(90-mod(4))) * uu - sin(pi/180*(90-mod(4))) * vv )/109*mod(3) ).^2 + ( (sin(pi/180*(90-mod(4))) * uu + cos(pi/180*(90-mod(4))) * vv )/109*mod(2) ).^2) ) .* exp( i * 2* pi / 206.26 * (mod(5)*uu + mod(6)*vv));