function [medDelay,SI, TauMap,uu,vv]=SI_estimator(mod, freq, rfreq)
% Usage: [medDelay,SI,TauMap,uu,vv]=SI_estimator(mod, freq, rfreq)
% This function returns medDelay, an estimated median delay in ps
% SI which is the structure index (1 + 2*log10(medDelay)). 
% It optionally returns a map of the delays and the corresponding uu,vv matrix it is calculated
% from. uu,vv are caluculated as 512x512 grids with a maximum corresponding
% to the diamater of the Earth divided by the reference wavelength

% freq is a vector of frequencies to be used in the delay estimation. It
% should be in MHz. freq=[8217 8256 8357 8517 8737 8857 8917 8937]; is
% appropriate for r1/4 type experiments.
% rfreq is the reference frequency - generally it is best to pick a
% frequency in the middle of the freq range. 

maxuv=12742E3 / (299792458 / rfreq/1E6)/1E6;
[uu,vv]=meshgrid(-maxuv:(2*maxuv)/511:maxuv);
% clear uu vv
% load('uu_cont11-1');
% load('vv_cont11-1');
TauMap=zeros(size(uu));
a=find(sqrt(uu.^2 +vv.^2 )<maxuv);
for i = 1:512
    for j=1:512
        for k=1:length(freq)
            A(k)=modVis(mod, uu(i,j)*freq(k)/rfreq, vv(i,j)*freq(k)/rfreq);
        end    
        P=polyfit( ((freq*1E6)-mean((freq*1E6))), angle(A), 1);
        TauMap(i,j)=P(1) / 2 / pi / 1E-12;
    end
    %i/512
end

medDelay=median(abs(TauMap(a)));
SI=1+2*log10(medDelay);