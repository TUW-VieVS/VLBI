% PSD taken from a time series - as it is in DTRF2020


function cpsd = cPostSeismDeformTS(mjd,antenna)

nEp=length(mjd);
nStat=length(antenna);

cpsd=zeros(3,nEp,nStat);

for iStat = 1:nStat
    if ~isempty(antenna(iStat).psdTS)
        cpsd(1,1:end,iStat) = antenna(iStat).psdTS(1);
        cpsd(2,1:end,iStat) = antenna(iStat).psdTS(2);
        cpsd(3,1:end,iStat) = antenna(iStat).psdTS(3);
    end

end

