% LOG
%  04.05.2015 M. Madzak: index i -> k
%  04.05.2015 M. Madzak: Added post-seismic deformation
%  30.05.2015 A. Girdiuk: The loop over the stations is replaced with logical ind_stat

function [x,y,z] = corap(antenna,stat,mjd)

ind_stat = ismember({antenna.name},stat);
ep = antenna(ind_stat).epoch;
x = antenna(ind_stat).x + (mjd-ep)/365.25* antenna(ind_stat).vx;
y = antenna(ind_stat).y + (mjd-ep)/365.25* antenna(ind_stat).vy;
z = antenna(ind_stat).z + (mjd-ep)/365.25* antenna(ind_stat).vz;

% add also post-seismic deformation (pdf)
if isfield(antenna(ind_stat),'psd')
    if ~isempty(antenna(ind_stat).psd)
        % get psd
        tempX=cPostSeismDeform(mjd,antenna(ind_stat)); % [m]

        x=x+tempX(1);
        y=y+tempX(2);
        z=z+tempX(3);

    end
end
