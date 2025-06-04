
function [rntw, vntw, transmat] = rv2ntw(reci, veci)

        % each of the components must be unit vectors
        % tangential component   
        tvec = veci/norm(veci);

        %cross-track component
        wvec = cross(reci, veci);
        wvec = wvec/norm(wvec);

        %normal component
        nvec = cross(tvec, wvec);
        nvec = nvec/norm(nvec);

        % assemble transformation matrix from to rsw frame (individual
        % components arranged in row vectors)
        transmat(1,1) = nvec(1);
        transmat(1,2) = nvec(2);
        transmat(1,3) = nvec(3);
        transmat(2,1) = tvec(1);
        transmat(2,2) = tvec(2);
        transmat(2,3) = tvec(3);
        transmat(3,1) = wvec(1);
        transmat(3,2) = wvec(2);
        transmat(3,3) = wvec(3);

        rntw = transmat*reci;
        vntw = transmat*veci;
end