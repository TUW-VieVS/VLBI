function [sat_c] = manipulateSourcesNTW(sat ,i_ntw )
    sat_c = sat; 
    for i=1:length(sat.x_crf)
        crfScPos = [sat.x_crf(i) sat.y_crf(i) sat.z_crf(i)];
        crfScVel = [sat.vx_crf(i) sat.vy_crf(i) sat.vz_crf(i)];
        [~, ~, transmatNTW] = rv2ntw(crfScPos', crfScVel');
        crfScPosNTW = transmatNTW*crfScPos';
        shift = 0.01; % in m --> 0.01 = 1 cm
        if i_ntw ==1 
            crfScPosNTW = [crfScPosNTW(1) + shift, crfScPosNTW(2), crfScPosNTW(3)]; % add to N
        elseif i_ntw==2
            crfScPosNTW = [crfScPosNTW(1), crfScPosNTW(2) + shift, crfScPosNTW(3)]; % add to T
        elseif i_ntw==3
            crfScPosNTW = [crfScPosNTW(1), crfScPosNTW(2), crfScPosNTW(3)+ shift]; % add to W
        end
        crfScPos_changed = transmatNTW'*crfScPosNTW';
        sat_c.x_crf(i) = crfScPos_changed(1);
        sat_c.y_crf(i) = crfScPos_changed(2);
        sat_c.z_crf(i) = crfScPos_changed(3);
    end
end