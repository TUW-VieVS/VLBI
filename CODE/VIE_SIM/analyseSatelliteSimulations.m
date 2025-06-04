function [] = analyseSatelliteSimulations(sources, name, antenna, parameter)
    if isempty(name)
        name = 'DUMMY';
    end

    path = pwd;
    if strcmp(path(end-3:end),'WORK')
        path = '../DATA/LEVEL3/';
    elseif strcmp(path(end-3:end),'MISC')
        path = '../../DATA/LEVEL3/';
    end
    
    if ~isfolder([path name])
        error('check folder name and path to LEVEL3 directory!')
    end
    
    load([path name '/x_' parameter.session_name '.mat'])
    load([path name '/atpa_' parameter.session_name '.mat'])
    nsat = length(sources.s);
    
    for isat=1:nsat
        % write Satellite Track to file (10-min interval)
        allLat={};
        allLon={};
        obsLat={};
        obsLon={};
        satName = sources.s(isat).name;
        if strcmp(parameter.vie_init.sc_orbit_file_type, 'tle')
            i = 1;
            u = length(sources.s.x_trf);
        else
            i = 289; %because there are 3 days in sources.s.x_trf -> one before and one after the day of the analysis
            u = 577;
        end
        
        while i <= u
            [lat, lon, ~] = xyz2ell([sources.s(isat).x_trf(i), sources.s(isat).y_trf(i), sources.s(isat).z_trf(i)]);	
            latDegtemp = rad2deg(lat);
            lonDegtemp = rad2deg(lon);      
            if latDegtemp < 0
               latDegString = strcat(num2str((-1)*latDegtemp) , 'S'); 
            else
               latDegString = strcat(num2str(latDegtemp) , 'N'); 
            end
            if lonDegtemp <= 0
                lonDegString = strcat(num2str((-1)*lonDegtemp) , 'W'); 
            else
               lonDegString = strcat(num2str(lonDegtemp) , 'E'); 
            end
            if sources.s(isat).obs(i) == 1 
                obsLat(end+1,1) = {latDegString}; 
                obsLon(end+1,1) = {lonDegString};
            else
                allLat(end+1,1) = {latDegString}; 
                allLon(end+1,1) = {lonDegString};
            end
            i = i + 2;
        end
        Tall = table(allLon, allLat);
        Tall.Properties.VariableNames = {'lat','lon'};
        Tobs = table(obsLon, obsLat);
        Tobs.Properties.VariableNames = {'lat','lon'};
        if nsat ==1 
            writetable(Tall, [path name '/' parameter.session_name '_SatelliteTrack_10min' '.txt'], 'Delimiter', ' ');
            writetable(Tobs, [path name '/' parameter.session_name '_SatelliteTrackObs_10min' '.txt'], 'Delimiter', ' ');
        else
            writetable(Tall, [path name '/' parameter.session_name '_SatelliteTrack_10min_' satName '.txt'], 'Delimiter', ' ');
            writetable(Tobs, [path name '/' parameter.session_name '_SatelliteTrackObs_10min_' satName '.txt'], 'Delimiter', ' ');
        end
    
        
        % write satellite track to file (15-min interval)
        allLat15={};
        allLon15={};
        obsLat15={};
        obsLon15={};
        if strcmp(parameter.vie_init.sc_orbit_file_type, 'tle')
            i = 1;
            u = length(sources.s.x_trf);
        else
            i = 289; %because there are 3 days in sources.s.x_trf -> one before and one after the day of the analysis
            u = 577;
        end
        while i <= u
            [lat, lon, ~] = xyz2ell([sources.s(isat).x_trf(i), sources.s(isat).y_trf(i), sources.s(isat).z_trf(i)]);	
            latDegtemp = rad2deg(lat);
            lonDegtemp = rad2deg(lon);      
            if latDegtemp < 0
               latDegString = strcat(num2str((-1)*latDegtemp) , 'S'); 
            else
               latDegString = strcat(num2str(latDegtemp) , 'N'); 
            end
            if lonDegtemp <= 0
                lonDegString = strcat(num2str((-1)*lonDegtemp) , 'W'); 
            else
               lonDegString = strcat(num2str(lonDegtemp) , 'E'); 
            end
            if sources.s(isat).obs(i) == 1 
                obsLat15(end+1,1) = {latDegString}; 
                obsLon15(end+1,1) = {lonDegString};
            else
                allLat15(end+1,1) = {latDegString}; 
                allLon15(end+1,1) = {lonDegString};
            end
            i = i + 3;
        end
        Tall15 = table(allLon15, allLat15);
        Tall15.Properties.VariableNames = {'lat','lon'};
        Tobs15 = table(obsLon15, obsLat15);
        Tobs15.Properties.VariableNames = {'lat','lon'};
    
        if nsat ==1 
            writetable(Tall15, [path name '/' parameter.session_name '_SatelliteTrack_15min' '.txt'], 'Delimiter', ' ');
            writetable(Tobs15, [path name '/' parameter.session_name '_SatelliteTrackObs_15min' '.txt'], 'Delimiter', ' ');
        else
            writetable(Tall15, [path name '/' parameter.session_name '_SatelliteTrack_15min_' satName '.txt'], 'Delimiter', ' ');
            writetable(Tobs15, [path name '/' parameter.session_name '_SatelliteTrackObs_15min_' satName '.txt'], 'Delimiter', ' ');
        end  
    end
    
    % write txt file with Station coordinates
    antennaLon ={};
    antennaLat = {};
    antennaName = {}; 
    for k=1:length([antenna.x])
        [lat, lon, ~] = xyz2ell([antenna(k).x, antenna(k).y, antenna(k).z]);
        latDeg = rad2deg(lat);
        lonDeg = rad2deg(lon);
        if latDeg < 0
           latDegString = strcat(num2str((-1)*latDeg) , 'S'); 
        else
           latDegString = strcat(num2str(latDeg) , 'N'); 
        end
        if lonDeg <= 0
            lonDegString = strcat(num2str((-1)*lonDeg) , 'W'); 
        else
           lonDegString = strcat(num2str(lonDeg) , 'E'); 
        end
        antennaLon(end+1,1) = {lonDegString};
        antennaLat(end+1,1) = {latDegString};
        antennaName(k) = {strcat(antenna(k).name)};
    end
    
    Tantenna = table(antennaLon, antennaLat, [antennaName(:,:)]');
    writetable(Tantenna, [path name '/' parameter.session_name '_Network.txt'], 'Delimiter', ' ', 'WriteVariableNames',0);
    
    
    if parameter.lsmopt.SatPos.pw_sat == 1
        SatPos =struct();
        for isat = 1:nsat

            clear T dateTime latTable lonTable pos1_val pos1_mx pos2_val pos2_mx pos3_val pos3_mx
        
            nEstTimes = size(x_.sat_pos1(isat).val,1);
            nsim = size(x_.sat_pos1(isat).val,2);
            satName = x_.sat_pos1(isat).name;
            col_satpos1(1,:) = x_.sat_pos1(isat).col;
            col_satpos2(1,:) = x_.sat_pos2(isat).col;
            col_satpos3(1,:) = x_.sat_pos3(isat).col;
            
            pos1_mx = mean(x_.sat_pos1(isat).mx, 2);
            pos2_mx = mean(x_.sat_pos2(isat).mx, 2);
            pos3_mx = mean(x_.sat_pos3(isat).mx, 2);
        
            pos1_val = mean(x_.sat_pos1(isat).val, 2);
            pos2_val = mean(x_.sat_pos2(isat).val, 2);
            pos3_val = mean(x_.sat_pos3(isat).val, 2);
       
            pos1_rep =  sqrt(sum( (x_.sat_pos1(isat).val - pos1_val).^2 , 2) * 1/(nsim-1));
            pos2_rep =  sqrt(sum( (x_.sat_pos2(isat).val - pos2_val).^2 , 2) * 1/(nsim-1));
            pos3_rep =  sqrt(sum( (x_.sat_pos3(isat).val - pos3_val).^2 , 2) * 1/(nsim-1));
       
            pos1_wrms =  sqrt(sum( (x_.sat_pos1(isat).val - pos1_val).^2 .* x_.sat_pos1(isat).mx , 2) .* 1./sum(x_.sat_pos1(isat).mx,2));
            pos2_wrms =  sqrt(sum( (x_.sat_pos2(isat).val - pos2_val).^2 .* x_.sat_pos2(isat).mx , 2) .* 1./sum(x_.sat_pos2(isat).mx,2));
            pos3_wrms =  sqrt(sum( (x_.sat_pos3(isat).val - pos3_val).^2 .* x_.sat_pos3(isat).mx , 2) .* 1./sum(x_.sat_pos3(isat).mx,2)) ;
	        
            for i=1:nEstTimes
                [lat, lon, ~] = xyz2ell(sources.s(isat).posEstIntXtrf(i,:));	
		        latDegtemp = rad2deg(lat);
                lonDegtemp = rad2deg(lon);
                [year, month, day, hour, minu, sec] = mjd2date(sources.s(isat).posEstIntMjd(i));
                dateTime(i,1) = string([num2str(day,'%02d') '.' num2str(month,'%02d') '.' num2str(year) ' ' num2str(hour,'%02d') ':' num2str(minu,'%02d') ':' num2str(sec,'%02d')]);
                
                if latDegtemp < 0
                   latDegString = strcat(num2str((-1)*latDegtemp) , 'S'); 
                else
                   latDegString = strcat(num2str(latDegtemp) , 'N'); 
                end
                if lonDegtemp <= 0
                    lonDegString = strcat(num2str((-1)*lonDegtemp) , 'W'); 
                else
                   lonDegString = strcat(num2str(lonDegtemp) , 'E'); 
                end
		        PosString(i,:) = {latDegString, lonDegString};
            end
        
            latTable = PosString(:,1);
            lonTable = PosString(:,2);
	               
            writematrix(x_.sat_pos1(isat).val,  [path name '/' parameter.session_name '_SatPos1_' satName '.txt'], 'Delimiter', ';');
            writematrix(x_.sat_pos2(isat).val,  [path name '/' parameter.session_name '_SatPos2_' satName '.txt'], 'Delimiter', ';');
            writematrix(x_.sat_pos3(isat).val,  [path name '/' parameter.session_name '_SatPos3_' satName '.txt'], 'Delimiter', ';');
           
            T = table(string(dateTime), pos1_val, pos1_mx, pos1_rep, pos1_wrms, pos2_val, pos2_mx, pos2_rep, pos2_wrms, pos3_val, pos3_mx, pos3_rep, pos3_wrms);
            
            if x_.units.sat_pos1_val(end-2:end)== 'RSW'
                T.Properties.VariableNames = {'Time', 'R_value','R_mx', 'R_rep', 'R_wrms','S_value','S_mx', 'S_rep', 'S_wrms', 'W_value','W_mx', 'W_rep', 'W_wrms'};
            elseif x_.units.sat_pos1_val(end-2:end) == 'NTW'
                T.Properties.VariableNames = {'Time', 'N_value','N_mx', 'N_rep', 'N_wrms','T_value','T_mx', 'T_rep', 'T_wrms', 'W_value','W_mx', 'W_rep', 'W_wrms'};
            else
                T.Properties.VariableNames = {'Time', 'SatPos1_val','SatPos1_mx', 'SatPos1_rep', 'SatPos1_wrms', 'SatPos2_val','SatPos2_mx', 'SatPos2_rep', 'SatPos2_wrms', 'SatPos3_val','SatPos3_mx', 'SatPos3_rep', 'SatPos3_wrms'};
            end
            idx = find(x_.units.sat_pos1_val == ' ', 1, 'last');
            fprintf(strcat('\n', satName,'\n'))
            fprintf(strcat('Number of simulations:  ', num2str(nsim) , '\n'))
            fprintf(strcat('Reference Frame:  ' , x_.units.sat_pos1_val(idx:end) , '\n\n'))
            disp(T)
            
            T2 = table(dateTime(1:end,1),latTable(1:end,1), lonTable(1:end,1), pos1_val(1:end,1), pos1_mx(1:end,1), pos1_rep, pos1_wrms, pos2_val(1:end,1), pos2_mx(1:end,1), pos2_rep, pos2_wrms, pos3_val(1:end,1), pos3_mx(1:end,1), pos3_rep, pos3_wrms);
            T2.Properties.VariableNames = {'date', 'lat','lon','SatPos1_val','SatPos1_mx', 'SatPos1_rep', 'SatPos1_wrms', 'SatPos2_val','SatPos2_mx', 'SatPos2_rep', 'SatPos2_wrms', 'SatPos3_val','SatPos3_mx', 'SatPos3_rep', 'SatPos3_wrms'};
            writetable(T2, [path name '/' parameter.session_name '_SatelliteEstimates_' satName '.txt'], 'Delimiter', ';');
            writetable(T2, [path name '/' parameter.session_name '_SatelliteEstimates_' satName '.xlsx']);   
        
            %Correlations
            N = atpa_.mat;
            Qxx = inv(N);
            Qxx=full(Qxx);
            QxxSatPos = Qxx([col_satpos1,col_satpos2, col_satpos3], [col_satpos1,col_satpos2, col_satpos3]);

            KorrMatrix = zeros(length(col_satpos1)*3,length(col_satpos1)*3);
            for i=1:length(col_satpos1)*3
                    for j=i+1:length(col_satpos1)*3
                            KorrMatrix(i,j) = QxxSatPos(i,j)./(sqrt(QxxSatPos(i,i).*QxxSatPos(j,j)));
                            KorrMatrix(j,i)=KorrMatrix(i,j);
                    end
            end
        
            writematrix(QxxSatPos, [path name '/' parameter.session_name '_QxxSatPos_' satName '.txt'], 'Delimiter','tab');
            writematrix(KorrMatrix, [path name '/' parameter.session_name '_CorrCoef_' satName '.txt'], 'Delimiter','tab');
        end
    end
    

    %Orbital elements
    if parameter.lsmopt.KepEle.estKepEle==1 

        
        for isat = 1:nsat
            ALL_rep_mx = [];
            k=1;
            wall=0;
            col =[];
            for iKep = 1:6

                clear T dateTime latTable lonTable omega_val omega_mx omega_rep omega_wrms PosString
                if parameter.lsmopt.KepEle.('estKepEle'+ string(iKep))  == 1
                    fprintf('Keplerian Element ' + string(iKep) + '\n');
                    if iKep == 1
                        fprintf('Units: [cm] \n');
                    elseif iKep==2
                        fprintf('Units: [ ] \n');
                    else
                        fprintf('Units: [mas]\n');
                    end
            
                    KepEleNum = 'KepEle' + string(iKep);
                    EstIntMjdText = 'EstIntMjdKepEle' + string(iKep);
                    EstIntXtrfText= 'EstIntXtrfKepEle' + string(iKep);

                    nEstTimes = size(x_.(KepEleNum)(isat).val,1);
                    nsim = size(x_.(KepEleNum)(isat).val,2);
                    satName = x_.(KepEleNum)(isat).name;
                    
                    KepEle_mx = mean(x_.(KepEleNum)(isat).mx, 2);
                    KepEle_val = mean(x_.(KepEleNum)(isat).val, 2);
                    KepEle_rep =  sqrt(sum( (x_.(KepEleNum)(isat).val - KepEle_val).^2 , 2) * 1/(nsim-1));
                    KepEle_wrms =  sqrt(sum( (x_.(KepEleNum)(isat).val - KepEle_val).^2 .* x_.(KepEleNum)(isat).mx , 2) .* 1./sum(x_.(KepEleNum)(isat).mx,2));
                    
                    for i=1:nEstTimes
                        [year, month, day, hour, minu, sec] = mjd2date(sources.s(isat).(EstIntMjdText)(i));
                        dateTime(i,1) = string([num2str(day,'%02d') '.' num2str(month,'%02d') '.' num2str(year) ' ' num2str(hour,'%02d') ':' num2str(minu,'%02d') ':' num2str(sec,'%02d')]);
                        [lat, lon, ~] = xyz2ell(sources.s(isat).(EstIntXtrfText)(i,:));	
		                latDegtemp = rad2deg(lat);
                        lonDegtemp = rad2deg(lon);           
                        if latDegtemp < 0
                           latDegString = strcat(num2str((-1)*latDegtemp) , 'S'); 
                        else
                           latDegString = strcat(num2str(latDegtemp) , 'N'); 
                        end
                        if lonDegtemp <= 0
                            lonDegString = strcat(num2str((-1)*lonDegtemp) , 'W'); 
                        else
                           lonDegString = strcat(num2str(lonDegtemp) , 'E'); 
                        end
		                PosString(i,:) = {latDegString, lonDegString};
                    end
                
                    latTable = PosString(:,1);
                    lonTable = PosString(:,2);
                    writematrix(x_.(KepEleNum)(isat).val,  [path name '/' parameter.session_name '_KepEle' num2str(iKep) '_' satName '.txt'], 'Delimiter', ';');
                    
                    fprintf(satName);
                    fprintf('\n')
                    T = table(string(dateTime), KepEle_val, KepEle_mx, KepEle_rep, KepEle_wrms);
                    T.Properties.VariableNames = {'Time', 'value','mx', 'rep ', 'wrms'};     
                    fprintf(strcat('Number of simulations:  ', num2str(nsim) , '\n'))
                    disp(T)
                    fprintf('\n')
            
                    Twrite = table(string(dateTime),latTable, lonTable, KepEle_val, KepEle_mx, KepEle_rep, KepEle_wrms);
                    Twrite.Properties.VariableNames = {'Time','lat','lon', 'value ','mx ', 'rep ', 'wrms '};
                    writetable(Twrite, [path name '/' parameter.session_name '_KepEleEstimates' num2str(iKep) '_' satName '.txt'], 'Delimiter', ';');
    
                    if parameter.lsmopt.KepEle.('estKepEle'+ string(iKep))  && isscalar(x_.(KepEleNum)(isat).col)
                        col = [col; x_.(KepEleNum)(isat).col];
                    end
                    if isscalar(x_.(KepEleNum)(isat).col) %if all orbital elements are estimated for one satellite
                        ALL_rep_mx(k,1) = iKep;
                        ALL_rep_mx(k,2) = KepEle_rep;
                        ALL_rep_mx(k,3) = KepEle_mx;
                        ALL_rep_mx(k,4) = KepEle_val;
                        k=k+1;
                        wall =1;
                    end
                    nEst(i)= nEstTimes;
                end
            end
            %Correlations only if every orbital element is estimated once
            %and all orbital elements are estimated
            if ~isempty(col) && max(nEst) ==1 && length(col) == 6
                N = atpa_.mat;
                Qxx = inv(N);
                QxxKepEle = Qxx(col, col);
                
                iQ = 1; 
                jQ = 2;
                k = 1;
                KorrKoef = zeros((size(QxxKepEle,1)*(size(QxxKepEle,1)-1)/2),2);
                KorrMatrix = zeros(6,6);
                for i=1:6
                    jQ = i+1;
                    if parameter.lsmopt.KepEle.('estKepEle' + string(i)) == 1
                        for j=i+1:6
                            if parameter.lsmopt.KepEle.('estKepEle' + string(j)) == 1
                                KorrMatrix(i,j) = QxxKepEle(iQ,jQ)./(sqrt(QxxKepEle(iQ,iQ).*QxxKepEle(jQ,jQ)));
                                KorrMatrix(j,i)=KorrMatrix(i,j);

                                KorrKoef(k,1) = i*10 + j;
                                KorrKoef(k,2) = KorrMatrix(i,j);
                                jQ = jQ +1;
                                k = k +1;
                            end
                        end
                        iQ = iQ +1;
                    end
                end
              
                writematrix(KorrKoef, [path name '/' parameter.session_name '_KorrKepEle_' satName '.txt'], 'Delimiter','tab');
            
                hearder_num = [0,1,2,3,4,5,6];
                output_matrix=[hearder_num; hearder_num(2:end)' KorrMatrix];
                writematrix(output_matrix, [path name '/' parameter.session_name '_KorrMatrix' satName '.xls']);
                output_matrix_small = output_matrix;
                
                k=2;
                for j=1:size(KorrMatrix,1)
                    if sum(KorrMatrix(j,1:6)) == 0
                        output_matrix_small(k,:) = [];
                    else
                        k=k+1;
                    end
                end
                k=2;
                for j=1:size(KorrMatrix,2)
                    if sum(KorrMatrix(1:6,j)) == 0
                        output_matrix_small(:,k) = [];
                    else
                        k=k+1;
                    end
                end
                writematrix(output_matrix_small, [path name '/' parameter.session_name '_KorrMatrixSmall' satName '.xls']);
                writematrix(output_matrix_small, [path name '/' parameter.session_name '_KepEleKorrMs' satName '.txt'], 'Delimiter',';');
                writematrix(output_matrix, [path name '/' parameter.session_name '_KepEleKorrM' satName '.txt'], 'Delimiter',';');
            end

            if wall
                row_header = {'iKep ', 'rep ', 'mx ', 'val '};
                ALL_rep_mx=[row_header; num2cell(ALL_rep_mx)];
                writecell(ALL_rep_mx, [path name '/' parameter.session_name '_KepEleRepMxALL_' satName '.txt'], 'Delimiter',';');
            end
        end
    end
end