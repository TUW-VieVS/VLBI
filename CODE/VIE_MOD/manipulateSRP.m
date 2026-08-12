function [sourcesChanged_list] = manipulateSRP(sources, FSO_prev, FSO_after, FRP_prev, FRP_after, FSO, FRP, T2C_s)
    sourceChanged_D0 = sources;
    sourceChanged_Y0 = sources;
    sourceChanged_B0 = sources;
    sourceChanged_DC = sources;
    sourceChanged_YC = sources;
    sourceChanged_BC = sources;
    sourceChanged_DS = sources;
    sourceChanged_YS = sources;
    sourceChanged_BS = sources;    

    dSec = 0.1;

    [coerpr1, trpr1, hrpr1, ~, ~, ~, tbnd1, ~, ~, sclpar1, ~, satnum1, ~, ~, ~] = readrpr(string('../ORBIT/FRP/') + string(FRP_prev));
    [coerpr, trpr, hrpr, ~, ~, ~, tbnd, ~, ~, sclpar, ~, satnum, ~, ~, ~] = readrpr(FRP);
    [coerpr3, trpr3, hrpr3, ~, ~, ~, tbnd3, ~, ~, sclpar3, ~, satnum3, ~, ~, ~] = readrpr(string('../ORBIT/FRP/') + string(FRP_after));

    [coeff1,t01,hstep1,~,~,~,tbound1,tosc1,oscele1,satnum1,narc1,descr1,source1] = readstd(string('../ORBIT/FSO/') + string(FSO_prev));
    [coeff,t0,hstep,~,~,~,tbound,tosc,oscele,satnum,narc,descr,source] = readstd(FSO);
    [coeff3,t03,hstep3,~,~,~,tbound3,tosc3,oscele3,satnum3,narc3,descr3,source3] = readstd(string('../ORBIT/FSO/') + string(FSO_after));

    for isat=1:length(sources.s)
        sat = sources.s(isat);

        for j=1:length(sources.s(isat).x_crf)
            sourceTimeUTC = cell2mat({sources.s(isat).mjd(j)}); 
            dt = datetime(sourceTimeUTC, 'ConvertFrom', 'mjd');
            [year, month, day, hour, minute, sec] = datevec(dt);
            leap_sec_tai_utc = tai_utc(sourceTimeUTC); 
            leap_sec_utc_gps = leap_sec_tai_utc - 19;
            sec = sec + leap_sec_utc_gps;
            time_gps = modjuldat(year, month, day, hour, minute, sec)';
            % time_gps_cells = num2cell(time_gps);
            % time_gps = sources.s(isat).mjd_gps(j);
            time = sourceTimeUTC;
            if j<289 
               [drdpar_FRP,~] = getrpr(sat.fso_name,time,coerpr1,trpr1,hrpr1,tbnd1,sclpar1,satnum1);
            elseif j>=289 && j<577
                [drdpar_FRP,~] = getrpr(sat.fso_name,time,coerpr,trpr,hrpr,tbnd,sclpar,satnum);
            else
               [drdpar_FRP,~] = getrpr(sat.fso_name,time,coerpr3,trpr3,hrpr3,tbnd3,sclpar3,satnum3);
            end
            drdD0 = drdpar_FRP(7,:); %m/m/s2
            drdY0 = drdpar_FRP(8,:); %m/m/s2
            drdB0 = drdpar_FRP(9,:); %m/m/s2
            drdDC = drdpar_FRP(10,:); %m/m/s2
            drdYC = drdpar_FRP(11,:); %m/m/s2
            drdBC = drdpar_FRP(12,:); %m/m/s2
            drdDS = drdpar_FRP(13,:); %m/m/s2
            drdYS = drdpar_FRP(14,:); %m/m/s2
            drdBS = drdpar_FRP(15,:); %m/m/s2
        
            faktor = 3*10^(-9); %m/s2 (3nm/s2)
            delta_rD0 = (faktor * drdD0); %m/s2 * m/m/s2 -> m
            delta_rY0 = (faktor * drdY0); %m
            delta_rB0 = (faktor * drdB0); %m
            delta_rDC = (faktor * drdDC); %m
            delta_rYC = (faktor * drdYC); %m
            delta_rBC = (faktor * drdBC); %m
            delta_rDS = (faktor * drdDS); %m
            delta_rYS = (faktor * drdYS); %m
            delta_rBS = (faktor * drdBS); %m

            sourceChanged_D0.s(isat).x_crf(j) = sources.s(isat).x_crf(j) + delta_rD0(1);
            sourceChanged_D0.s(isat).y_crf(j) = sources.s(isat).y_crf(j) + delta_rD0(2);
            sourceChanged_D0.s(isat).z_crf(j) = sources.s(isat).z_crf(j) + delta_rD0(3);

            sourceChanged_Y0.s(isat).x_crf(j) = sources.s(isat).x_crf(j) + delta_rY0(1);
            sourceChanged_Y0.s(isat).y_crf(j) = sources.s(isat).y_crf(j) + delta_rY0(2);
            sourceChanged_Y0.s(isat).z_crf(j) = sources.s(isat).z_crf(j) + delta_rY0(3);

            sourceChanged_B0.s(isat).x_crf(j) = sources.s(isat).x_crf(j) + delta_rB0(1);
            sourceChanged_B0.s(isat).y_crf(j) = sources.s(isat).y_crf(j) + delta_rB0(2);
            sourceChanged_B0.s(isat).z_crf(j) = sources.s(isat).z_crf(j) + delta_rB0(3);

            sourceChanged_DC.s(isat).x_crf(j) = sources.s(isat).x_crf(j) + delta_rDC(1);
            sourceChanged_DC.s(isat).y_crf(j) = sources.s(isat).y_crf(j) + delta_rDC(2);
            sourceChanged_DC.s(isat).z_crf(j) = sources.s(isat).z_crf(j) + delta_rDC(3);

            sourceChanged_YC.s(isat).x_crf(j) = sources.s(isat).x_crf(j) + delta_rYC(1);
            sourceChanged_YC.s(isat).y_crf(j) = sources.s(isat).y_crf(j) + delta_rYC(2);
            sourceChanged_YC.s(isat).z_crf(j) = sources.s(isat).z_crf(j) + delta_rYC(3);

            sourceChanged_BC.s(isat).x_crf(j) = sources.s(isat).x_crf(j) + delta_rBC(1);
            sourceChanged_BC.s(isat).y_crf(j) = sources.s(isat).y_crf(j) + delta_rBC(2);
            sourceChanged_BC.s(isat).z_crf(j) = sources.s(isat).z_crf(j) + delta_rBC(3);

            sourceChanged_DS.s(isat).x_crf(j) = sources.s(isat).x_crf(j) + delta_rDS(1);
            sourceChanged_DS.s(isat).y_crf(j) = sources.s(isat).y_crf(j) + delta_rDS(2);
            sourceChanged_DS.s(isat).z_crf(j) = sources.s(isat).z_crf(j) + delta_rDS(3);

            sourceChanged_YS.s(isat).x_crf(j) = sources.s(isat).x_crf(j) + delta_rYS(1);
            sourceChanged_YS.s(isat).y_crf(j) = sources.s(isat).y_crf(j) + delta_rYS(2);
            sourceChanged_YS.s(isat).z_crf(j) = sources.s(isat).z_crf(j) + delta_rYS(3);

            sourceChanged_BS.s(isat).x_crf(j) = sources.s(isat).x_crf(j) + delta_rBS(1);
            sourceChanged_BS.s(isat).y_crf(j) = sources.s(isat).y_crf(j) + delta_rBS(2);
            sourceChanged_BS.s(isat).z_crf(j) = sources.s(isat).z_crf(j) + delta_rBS(3);

            if j> 2 && j < 700
                if j<289 
                    [drdpar_FRPm,~] = getrpr(sat.fso_name,time-dSec/86400,coerpr1,trpr1,hrpr1,tbnd1,sclpar1,satnum1); %1sec before
                    [drdpar_FRPp,~] = getrpr(sat.fso_name,time+dSec/86400,coerpr1,trpr1,hrpr1,tbnd1,sclpar1,satnum1); %1sec after
                elseif j==289
                    [drdpar_FRPm,~] = getrpr(sat.fso_name,time-dSec/86400,coerpr1,trpr1,hrpr1,tbnd1,sclpar1,satnum1); %1sec before
                    [drdpar_FRPp,~] = getrpr(sat.fso_name,time+dSec/86400,coerpr,trpr,hrpr,tbnd,sclpar,satnum); %1sec after
                elseif j>289 && j<577
                    [drdpar_FRPm,~] = getrpr(sat.fso_name,time-dSec/86400,coerpr,trpr,hrpr,tbnd,sclpar,satnum); %1sec before
                    [drdpar_FRPp,~] = getrpr(sat.fso_name,time+dSec/86400,coerpr,trpr,hrpr,tbnd,sclpar,satnum); %1sec after
                elseif j == 577
                    [drdpar_FRPm,~] = getrpr(sat.fso_name,time-dSec/86400,coerpr,trpr,hrpr,tbnd,sclpar,satnum); %1sec before
                    [drdpar_FRPp,~] = getrpr(sat.fso_name,time+dSec/86400,coerpr3,trpr3,hrpr3,tbnd3,sclpar3,satnum3); %1sec after
                else
                    [drdpar_FRPm,~] = getrpr(sat.fso_name,time-dSec/86400,coerpr3,trpr3,hrpr3,tbnd3,sclpar3,satnum3); %1sec before
                    [drdpar_FRPp,~] = getrpr(sat.fso_name,time+dSec/86400,coerpr3,trpr3,hrpr3,tbnd3,sclpar3,satnum3); %1sec after
                end

                drdD0m = drdpar_FRPm(7,:);
                drdY0m = drdpar_FRPm(8,:);
                drdB0m = drdpar_FRPm(9,:);
                drdDCm = drdpar_FRPm(10,:);
                drdYCm = drdpar_FRPm(11,:);
                drdBCm = drdpar_FRPm(12,:);
                drdDSm = drdpar_FRPm(13,:);
                drdYSm = drdpar_FRPm(14,:);
                drdBSm = drdpar_FRPm(15,:);
    
                drdD0p = drdpar_FRPp(7,:);
                drdY0p = drdpar_FRPp(8,:);
                drdB0p = drdpar_FRPp(9,:);
                drdDCp = drdpar_FRPp(10,:);
                drdYCp = drdpar_FRPp(11,:);
                drdBCp = drdpar_FRPp(12,:);
                drdDSp = drdpar_FRPp(13,:);
                drdYSp = drdpar_FRPp(14,:);
                drdBSp = drdpar_FRPp(15,:);
                
                delta_rD0m = (faktor * drdD0m);
                delta_rY0m = (faktor * drdY0m);
                delta_rB0m = (faktor * drdB0m);
                delta_rDCm = (faktor * drdDCm);
                delta_rYCm = (faktor * drdYCm);
                delta_rBCm = (faktor * drdBCm);
                delta_rDSm = (faktor * drdDSm);
                delta_rYSm = (faktor * drdYSm);
                delta_rBSm = (faktor * drdBSm);
    
                delta_rD0p = (faktor * drdD0p);
                delta_rY0p = (faktor * drdY0p);
                delta_rB0p = (faktor * drdB0p);
                delta_rDCp = (faktor * drdDCp);
                delta_rYCp = (faktor * drdYCp);
                delta_rBCp = (faktor * drdBCp);
                delta_rDSp = (faktor * drdDSp);
                delta_rYSp = (faktor * drdYSp);
                delta_rBSp = (faktor * drdBSp);
                                
                if j<289 
                     [posp,~] = getorb_par(sat.fso_name,time+dSec/86400,coeff1,t01,hstep1,tbound1,satnum1);
                     [posm,~] = getorb_par(sat.fso_name,time-dSec/86400,coeff1,t01,hstep1,tbound1,satnum1);
                elseif j==289
                     [posp,~] = getorb_par(sat.fso_name,time+dSec/86400,coeff,t0,hstep,tbound,satnum);
                     [posm,~] = getorb_par(sat.fso_name,time-dSec/86400,coeff1,t01,hstep1,tbound1,satnum1);
                elseif j>289 && j<577
                     [posp,~] = getorb_par(sat.fso_name,time+dSec/86400,coeff,t0,hstep,tbound,satnum);
                     [posm,~] = getorb_par(sat.fso_name,time-dSec/86400,coeff,t0,hstep,tbound,satnum);
                elseif j == 577
                     [posp,~] = getorb_par(sat.fso_name,time+dSec/86400,coeff3,t03,hstep3,tbound3,satnum3);
                     [posm,~] = getorb_par(sat.fso_name,time-dSec/86400,coeff,t0,hstep,tbound,satnum);
                else
                     [posp,~] = getorb_par(sat.fso_name,time+dSec/86400,coeff3,t03,hstep3,tbound3,satnum3);
                     [posm,~] = getorb_par(sat.fso_name,time-dSec/86400,coeff3,t03,hstep3,tbound3,satnum3);
                end
    
                vel_D0 = ((posp + delta_rD0p) - (posm + delta_rD0m))/(dSec*2);
                vel_Y0 = ((posp + delta_rY0p) - (posm + delta_rY0m))/(dSec*2);
                vel_B0 = ((posp + delta_rB0p) - (posm + delta_rB0m))/(dSec*2);
                vel_DC = ((posp + delta_rDCp) - (posm + delta_rDCm))/(dSec*2);
                vel_YC = ((posp + delta_rYCp) - (posm + delta_rYCm))/(dSec*2);
                vel_BC = ((posp + delta_rBCp) - (posm + delta_rBCm))/(dSec*2);
                vel_DS = ((posp + delta_rDSp) - (posm + delta_rDSm))/(dSec*2);
                vel_YS = ((posp + delta_rYSp) - (posm + delta_rYSm))/(dSec*2);
                vel_BS = ((posp + delta_rBSp) - (posm + delta_rBSm))/(dSec*2);
    
                sourceChanged_D0.s(isat).vx_crf(j) = vel_D0(1);
                sourceChanged_D0.s(isat).vy_crf(j) = vel_D0(2);
                sourceChanged_D0.s(isat).vz_crf(j) = vel_D0(3);
    
                sourceChanged_Y0.s(isat).vx_crf(j) = vel_Y0(1);
                sourceChanged_Y0.s(isat).vy_crf(j) = vel_Y0(2);
                sourceChanged_Y0.s(isat).vz_crf(j) = vel_Y0(3);
    
                sourceChanged_B0.s(isat).vx_crf(j) = vel_B0(1);
                sourceChanged_B0.s(isat).vy_crf(j) = vel_B0(2);
                sourceChanged_B0.s(isat).vz_crf(j) = vel_B0(3);
    
                sourceChanged_DC.s(isat).vx_crf(j) = vel_DC(1);
                sourceChanged_DC.s(isat).vy_crf(j) = vel_DC(2);
                sourceChanged_DC.s(isat).vz_crf(j) = vel_DC(3);
    
                sourceChanged_YC.s(isat).vx_crf(j) = vel_YC(1);
                sourceChanged_YC.s(isat).vy_crf(j) = vel_YC(2);
                sourceChanged_YC.s(isat).vz_crf(j) = vel_YC(3);
    
                sourceChanged_BC.s(isat).vx_crf(j) = vel_BC(1);
                sourceChanged_BC.s(isat).vy_crf(j) = vel_BC(2);
                sourceChanged_BC.s(isat).vz_crf(j) = vel_BC(3);
    
                sourceChanged_DS.s(isat).vx_crf(j) = vel_DS(1);
                sourceChanged_DS.s(isat).vy_crf(j) = vel_DS(2);
                sourceChanged_DS.s(isat).vz_crf(j) = vel_DS(3);
    
                sourceChanged_YS.s(isat).vx_crf(j) = vel_YS(1);
                sourceChanged_YS.s(isat).vy_crf(j) = vel_YS(2);
                sourceChanged_YS.s(isat).vz_crf(j) = vel_YS(3);
    
                sourceChanged_BS.s(isat).vx_crf(j) = vel_BS(1);
                sourceChanged_BS.s(isat).vy_crf(j) = vel_BS(2);
                sourceChanged_BS.s(isat).vz_crf(j) = vel_BS(3);
            end    
        end
    end

    [sourceChanged_D0] = getTRF_PosVelSatellites(sourceChanged_D0, T2C_s);
    [sourceChanged_Y0] = getTRF_PosVelSatellites(sourceChanged_Y0, T2C_s);
    [sourceChanged_B0] = getTRF_PosVelSatellites(sourceChanged_B0, T2C_s);
    [sourceChanged_DC] = getTRF_PosVelSatellites(sourceChanged_DC, T2C_s);
    [sourceChanged_YC] = getTRF_PosVelSatellites(sourceChanged_YC, T2C_s);
    [sourceChanged_BC] = getTRF_PosVelSatellites(sourceChanged_BC, T2C_s);
    [sourceChanged_DS] = getTRF_PosVelSatellites(sourceChanged_DS, T2C_s);
    [sourceChanged_YS] = getTRF_PosVelSatellites(sourceChanged_YS, T2C_s);
    [sourceChanged_BS] = getTRF_PosVelSatellites(sourceChanged_BS, T2C_s);

    sourcesChanged_list = {sourceChanged_D0, sourceChanged_Y0, sourceChanged_B0, sourceChanged_DC, sourceChanged_YC, sourceChanged_BC, sourceChanged_DS, sourceChanged_YS, sourceChanged_BS};
end