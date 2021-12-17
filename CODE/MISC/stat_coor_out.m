% Hana Krasna
% compute total XYZ and dHEN w.r.t. a priori from LEVEL3 data

function [ ] = stat_coor_out( folder )

    %folder = runp.lsm_path;
    pathl3 = ['../DATA/LEVEL3/' folder '/'];
    files = dir([pathl3 'x_*.mat']);

    
    fileNames = {files.name};
    
    fid=fopen(['../OUT/XYZ_dHEN_' folder '.txt'],'wt');
    fprintf(fid,'#col 1: station \n');
    fprintf(fid,'#col 2: mjd of the session  \n');
    fprintf(fid,'#col 3-5: a priori XYZ [m]\n');
    fprintf(fid,'#col 6-8: estimated XYZ [m] \n');
    fprintf(fid,'#col 9-11: formal errors of estimated XYZ [m] \n');
    fprintf(fid,'#col 12-14: dHEN w.r.t. a priori [m] \n'); 
    fprintf(fid,'#col 15-17: formal errors of dHEN [m] \n');
    fprintf(fid,'#col 18: session \n');
    fprintf(fid,'#\n');
    
    for i = 1:length(files)
        load([pathl3 fileNames{i}])
        load([pathl3 fileNames{i}(3:end-4) '_antenna.mat'])
        load([pathl3 fileNames{i}(3:end-4) '_parameter.mat'])
        
        mjda = [x_.coorx.mjd];
        mjd=mjda(1);
        
        % check for fixed stations
        clear idest
        for j=1:length(x_.coorx)
            idest(j) = ~isempty(x_.coorx(j).mjd);
        end
        
        cpsd_all = cPostSeismDeform(mjd,antenna);
          
        % calculate apriori values
        cpsdX=[]; cpsdY=[]; cpsdZ=[];aprX=[];aprY=[];aprZ=[];
        cpsdX(1,:)=cpsd_all(1,1,idest);
        aprX(:,1)=[antenna(idest).x]+[antenna(idest).vx].*(mjd-[antenna(idest).epoch])./365.25 + cpsdX; %m
        cpsdY(1,:)=cpsd_all(2,1,idest);
        aprY(:,1)=[antenna(idest).y]+[antenna(idest).vy].*(mjd-[antenna(idest).epoch])./365.25 + cpsdY; %m
        cpsdZ(1,:)=cpsd_all(3,1,idest);
        aprZ(:,1)=[antenna(idest).z]+[antenna(idest).vz].*(mjd-[antenna(idest).epoch])./365.25 + cpsdZ; %m                 




        totXYZ=[];
       % calculate total estimated values
        totXYZ(:,1)=aprX + [x_.coorx.val]'./100;
        totXYZ(:,2)=aprY + [x_.coory.val]'./100;
        totXYZ(:,3)=aprZ + [x_.coorz.val]'./100;

        
       % formal errors
       mxyzT = [[x_.coorx.mx]./100; [x_.coory.mx]./100; [x_.coorz.mx]./100];
       mxyz=mxyzT';
        
       
        % transform xyz->hen
        [lat,lon,~]=xyz2ell(totXYZ);
        dxyz = totXYZ-[aprX,aprY,aprZ];
        [dhen]=xyz2ren(dxyz,lat,lon);
        [mhen,~]=xyz2ren_sigma(mxyz,lat,lon);
       
       
       
       anames = {antenna(idest).name};
       for j = 1:length(antenna(idest))
            fprintf(fid,'%s  %8.0f %15.4f %15.4f %15.4f       %15.4f %15.4f %15.4f   %10.4f %10.4f %10.4f    %10.4f %10.4f %10.4f   %10.4f %10.4f %10.4f   %s\n', anames{j}, mjd, aprX(j), aprY(j), aprZ(j),...
                totXYZ(j,:), mxyz(j,:),...
                dhen(j,:), mhen(j,:), parameter.session_name);
       end
    end
    fclose(fid);
end