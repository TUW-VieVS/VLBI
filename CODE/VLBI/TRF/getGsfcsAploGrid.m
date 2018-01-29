% stupid long - writes a txt file (.dat) which is then read to get
% values... but anyway:
% This function reads the netCDF file consisting of atmo pressure loading
% values from Leonid Petrovs model. The output are grids (lon,lat) and the
% corresponding data values (in the same grid).

function [mesh_lon, mesh_lat, cos_s1_dr, sin_s1_dr, cos_s2_dr, ...
    sin_s2_dr, cos_s1_dn, sin_s1_dn, cos_s2_dn, sin_s2_dn, cos_s1_de, ...
    sin_s1_de, cos_s2_de, sin_s2_de]=getGsfcsAploGrid(aploNetcdfFile)
%clc
%aploNetcdfFile='../neededFiles/aplo_s1_s2_noib_1.0x1.0deg.nc';
%ncid=netcdf.open(aploNetcdfFile);
S = netcdf(aploNetcdfFile);

d = S.VarArray(1,6).Data;

% preallocateing and creating running index
sin_s1_dr=zeros(181,360); cos_s1_dr=sin_s1_dr; sin_s2_dr=sin_s1_dr; 
cos_s2_dr=sin_s1_dr; sin_s1_de=sin_s1_dr; cos_s1_de=sin_s1_dr;
sin_s2_de=sin_s1_dr; cos_s2_de=sin_s1_dr; sin_s1_dn=sin_s1_dr; 
cos_s1_dn=sin_s1_dr; sin_s2_dn=sin_s1_dr; cos_s2_dn=sin_s1_dr; 
k = 0;
for i = 1:360
    for j = 1:181
        k = k + 1;
        sin_s1_dr(182-j,i) = double(d(j,i,1,2,1))/100;
        cos_s1_dr(182-j,i) = double(d(j,i,1,1,1))/100;
        sin_s2_dr(182-j,i) = double(d(j,i,2,2,1))/100;
        cos_s2_dr(182-j,i) = double(d(j,i,2,1,1))/100;
        sin_s1_de(182-j,i) = double(d(j,i,1,2,2))/100;
        cos_s1_de(182-j,i) = double(d(j,i,1,1,2))/100;
        sin_s2_de(182-j,i) = double(d(j,i,2,2,2))/100;
        cos_s2_de(182-j,i) = double(d(j,i,2,1,2))/100;
        sin_s1_dn(182-j,i) = double(d(j,i,1,2,3))/100;
        cos_s1_dn(182-j,i) = double(d(j,i,1,1,3))/100;
        sin_s2_dn(182-j,i) = double(d(j,i,2,2,3))/100;
        cos_s2_dn(182-j,i) = double(d(j,i,2,1,3))/100;
    end
end

%fid = fopen('s1_s2_def_cm_noib_leonid.dat','w');

% preallocating and create starting index
s1_s2=zeros(181*360+181, 14);
rowInd=1;

for i = 1:360
    for j = 1:181
%         fprintf(fid,'%8.3f %8.3f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n', i-1,91-j, ...
%         -cos_s1_dr(j,i),-sin_s1_dr(j,i),cos_s2_dr(j,i),sin_s2_dr(j,i), ...
%         -cos_s1_dn(j,i),-sin_s1_dn(j,i),cos_s2_dn(j,i),sin_s2_dn(j,i), ...
%         -cos_s1_de(j,i),-sin_s1_de(j,i),cos_s2_de(j,i),sin_s2_de(j,i));
        s1_s2(rowInd,:)=[i-1,91-j, ...
            -cos_s1_dr(j,i),-sin_s1_dr(j,i),cos_s2_dr(j,i),sin_s2_dr(j,i), ...
            -cos_s1_dn(j,i),-sin_s1_dn(j,i),cos_s2_dn(j,i),sin_s2_dn(j,i), ...
            -cos_s1_de(j,i),-sin_s1_de(j,i),cos_s2_de(j,i),sin_s2_de(j,i)];
        rowInd=rowInd+1;
    end
end
for i = 1
    for j = 1:181
%         fprintf(fid,'%8.3f %8.3f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n', 360,91-j, ...
%             -cos_s1_dr(j,i),-sin_s1_dr(j,i),cos_s2_dr(j,i),sin_s2_dr(j,i), ...
%             -cos_s1_dn(j,i),-sin_s1_dn(j,i),cos_s2_dn(j,i),sin_s2_dn(j,i), ...
%             -cos_s1_de(j,i),-sin_s1_de(j,i),cos_s2_de(j,i),sin_s2_de(j,i));
        s1_s2(rowInd,:)=[360,91-j, ...
            -cos_s1_dr(j,i),-sin_s1_dr(j,i),cos_s2_dr(j,i),sin_s2_dr(j,i), ...
            -cos_s1_dn(j,i),-sin_s1_dn(j,i),cos_s2_dn(j,i),sin_s2_dn(j,i), ...
            -cos_s1_de(j,i),-sin_s1_de(j,i),cos_s2_de(j,i),sin_s2_de(j,i)];
        rowInd=rowInd+1;
    end
end

%fclose(fid);

%% My try of Step 2
%s1_s2 = load('s1_s2_def_cm_noib_leonid.dat');

% get all longitudes and latitudes
allLons=unique(s1_s2(:,1));
allLats=unique(s1_s2(:,2));

[mesh_lon, mesh_lat]=meshgrid(allLons, allLats);

% preallocate data matrices
emptyMesh=zeros(size(mesh_lon));
cos_s1_dr=emptyMesh;
sin_s1_dr = emptyMesh;
cos_s2_dr = emptyMesh;
sin_s2_dr = emptyMesh;
cos_s1_dn = emptyMesh;
sin_s1_dn = emptyMesh;
cos_s2_dn = emptyMesh;
sin_s2_dn = emptyMesh;
cos_s1_de = emptyMesh;
sin_s1_de = emptyMesh;
cos_s2_de = emptyMesh;
sin_s2_de = emptyMesh;

% fill these empty matrices with data
cos_s1_dr(:)=s1_s2(:,3);
sin_s1_dr(:)=s1_s2(:,4);
cos_s2_dr(:)=s1_s2(:,5);
sin_s2_dr(:)=s1_s2(:,6);
cos_s1_dn(:)=s1_s2(:,7);
sin_s1_dn(:)=s1_s2(:,8);
cos_s2_dn(:)=s1_s2(:,9);
sin_s2_dn(:)=s1_s2(:,10);
cos_s1_de(:)=s1_s2(:,11);
sin_s1_de(:)=s1_s2(:,12);
cos_s2_de(:)=s1_s2(:,13);
sin_s2_de(:)=s1_s2(:,14);


% 
% % interpolate to 'our' coordinates
% %% Step 2 - load the data file and
% s1_s2 = load('s1_s2_def_cm_noib_leonid.dat');
% 
% outname = 's12_cm_noib_leonid.dat';
% 
% klm = size(s1_s2,1);
% 
% tic
% 
% dlatdeg = zeros(361,181);
% dlondeg = zeros(361,181);
% cos_s1_dr = zeros(361,181);
% sin_s1_dr = zeros(361,181);
% cos_s2_dr = zeros(361,181);
% sin_s2_dr = zeros(361,181);
% cos_s1_dn = zeros(361,181);
% sin_s1_dn = zeros(361,181);
% cos_s2_dn = zeros(361,181);
% sin_s2_dn = zeros(361,181);
% cos_s1_de = zeros(361,181);
% sin_s1_de = zeros(361,181);
% cos_s2_de = zeros(361,181);
% sin_s2_de = zeros(361,181);
% 
% k = 0;
% for ilon = 1:361
%     for ilat = 1:181
%         k = k + 1;
%         dlondeg(ilon,ilat) = s1_s2(k,1);
%         dlatdeg(ilon,ilat) = s1_s2(k,2);
%         cos_s1_dr(ilon,ilat) = s1_s2(k,3);
%         sin_s1_dr(ilon,ilat) = s1_s2(k,4);
%         cos_s2_dr(ilon,ilat) = s1_s2(k,5);
%         sin_s2_dr(ilon,ilat) = s1_s2(k,6);
%         cos_s1_dn(ilon,ilat) = s1_s2(k,7);
%         sin_s1_dn(ilon,ilat) = s1_s2(k,8);
%         cos_s2_dn(ilon,ilat) = s1_s2(k,9);
%         sin_s2_dn(ilon,ilat) = s1_s2(k,10);
%         cos_s1_de(ilon,ilat) = s1_s2(k,11);
%         sin_s1_de(ilon,ilat) = s1_s2(k,12);
%         cos_s2_de(ilon,ilat) = s1_s2(k,13);
%         sin_s2_de(ilon,ilat) = s1_s2(k,14);
%         
%     end
% end
% toc
% 
% 
% %% Interpolate to "our" coordinates
% keyboard; % i don't have that file...
% % 'Parameter des GRS80'
% a = 6378137;
% b = 6356752.3141;
% 
% fid = fopen('statlist.dat');
% k = 0;
% while ~feof(fid)
%     
%     line = fgetl(fid);
%     k = k + 1;
%     stat8(k,1:8) = line(1:8);
%     P = sscanf(line(9:length(line)),'%f',[3 1]);
%     [phi,lam,H] = kart2ell(P,a,b);
%     lat = phi*180/pi;
%     lon = lam*180/pi;
%     if lon < 0
%         lon = lon + 360;
%     end
%     dlat(k) = lat;
%     dlon(k) = lon;
%     hght(k) = H;
% 
% end
% 
% fclose(fid);
% 
% cos_s1_dr1 = griddata(dlondeg,dlatdeg,cos_s1_dr,dlon,dlat);
% sin_s1_dr1 = griddata(dlondeg,dlatdeg,sin_s1_dr,dlon,dlat);
% cos_s2_dr1 = griddata(dlondeg,dlatdeg,cos_s2_dr,dlon,dlat);
% sin_s2_dr1 = griddata(dlondeg,dlatdeg,sin_s2_dr,dlon,dlat);
% cos_s1_dn1 = griddata(dlondeg,dlatdeg,cos_s1_dn,dlon,dlat);
% sin_s1_dn1 = griddata(dlondeg,dlatdeg,sin_s1_dn,dlon,dlat);
% cos_s2_dn1 = griddata(dlondeg,dlatdeg,cos_s2_dn,dlon,dlat);
% sin_s2_dn1 = griddata(dlondeg,dlatdeg,sin_s2_dn,dlon,dlat);
% cos_s1_de1 = griddata(dlondeg,dlatdeg,cos_s1_de,dlon,dlat);
% sin_s1_de1 = griddata(dlondeg,dlatdeg,sin_s1_de,dlon,dlat);
% cos_s2_de1 = griddata(dlondeg,dlatdeg,cos_s2_de,dlon,dlat);
% sin_s2_de1 = griddata(dlondeg,dlatdeg,sin_s2_de,dlon,dlat);
% 
% fid = fopen(outname,'w');
% for i = 1:length(dlon)
%     fprintf(fid,'%s%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f \n', ...
%         stat8(i,1:8),dlat(i),dlon(i),hght(i), ...
%         cos_s1_dr1(i),sin_s1_dr1(i),cos_s2_dr1(i),sin_s2_dr1(i), ...
%         cos_s1_dn1(i),sin_s1_dn1(i),cos_s2_dn1(i),sin_s2_dn1(i), ...
%         cos_s1_de1(i),sin_s1_de1(i),cos_s2_de1(i),sin_s2_de1(i));
% end
% fclose(fid);
% 
% %% Step 3
% 
% close all;
% clear all;
% 
% fid = fopen('s12_cm_noib_leonid.dat');
% i=0;
% while ~feof(fid)
%     line = [fgetl(fid),'                                                             '];
%     if line(1) ~= '#'
%         i = 1+i;
%         
%         atide(i).ivsname=line(1:8); 
%         atide(i).mat = str2num(line(9:length(line)));
%     end
% end 
% fclose(fid);
% 
% save s12_cm_noib_leonid atide
% 
% 
% 
% 
