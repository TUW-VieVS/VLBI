% Script which removes offset and rate from the hydrology loading displacements provided 
% by the NASA GSFC VLBI group which are avaliable at http://lacerta.gsfc.nasa.gov/hydlo
% (Eriksson & MacMillan). The displacement series are computed from the
% GLDAS NOAH model provided by the NASA GSFC GLDAS team" 

% Created by Hana Krasna, March 17, 2017


clear all
clc
curDate=clock;          % current date and time

fid = fopen(['stations_solve.txt']);
k = 0;
while ~feof(fid)
    line = fgetl(fid);
    k = k + 1;
    statlist(k,1:8) = line(3:10);
end
fclose(fid);


% chech which files we have (Download from http://lacerta.gsfc.nasa.gov/hydlo)
files_all = dir(['hydlo/loadingfiles/cmte_series/*.*']);
files_all([1 2])=[];
lfi = length(files_all);


for i=1:length(files_all)
    
    up_detrend=[]; e_detrend=[]; n_detrend=[];
    
    fileID = fopen([files_all(i).folder '/' files_all(i).name]);
    GSFC = textscan(fileID,'%f %f %f %f','CommentStyle','#');
    fclose(fileID);
    
    % up
    ii=2;
    l=GSFC{ii};
    A=zeros(length(GSFC{ii}),2);
    A(:,1)=1;
    A(:,2)=(GSFC{1}-51544)/365.25;
    
    N=A'*A;
    b=A'*l;
    x=inv(N)*b;
    
    up_detrend= GSFC{ii} - (x(1)+ ((GSFC{1}-51544)/365.25)*x(2));
    
    
    figure(1)
    subplot(3,1,ii-1)
    plot((GSFC{1}-51544)/365,GSFC{ii},'b')
    hold on
    plot((GSFC{1}-51544)/365,x(1)+ ((GSFC{1}-51544)/365.25)*x(2),'g')
    hold on
    plot((GSFC{1}-51544)/365,up_detrend,'r')
    hold off
    ylabel('up [m]' )
    title( files_all(i).name(1:end-9))

    
    % e
    ii=3;
    l=GSFC{ii};
    A=zeros(length(GSFC{ii}),2);
    A(:,1)=1;
    A(:,2)=(GSFC{1}-51544)/365.25;
    
    N=A'*A;
    b=A'*l;
    x=inv(N)*b;
    
    e_detrend= GSFC{ii} - (x(1)+ ((GSFC{1}-51544)/365.25)*x(2));
    
    
    subplot(3,1,ii-1)
    plot((GSFC{1}-51544)/365,GSFC{ii},'b')
    hold on
    plot((GSFC{1}-51544)/365,x(1)+ ((GSFC{1}-51544)/365.25)*x(2),'g')
    hold on
    plot((GSFC{1}-51544)/365,e_detrend,'r')
    hold off
    ylabel('e [m]' )
    
    
    % n
    ii=4;
    l=GSFC{ii};
    A=zeros(length(GSFC{ii}),2);
    A(:,1)=1;
    A(:,2)=(GSFC{1}-51544)/365.25;
    
    N=A'*A;
    b=A'*l;
    x=inv(N)*b;
    
    n_detrend= GSFC{ii} - (x(1)+ ((GSFC{1}-51544)/365.25)*x(2));
    
    
    subplot(3,1,ii-1)
    plot((GSFC{1}-51544)/365,GSFC{ii},'b')
    hold on
    plot((GSFC{1}-51544)/365,x(1)+ ((GSFC{1}-51544)/365.25)*x(2),'g')
    hold on
    plot((GSFC{1}-51544)/365,n_detrend,'r')
    hold off
    ylabel('n [m]' )
    xlabel('years since 2000')

    % create txt file
    fid=fopen(['hydlo_detrend/' files_all(i).name],'wt');
    
    
    fprintf(fid, '# Hydrology loading displacements computed from the monthly GLDAS NOAH model\n');
    fprintf(fid, '# This is a product of the Hydrology loading service: http://lacerta.gsfc.nasa.gov/hydlo/\n');
    fprintf(fid, '# Displacements are with respect to the center of mass of the total earth\n');
    fprintf(fid, '#\n');
    fprintf(fid, '# Created for VieVS on %02.0f.%02.0f.%04.0f at %02.0f:%02.0f:%02.0f local time\n', curDate(3), curDate(2), curDate(1), curDate(4), curDate(5), curDate(6));
    fprintf(fid, '# OFFSET + TREND REMOVED!!! \n');
    fprintf(fid, '# mjd up (meters) eastwest (meters) northsouth (meters)\n');
    
    mjduen=[GSFC{1} up_detrend e_detrend n_detrend];
    fprintf(fid,'%5.0f  %10.5f  %10.5f  %10.5f\n',mjduen');
    fclose(fid);
    
%     print('-dpng','-r600',[pth 'hydlo_detrend/plots/' files_all(i).name(1:end-4) '.png'])

end