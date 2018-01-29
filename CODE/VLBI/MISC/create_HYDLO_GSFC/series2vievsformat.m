% Script which converts the station-wise .txt files 
% (hydrology loading displacement)to the yearly matlab structures used in VieVS
% It is recommended to use the detrended time series.
% Created on March 17, 2017 by Hana Krasna



clear all
clc
pth1= 'hydlo_detrend/';

DELIMITER = ' ';
HEADERLINES = 7



fid = fopen(['stations_solve.txt']);
k = 0;
while ~feof(fid)
    line = fgetl(fid);
    k = k + 1;
    statlist(k,1:8) = line(3:10);
end
fclose(fid);


% chech which files we have
files_all = dir([ pth1 '*.*']);
files_all([1 2])=[];
lfi = length(files_all);


for y = 1979:2016
    
    % function modjuldat converts year to mjd
    mjd = modjuldat(y);
    mjd1 = mjd-60; % 60 days before beginning of the year (1.Dec)
    mjd2 = mjd+425; % 60 days after end of the year (1.Jan)

    % loop for stations from statlist
    for ist=1:size(statlist,1)
        A=[];
        hydlo(ist).ivsname = statlist(ist,:);
        
        % get rid off of the blank spaces for comparison
        stcel=cellstr(statlist(ist,:));
        aname=char(stcel);
        
        % loop for stations with the corrections
        for j=1:length(files_all)
            namecorr1=files_all(j).name(1:end-9);
            
            % instead of a blank space write '_'
            idblank=find(namecorr1(1:max(find(namecorr1(:)~=' ')))==' ');
            namecorr=namecorr1;
            namecorr(idblank)='_';

            % check if we have corrections for the station from statlist
            if strcmp(namecorr,aname)
                A = importdata([pth1 files_all(j).name],DELIMITER,HEADERLINES);
            end
        end
        
        if ~isempty(A)
            id=find(A.data(:,1)>mjd1 & A.data(:,1)<mjd2);
            hydlo(ist).mjd = A.data(id,:);
        end
        clear aname stcel A id

    end
    
    
    % delete stations taken from the statlist, for which we don't have
    % corrections
    k=[];
    for i=1:length(hydlo)
        if isempty(hydlo(i).mjd)
            k=[k i];
        end
    end
    hydlo(k)=[];
   

    save([num2str(y),'_CMTE_HYDLO'],'hydlo');
    
    clear hydlo
    
end
    

