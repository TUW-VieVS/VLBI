% update CATALOG files
% to be stored in VieVS/WORK/
% lucia, 25-3-2014
% M. Schartner, 26-4-2016: Minor bug fix
% Lucia, 10-5-18: new server for catalogs, replace 'urlwrite' with
% 'websave'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% files to be updated:
f.antenna  = 1;
f.position = 1;
f.equip    = 1;
f.mask     = 1;
f.modes    = 1;
f.source   = 1;
f.flux     = 1;
f.freq     = 1;
f.rx       = 1;
f.loif     = 1;
f.rec      = 1;
f.hdpos    = 1;
f.tracks   = 1;
f.stations = 1;
% delete Level5 .mat files
f.delL5    = 1;

% +++++++++++++++++++++++++++++

outfolder='../CATALOGS/';
% gsfc='ftp://gemini.gsfc.nasa.gov/pub/sked/catalogs/';
gsfc='https://vlbi.gsfc.nasa.gov/software/sked/catalogs/';

%disp(['download new files from ',gsfc]);
h = waitbar(0,'Please wait...');
i=1;


if f.antenna==1
    filename='antenna.cat';
    url=[gsfc, filename,];
    % download
%     urlwrite(url, [outfolder, filename]);  
    websave([outfolder, filename],url);  
    waitbar(i/15,h,['updated file: ',filename]);
    i=i+1;
end
if f.position==1
    filename='position.cat';
    url=[gsfc, filename,];
    % download
%     urlwrite(url, [outfolder, filename]);
    websave([outfolder, filename],url);  
    waitbar(i/15,h,['updated file: ',filename]);
    i=i+1;
end
if f.equip==1
    filename='equip.cat';
    url=[gsfc, filename,];
    % download
%     urlwrite(url, [outfolder, filename]);
    websave([outfolder, filename],url);  
    waitbar(i/15,h,['updated file: ',filename]);
    i=i+1;
end
if f.mask==1
    filename='mask.cat';
    url=[gsfc, filename,];
    % download
%     urlwrite(url, [outfolder, filename]); 
    websave([outfolder, filename],url);  
    waitbar(i/15,h,['updated file: ',filename]);
    i=i+1;
end
if f.modes==1
    filename='modes.cat';
    url=[gsfc, filename,];
    % download
%     urlwrite(url, [outfolder, filename]); 
    websave([outfolder, filename],url);  
    waitbar(i/15,h,['updated file: ',filename]);
    i=i+1;
end
if f.source==1
    filename='source.cat';
    url=[gsfc, filename,];
    % download
%     urlwrite(url, [outfolder, filename]);  
    websave([outfolder, filename],url);  
    waitbar(i/15,h,['updated file: ',filename]);
    i=i+1;
end
if f.flux==1
    filename='flux.cat';
    url=[gsfc, filename,];
    % download
%     urlwrite(url, [outfolder, filename]);  
    websave([outfolder, filename],url);  
    waitbar(i/15,h,['updated file: ',filename]);
    i=i+1;
end
if f.freq==1
    filename='freq.cat';
    url=[gsfc, filename,];
    % download
%     urlwrite(url, [outfolder, filename]);  
    websave([outfolder, filename],url);  
    waitbar(i/15,h,['updated file: ',filename]);
    i=i+1;
end
if f.rx==1
    filename='rx.cat';
    url=[gsfc, filename,];
    % download
%     urlwrite(url, [outfolder, filename]); 
    websave([outfolder, filename],url);  
    waitbar(i/15,h,['updated file: ',filename]);
    i=i+1;
end
if f.loif==1
    filename='loif.cat';
    url=[gsfc, filename,];
    % download
%     urlwrite(url, [outfolder, filename]);  
    websave([outfolder, filename],url);  
    waitbar(i/15,h,['updated file: ',filename]);
    i=i+1;
end
if f.rec==1
    filename='rec.cat';
    url=[gsfc, filename,];
    % download
%     urlwrite(url, [outfolder, filename]);  
    websave([outfolder, filename],url);  
    waitbar(i/15,h,['updated file: ',filename]);
    i=i+1;
end
if f.hdpos==1
    filename='hdpos.cat';
    url=[gsfc, filename,];
    % download
%     urlwrite(url, [outfolder, filename]);
    websave([outfolder, filename],url);  
    waitbar(i/15,h,['updated file: ',filename]);
    i=i+1;
end
if f.tracks==1
    filename='tracks.cat';
    url=[gsfc, filename,];
    % download
%     urlwrite(url, [outfolder, filename]); 
    websave([outfolder, filename],url);  
    waitbar(i/15,h,['updated file: ',filename]);
    i=i+1;
end
if f.stations==1
    filename='stations.tmp';
    url=[gsfc, filename,];
    % download
%     urlwrite(url, [outfolder, 'stations.cat']);  
    websave([outfolder, filename],url);  
    waitbar(i/15,h,['updated file: ',filename]);
    i=i+1;
end
if f.delL5==1
    if exist('../DATA/LEVEL5/source.mat', 'file')
        delete ../DATA/LEVEL5/source.mat
        waitbar(i/15,h,'source file in ../DATA/LEVEL5/ deleted');
        i=i+1;
    end
end
close(h);
  
  