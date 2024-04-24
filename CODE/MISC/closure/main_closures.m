% Script for group delay closure computation
% Written by Matthias Schartner
% - baseline delays only
% - order of stations: ab, ac, bc
% Updated by Hana Krasna
% - updated for the order of stations: ab, ac, cb
% - geocentric delays added
%
% typ: 1 - baseline delay, 2 - geocentric delay

function main_closures(process_list,typ)

process_list=char(process_list);

format compact;
% clc, clear, close all;

%% generate list of closure delays from CONT17 XA sessions:
% XA = dir(fullfile('Data','17*XA'))
% 
% %for i = 1:length(XA) % <-- run over all XA sessions
% for i = 1:1 % <-- run only one XA session
%     fprintf('Read: %s \n',XA(i).name);
%     scans = getScanList(fullfile(XA(i).folder,XA(i).name));
%     fprintf('Find closures: %s \n',XA(i).name);
%     scans2closures(XA(i).name, scans);
% end


% % sessions from process list
pth = '../DATA/vgosDB/';
pthOUT = '../OUT/CLOSURES/';
% process_list=['2020/20OCT29PI [vgosDB]'
%               '2020/20OCT29VE [vgosDB]'];
% typ = 1; % 1 - baseline delay, 2 - geocentric delay



%% uncompress files vgosDB *.tar.gz or *.tgz file
for i = 1:size(process_list,1)
%     curNcFolder = [pth ,process_list(i,1:14),'/'];
    
    session_pl = process_list(i,:);
    session_name  = session_pl(1 : (strfind(session_pl, ' [vgosDB]')-1));

    curNcFolder = [pth ,session_name '/'];

    % uncompress vgosDB *.tar.gz or *.tgz file
    wasCompressed = false;
    if ~isfolder(curNcFolder)
        wasCompressed = true;

        % uncompress vgosDB *.tar.gz or *.tgz file
        vgosTargz = [curNcFolder(1:end-1),'.tar.gz'];
        vgosTgz = [curNcFolder(1:end-1),'.tgz'];
        curSlash = sort([strfind(curNcFolder,'/'), strfind(curNcFolder,'\')]);
        vgosTgzFolder = curNcFolder(1:curSlash(end-1));

        if exist(vgosTgz,'file')
            untar(vgosTgz,vgosTgzFolder);
        elseif exist(vgosTargz,'file')
            untar(vgosTargz,vgosTgzFolder);
        else
            fprintf('ERROR: %s does not exist!\n',vgosTgz);
        end       
    end

    PATH_in = [pth session_name];   
    fprintf('Read: %s \n',PATH_in);    
    scans = getScanList(PATH_in, typ);

    fprintf('Find closures: %s \n',PATH_in);
%     scans2closures(pthOUT, process_list(i,6:14), scans, typ);
    scans2closures(pthOUT,session_name(6:end), scans, typ);
    if wasCompressed
        rmdir(PATH_in, 's')
    end
end

