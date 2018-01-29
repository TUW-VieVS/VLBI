% saves statistics in excel file
%
% created: Matthias Schartner
% Changelog: 2017-05-29: M. Schartner: Bugfix 

function [ ] = writeLEVEL3toXLSX( runp )

    folder = runp.lsm_path;
    path = ['../DATA/LEVEL3/' folder '/'];
    files = dir([path 'x_*.mat']);
    path_xlsx = ['../DATA/SCHED/' folder '/'];

    saveStatAsXLSX( [],path_xlsx,[],0,1 )
    fileNames = {files.name};

    for isched = 1:str2double(files(end).name(14:16))
        idxsim = strfind(fileNames,sprintf('V%03d',isched));
        idxsim = cellfun(@isempty,idxsim);
        idxsim = ~idxsim;

        thisFiles = files(idxsim);
        c = struct('coorx',[],'coorxm',[],...
                      'coory',[],'coorym',[],...
                      'coorz',[],'coorzm',[],...
                      'coor3d',[],'coor3dm',[],...
                      'xpol',[],'xpolm',[],...
                      'ypol',[],'ypolm',[],...
                      'dut1',[],'dut1m',[],...
                      'nutdx',[],'nutdxm',[],...
                      'nutdy',[],'nutdym',[],...
                      'soura',[],'souram',[],...
                      'sourainNNR',[],'sourainNNRm',[],...
                      'souranotinNNR',[],'souranotinNNRm',[],...
                      'soude',[],'soudem',[],...
                      'soudeinNNR',[],'soudeinNNRm',[],...
                      'soudenotinNNR',[],'soudenotinNNRm',[]);

    %     allNames = {files.name};
        for isim = 1:length(thisFiles)

            line = str2double(thisFiles(isim).name(14:16));
            load([path thisFiles(isim).name])

            % Station coordinates
            c.coorx(isim,:) = [x_.coorx.val];
            c.coorxm(isim,:) = [x_.coorx.mx];

            c.coory(isim,:) = [x_.coory.val];
            c.coorym(isim,:) = [x_.coory.mx];

            c.coorz(isim,:) = [x_.coorz.val];
            c.coorzm(isim,:) = [x_.coorz.mx];

            c.coor3d(isim,:) = sqrt([x_.coorx.val].^2+[x_.coory.val].^2+[x_.coorz.val].^2);
            c.coor3dm(isim,:) = sqrt([x_.coorx.mx].^2+[x_.coory.mx].^2+[x_.coorz.mx].^2); % Fehlerfortpflanzung

            % EOPs
            c.xpol(isim,:) = [x_.xpol.val];
            c.xpolm(isim,:) = [x_.xpol.mx];

            c.ypol(isim,:) = [x_.ypol.val];
            c.ypolm(isim,:) = [x_.ypol.mx];

            c.dut1(isim,:) = [x_.dut1.val];
            c.dut1m(isim,:) = [x_.dut1.mx];

            c.nutdx(isim,:) = [x_.nutdx.mx];
            c.nutdxm(isim,:) = [x_.nutdx.val];

            c.nutdy(isim,:) = [x_.nutdx.val];
            c.nutdym(isim,:) = [x_.nutdx.mx];

            % Source coordinates
            souinNNR = logical([x_.soura.inNNR]);

            soude = [x_.soude.val];
            soura = [x_.soura.val]; % cos(soude)

            souram = [x_.soura.mx];
            soudem = [x_.soude.mx];

            if ~isempty(soura)
                c.soura(isim,:) = soura;
                c.souram(isim,:) = souram;

                c.soude(isim,:) = soude;
                c.soudem(isim,:) = soudem;

                c.sou2d(isim,:) =  sqrt(soura.^2+soude.^2);
                c.sou2dm(isim,:) = sqrt(souram.^2+soudem.^2); % Fehlerfortpflanzung

                c.sourainNNR(isim,:) = soura(souinNNR);
                c.sourainNNRm(isim,:) = souram(souinNNR);

                c.soudeinNNR(isim,:) = soude(souinNNR);
                c.soudeinNNRm(isim,:) = soudem(souinNNR);

                c.souranotinNNR(isim,:) = soura(~souinNNR);
                c.souranotinNNRm(isim,:) = souram(~souinNNR);

                c.soudenotinNNR(isim,:) = soude(~souinNNR);
                c.soudenotinNNRm(isim,:) = soudem(~souinNNR);
            else
                c.soura(isim,:) = NaN;
                c.souram(isim,:) = NaN;

                c.soude(isim,:) = NaN;
                c.soudem(isim,:) = NaN;

                c.sou2d(isim,:) =  NaN;
                c.sou2dm(isim,:) = NaN;

                c.sourainNNR(isim,:) = NaN;
                c.sourainNNRm(isim,:) = NaN;

                c.soudeinNNR(isim,:) = NaN;
                c.soudeinNNRm(isim,:) = NaN;

                c.souranotinNNR(isim,:) = NaN;
                c.souranotinNNRm(isim,:) = NaN;

                c.soudenotinNNR(isim,:) = NaN;
                c.soudenotinNNRm(isim,:) = NaN;
            end
        end

        % Station coordinates
        stat.coorx = mean(std(c.coorx));
        stat.coorxm = mean(mean(c.coorxm));

        stat.coory = mean(std(c.coory));
        stat.coorym = mean(mean(c.coorym));

        stat.coorz = mean(std(c.coorz));
        stat.coorzm = mean(mean(c.coorzm));

        stat.coor3d = mean(std(c.coor3d));
        stat.coor3dm = mean(mean(c.coor3dm));

        % EOP
        stat.xpol = mean(std(c.xpol));
        stat.xpolm = mean(mean(c.xpolm));

        stat.ypol = mean(std(c.ypol));
        stat.ypolm = mean(mean(c.ypolm));

        stat.dut1 = mean(std(c.dut1));
        stat.dut1m = mean(mean(c.dut1m));

        stat.nutdx = mean(std(c.nutdx));
        stat.nutdxm = mean(mean(c.nutdxm));

        stat.nutdy = mean(std(c.nutdy));
        stat.nutdym = mean(mean(c.nutdym));

        % Source coordinates
        stat.soura = mean(std(c.soura));
        stat.souram = mean(mean(c.souram));

        stat.soude = mean(std(c.soude));
        stat.soudem = mean(mean(c.soudem));

        stat.sou2d = mean(std(c.sou2d));
        stat.sou2dm = mean(mean(c.sou2dm));

        stat.sourainNNR = mean(std(c.sourainNNR));
        stat.sourainNNRm = mean(mean(c.sourainNNRm));

        stat.soudeinNNR = mean(std(c.soudeinNNR));
        stat.soudeinNNRm = mean(mean(c.soudeinNNRm));

        stat.souranotinNNR = mean(std(c.souranotinNNR));
        stat.souranotinNNRm = mean(mean(c.souranotinNNRm));

        stat.soudenotinNNR = mean(std(c.soudenotinNNR));
        stat.soudenotinNNRm = mean(mean(c.soudenotinNNRm));

        stat.isim = length(thisFiles);
        stat.ista = length([x_.coorx.val]);
        stat.isrc = length([x_.soude.val]);

        saveStatAsXLSX( stat,path_xlsx,[],line,1)

    end
end
