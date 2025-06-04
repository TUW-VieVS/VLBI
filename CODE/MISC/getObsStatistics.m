clear all
t = "1GalA_2GalBn";
t2 = "1GalA_2GalB";
%folder = ["cVGOS12_" + t];
folder2 = ["cVGOS12_" + t2];
folder = folder2;
path = "\\project5\data-write\HG\hwolf\VieVSSatellites\VLBI\DATA";
%path = "\data-write\HG\hwolf\VieVSSatellites\VLBI\DATA";

v = [10, 20, 30, 40, 50, 60];
%v = [10, 20, 30, 40];
for vi =1:length(v)
    fname = [folder2 + (v(vi)) + "n"];
    %fname = [folder + (v(vi))];

    load([path + "\LEVEL1\cVGOS12_" + t + "\22AUG27VS_" + fname +  "_antenna.mat"])
    load([path + "\LEVEL1\cVGOS12_" + t + "\22AUG27VS_" + fname +  "_scan.mat"])
    n_scan = length(scan);
    na = length(antenna);
    nscanSat = zeros(na, 1); nobsSat =  zeros(na, 1);
    nscanQu = zeros(na, 1); nobsQu =  zeros(na, 1);
    nobs = zeros(na, 1);
    names = string (extractfield(antenna, 'name'))';

    for itim = 1:n_scan 
            i1 = extractfield(scan(itim).obs(:), 'i1');
            i2 = extractfield(scan(itim).obs(:), 'i2');
            i = [i1, i2];
            iu = unique([i1, i2])';
        
            if strcmp (scan(itim).obs_type,"s")
               nscanSat(iu,1) = nscanSat(iu,1) + 1;
            else
               nscanQu(iu,1) = nscanQu(iu,1) + 1;
            end
    
            for k=1:length(i)
                if strcmp (scan(itim).obs_type,"s")
                   nobsSat(i(k),1) = nobsSat(i(k),1) + 1;
                else
                   nobsQu(i(k),1) = nobsQu(i(k),1) + 1;
                end
            end
    end   
    nobsTot = nobsSat + nobsQu;
    nscanTot = nscanSat + nscanQu;

    pobsSat = nobsSat.*100./nobsTot;
    pscanSat = nscanSat.*100./nscanTot;

    pdiffobsSat = pobsSat - v(vi);
    pdiffscanSat = pscanSat - v(vi);

    [names, order] = sort(names);
    nscanSat = nscanSat(order,:);
    nscanQu = nscanQu(order,:);
    nscanTot = nscanTot(order,:);
    nobsSat = nobsSat(order,:);
    nobsQu = nobsQu(order,:);
    nobsTot = nobsTot(order,:);
    pscanSat = pscanSat(order,:);
    pobsSat = pobsSat(order,:);
    pdiffscanSat = pdiffscanSat(order,:);
    pdiffobsSat = pdiffobsSat(order,:);
    

    clear E
    E = table(names, nscanSat, nscanQu, nscanTot,  nobsSat, nobsQu, nobsTot, pscanSat, pobsSat, pdiffscanSat, pdiffobsSat);
    E.Properties.VariableNames = {'Name', '#scanSat','#scanQu','#scanTot', '#obsSat','#obsQu','#obsTot', 'pscanSat', 'pobsSat', 'pdiffscanSat', 'pdiffobsSat'};
    writetable(E , [path + "\LEVEL3\" + folder + "\" + fname + "_ObsStatistics" + ".xlsx"]);
    writetable(E , [path + "\LEVEL3\" + folder + "\" + fname + "_ObsStatistics" + ".txt"], 'Delimiter', 'tab');
    fprintf('\n statistics\n')
    disp(E)
    clear E

end
