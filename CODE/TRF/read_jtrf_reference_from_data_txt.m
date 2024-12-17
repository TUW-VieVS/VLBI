

function [ns_codes] = read_jtrf_reference_from_data_txt(ns_codes,datafile)

infls = dir(datafile);
for i = 3:length(infls)
    idcdp = strcmp(cell({ns_codes.CDP}),{infls(i).name(50:53)});
    iddom = strcmp(cell({ns_codes.domes}),{infls(i).name(40:48)});
    idnsc = iddom & idcdp;
    if sum(idnsc)>1
        if sum(strcmp(deblank({ns_codes(idnsc).name}),'NOTO'))>0
            idnsc=idnsc & strcmp(deblank({ns_codes.name}),'NOTO'); %NOTO, NOTOX
        elseif sum(strcmp(deblank({ns_codes(idnsc).name}),'NRAO85_3'))
            idnsc=idnsc & strcmp(deblank({ns_codes.name}),'NRAO85_3'); %NRAO85_3, VLBA085_3, WIDE85_3
        end
    
    end   
    ns_codes = read_jtrf2020_station_position_xyz_vlbi_data(ns_codes, idnsc, datafile, infls(i).name);
end

% 12717S001 NOTO, NOTOX
% 40441S004_7214 NRAO85_3, VLBA085_3, WIDE85_3