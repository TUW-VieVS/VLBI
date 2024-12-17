


function ns_codes = read_jtrf2020_station_position_xyz_vlbi_data(ns_codes, idnsc, indire, infile0)

infile = [indire '/' infile0];

fid=fopen(infile);


while ~feof(fid)
    curl=fgetl(fid);
    if isempty(curl)
        curl='                                   ';
    end
    
    if strcmpi(curl(1:20),'# REFERENCE EPOCH T0')
        y = str2double(split(curl(52:70)));
        ep=date2mjd(y');
        ns_codes(idnsc).jtrf2020.break(1).epoch = ep;
        ns_codes(idnsc).jtrf2020.break(1).start = ep;
    elseif strcmpi(curl(1:21),'# EPOCH LAST FORECAST')
        y = str2double(split(curl(52:70)));
        ep=date2mjd(y');
        ns_codes(idnsc).jtrf2020.break(1).end = ep;
    elseif strcmpi(curl(1:24),'# REFERENCE GEOCENTRIC X')
        ns_codes(idnsc).jtrf2020.break(1).x = str2double(curl(52:70));
    elseif strcmpi(curl(1:24),'# REFERENCE GEOCENTRIC Y')
        ns_codes(idnsc).jtrf2020.break(1).y = str2double(curl(52:70));
    elseif strcmpi(curl(1:24),'# REFERENCE GEOCENTRIC Z')
        ns_codes(idnsc).jtrf2020.break(1).z = str2double(curl(52:70));
    end
end
fclose(fid);