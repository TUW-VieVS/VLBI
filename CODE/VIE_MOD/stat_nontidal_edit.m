


function stcorr_xyz =  stat_nontidal_edit(data, aname, ant, mjd1, mjd2 )

    % reduce data to lines that contain the respective station
    wantedLines = ismember(data{1,1}, strtrim(aname));
    for k=1:length(data)
        data_short{k} = data{k}(wantedLines);
    end

    % make data smaller by extracting only data on the respective day + 1 day before and after
    wantedLines = data_short{2}>mjd1-1.25 & data_short{2}<mjd2+1;
    for k=1:length(data_short)
        data_short{k} = data_short{k}(wantedLines);
    end

    if ~isempty(data_short)
        % ren -> xyz
        matm2 = [data_short{2} data_short{3} data_short{4} data_short{5}];
        [stcorr_xyz] = call_ren2xyz(matm2,ant);
    end

