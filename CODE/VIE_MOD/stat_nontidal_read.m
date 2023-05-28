
function data = stat_nontidal_read(idoy, numyrs,iye, ntsl_path,ntsl_model,ntsl_suffix)

data = cell(1,5);
for idyr=0:numyrs-1
    f = dir([ntsl_path ntsl_model '/*' ntsl_suffix]); 
    if idoy>360
        id = contains({f.name},num2str(iye+idyr)); % session in Dec
    else
        id = contains({f.name},num2str(iye-idyr)); % session in Jan
    end
    fil = [ntsl_path ntsl_model '/' f(id).name];

    data_y = cell(1,5);
    if sum(id)>0
        if exist(fil,'file')
            fid = fopen(fil);
            data_y = textscan(fid,'%s%f%f%f%f','CommentStyle','!');
            fclose(fid);
        else
            data =cell(1,5);
        end
    else
        data =cell(1,5);
    end
    data = cellfun(@(a,b) [a;b], data, data_y, 'uni',0 );
end
