% This function is mainly a wrapper to call all necessary subfunctions of
% vie_optimize_conditions
%
% CREATED: 24.02.17 - Matthias Schartner
%
% CHANGELOG: 
function [ flag_break, sched, source ] = optimizer_main( sched, station, twin, source, obsmode, PARA, iOpt,isched )

% if you have one check your schedule
[ flag, goodSources, badSources ] = optimizer_get_good_sources( sched, PARA.optimization_condition );
flag_break = 0;

if flag
    % break if everything is ok
    flag_break = 1;
elseif length(goodSources) < 15
    % if there are less than 15 sources left, also break
    fprintf(PARA.fid_body,'optimization condition was not fulfilled but there are only %s sources left!\n iterative scheduling was stopped',length(goodSources));
    warning('optimization condition was not fulfilled but there are only %s sources left!\n iterative scheduling was stopped!',length(goodSources));
    flag_break = 1;
elseif strcmp(PARA.optimization_type,'quick fillin')
    [ sched ] = optimizer_quick_fillin( sched, station, twin, source, obsmode, PARA, goodSources, badSources );
    flag_break = 1;
else
    % select only the good sources and redo everything
    iOpt = iOpt+1;
    fprintf(PARA.fid_body,'\n\n******************************************************************\n');
    fprintf(PARA.fid_body,'** optimization condition was not fulfilled, start new schedule ** \n**             schedule %3d iteration %3d sources %3d           **\n',isched,iOpt,length(goodSources));
    fprintf(PARA.fid_body,'******************************************************************\n');

    % get all source names which are not used:
    allSource = 1:length(source);
    usedSources = [goodSources,badSources];
    bool = false(size(allSource));
    bool(usedSources)=1;
    notUsedNames = {source(~bool).name};
    fprintf(PARA.fid_body,'Removed all sources which were not observed in the previouse schedule\n');
    for i = 1:length(notUsedNames)
        fprintf(PARA.fid_body,'%s ',notUsedNames{i});
        if mod(i,7)==0
            fprintf(PARA.fid_body,'%s\n',notUsedNames{i});
        end
    end
    fprintf(PARA.fid_body,'\n');

    % get all source names which are bad
    fprintf(PARA.fid_body,'Removed all sources which do not fulfill the condition:\n');
    allBadNames = {source(badSources).name};
    for i = 1:length(allBadNames)
        fprintf(PARA.fid_body,'%s ',allBadNames{i});
        if mod(i,7)==0
            fprintf(PARA.fid_body,'%s\n',allBadNames{i});
        end
    end
    fprintf(PARA.fid_body,'\n\n');
    source = source(goodSources);
end

end

