function [] = cleanVgosDB()

    folders = dir('../DATA/vgosDB');
    folders = folders(3:end);
     
    % loop over all directories
   for i = 1:length(folders)
        if strcmp(folders(i).name,'.gitignore')
            continue;
        end
        
        files = dir([folders(i).folder '/' folders(i).name]);
        files = files(3:end);
            
        % get names without file ending
        names = {files.name};
        for j = 1:length(names)
            tmp = names{j};
            idx = strfind(tmp,'.');
            if ~isempty(idx)
                tmp = extractBefore(tmp,min(idx));
            end
            names{j} = tmp;
        end
        
        % loop over all files
        for j = 1:length(files)
            if strcmp(files(j).name,'.gitignore')
                continue;
            end
            
            if files(j).isdir
                try
                    % check if file name is defined multiple times (then one
                    % file should be the .tar or .tar.gz)
                    if sum(strcmp(files(j).name,names)) > 1
                        rmdir([files(j).folder '/' files(j).name],'s')
                    end
                catch
                    warning('%s could not be removed from VgosDB folder!', files(j).name)
                end
            end
        end
    end
end

