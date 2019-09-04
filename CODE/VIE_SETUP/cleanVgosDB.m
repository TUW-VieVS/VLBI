function [] = cleanVgosDB()

    folders = dir('../DATA/vgosDB');
    folders = folders(3:end);
    for i = 1:length(folders)
        if strcmp(folders(i).name,'.gitignore')
            continue;
        end
        
        files = dir([folders(i).folder '/' folders(i).name]);
        files = files(3:end);
        for j = 1:length(files)
            if strcmp(files(j).name,'.gitignore')
                continue;
            end
            
            if files(j).isdir
                try
                    rmdir([files(j).folder '/' files(j).name],'s')
                catch
                    warning('%s could not be removed from VgosDB folder!', files(j).name)
                end
            end
        end
    end
end

