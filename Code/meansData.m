function meansData() 
    outputDir = sprintf('%s\\output',pwd);
    csvfiles = dir(fullfile(outputDir,'*.csv'))
    for k = 1:length(csvfiles)
        baseFileName = csvfiles(k).name;
        fullFileName = fullfile(outputDir, baseFileName);
        [folder, baseFileNameNoExt, extension] = fileparts(fullFileName);
        
        output = csvread(fullFileName);
        outputFolder = sprintf('%s\\outputmeans',pwd)
        if not(isfolder(outputFolder))
            mkdir (outputFolder);
        end
        outputfile = sprintf('%s\\outputmeans\\%s_means%s', pwd, baseFileNameNoExt, extension);
        
        outputmeans = mean(output);
        csvwrite(outputfile, outputmeans);       
    end
end