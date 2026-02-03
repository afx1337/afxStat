function destFolder = afxStatXlsx(designFile, contrasts, maskFile, nPerms, inference, FWE, threshVox, threshClust)
    % Read design file
    [imgFiles,colLabels,X] = afxReadDesign(designFile);
    comment.designFile = designFile;
    comment.colLabels = colLabels;
    
    % output dir
    if FWE, FweStr = 'FWE'; else, FweStr = 'uncorr'; end
    [~,designName,~] = fileparts(designFile);
    destFolder = fullfile('results','Stat',designName,[inference '-' FweStr ]);

    % pass to afxStatFiles()
    destFolder = afxStatFiles(imgFiles, [], [], [], X, contrasts, maskFile, nPerms, inference, FWE, threshVox, threshClust, destFolder, comment);
end
