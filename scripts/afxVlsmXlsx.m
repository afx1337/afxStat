function destFolder = afxVlsmXlsx(designFile, minOverlap, regressLesion, nPerms, inference, FWE, threshVox, threshClust)
    % Read design file
    [imgFiles,colLabels,X] = afxReadDesign(designFile);
    comment.designFile = designFile;
    comment.colLabels = colLabels;
    
    % output dir
    if FWE, FweStr = 'FWE'; else, FweStr = 'uncorr'; end
    [~,designName,~] = fileparts(designFile);
    destFolder = fullfile('results','VLSM',designName,[inference '-' FweStr '-' regressLesion 'Lesion']);
    
    % pass to afxVlsmFiles()
    destFolder = afxVlsmFiles(imgFiles, X, minOverlap, regressLesion, nPerms, inference, FWE, threshVox, threshClust, destFolder, comment);
end

