function [destFolder,tCrit, kCrit] = afxStat(designFile, contrasts, maskFile, nPerms, inference, FWE, threshVox, threshClust)
    [Y,X,rowLabels,colLabels,XYZ] = afxReadDesign(designFile);
    % load mask
    Y.mask   = afxVolumeResample(maskFile,XYZ,0)' > .5;
    Y.dat = Y.dat(:,Y.mask);
    % design
    if ~any(var(X) == 0 & mean(X) ~= 0) % constant term?
        X = [X ones(size(Y.dat,1),1) ];
    end
    
    % output dir
    tmp = designFile(end:-1:1); tmp = fileparts(tmp);
    tmp = tmp(end:-1:1); tmp = fileparts(tmp); tmp = fileparts(tmp);
    if FWE, tmp2 = 'FWE'; else, tmp2 = 'uncorr'; end
    [~,designName,~] = fileparts(designFile);
    destFolder = fullfile('results','stats',tmp,designName,[inference '-' tmp2 ]);
    
    for i = 1:length(contrasts)
        contrast = contrasts{i};
        contrast = [contrast zeros(1,size(X,2)-size(contrast,2))];
        % stat
        [t, tCrit, kCrit, pVal, k] = afxGlmPerm(Y, X, contrast, nPerms, inference, FWE, threshVox, threshClust);

        % save results
        info = struct('designFile',designFile,'rowLabels',{rowLabels},'colLabels',{colLabels},'design',X,'contrast',contrast,'inference',inference,'correction',tmp2,'mask',maskFile,'nPerms',nPerms,'threshVox',threshVox,'threshClust',threshClust,'tCrit',tCrit,'kCrit',kCrit,'pValues',pVal,'clusterSizes',k);
        afxGlmWrite(destFolder,i,Y.dim,Y.mat,Y.mask,t,tCrit,kCrit,info);
    end
end