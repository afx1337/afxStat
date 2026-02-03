function destFolder = afxStat(Y, X, contrasts, nPerms, inference, FWE, threshVox, threshClust, destFolder, comment)
    % design -> add constant term
    if ~any(var(X) <= eps & mean(X) ~= 0) % constant term?
        X = [X ones(size(Y.dat,1),1) ];
    end
    
    % FWE -> string
    if FWE, FweStr = 'FWE'; else, FweStr = 'uncorr'; end
    
    % for each contrast
    for i = 1:length(contrasts)
        contrast = contrasts{i};
        contrast = [contrast zeros(1,size(X,2)-size(contrast,2))]; % pad contrast vector with zeros
        % stat
        [t, tCrit, kCrit, pVal, k] = afxGlmPerm(Y, X, contrast, nPerms, inference, FWE, threshVox, threshClust);

        % save results
        comment.design = X;
        comment.contrast = contrast;
        comment.inference = inference;
        comment.FWE = FweStr;
        comment.nPerms = nPerms;
        comment.threshVox = threshVox;
        comment.threshClust = threshClust;
        comment.tCrit = tCrit;
        comment.kCrit = kCrit;
        comment.pValues = pVal;
        comment.clusterSizes = k;
        afxGlmWrite(destFolder,i,Y.dim,Y.mat,Y.mask,t,tCrit,kCrit,comment);
    end
end
