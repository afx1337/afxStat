function [destFolder, tCrit, kCrit] = afxVLSMExternalData(Y, rowLabels, minOverlap, X, contrast, regressLesion, nPerms, inference, FWE, threshVox, threshClust, destFolder, comment, flipLR, thr)
    % defaults
    if ~exist('thr','var'),     thr = .5; end
    if ~exist('flipLR','var'),  flipLR = false(1,length(imgFiles)); end

    % prepare VLSM
    [Y,X,contrast,lesionOverlay,minOverlapAbs,minOverlapPct] = afxPrepareVLSM(Y,minOverlap,X,contrast,regressLesion,flipLR,thr);
    
    % perform GLM permutations
    [t, tCrit, kCrit, pVal, k] = afxGlmPerm(Y, X, contrast, nPerms, inference, FWE, threshVox, threshClust);
    
    % save results
    if FWE, FWE = 'FWE'; else, FWE = 'uncorr'; end
    info = struct('images',{imgFiles},'rowLabels',{rowLabels},'design',X,'contrast',contrast,'comment',comment,'inference',inference,'correction',FWE,'controlForLesionVolume',regressLesion,'minOverlapPct',minOverlapPct,'minOverlapAbs',minOverlapAbs,'nPerms',nPerms,'threshVox',threshVox,'threshClust',threshClust,'tCrit',tCrit,'kCrit',kCrit,'pValues',pVal,'clusterSizes',k);
    afxGlmWrite(destFolder,1,Y.dim,Y.mat,Y.mask,t,tCrit,kCrit,info);
    afxVolumeWrite(fullfile(destFolder,'sumMap.nii'),lesionOverlay,Y.dim,'uint8',Y.mat,'lesion overlay',false);
end