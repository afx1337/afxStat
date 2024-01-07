function [destFolder, tCrit, kCrit] = afxSDSM(imgDisco, imgLesion, minOverlap, X, contrast, nPerms, inference, FWE, threshVox, threshClust, destFolder, comment, flipLR, thr)
    % defaults
    if ~exist('thr','var'),     thr = .6; end
    if ~exist('flipLR','var'),  flipLR = false(1,size(X,1)); end

    % load disconnection images
    if ~isstruct(imgDisco)
        [Y.dat,~,Y.dim,Y.mat] = afxLoadFunc(imgDisco);
        rowLabels = imgDisco;
    else
        Y = imgDisco;
        rowLabels = 'unknown';
    end
    
    % add regressor for lesion size
    [lesions.dat,~,~,lesions.mat] = afxLoadFunc(imgLesion);
    lesionSizes = sum(lesions.dat,2)*prod(sqrt(sum(lesions.mat(1:3,1:3).^2)))/1000;
    X = [ X lesionSizes./max(lesionSizes)];
    
    % prepare VLSM
    [Y,X,contrast,lesionOverlay,minOverlapAbs,minOverlapPct] = afxPrepareVLSM(Y,minOverlap,X,contrast,'sdsm',flipLR,1);
    % thresholding (greater or equal, see Wawrzyniak et al., 2022)
    Y.dat = Y.dat >= thr;
    
    % perform GLM permutations
    [t, tCrit, kCrit, pVal, k] = afxGlmPerm(Y, X, contrast, nPerms, inference, FWE, threshVox, threshClust);
    
    % save results
    if FWE, FWE = 'FWE'; else, FWE = 'uncorr'; end
    info = struct('images',{imgDisco},'rowLabels',{rowLabels},'design',X,'contrast',contrast,'comment',comment,'inference',inference,'correction',FWE,'controlForLesionVolume',regressLesion,'minOverlapPct',minOverlapPct,'minOverlapAbs',minOverlapAbs,'nPerms',nPerms,'threshVox',threshVox,'threshClust',threshClust,'flipLR',flipLR,'thr',thr,'tCrit',tCrit,'kCrit',kCrit,'pValues',pVal,'clusterSizes',k);
    afxGlmWrite(destFolder,1,Y.dim,Y.mat,Y.mask,t,tCrit,kCrit,info);
    afxVolumeWrite(fullfile(destFolder,'sumMap.nii'),lesionOverlay,Y.dim,'uint8',Y.mat,'lesion overlay',false);
end