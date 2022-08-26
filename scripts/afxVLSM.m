function destFolder = afxVLSM(designFile, minOverlap, regressLesion, nPerms, inference, FWE, threshVox, threshClust)
    % load design and data
    [pth,designName,~] = fileparts(designFile);
    tmpFile = fullfile(pth,[designName '-tmp.mat']);
    if ~exist(tmpFile,'file')
        [Y,X,rowLabels,colLabels] = afxReadDesign(designFile);
        save(tmpFile,'Y','X','rowLabels','colLabels');
    else
        load(tmpFile);
    end
    
    % contrast
    if afxIsBinomial(X(:,1))
		contrast = [1 -1 zeros(1,size(X,2)-2)];
	else
		contrast = [1 zeros(1,size(X,2)-1)];
	end
    
    % prepare VLSM
    thr = .5;
    flipLR = false(1,size(Y.dat,1));
    [Y,X,contrast,lesionOverlay,minOverlapAbs,minOverlapPct] = afxPrepareVLSM(Y,minOverlap,X,contrast,regressLesion,flipLR,thr);

    % perform GLM permutations
    [t, tCrit, kCrit, pVal, k] = afxGlmPerm(Y, X, contrast, nPerms, inference, FWE, threshVox, threshClust);
    
    % save results
    if FWE, FWE = 'FWE'; else, FWE = 'uncorr'; end
    nowstr = datestr(datetime('now'),'yyyymmdd_HHMMSS');
    destFolder = fullfile('results','VLSM',designName,[inference '-' FWE '-' regressLesion 'Lesion'],nowstr);
    info = struct('designFile',designFile,'rowLabels',{rowLabels},'colLabels',{colLabels},'design',X,'contrast',contrast,'inference',inference,'correction',FWE,'controlForLesionVolume',regressLesion,'minOverlapPct',minOverlapPct,'minOverlapAbs',minOverlapAbs,'nPerms',nPerms,'threshVox',threshVox,'threshClust',threshClust,'tCrit',tCrit,'kCrit',kCrit,'pValues',pVal,'clusterSizes',k);
    afxGlmWrite(destFolder,1,Y.dim,Y.mat,Y.mask,t,tCrit,kCrit,info);
    afxVolumeWrite(fullfile(destFolder,'sumMap.nii'),lesionOverlay,Y.dim,'uint8',Y.mat,'lesion overlay',false);
end

function isBinom = afxIsBinomial(dat)
    % check if data is binominal
    mn = min(dat(:));
    mx = max(dat(:));
    nMin = sum(dat(:)==mn);
    nMax = sum(dat(:)==mx);
    if (nMin+nMax) ~= length(dat(:))
        isBinom = false;
    else
        isBinom = true;
    end
end