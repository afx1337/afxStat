function destFolder = afxVlsmFiles(imgFiles, flipLR, FWHM, X, minOverlap, regressLesion, nPerms, inference, FWE, threshVox, threshClust, destFolder, comment)
    if ~exist('comment','var'), comment = struct([]); end

    % pass to afxStatPrepare()
    [Y, comment] = afxStatPrepare(imgFiles, flipLR, FWHM, 0.5, [], comment);

    % pass to afxVlsmPrepare()
    [Y,X,lesionOverlay,comment] = afxVlsmPrepare(Y,minOverlap,X,regressLesion,comment);

    % contrast
    if afxIsBinomial(X(:,1))
		contrast = [1 -1 zeros(1,size(X,2)-2)];
	else
		contrast = [1 zeros(1,size(X,2)-1)];
    end
    
    % pass to afxStat()
    destFolder = afxStat(Y, X, {contrast}, nPerms, inference, FWE, threshVox, threshClust, destFolder, comment);
    
    % save lesion overlay
    afxVolumeWrite(fullfile(destFolder,'sumMap.nii'),lesionOverlay,Y.dim,'uint8',Y.mat,'lesion overlay',false);
end

function isBinom = afxIsBinomial(dat)
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