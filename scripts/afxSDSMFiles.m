function destFolder = afxSDSMFiles(imgNetwork, flipLR, FWHM, thr, imgLesion, X, minOverlap, contrasts, nPerms, inference, FWE, threshVox, threshClust, destFolder, comment)
    if ~exist('comment','var'), comment = struct([]); end

    % add regressor for lesion size
    [lesions.dat,~,~,lesions.mat] = afxVolumeRead(imgLesion);
    lesionSizes = sum(lesions.dat,2)*prod(sqrt(sum(lesions.mat(1:3,1:3).^2)))/1000;
    clear lesions;
    X = [ X lesionSizes./max(lesionSizes)];

    % pass to afxStatPrepare()
    [Y, comment] = afxStatPrepare(imgNetwork, flipLR, FWHM, thr, [], comment);

    % pass to afxVlsmPrepare()
    [Y,X,lesionOverlay,comment] = afxVlsmPrepare(Y,minOverlap,X,'sdsm',comment);

    % pass to afxStat()
    destFolder = afxStat(Y, X, contrasts, nPerms, inference, FWE, threshVox, threshClust, destFolder, comment);
    
    % save lesion overlay
    afxVolumeWrite(fullfile(destFolder,'sumMap.nii'),lesionOverlay,Y.dim,'uint8',Y.mat,'lesion overlay',false);
end