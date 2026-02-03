function destFolder = afxStatFiles(imgFiles, flipLR, FWHM, thr, X, contrasts, maskFile, nPerms, inference, FWE, threshVox, threshClust, destFolder, comment)
    if ~exist('comment','var'), comment = struct([]); end

	% pass to afxStatPrepare()
    [Y, comment] = afxStatPrepare(imgFiles, flipLR, FWHM, thr, maskFile, comment);
    
    % pass to afxStat()
    destFolder = afxStat(Y, X, contrasts, nPerms, inference, FWE, threshVox, threshClust, destFolder, comment);
end