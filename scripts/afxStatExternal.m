function [destFolder,tCrit, kCrit] = afxStatExternal(imgFiles, FWHM, X, contrasts, maskFile, nPerms, inference, FWE, threshVox, threshClust, destFolder, comment, flipLR)
    % load images
    [Y.dat,XYZ,Y.dim,Y.mat] = afxLoadFunc(imgFiles);
    
    % flipLR
    if exist('flipLR','var')
        for iImg = 1:length(imgFiles)
            if flipLR(iImg)
                tmp = reshape(Y.dat(iImg,:),Y.dim);
                tmp = tmp(end:-1:1,:,:);
                Y.dat(iImg,:) = tmp(:);
            end
        end
    end
    
    % smooth
    if ~isempty(FWHM)
        Y.dat = afxSmoooth(Y.dat,FWHM,Y.dim,Y.mat);
    end
    rowLabels = imgFiles;
    
    % load mask
    Y.mask = afxVolumeResample(maskFile,XYZ,0)' > .5;
    Y.dat = Y.dat(:,Y.mask);
    
    % design
    if ~any(var(X) == 0 & mean(X) ~= 0) % constant term?
        X = [X ones(size(Y.dat,1),1) ];
    end
    if FWE, tmp2 = 'FWE'; else, tmp2 = 'uncorr'; end
    
    for i = 1:length(contrasts)
        contrast = contrasts{i};
        contrast = [contrast zeros(1,size(X,2)-size(contrast,2))];
        % stat
        [t, tCrit, kCrit, pVal, k] = afxGlmPerm(Y, X, contrast, nPerms, inference, FWE, threshVox, threshClust);

        % save results
        info = struct('images',{imgFiles},'rowLabels',{rowLabels},'design',X,'contrast',contrast,'comment',comment,'inference',inference,'correction',tmp2,'mask',maskFile,'nPerms',nPerms,'threshVox',threshVox,'threshClust',threshClust,'tCrit',tCrit,'kCrit',kCrit,'pValues',pVal,'clusterSizes',k);
        afxGlmWrite(destFolder,i,Y.dim,Y.mat,Y.mask,t,tCrit,kCrit,info);
    end
end