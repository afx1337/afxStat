function [Y, comment] = afxStatPrepare(imgFiles, flipLR, FWHM, thr, maskFile, comment)
    % load images
    [Y.dat,XYZ,Y.dim,Y.mat] = afxVolumeRead(imgFiles);
    comment(1).imgFiles = imgFiles;
    
    % load mask
    if ~isempty(maskFile)
        Y.mask = afxVolumeResample(maskFile,XYZ,0)' > .5;
        clear XYZ;
        Y.dat = Y.dat(:,Y.mask);
        comment.maskFile = maskFile;
    end
    
    % flipLR
    if ~isempty(flipLR)
        for iImg = 1:numel(imgFiles)
            if flipLR(iImg)
                tmp = reshape(Y.dat(iImg,:),Y.dim);
                tmp = tmp(end:-1:1,:,:);
                Y.dat(iImg,:) = tmp(:);
            end
        end
        comment.flipLR = flipLR;
    end
    
    % smooth
    if ~isempty(FWHM) && FWHM ~= 0
        Y.dat = afxSmoooth(Y.dat,FWHM,Y.dim,Y.mat);
        comment.FWHM = FWHM;
    end
    
    % threshold
    if ~isempty(thr)
        Y.dat = Y.dat > thr;
        comment.binarizationThreshold = thr;
    end
end
