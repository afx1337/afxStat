function [Y,X,contrast,lesionOverlay,minOverlapAbs,minOverlapPct] = afxPrepareVLSM(Y,minOverlap,X,contrast,regressLesion,flipLR,thr)
    % thresholding
    Y.dat = Y.dat > thr;
    
    % lesion volumes and overlay
    lesionVolumes = sum(Y.dat,2);
    lesionOverlay = sum(Y.dat,1);
    
    % flipLR
    for iImg = 1:size(Y.dat,1)
        if flipLR(iImg)
            tmp = reshape(Y.dat(iImg,:),Y.dim);
            tmp = tmp(end:-1:1,:,:);
            Y.dat(iImg,:) = tmp(:);
        end
    end
    
    % calculate min overlap
    if minOverlap < 0
        minOverlapAbs = -minOverlap;
        minOverlapPct = round(minOverlapAbs/size(Y.dat,1)*100,2);
        Y.mask = lesionOverlay >= minOverlapAbs;
    else
        minOverlapAbs = ceil(size(Y.dat,1)*minOverlap/100);
        minOverlapPct = minOverlap;
        Y.mask = lesionOverlay >= minOverlapAbs;
    end
    if ~any(Y.mask), error('No voxels with enough lesion coverage.'); end
    % mask data
    Y.dat = Y.dat(:,Y.mask);
    
    % adress lesion volume
    switch regressLesion
        case 'sdsm'
            ; % do nothing
        case 'none'
            warning('Lesion volume will not be addresed in the current analysis.');
        case 'regress' % niiStat strategy
            Xtmp = [lesionVolumes./max(lesionVolumes) ones(size(lesionVolumes,1),1)];
            for i = 1:size(X,2)
                beta = Xtmp\X(:,i);
                X(:,i) = X(:,i) - Xtmp*beta;
            end
        case 'covariate' % e.g. https://www.nature.com/articles/s41598-017-08728-x (also uses cluster based inference) (https://doi.org/10.1038/s41598-017-08728-x)
            X = [X lesionVolumes./max(lesionVolumes)];
        case 'direct' % "dTLVC", see http://europepmc.org/backend/ptpmcrender.fcgi?accid=PMC4400840&blobtype=pdf (doi:10.1038/ncomms7762)
            Y.dat = Y.dat./repmat(sqrt(lesionVolumes),1,size(Y.dat,2));
        otherwise
            error(['Unkonwn option for lesion volume treatmend: ' regressLesion ]);
    end

    % add constant term to design matrix
    if ~any(var(X) == 0 & mean(X) ~= 0) % constant term?
        X = [X ones(size(Y.dat,1),1) ];
    end
    
    % pad contrast vector with zeros
    contrast = [contrast zeros(1,size(X,2)-size(contrast,2))];
end