function [Y,X,lesionOverlay,comment] = afxVlsmPrepare(Y,minOverlap,X,regressLesion,comment)
    % binarize data
    Y.dat = Y.dat ~= 0;

    % calculate lesion volumes and overlay
    lesionVolumes = sum(Y.dat,2);
    lesionOverlay = sum(Y.dat,1);

    % calculate minimum lesion overlap (neg. val. -> abs; pos. val. -> pct)
    if minOverlap < 0
        comment(1).minOverlapAbs = -minOverlap;
        comment.minOverlapPct = round(-minOverlap/size(Y.dat,1)*100,2);
    else
        comment.minOverlapAbs = ceil(size(Y.dat,1)*minOverlap/100);
        comment.minOverlapPct = minOverlap;
    end
    Y.mask = lesionOverlay >= comment.minOverlapAbs;
    if ~any(Y.mask), error('No voxels with enough lesion coverage.'); end
    % mask data
    Y.dat = Y.dat(:,Y.mask);
    
    % adress lesion volume
    switch regressLesion
        case 'sdsm'
            % do nothing
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
            Y.dat = Y.dat ./ sqrt(lesionVolumes);
        otherwise
            error(['Unkonwn option for lesion volume treatment: ' regressLesion ]);
    end
    comment.regressLesion = regressLesion;
end