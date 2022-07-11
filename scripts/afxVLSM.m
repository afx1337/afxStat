function [destFolder,tCrit, kCrit] = afxVLSM(designFile, minOverlap, regressLesion, nPerms, inference, FWE, threshVox, threshClust)
    [pth,designName,~] = fileparts(designFile);
    tmpFile = fullfile(pth,[designName '-tmp.mat']);
    if ~exist(tmpFile,'file')
        [Y,X,rowLabels,colLabels] = afxReadDesign(designFile);
        Y.dat = Y.dat > .5;
        lesionVolumes = sum(Y.dat,2);
        lesionOverlay = sum(Y.dat,1);
        save(tmpFile,'Y','X','rowLabels','colLabels','lesionVolumes','lesionOverlay');
    else
        load(tmpFile);
    end
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
    Y.dat = Y.dat(:,Y.mask);
    
    switch regressLesion
        case 'none'
            warning('Lesion volume will not be addresed in the current analysis.');
        case 'regress' % niiStat strategy
            Xtmp = [lesionVolumes ones(size(lesionVolumes,1),1)];
            for i = 1:size(X,2)
                beta = Xtmp\X(:,i);
                X(:,i) = X(:,i) - Xtmp*beta;
            end
        case 'covariate' % e.g. https://www.nature.com/articles/s41598-017-08728-x (also uses cluster based inference) (https://doi.org/10.1038/s41598-017-08728-x)
            X = [X lesionVolumes];
        case 'direct' % "dTLVC", see http://europepmc.org/backend/ptpmcrender.fcgi?accid=PMC4400840&blobtype=pdf (doi:10.1038/ncomms7762)
            Y.dat = Y.dat./repmat(sqrt(lesionVolumes),1,size(Y.dat,2));
        otherwise
            error(['Unkonwn option for lesion volume treatmend: ' regressLesion ]);
    end
        
    if ~any(var(X) < eps & mean(X) ~= 0) % constant term?
        X = [X ones(size(Y.dat,1),1) ];
    end
    contrast = [1 zeros(1,size(X,2)-1)];

    [t, tCrit, kCrit, pVal, k] = afxGlmPerm(Y, X, contrast, nPerms, inference, FWE, threshVox, threshClust);
    if FWE, FWE = 'FWE'; else, FWE = 'uncorr'; end
    nowstr = datestr(datetime('now'),'yyyymmdd_HHMMSS');
    destFolder = fullfile('results','VLSM',designName,[inference '-' FWE '-' regressLesion 'Lesion'],nowstr);
    
    info = struct('designFile',designFile,'rowLabels',{rowLabels},'colLabels',{colLabels},'design',X,'inference',inference,'correction',FWE,'controlForLesionVolume',regressLesion,'minOverlapPct',minOverlap,'minOverlapAbs',minOverlapAbs,'nPerms',nPerms,'threshVox',threshVox,'threshClust',threshClust,'tCrit',tCrit,'kCrit',kCrit,'pValues',pVal,'clusterSizes',k);
    afxGlmWrite(destFolder,1,Y.dim,Y.mat,Y.mask,t,tCrit,kCrit,info);
    afxVolumeWrite(fullfile(destFolder,'sumMap.nii'),lesionOverlay,Y.dim,'uint8',Y.mat,'lesion overlay',false);
end