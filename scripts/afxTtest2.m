function afxTtest2()
    % user input
    images1 = cellstr(spm_select([1 Inf],'image','Select images 1'));
    images2 = cellstr(spm_select([1 Inf],'image','Select images 2'));
    maskFile = cellstr(spm_select([1 1],'image','Select explicit mask'));
    maskFile = maskFile{1};
    title = spm_input('Name',1,'s',[date '-']);
    nPerms = spm_input('Permutations',2,'e','5000',1);
    inference = spm_input('Inference',3,'b',{'voxel' 'cluster'},[true false],1);
    if inference, inference = 'voxel'; else, inference = 'cluster'; end
    FWE = spm_input('correction',4,'b',{'FWE' 'uncorr'},[true false],1);
    if strcmp(inference,'cluster')
        threshVox = spm_input('p(vox)',5,'e','0.001',1);
        threshClust = spm_input('p(clust)',6,'e','0.05',1);
    else
        threshVox = spm_input('p(vox)',5,'e','0.05',1);
        threshClust = [];
    end
    if FWE, FWEstr = 'FWE'; else, FWEstr = 'uncorr'; end
    destFolder = fullfile('results','ttest2',title,[inference '-' FWEstr ]);
    % load data
    [Y.dat,XYZ,Y.dim,Y.mat] = afxVolumeRead([images1 ; images2]);
    % load mask
    Y.mask = afxVolumeResample(maskFile,XYZ,0)' > .5;
    Y.dat = Y.dat(:,Y.mask);
    % design
    X = [[ones(length(images1),1) zeros(length(images1),1); zeros(length(images2),1) ones(length(images2),1)] ones(size(Y.dat,1),1)];
    % contrast (Group1 > Group2?)
    contrasts{1} = [1 -1 0];
    contrasts{2} = [-1 1 0];
    for i = 1:length(contrasts)
        % stat
        [t, tCrit, kCrit, pVal, k] = afxGlmPerm(Y, X, contrasts{i}, nPerms, inference, FWE, threshVox, threshClust);
        % save everything
        info = struct('images1',{images1},'images2',{images2},'design',X,'contrast',contrasts{i},'inference',inference,'correction',FWEstr,'mask',maskFile,'nPerms',nPerms,'threshVox',threshVox,'threshClust',threshClust,'tCrit',tCrit,'kCrit',kCrit,'pValues',pVal,'clusterSizes',k);
        afxGlmWrite(destFolder,i,Y.dim,Y.mat,Y.mask,t,tCrit,kCrit,info);
    end
end
