function afxTtest()
    % user input
    
    images = cellstr(spm_select([1 Inf],'image','Select images'));
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
    destFolder = fullfile('results','ttest',title,[inference '-' FWEstr ]);
    
    % preare
    imgFiles = images;
    % design
    X = ones(length(images),1);
    % contrast
    contrasts{1} = [1];
    
    % pass to afxStatFiles()
    afxStatFiles(imgFiles, [], [], [], X, contrasts, maskFile, nPerms, inference, FWE, threshVox, threshClust, destFolder);
end
