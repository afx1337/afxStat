function afxGlmWrite(destFolder, c, dim, mat, mask, t, tCrit, kCrit, info)
    info.df = size(info.design,1)-rank(info.design);
	info.matlabversion = version();
	
    if exist(destFolder,'dir')
        warning(['Output folder ' destFolder ' already exists. Files will be overwritten.']);
    else
        mkdir(destFolder);
    end
    cstr = sprintf('%03d',c);
    fileMask = fullfile(destFolder,'mask.nii');
    fileTraw = fullfile(destFolder,['TMap_' cstr '.nii']);
    fileTthresh = fullfile(destFolder,['TMap_' cstr '_filtered.nii']);
    fileInfoMat = fullfile(destFolder,['info_' cstr '.mat']);
    fileInfoTxt = fullfile(destFolder,['info_' cstr '.txt']);
    
    afxVolumeWrite(fileMask,mask,dim,'uint8',mat,'mask');
    
    img = zeros(size(mask));
    img(mask) = t;
    afxVolumeWrite(fileTraw,img,dim,'int16',mat,'untresholded t-statistic');
    
    if kCrit > 0
        [L,numClust] = spm_bwlabel(reshape(double(img > tCrit),dim));
        k = histc(L(:),1:numClust);
        sigClust = find(k >= kCrit);
        sigVox = ismember(L(:),sigClust);
        afxVolumeWrite(fileTthresh,sigVox'.*img,dim,'int16',mat,'tresholded t-statistic');
    else
        sigClust = [];
        afxVolumeWrite(fileTthresh,(img > tCrit).*img,dim,'int16',mat,'tresholded t-statistic');
    end
    
    save(fileInfoMat, 'info');
    if exist(fileInfoTxt,'file'), delete(fileInfoTxt); end
    diary(fileInfoTxt);
    disp(info);
    fprintf('\nSuprathreshold voxels: %i; suprathreshold clusters: %i\n',nnz(t>tCrit),length(sigClust));
    
    % print cluster table
    if strcmp(info.inference,'cluster')
        fprintf('\n cluster | p      | k\n');
        for i=1:length(sigClust)
            fprintf('      %2d | %.4f | %5d\n',i,info.pValues(sigClust(i)),info.clusterSizes(sigClust(i)));
        end
    end
    diary off;
end