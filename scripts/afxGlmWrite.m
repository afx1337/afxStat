function afxGlmWrite(destFolder, c, dim, mat, mask, t, tCrit, kCrit, info)
    info.df = size(info.design,1)-rank(info.design);
    info.metaDate = datestr(now);
    info.metaVersion = afxVersion();
    info.metaMatlab = version();
    info.metaToolboxes = ver();

    cstr = sprintf('%03d',c);
    fileMask = fullfile(destFolder,'mask.nii');
    fileTraw = fullfile(destFolder,['TMap_' cstr '.nii']);
    fileTthresh = fullfile(destFolder,['TMap_' cstr '_filtered.nii']);
    fileInfoMat = fullfile(destFolder,['info_' cstr '.mat']);
    fileInfoTxt = fullfile(destFolder,['info_' cstr '.txt']);

    if exist(destFolder,'dir')
        if c == 1
            warning(['Output folder ' destFolder ' already exists. Files will be overwritten.']);
        end
    else
        mkdir(destFolder);
    end
    
    if c == 1, afxVolumeWrite(fileMask,mask,dim,'uint8',mat,'mask',false); end
    
    nSigVox = nnz(t>tCrit);
    img = zeros(size(mask),'like',t);
    img(mask) = t;
    clear t;
    afxVolumeWrite(fileTraw,img,dim,'int16',mat,'untresholded t-statistic');
    
    if kCrit > 0
        [L,numClust] = spm_bwlabel(reshape(double(img > tCrit),dim));
        k = accumarray(L(L>0),1,[numClust 1]);
        sigClust = find(k >= kCrit);
        sigVox = ismember(L(:),sigClust);
        clear L k;
        afxVolumeWrite(fileTthresh,sigVox.*img,dim,'int16',mat,'tresholded t-statistic');
        clear sigVox;
    else
        sigClust = [];
        afxVolumeWrite(fileTthresh,(img > tCrit).*img,dim,'int16',mat,'tresholded t-statistic');
    end
    
    save(fileInfoMat, 'info');
    if exist(fileInfoTxt,'file'), delete(fileInfoTxt); end
    diary(fileInfoTxt);
    disp(info);
    fprintf('\nSuprathreshold voxels: %i; suprathreshold clusters: %i\n',nSigVox,length(sigClust));
    
    % print cluster table
    if strcmp(info.inference,'cluster')
        fprintf('\n cluster | p      | k\n');
        for i=1:length(sigClust)
            fprintf('      %2d | %.4f | %5d\n',i,info.pValues(sigClust(i)),info.clusterSizes(sigClust(i)));
        end
    end
    diary off;
end

function version = afxVersion()
    thisFile = mfilename('fullpath');
    thisDir = fileparts(thisFile);
    versionFile = fullfile(thisDir,'..', 'VERSION');
    if exist(versionFile, 'file')
        version = strtrim(fileread(versionFile));
    else
        version = 'unknown';
        warning('VERSION file not found.');
    end
end
