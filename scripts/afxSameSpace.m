function afxSameSpace()
    spaceImg = cellstr(spm_select([1 1],'image','Select space defining image'));
    images = cellstr(spm_select([1 Inf],'image','Select other images'));
    hold = spm_input('Interpolation method',1,'m',{'nearest neighbour' 'trilinear'},{ 0 1 },1);
    dt = spm_input('Data type',2,'m',{'uint8', 'int16', 'int32', 'float32', 'float64', 'int8', 'uint16', 'uint32'},{'uint8', 'int16', 'int32', 'float32', 'float64', 'int8', 'uint16', 'uint32'},2);
    gz = spm_input('Data type',3,'m',{'nii' 'nii.gz'},{false true},1);
    ss = spm_input('Slope scaling',4,'m',{'yes' 'no'},{true false},1);
    h = waitbar(0,'Progress ...');
    
    V = spm_vol(spaceImg{1});
    [~,XYZmm] = spm_read_vols(spm_vol(V));
    XYZmm = [XYZmm; ones(1,size(XYZmm,2))];
    
    for i = 1:length(images)
        dat = afxVolumeResample(images{i},XYZmm,hold{1});
        afxSave(images{i},dat,V.dim,dt{1},V.mat,gz{1},ss{1});
        waitbar(i/length(images),h)
    end
end

function afxSave(fnameOrig,dat,dim,dt,mat,gz,ss)
    [p,f,~] = fileparts(fnameOrig);
    np = fullfile(p,'SS');
    if ~exist(np,'dir')
        mkdir(p,'SS');
    end
    if strcmp(f(1:3),'SS_')
        f = f(4:end);
    end
    fname = fullfile(np,strcat('SS_',f,'.nii'));
    afxVolumeWrite(fname,dat,dim,dt,mat,'',ss);
    if gz
        gzip(fname);
        delete(fname);
    end
end
