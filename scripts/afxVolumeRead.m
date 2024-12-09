function [y,XYZmm,dim,mat] = afxVolumeRead(func)
    fprintf('Loading imaging data ')
    % read imaging data and convert to one dimensional vector
    
    % check cache
    if ismember('afxCache',evalin('base','who'))
        afxCache = evalin('base','afxCache');
        if strcmp(jsonencode(func),jsonencode(afxCache.func))
            fprintf('... using cached data ... done.\n');
            y = afxCache.y;
            XYZmm = afxCache.XYZmm;
            dim = afxCache.dim;
            mat = afxCache.mat;
            return
        end
    end

    % new (~4 seconds for a NKI data set)
    for  i = 1:length(func)
        fprintf('.');
        tmp = nifti(func{i});
        if size(tmp.dat,4) > 1
            parfor j = 1:size(tmp.dat,4)
                y(j,:) = reshape(tmp.dat(:,:,:,j),1,[]);
            end
            break
        else
            y(i,:) = tmp.dat(:);
        end
    end
    
    % load first image for MNI coordinates in mm and for
    % dimensions etc.
    Vfunc = spm_vol(func{1});
    [~,XYZmm] = spm_read_vols(spm_vol(Vfunc));
    XYZmm = [XYZmm; ones(1,size(XYZmm,2))];
    dim = Vfunc.dim;
    mat = Vfunc.mat;
    
    % save data to cache
    afxCache = struct([]);
    afxCache(1).y = y;
    afxCache.XYZmm = XYZmm;
    afxCache.dim = dim;
    afxCache.mat = mat;
    afxCache.func = func;
    assignin('base','afxCache',afxCache);
    
    fprintf(' done.\n');
end
