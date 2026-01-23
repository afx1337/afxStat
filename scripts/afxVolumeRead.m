function [y,XYZmm,dim,mat] = afxVolumeRead(func, varargin)
    % function [y,XYZmm,dim,mat] = afxVolumeRead(filenames, varargin)
    %
    % Load NIfTI images with optional caching and dynamic cutoff
    %
    % Optional parameters:
    %   'precision' : 'single' (default) or 'logical'
    %   'cacheCutoffGB' : maximum RAM for cache (default: 1 GB)

    p = inputParser;
    addParameter(p,'precision','single',@(x) ismember(x,{'single','logical'}));
    addParameter(p,'cacheCutoffGB',1,@(x) isnumeric(x) && x>0);
    parse(p,varargin{:});

    precision = p.Results.precision;
    cacheCutoffGB = p.Results.cacheCutoffGB;

    fprintf('Loading imaging data ');

    %% --- Estimate memory requirement ---
    tmp = nifti(func{1});
    nVolumes = size(tmp.dat,4);
    nVox = prod(size(tmp.dat,1:3));

    x = zeros(1,1,precision);
    bytesPerVoxel = whos('x');
    bytesPerVoxel = bytesPerVoxel.bytes;

    estRAMBytes = nVox * max(nVolumes, numel(func)) * bytesPerVoxel;
    useCache = estRAMBytes <= cacheCutoffGB*1e9;

    %% --- Check cache if allowed ---
    if useCache && evalin('base','exist(''AFX_CACHE'',''var'')')
        afxCache = evalin('base','AFX_CACHE');
        if isstruct(afxCache) && isfield(afxCache,'func') && isequal(func,{afxCache.func}')
            fprintf('... using cached data ... done.\n');
            y = afxCache.y;
            XYZmm = afxCache.XYZmm;
            dim = afxCache.dim;
            mat = afxCache.mat;
            return
        end
    end

    %% --- Preallocate y ---
    if nVolumes > 1
        assert(numel(func)==1,'4D data only supported for single file');
        y = zeros(nVolumes,nVox,precision);
    else
        y = zeros(numel(func),nVox,precision);
    end

    %% --- Load in batches ---
    nFiles = numel(func);
    for i = 1:nFiles
        fprintf('.');
        if i > 1, tmp = nifti(func{i}); end

        if nVolumes > 1
            parfor j = 1:nVolumes
                if strcmp(precision,'logical')
                    y(j,:) = reshape(tmp.dat(:,:,:,j),[1,nVox])~=0;
                else
                    y(j,:) = single(reshape(tmp.dat(:,:,:,j),[1,nVox]));
                end
            end
            break
        else
            if strcmp(precision,'logical')
                y(i,:) = tmp.dat(:)~=0;
            else
                y(i,:) = single(tmp.dat(:));
            end
        end
    end


    %% --- Load coordinates and affine ---
    Vfunc = spm_vol(func{1});
    [~,XYZmm] = spm_read_vols(Vfunc);
    XYZmm = [XYZmm; ones(1,size(XYZmm,2))];
    dim = Vfunc.dim;
    mat = Vfunc.mat;

    %% --- Save to cache if allowed ---
    if useCache
        afxCache = struct( ...
            'y',y, ...
            'XYZmm',XYZmm, ...
            'dim',dim, ...
            'mat',mat, ...
            'func',func );
        assignin('base','AFX_CACHE',afxCache);
    end

    fprintf(' done.\n');
end
