function dat = afxVolumeResample(volume,XYZorig,hold)
    % dat = afxVolumeResample(fname,XYZ,hold)
    % resample fname to space defined by XYZ and return data as (linearly
    % packed) vector
    % hold refers to interpolation method for the resampling (see
    % spm_sample_vol.m), (default: nearest neighbour)
    
    if nargin < 3, hold = 0; end
    % load volume
    Vresample = spm_vol(volume);
    
    % preallocate
    nImg = numel(Vresample);
    nVox = size(XYZorig,2);
    dat  = nan(nImg, nVox, 'single');
        
    for i = 1:nImg
        % compute voxel coordinates in voxel space of volume wich correspond to
        % the space given by XYZ
        RCPresample = Vresample(i).mat\XYZorig;
        % resample volume
        datTmp  = spm_sample_vol(Vresample(i),RCPresample(1,:),RCPresample(2,:),RCPresample(3,:),hold);
        % pack to vector
        dat(i,:) = datTmp(:);
    end
end
