function y = afxSmoooth(y,FWHM,dim,mat)
    % y = afxSmoooth(y,FWHM,dim,mat)
    % convolve data in y (linearly packed in the spatial dimension with
    % original dimensions dim and q-form matrix mat) with gausian smoothing
    % kernel with full width at half maximum of FWHM

    fprintf('   Smoothing ...')
    if length(FWHM) == 1, FWHM = [FWHM FWHM FWHM]; end
    if FWHM ~= [0 0 0]
        % adapt fwhm to voxel size
        FWHM = FWHM./sqrt(sum(mat(1:3,1:3).^2));
        % generate dummy matrix for smoothed data
        sdat = nan(dim);
        % iterate over all timepoints
        for i = 1:size(y,1)
            % perform smoothing
            spm_smooth(reshape(y(i,:),dim),sdat,FWHM);
            % copy data back to y
            y(i,:) = sdat(:);
        end
    end
    fprintf(' done\n')
end
