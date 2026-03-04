function alpha = afxGlmLesionContributionMaps(imgFiles, imgOptions, X, contrast, sigMask, lesFiles, lesOptions, destFolder)
    % alpha = afxGlmLesionContributionMaps(imgFiles, imgOptions, X, contrast, sigMask, lesFiles, lesOptions, destFolder)
    %
    % Estimates individual contribution/alignment with Lesion Network Mapping
    % and projects these estimates back to lesion space.
    %
    % Input:
    %   imgFiles    cell array  Lesion network maps (filenames)
    %   imgOptions  struct      Options (flipLR FWHM thr maskFile)
    %   X           matrix      Design matrix
    %   contrast    matrix      Contrast vector
    %   sigMask     char        Mask of significant voxels (J)
    %   lesFiles    cell array  Lesion masks (filenames)
    %   lesOptions  struct      Options (flipLR FWHM thr)
    %   destFolder  char        Output folder
    %
    % Output:
    %   alpha       vector      Individual contributions to LNM effect
    %
    % Example:
    %   lnmFiles = {'lnm1.nii' 'lnm2.nii' ... };
    %   imgOptions.flipLR = [];
    %   imgOptions.FWHM = [];
    %   imgOptions.thr = [];
    %   imgOptions.maskFile = 'lnmMask.nii';
    %   X = [0 1; 0 1; 1 0; ...];
    %   contrast = [1 -1];
    %   sigMask = 'significant_voxels.nii';
    %   lesionFiles = {'lesion1.nii' 'lesion2.nii' ... };
    %   lesOptions.flipLR = [];
    %   lesOptions.FWHM = [];
    %   lesOptions.thr = [];
    %   destDir = 'backProjection';
    %
    %   afxGlmLesionContributionMaps(lnmFiles, imgOptions, X, contrast, sigMask, lesionFiles, lesOptions, destDir);

    % pass to afxStatPrepare() to load LNM data
    [Y, ~] = afxStatPrepare(imgFiles, imgOptions.flipLR, imgOptions.FWHM, imgOptions.thr, imgOptions.maskFile, struct([]));
    
    % load significant voxels (J)
    J = afxStatPrepare(imgFiles(1), [], [], [], sigMask);
    J = J.mask(Y.mask);
    
    % OLS
    pXX = pinv(X'*X);
    pX  = pXX * X';
    b   = pX * Y.dat; % beta

    % contrast specific weights (w)
    %W = (X / (X' * X)) * contrast';
    W = pX' * contrast';

    % contributions of individual observations (Y scaled)
    Ys = W .* Y.dat; % contribution
    
    % effect size and t-statistic (for reference)
    df = size(X,1) - rank(X); 
    R = Y.dat - X * b;                  % residuals
    s = sum(R .* R);                    % sum squared error
    s = sqrt(s .* (contrast * pXX * contrast') / df); % standard error
    theta = contrast * b;               % effect size estimate
    t = theta ./ s;                     % t-statistic
    
    % weight contributions in significant voxels (J) by effect size (alpha)
    alpha = Ys(:,J).*theta(J);
    
    % calculate sum and scale to std
    alpha = sum(alpha,2);
    alpha = alpha/std(alpha);
        
    % load and binarize lesions (L)
    [L, ~] = afxStatPrepare(lesFiles, lesOptions.flipLR, lesOptions.FWHM, lesOptions.thr, [], struct([]));
    L.dat = L.dat ~= 0;
    
    % make output folder if non-existant
    if ~exist(destFolder','dir'), mkdir(destFolder); end
   
    % save general data
    save(fullfile(destFolder,'info.mat'),'imgFiles','imgOptions','X','contrast','sigMask','lesFiles','lesOptions','alpha');

    % save t-map (for reference)
    tMap = nan(size(Y.mask), 'like', t);
    tMap(Y.mask) = t;
    afxVolumeWrite(fullfile(destFolder,'tMap.nii'),tMap,Y.dim,'int16',Y.mat);
   
    % save mean weighted lesion map
    L.dat = L.dat.*alpha;
    afxVolumeWrite(fullfile(destFolder,'meanMap.nii'),mean(L.dat),L.dim,'int16',L.mat);
    
    % save group specific maps
    for iCol = 1:size(X,2)
        regressor = X(:,iCol);
        if numel(unique(regressor)) == 2
            regressor = regressor > (min(regressor)+max(regressor))/2;
            afxVolumeWrite(fullfile(destFolder,sprintf('meanMap_G%02d.nii',iCol)),sum(L.dat(regressor,:))./numel(regressor),L.dim,'int16',L.mat);
        end
    end
    
end