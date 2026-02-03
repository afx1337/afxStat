function [t, tCrit, kCrit, pVal, k] = afxGlmPerm(Y, X, contrast, nPerms, inference, FWE, threshVox, threshClust)
    % [t, tCrit, kCrit] = afxGlmPerm(Y, X, contrast, nPerms, inference, threshVox, threshClust)
    % returns uncorrected t-score for all voxels/rois in Y given design matrix in X
    %   if either X or Y is binary, results are pooled variance t-test
    %   (=two-sample), else correlation coefficients
    %
    % Y.dat:       Volume/roi data
    % Y.dim:       Necessary for cluster level inference
    % Y.mask:      Inform cluster level inference about mask; note, that Y.dat needs to be already masked
    % X:           Design matrix (without constant term) (matrix)
    % contrast:    Contrast vecror
    % nPerms:      Number of permutations to compute (integer)
    % inference:   Inference based on voxel value or cluster size ('voxel'/'cluster')
    % FWE:         Multiple comparisons correction? (true/false)
    % threshVox:   Critical p
    % threshClust: Critical p for cluster level inference
    %
    % inspired by niiStat by Chris Rorden
    
    % defaults
    if ~isfield(Y,'mask'), Y.mask = []; end
    if ~exist('nPerms','var'), nPerms = 5000; end
    if ~exist('inference','var'), inference = 'cluster'; end
    if ~exist('FWE','var'), FWE = true; end
    if ~exist('threshVox','var'), threshVox = 0.001; end
    if ~exist('threshClust','var'), threshClust = 0.05; end
    if strcmp(inference,'cluster'),   clusterInference = true;
    elseif strcmp(inference,'voxel'), clusterInference = false;
    else, error(['afxGLMPerm: unknown inference mode: ' inference]);
    end
    pVal = []; k = [];
    
    % Basics and reusable components
    [n, p] = size(X);  % rows = observations, columns = factors
    [nY, v] = size(Y.dat); % v is number of tests/voxels
    if nY ~= n, error('afxGlmPerm X and Y data sizes are inconsistent'); end
    df = n - rank(X);  % degrees of freedom

    [isBinomial, Y.dat, X] = afxPrepareLiebermeister(Y.dat,X);
   
    % Freedman-Lane procedure (see Winkler et al., 2014)    
    if nPerms < 0
        fprintf('Performing Freedman-Lane procedure ...\n');
        % niiStat strategy
        c0 = eye(p) - contrast' * pinv(contrast');
        X0 = X * c0;
        R0 = eye(n) - X0 * pinv(X0);
        Y.dat = R0 * Y.dat;
    end
    
    % Original design (identity permutation)
    t = afxGlmCompute(Y.dat,X,contrast,df,isBinomial);    
    if any(~isfinite(t(:)))
        warning('afxGlmPerm zeroed NaN t-scores');
        t(~isfinite(t))= 0;
    end
    
    if clusterInference
        % critical t value
        if threshVox > 1
            tCrit = threshVox;
        else
            % TODO: chage tinv some spm_builtin
            if isBinomial
                tCrit = sqrt(2) * erfcinv(threshVox*2);
            else
                tCrit = tinv(1-threshVox,df);
            end
        end
        % build voxel space
        [R,C,P]  = ndgrid(1:Y.dim(1),1:Y.dim(2),1:Y.dim(3));
        RCP = [R(:)';C(:)';P(:)'];
        clear R C P;
        if ~isempty(Y.mask), RCP = RCP(:,Y.mask); end
        if abs(nPerms) < 2
            kCrit = Inf;
            warning('No permutations performed. Critical cluster threshold unknown and set to Inf.');
            return
        end
    else
        if abs(nPerms) < 2
            kCrit = 0;
            if FWE
                tCrit = tinv(1-threshVox/v,df);
                warning('No permutations performed. Critical voxel threshold is Bonferroni corrected.');
            else
                tCrit = tinv(1-threshVox,df);
                warning('No permutations performed. Critical voxel threshold is uncorrected.');
            end
            return
        end
    end

    if FWE
        % initialize distribution of largest cluster/voxel value
        nullDist = nan(1,abs(nPerms));
    else
        nullDist = [];
    end
    
    signPerm = false;
    if ~isBinomial
        c0 = eye(p) - contrast' * pinv(contrast');
        v0 = sum(var(X - X * c0));
        if v0 == 0
            fprintf('Perform sign permutation ...\n');
            signPerm = true;
        end
    end
    
    perm5pct = abs(nPerms) * .05;
    perm2pct = ceil(abs(nPerms) * .02);

    fprintf('Performing permutations ...\n')
    % actual permutations to obtain null distribution
    tic
    for p = 1:abs(nPerms)
        if p == perm2pct
            fprintf('  Estimated time: %.1f minutes\n  |                  |\n  ',round(toc*50/60,1));
        end
        if mod(p,perm5pct) < 1
            fprintf('*');
        end
        
        % permutation of the design matrix
        if signPerm
            Xp  = X.*(randi(2,nY,1)*2-3);
        else
            Xp  = X(randperm(n), :);
        end
        
        % permutation distribution "should" include identity permutation
        % (Winkler et. al., 2014)
        if p == 1
            tp = t;
        else
            tp = afxGlmCompute(Y.dat,Xp,contrast,df,isBinomial);
        end
        
        if FWE
            if clusterInference
                idx = tp > tCrit;
                if any(idx)
                    clust = spm_clusters(RCP(:,idx));
                    nullDist(p) = nnz(clust == mode(clust));
                else
                    nullDist(p) = 0;
                end
            else
                nullDist(p) = max(tp);
            end
        else
            if clusterInference
                idx = tp > tCrit;
                if any(idx)
                    clust = spm_clusters(RCP(:,idx));
                    nullDist = [nullDist histcounts(clust,unique(clust)) ];
                else
                    nullDist(end+1) = 0;
                end
            else
                nullDist = [nullDist tp];
            end
        end
    end
    fprintf(' done.\n');
     
    if clusterInference
        kCrit = afxPermThreshHigh(nullDist,threshClust);
    else
        kCrit = 0;
        tCrit = afxPermThreshHigh(nullDist,threshVox);
    end
    
    if nargout > 3
        pVal = [];
        if clusterInference
            idx = t > tCrit;
            if any(idx)
                clust = spm_clusters(RCP(:,idx));
                numClust = max(clust);
                k = histc(clust,1:numClust);
                for i = 1:length(k)
                    pVal(i) = nnz(k(i) <= nullDist)/abs(nPerms);
                end
            end
        else
            for i = 1:length(t)
                pVal(i) = nnz(t(i) <= nullDist)/abs(nPerms);
            end
        end
    end
end