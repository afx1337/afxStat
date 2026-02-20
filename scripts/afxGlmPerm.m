function [t, tCrit, kCrit, pVal, k] = afxGlmPerm(Y, X, contrast, nPerms, inference, FWE, threshVox, threshClust)
    % [t, tCrit, kCrit] = afxGlmPerm(Y, X, contrast, nPerms, inference, FWE, threshVox, threshClust)
    % returns uncorrected t-scores for all voxels/rois in Y given design matrix in X as well as
    % significance thresholds tCrit and kCrit
    %
    % Y.dat:       Volume/roi data
    % Y.dim:       Necessary for cluster level inference
    % Y.mask:      Inform cluster level inference about mask; note, that Y.dat needs to be already masked
    % X:           Design matrix
    % contrast:    Contrast vecror
    % nPerms:      Number of permutations to compute, negative values for Freedman-Lane procedure
    % inference:   Inference based on voxel value or cluster size ('voxel'/'cluster')
    % FWE:         Multiple comparisons correction? (true/false)
    % threshVox:   Critical p on the voxel level (if negativ value tCrit = abs(val))
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
        fprintf('Will perform Freedman-Lane procedure.\n');
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
        error('Permutations without FWE correction are not supported.')
    end
    
    signPerm = false;
    if ~isBinomial
        c0 = eye(p) - contrast' * pinv(contrast');
        v0 = sum(var(X - X * c0));
        if v0 == 0
            fprintf('Will perform sign permutations.\n');
            signPerm = true;
        end
    end
    
    if nnz(var(X) > eps & abs(contrast) < eps) > 0 && nPerms > 0
        % forces user to think about the design matrix or to give up and
        % use other software
        error('Design contains nuissance variables but Freedman-Lane has not been enabled.');
    end
    
    fprintf('Performing permutations ...\n')
    perm1pct = abs(nPerms) * .01;
    lastLen = printProgress(20, 0, 0, 0);

    % actual permutations to obtain null distribution
    tic
    for iPerm = 1:abs(nPerms)
        if mod(iPerm,perm1pct) < 1
            lastLen = printProgress(20, round(iPerm/perm1pct), toc/60, lastLen);
        end
        
        % permutation of the design matrix
        if signPerm
            Xp  = X.*sign(randn(n,1));
        else
            Xp  = X(randperm(n), :);
        end
        
        % permutation distribution "should" include identity permutation
        % (Winkler et. al., 2014)
        if iPerm == 1
            tp = t;
        else
            tp = afxGlmCompute(Y.dat,Xp,contrast,df,isBinomial);
        end
        
        if clusterInference
            idx = tp > tCrit;
            if any(idx)
                clust = spm_clusters(RCP(:,idx));
                nullDist(iPerm) = nnz(clust == mode(clust));
            else
                nullDist(iPerm) = 0;
            end
        else
            nullDist(iPerm) = max(tp);
        end
    end
     
    if clusterInference
        idx = t > tCrit;
        k = [];
        if any(idx)
            clust = spm_clusters(RCP(:,idx));
            numClust = max(clust);
            k = histc(clust,1:numClust);
        end
        [kCrit, pVal] = afxPermThreshHigh(nullDist,threshClust,k);
        kCrit = ceil(kCrit);
    else
        kCrit = 0;
        [tCrit, pVal] = afxPermThreshHigh(nullDist,threshVox,t);
    end
end

function lastLen = printProgress(barWidth, pct, elapsedMin, lastLen)
    %PRINTPROGRESS Print single-line console progress bar with ETA.
    %   lastLen = PRINTPROGRESS(barWidth, pct, elapsedMin, lastLen)
    %   Pass lastLen from the previous call (use 0 initially).

    nEqual = round(barWidth * pct / 100);
    backStr = repmat('\b',1,lastLen);
    barStr = ['[', repmat('=',1,nEqual), repmat(' ',1,barWidth-nEqual), ']'];
    outStr = sprintf('  %3d %% %s', pct, barStr);
    
    if pct > 0 && pct < 100
        totalEst   = elapsedMin / (pct/100);
        totalSeconds = round(max(0, totalEst - elapsedMin)*60);
        mm = floor(totalSeconds / 60);
        ss = mod(totalSeconds, 60);
        outStr = [outStr, sprintf(' (%d:%02d min left)', mm, ss)];
    end

    fprintf([backStr '%s'], outStr);
    lastLen = length(outStr);

    if pct == 100
        fprintf(' done.\n');
        lastLen = 0;
    end
end
