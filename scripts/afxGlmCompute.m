function s = afxGlmCompute(Y,X,c,df,isBinomial)
    if isBinomial
        s = lieberMeister(Y,X);
    else
        [n, reg] = size(X);
        [nY, ~] = size(Y);
        if nY ~= n, error('X and Y data sizes are inconsistent.'); end
        if reg ~= size(c,2), error('Contrast vector has wrong length.'); end

        %df = n - rank(X);                  % degrees of freedom
        %Xpinv = pinv(X);                   % pseudoinverse of X
        %pXX = Xpinv*Xpinv';                % = pinv(X'*X)
        pXX = pinv(X'*X);                   % a bit faster
        pX  = pXX * X';                     % = pinv(X'*X)*X'
        b = pX * Y;                         % parameters (betas)
        Y = Y - X * b;                      % residuals
        s = sum(Y .* Y);                    % sum squared error
        s = sqrt(s .* (c * pXX * c') / df); % standard error
        s = c * b ./ s;                     % t-statistic
    end
end

function z = lieberMeister(groupVector, obsVector)
    %Liebermeister measure
    %[z] = lieber([1 1 1 1 0 0 0 0; 0 1 0 0 1 1 1 1]',[1 1 1 1 0 0 1 0]')
    %[z] = lieber2([1 1 1 1 0 0 0 0]',[0 1 0 0 1 1 1 1]')
    if size(obsVector,2) ~= 1
        error('Error: obsVector must be a single column\n');
    end
    M = size(obsVector,1); %number of observations
    if (M ~= size(groupVector,1))
        error('%s error: group and observations matrices must have same number of rows\n',mfilename);
    end
    if M < 3
        error('%s error: Unabled to test voxels with less than 3 observations',mfilename);
    end
    N = sum(obsVector(:)); %behavior = 1
    a = sum(bsxfun(@min,groupVector,obsVector));
    K = sum(groupVector,1);
    pL = fexact(a,M,K,N,'test','l','tail','l');
    pR = fexact(a,M,K,N,'test','l','tail','r');
    indexL = find(pL(:) < pR(:));
    p2z = @(p) -sqrt(2) * erfcinv(p*2);
    z = pR; %vector Z now probability of right tail
    z(indexL) = pL(indexL); %vector Z now probability of extreme tail
    z = -p2z(z); %convert probabilities to Z-scores, 1/2014: minus so p0.05 is positive 1.64
    z(indexL) = -z(indexL); %negative correlations have negative Z scores
end