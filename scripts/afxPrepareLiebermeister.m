function [isBinomial, dat1, dat2] = afxPrepareLiebermeister(dat1,dat2)
    isBinom1 = afxIsBinomial(dat1);
    isBinom2 = afxIsBinomial(dat2);
    if isBinom1 && isBinom2 && size(dat2,2) == 2
        disp("Perform Liebermeister tests");
        isBinomial = true;
        dat2 = dat2(:,1);
        dat1 = afxScaleBinary(dat1);
        dat2 = afxScaleBinary(dat2);
    else
        isBinomial = false;
    end
end

function isBinom = afxIsBinomial(dat)
    % check if data is binominal
    mn = min(dat(:));
    mx = max(dat(:));
    nMin = sum(dat(:)==mn);
    nMax = sum(dat(:)==mx);
    if (nMin+nMax) ~= length(dat(:))
        isBinom = false;
    else
        isBinom = true;
    end
end

function dat = afxScaleBinary(dat)
    thr = mean(dat(:));
    dat(dat(:) < thr) = 0;
    dat(dat(:) >= thr) = 1;
end