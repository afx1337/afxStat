function [thresh, pVal] = afxPermThreshHigh(nullDist, pcrit, vals)
    % the old way
    thresh = prctile(nullDist(:),(1-pcrit)*100);
    pVal = sum(vals(:) <= nullDist(:)', 2)' / numel(nullDist);
end