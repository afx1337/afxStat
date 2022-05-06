function thresh = afxPermThreshHigh(permScores, kPcrit)
    permScores = sort(permScores(:),'descend');
    thresh = permScores(ceil(numel(permScores) * kPcrit));
    %report next most significant score in case of ties
    permScores = permScores(permScores > thresh);
    if ~isempty(permScores)
        thresh = min(permScores(:));
    end
end