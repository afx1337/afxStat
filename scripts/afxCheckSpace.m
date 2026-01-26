function sameSpace = afxCheckSpace(images, interactive)

    if nargin < 2
        interactive = false;
    end

    % get images interactively
    if ~exist('images','var')
        interactive = true;
        images = cellstr(spm_select([1 Inf],'image','Select images'));
    end
    
    sameSpace = true;

    % dim and q-form matrix of first image
    V = spm_vol(images{1});
    if interactive, spaceOrigStr = strcat('dim=',jsonencode(V.dim),', mat=',jsonencode(V.mat)); end
    spaceOrig.dim = V.dim;
    spaceOrig.mat = V.mat;
    if interactive, fprintf('[first] %s: %s\n',spaceOrigStr,images{1}(end-45:end)); end
    
    % compare against all other images
    fails  ={};
    for i = 2:length(images)
        V = spm_vol(images{i});
        if interactive, spaceCurStr = strcat('dim=',jsonencode(V.dim),', mat=',jsonencode(V.mat)); end
        spaceCur.dim = V.dim;
        spaceCur.mat = V.mat;
        if isequal(spaceOrig.dim,spaceCur.dim) && isequal(spaceOrig.mat,spaceCur.mat)
            if interactive, fprintf('   [ok] %s: %s\n',spaceCurStr,images{i}(end-45:end)); end
        else
            sameSpace = false;
            if interactive
                fprintf(' [fail] %s: %s\n',spaceCurStr,images{i}(end-45:end));
                fails{end+1} = images{i};
            else
                return
            end
        end
    end
    if interactive
        fprintf('\nDifferent space:\n');
        for i = 1:length(fails)
            fprintf(' [fail] %s\n',fails{i})
        end
    end
end