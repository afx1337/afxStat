function sameSpace = afxCheckSpace(images)

    % get images interactively
    if ~exist('images','var')
        interactive = true;
        images = cellstr(spm_select([1 Inf],'image','Select images'));
    else
        interactive = false;
    end
    
    sameSpace = true;

    % dim and q-form matrix of first image
    V = spm_vol(images{1});
    spaceOrig = strcat('dim=',jsonencode(V.dim),', mat=',jsonencode(V.mat));
    if interactive, fprintf('[first] %s: %s\n',spaceOrig,images{1}(end-45:end)); end
    
    % compare against all other images
    fails  ={};
    for i = 2:length(images)
        V = spm_vol(images{i});
        spaceCur = strcat('dim=',jsonencode(V.dim),', mat=',jsonencode(V.mat));
        if strcmp(spaceOrig,spaceCur)
            if interactive, fprintf('   [ok] %s: %s\n',spaceCur,images{i}(end-45:end)); end
        else
            sameSpace = false;
            if interactive
                fprintf(' [fail] %s: %s\n',spaceCur,images{i}(end-45:end));
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