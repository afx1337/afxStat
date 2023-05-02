function varargout = afxOrthviews(action,varargin)
    % afxOrthviews is a simple interface to SPMs secions viewer. It
    % supports contrasts from SPM.mat, niftis, masks, spheres and raw data.
    %
    % afxOrthviews('init' [, underlay])
    %    Initialize SPM graphics window and load underlay.
    %    underlay: anatomical underlay (default: ch2better_1mm.nii)
    %
    % afxOrthviews('spm', filename, contrast [, threshVox [, threshClust [, reposition ]]])
    %    Visualize a contrast from a SPM.mat.
    %    filename: path and filename of SPM.mat
    %    contrast: contrast number
    %    threshVox: high threshold (default: 0.001)
    %    threshClust: cluster-level threshold (default: 0.05)
    %    reposition: go to global maximum (default: true)
    %
    % afxOrthviews('nii', filename, threshVox [, threshClust [, reposition [, title]]])
    %    Visualize nifti.
    %    filename: name of nifti
    %    threshVox: voxel-level threshold
    %    threshClust: cluster extent (default: 0)
    %    reposition: go to global maximum (default: true)
    %    title: title string
    %
    % afxOrthviews('nii_resample', filename, threshVox [, threshClust [, size [, reposition [, title]]]])
    %    Visualize nifti.
    %    filename: name of nifti
    %    threshVox: voxel-level threshold
    %    threshClust: cluster extent (default: 0)
    %    size: resample to size mm (default: 1)
    %    reposition: go to global maximum (default: true)
    %    title: title string
    %
    % afxOrthviews('dat', data, dim, mat, threshVox [, threshClust [, reposition [, title]]])
    %    Visualize data.
    %    data: data (linear)
    %    dim: dimensions
    %    mat: q-form matrix
    %    threshVox: voxel-level threshold
    %    threshClust: cluster extent (default: 0)
    %    reposition: go to global maximum (default: true)
    %    title: title string
    %
    % afxOrthviews('mask', filename [, title])
    %    Visualize mask (=binary nifti).
    %    filename: name of nifti
    %    title: title string
    %
    % afxOrthviews('sphere', coord, radius [, title])
    %    Visualise coordinate with a sphere.
    %    coord: coordinate vector [x y z] in mm in MNI-space.
    %    radius: radius in mm
    %    title: title string
    %
    % afxOrthviews('reset')
    %    Clears all added visualisations.
    %
    % afxOrthviews('pos')
    %    Returns current position in MNI space.
    %
    % afxOrthviews('title',title)
    %    Sets new title.
    global hTitAx hResAx;
    if ~strcmpi(action,'init') && (isempty(spm_figure('FindWin','Graphics')) || ~exist('hTitAx','var') || ~exist('hResAx','var'))
        afxOrthviews('init');
    end
    
    switch lower(action)
        case 'init'
            if length(varargin) >= 1
                niiUnderlay = varargin{1};
                customBB = false;
            else
                niiUnderlay = fullfile('templates','MNI152_T1_0.5mm_masked.nii');
                customBB = true;
            end
            spm_figure('GetWin','Graphics');
            spm_orthviews('Reset');
            spm_figure('Clear','Graphics');
            spm_orthviews('Interp',0)
            spm_orthviews('Image', niiUnderlay);
            if customBB
                spm_orthviews('BB',[-74 -108 -60; 74 76 87]);
            else
                spm_orthviews('MaxBB');
            end
            FS = spm('FontSizes');
            hTitAx = axes('Parent',1,...
                'Position',[0.02 0.96 0.96 0.04],...
                'Visible','off');
            hResAx = axes('Parent',1,...
                'Position',[0.05 0.05 0.45 0.05],...
                'DefaultTextVerticalAlignment','baseline',...
                'DefaultTextFontSize',FS(9),...
                'DefaultTextColor',[1,1,1]*.5,...
                'Units','points',...
                'Visible','off');
        case 'spm'
            SPMfname = varargin{1};
            SPMcontrast = varargin{2};
            if length(varargin) < 3
                threshVox = .001;
            else
                threshVox = varargin{3};
            end
            if length(varargin) < 4
                threshClust = .05;
            else
                threshClust = varargin{4};
            end
            if length(varargin) < 5
                reposition = true;
            else
                reposition = varargin{5};
            end
            % load SPM
            load(SPMfname);
            % get search space voxels
            RCP = SPM.xVol.XYZ;
            % degrees of freedom
            df = [SPM.xCon(SPMcontrast).eidf SPM.xX.erdf];
            % statistic (T/F/...)
            STAT = SPM.xCon(SPMcontrast).STAT;
            % number of contrasts
            n = 1;
            % residuals
            R = SPM.xVol.R;
            % get get high threshold
            u = spm_u(threshVox^(1/n),df,STAT);
            % get cluster extent (gaussian rft)
            V2R = 1/prod(SPM.xVol.FWHM);
            [swd,~,~] = fileparts(SPMfname);
            V = spm_vol(fullfile(swd,SPM.xCon(SPMcontrast).Vspm.fname));
            Z = spm_data_read(V,'xyz',RCP);
            [~,~,ue] = spm_uc_clusterFDR(threshClust,df,STAT,R,n,Z,RCP,V2R,u);
            % apply high threshold
            ind = find(Z(:) > u);
            % get clusters
            clust = spm_clusters(RCP(:,ind));
            numClust = max(clust);
            N = histc(clust,1:numClust);
            % apply critical cluster extent
            sigClust = find(N >= ue);
            ind = ind(ismember(clust,sigClust));
            % add overlay
            spm_orthviews('RemoveBlobs',1);
            spm_figure('ColorMap','gray-hot');
            V = SPM.xCon(SPMcontrast).Vspm;
            spm_orthviews('AddBlobs', 1, RCP(:,ind), Z(ind), V.mat);
            afxTitle(SPM.xCon(SPMcontrast).name);
            afxText(sprintf('cluster-level: p(FWE) < %0.2f (k = %d)\nvoxel-level: p < %0.3f (%s > %0.5f)',threshClust,ue,threshVox,STAT,u));
            spm_orthviews('Redraw');
            % goto global max
            if ~isempty(ind) && reposition
                peakCoord = V.mat*[RCP(:,Z==max(Z(ind)));1];
                spm_orthviews('Reposition',peakCoord(1:3));
            end
        case 'dat' % Z, dim, mat, u, k, rep, title
            Z = varargin{1};
            V.dim = varargin{2};
            V.mat = varargin{3};
            threshVox = varargin{4};
            if length(varargin) < 5
                threshClust = 0;
            else
                threshClust = varargin{5};
            end
            if length(varargin) < 6
                reposition = true;
            else
                reposition = varargin{6};
            end
            if length(varargin) < 7
                myTitle = '';
            else
                myTitle = varargin{7};
            end
            % load overlay
            [R,C,P]  = ndgrid(1:V.dim(1),1:V.dim(2),1:V.dim(3));
            RCP = [R(:)';C(:)';P(:)'];
            clear R C P
            % apply high threshold
            if threshVox < 0
                Z = -Z;
                threshVox = -threshVox;
                spm_figure('ColorMap','gray-cool');
            else
                spm_figure('ColorMap','gray-hot');
            end
            ind = find(Z(:) > threshVox);
            % apply cluster threshold
            clust = spm_clusters(RCP(:,ind));
            numClust = max(clust);
            N = histc(clust,1:numClust);
            sigClust = find(N >= threshClust);
            ind = ind(ismember(clust,sigClust));
            % add overlay
            spm_orthviews('RemoveBlobs',1);
            spm_orthviews('AddBlobs', 1, RCP(:,ind), Z(ind), V.mat);
            afxTitle(myTitle);
            afxText(sprintf('threshold: %0.3f; cluster extent: k = %d',threshVox,threshClust));
            spm_orthviews('Redraw');
            % goto global max
            if ~isempty(ind) && reposition
                peakCoord = V.mat*[mean(RCP(:,Z==max(Z(ind))),2);1];
                spm_orthviews('Reposition',peakCoord(1:3));
            end
        case 'nii'
            niiOverlay = varargin{1};
            threshVox = varargin{2};
            if length(varargin) < 3
                threshClust = 0;
            else
                threshClust = varargin{3};
            end
            if length(varargin) < 4
                reposition = true;
            else
                reposition = varargin{4};
            end
            if length(varargin) < 5
                [~,myTitle,~] = fileparts(niiOverlay);
            else
                myTitle = varargin{5};
            end
            % load overlay
            V = spm_vol(niiOverlay);
            Z = spm_read_vols(V);
            [R,C,P]  = ndgrid(1:V(1).dim(1),1:V(1).dim(2),1:V(1).dim(3));
            RCP = [R(:)';C(:)';P(:)'];
            clear R C P
            % apply high threshold
            if threshVox < 0
                Z = -Z;
                threshVox = -threshVox;
                spm_figure('ColorMap','gray-cool');
            else
                spm_figure('ColorMap','gray-hot');
            end
            ind = find(Z(:) > threshVox);
            % apply cluster threshold
            clust = spm_clusters(RCP(:,ind));
            numClust = max(clust);
            N = histc(clust,1:numClust);
            sigClust = find(N >= threshClust);
            ind = ind(ismember(clust,sigClust));
            % add overlay
            spm_orthviews('RemoveBlobs',1);
            spm_orthviews('AddBlobs', 1, RCP(:,ind), Z(ind), V.mat);
            afxTitle(myTitle);
            afxText(sprintf('threshold: %0.3f; cluster extent: k = %d',threshVox,threshClust));
            spm_orthviews('Redraw');
            % goto global max
            if ~isempty(ind) && reposition
                peakCoord = V.mat*[mean(RCP(:,Z==max(Z(ind))),2);1];
                spm_orthviews('Reposition',peakCoord(1:3));
            end
        case 'nii_resample'
            niiOverlay = varargin{1};
            threshVox = varargin{2};
            if length(varargin) < 3
                threshClust = 0;
            else
                threshClust = varargin{3};
            end
            if length(varargin) < 4
                sz = 1;
            else
                sz = varargin{4};
            end
            if length(varargin) < 5
                reposition = true;
            else
                reposition = varargin{5};
            end
            if length(varargin) < 6
                [~,myTitle,~] = fileparts(niiOverlay);
            else
                myTitle = varargin{6};
            end
            % load overlay
            [V,Z,RCP,factor] = afxResample1mm(niiOverlay,sz);
            % apply high threshold
            if threshVox < 0
                Z = -Z;
                threshVox = -threshVox;
                spm_figure('ColorMap','gray-cool');
            else
                spm_figure('ColorMap','gray-hot');
            end
            ind = find(Z(:) > threshVox);
            % apply cluster threshold
            clust = spm_clusters(RCP(1:3,ind));
            numClust = max(clust);
            N = histc(clust,1:numClust);
            sigClust = find(N >= threshClust*factor);
            ind = ind(ismember(clust,sigClust));
            % add overlay
            spm_orthviews('RemoveBlobs',1);
            spm_orthviews('AddBlobs', 1, RCP(1:3,ind), Z(ind), V.mat);
            afxTitle(myTitle);
            afxText(sprintf('resampled to %gx%gx%g mm\nthreshold: %0.3f; cluster extent: k = %d*%d',sz,sz,sz,threshVox,threshClust,factor));
            spm_orthviews('Redraw');
            % goto global max
            if ~isempty(ind) && reposition
                peakCoord = V.mat*[mean(RCP(1:3,Z==max(Z(ind))),2);1];
                spm_orthviews('Reposition',peakCoord(1:3));
            end       
        case 'mask'
            niiOverlay = varargin{1};
            if length(varargin) < 2
                [~,myTitle,~] = fileparts(niiOverlay);
            else
                myTitle = varargin{2};
            end
            V = spm_vol(niiOverlay);
            Z = spm_read_vols(V);
            [R,C,P]  = ndgrid(1:V(1).dim(1),1:V(1).dim(2),1:V(1).dim(3));
            RCP = [R(:)';C(:)';P(:)'];
            ind = Z(:)>.5;
            clear R C P;
            % add overlay
            spm_orthviews('AddColouredBlobs',1,RCP(:,ind),1,V.mat,[.85 0 0]);
            afxTitle(myTitle);
            afxText('');
            spm_orthviews('Redraw');
            if ~isempty(ind)
                peakCoord = V.mat*[mean(RCP(:,ind),2);1];
                spm_orthviews('Reposition',peakCoord(1:3));
            end
        case 'sphere'
            coord = varargin{1};
            radius = varargin{2};
            if length(varargin) < 3
                myTitle = mat2str(coord);
            else
                myTitle = varargin{3};
            end
            V.mat = [1 0 0 -80; 0 1 0 -120; 0 0 1 -70; 0 0 0 1];
            V.dim = [160 210 160];
            [R,C,P]  = ndgrid(1:V(1).dim(1),1:V(1).dim(2),1:V(1).dim(3));
            RCP = [R(:)';C(:)';P(:)'];
            clear R C P;
            RCP(4,:) = 1;
            XYZ = V.mat(1:3,:)*RCP;
            D = sqrt((XYZ(1,:)-coord(1)).^2+(XYZ(2,:)-coord(2)).^2+(XYZ(3,:)-coord(3)).^2);
            ind = D < radius;
            Z = ind;
            % add overlay
            spm_orthviews('AddColouredBlobs',1,RCP(:,ind),1,V.mat,[.85 0 0]);
            afxTitle(myTitle);
            afxText('');
            spm_orthviews('Redraw');
            if ~isempty(ind)
                peakCoord = V.mat*mean(RCP(:,ind),2);
                spm_orthviews('Reposition',peakCoord(1:3));
            end
        case 'reset'
            spm_orthviews('RemoveBlobs',1);
            afxTitle('');
            afxText('');
        case 'pos'
            varargout{1} = spm_orthviews('Pos');
        case 'title'
            afxTitle(varargin{1});
    end
    
    % write nifti
    if exist('RCP','var') && nargout == 1
        imgDat.Z = Z(ind);
        imgDat.RCP = RCP(:,ind);
        imgDat.dim = V.dim';
        imgDat.mat = V.mat;
        varargout{1} = imgDat;
        %spm_write_filtered(Z(ind),RCP(:,ind),V.dim',V.mat);
    end
end

function [Vnew,Z,RCP,factor] = afxResample1mm(fName,sz)
    Vnew.dim = round([160/sz 210/sz 160/sz]);
    Vnew.mat = [sz 0 0 -80; 0 sz 0 -120; 0 0 sz -70; 0 0 0 1];
    [R,C,P]  = ndgrid(1:Vnew(1).dim(1),1:Vnew(1).dim(2),1:Vnew(1).dim(3));
    RCP = [R(:)';C(:)';P(:)'];
    clear R C P;
    RCP(4,:) = 1;
    XYZorig = Vnew.mat*RCP;
    Vold = spm_vol(fName);
    % compute voxel coordinates in voxel space of volume wich correspond to
    % the space given by XYZ
    XYZresample = Vold.mat\XYZorig;
    % resample volume
    Z = spm_sample_vol(Vold,XYZresample(1,:),XYZresample(2,:),XYZresample(3,:),1);
    factor = prod(sqrt(sum(Vold.mat(1:3,1:3).^2,2)))/sz^3;
end

function afxTitle(tit)
    FS = spm('FontSizes');
    global hTitAx;
    cla(hTitAx);
    text(0.5,0.5,tit,'Parent',hTitAx,...
        'HorizontalAlignment','center',...
        'VerticalAlignment','top',...
        'FontWeight','Bold','FontSize',FS(14))
end

function afxText(txt)
    global hResAx;
    cla(hResAx);
    AxPos = get(hResAx,'Position'); set(hResAx,'YLim',[0,AxPos(4)])
    text(0,12,txt,'Parent',hResAx)
end