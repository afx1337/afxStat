function afxStatGUI()
    % define toolbox directories
    afxRsDir = dir(fullfile('..','afxRs*'));
    [~,idx] = sort({afxRsDir.name});
    afxRsDir = afxRsDir(idx);
    afxRsDir = fullfile(afxRsDir(end).folder,afxRsDir(end).name);
    if ~exist(afxRsDir,'dir')
        error('No resting-state toolbox detected ...');
    else
        fprintf('afxRs: %s\n',afxRsDir)
    end
	afxStatDir = pwd();
    % datestr
    nowstr = datestr(datetime('now'),'yyyymmdd_HHMMSS');
    % prompt for mode
    stat = spm_input('Modus',1,'m',{'VLSM' 'VLSM (binär)' 'LNM (1stlevel)' 'LNSM (2nd-level)' 'LNSM (2nd-level, binär)' 'SDSM (binär)'},{'VLSM' 'VLSMbin' 'LNM' 'Perm' 'Permbin' 'SDSMbin'},1);
    stat = stat{1};
    % prompt for design file
    designFile = cellstr(spm_select(1 ,'^.*\.xlsx$','Select design file'));
    designFile = designFile{1};
    switch stat
        case {'VLSM' 'VLSMbin'}
            % get parameters
            minOverlapPct = spm_input('Minimum overlap (%)',2,'e','10',1);
            regressLesion = spm_input('Account for lesion volume?',3,'m',{'no' 'regress' 'covariate' 'direct'},{'none' 'regress' 'covariate' 'direct'},2);
            regressLesion = regressLesion{1};
            [nPerms,inference,FWE,threshVox,threshClust] = afxPermutationSetup(4);
            % perform vlsm
            addpath('scripts');
            if strcmp(stat,'VLSM')
                afxVLSM(designFile, minOverlapPct, regressLesion, nPerms, inference, FWE, threshVox, threshClust)
            else
                afxVLSMbin(designFile, minOverlapPct, regressLesion, nPerms, inference, FWE, threshVox, threshClust)
            end
            rmpath('scripts');
        case 'LNM'
            % setup
            [denoisingOptions,cohort,gmMask] = afxLNMSetup(2,afxRsDir);
            [~,denoisingOptionsN,~] = fileparts(denoisingOptions);
            % perform lesion network mapping
            [designFileP,designFileN,~] = fileparts(designFile);
            projectName = fullfile('LNM',designFileN,[denoisingOptionsN '_' nowstr]);
            [X,rowLabels,colLabels] = afxReadDesignLNSM(designFile);
            cd(afxRsDir);
            addpath('scripts');
            s = load(cohort);
            s.subjects = afxUpdatePaths(s.subjects);
            roisAll = createROIs(rowLabels,gmMask);
            [rois,delIdx] = afxCheckRois(roisAll,s.subjects);
            X = X(~delIdx,:);
            firstlevelDir = afxFirstlevel(s.subjects,denoisingOptions,rois,projectName,{'fc_wholebrain'});
            firstlevelInfo = fullfile(firstlevelDir,'firstlevel_info.mat');
            [filesLNM] = afxLnmSecondlevel(firstlevelInfo,[],true);
            filesLNM = fullfile(pwd,filesLNM);
            rmpath('scripts');
            % save LNSM design file
            [~,cohortN,~] = fileparts(cohort);
            xlswrite(fullfile(designFileP,strcat(designFileN,'_LNSM_',[cohortN '-' denoisingOptionsN],'.xlsx')),['lesionNetwork' colLabels; filesLNM num2cell(X)]);
            afxSaveCellArray(fullfile(designFileP,strcat(designFileN,'_LNSM_',[cohortN '-' denoisingOptionsN],'_exlusion.txt')),'These patients have been excluded from lesion network mapping because their lesions were entirely located in the white matter in at least on healthy subject:',rowLabels(delIdx))
        case {'Perm' 'Permbin'}
            [nPerms,inference,FWE,threshVox,threshClust] = afxPermutationSetup(4);
            allMasks = dir('masks');
            allMasks = {allMasks(~[allMasks.isdir]).name 'other'};
            maskFile = spm_input('Select mask',7,'m',allMasks,allMasks,1);
            maskFile = maskFile{1};
            if strcmp(maskFile,'other')
                maskFile = cellstr(spm_select(1 ,'image','Select mask file',{},'masks'));
                maskFile = maskFile{1};
            else
                maskFile = fullfile('masks',maskFile);
            end
            [X,rowLabels,~] = afxReadDesignLNSM(designFile);
            [designFileP,designFileN,~] = fileparts(designFile);
            if FWE, FWEstr = 'FWE'; else FWEstr = 'uncorr'; end
            destFolder = fullfile('results','Perm',designFileN,[inference '-' FWEstr],nowstr);
            addpath('scripts');
            if strcmp(stat,'Perm')
                afxStatExternal(rowLabels, [], X, {[1] [-1]}, maskFile, nPerms, inference, FWE, threshVox, threshClust, destFolder,'LNSM')
            else
                afxStatExternal(rowLabels, [], X, {[1 -1] [-1 1]}, maskFile, nPerms, inference, FWE, threshVox, threshClust, destFolder,'LNSM')
            end
            rmpath('scripts');
            copyfile(fullfile(designFileP,strcat(designFileN,'_exlusion.txt')),fullfile(destFolder,'excluded_patients.txt'));
        case {'SDSMbin'}
            % get parameters
            minOverlapPct = spm_input('Minimum overlap (%)',2,'e','-5',1);
            thr = spm_input('SDSM Threshold',3,'e','.6',1);
            [nPerms,inference,FWE,threshVox,threshClust] = afxPermutationSetup(4,true);
            [X, imgDisco, imgLesion] = afxReadDesignSDSM(designFile);
            [~,designFileN,~] = fileparts(designFile);
            if FWE, FWEstr = 'FWE'; else FWEstr = 'uncorr'; end
            destFolder = fullfile('results','SDSM',designFileN,[inference '-' FWEstr],nowstr);
            addpath('scripts');
            afxSDSM(imgDisco, imgLesion, minOverlapPct, X, [1 -1], nPerms, inference, FWE, threshVox, threshClust, destFolder, 'SDSM', thr);
            rmpath('scripts');
    end
end

function [X,rowLabels,colLabels] = afxReadDesignLNSM(designFname)
  if ~exist(designFname,'file'), error(['Could not read file ' designFname]); end
  [~,~,raw] = xlsread(designFname);
  fprintf('Read design file %s with %i rows and %i columns.\n',designFname,size(raw,1),size(raw,2));
  X = cell2mat(raw(2:end,2:end));
  colLabels = raw(1,2:end);
  rowLabels = raw(2:end,1);
end

function [X, imgDisco, imgLesion] = afxReadDesignSDSM(designFname)
  if ~exist(designFname,'file'), error(['Could not read file ' designFname]); end
  [~,~,raw] = xlsread(designFname);
  fprintf('Read design file %s with %i rows and %i columns.\n',designFname,size(raw,1),size(raw,2));
  X = cell2mat(raw(2:end,4:end));
  dirDisco = raw(2:end,1);
  dirLesion = raw(2:end,2);
  fnames = raw(2:end,3);
  imgDisco = strcat(dirDisco,fnames);
  imgLesion = strcat(dirLesion,fnames);
end

function rois = createROIs(fnames,gmMask)
    rois = struct([]);
    for iRoi = 1:length(fnames)
        [~,fname,~] = fileparts(fnames{iRoi});
        rois(end+1).name = [fname];
        rois(end).type = 'image';
        rois(end).file = fullfile(fnames{iRoi});
        rois(end).gmMask = gmMask;
    end
end

function [nPerms,inference,FWE,threshVox,threshClust] = afxPermutationSetup(start,fl)
    nPerms = spm_input('Permutations',start,'e','5000',1);
    if ~exist('fl','var')
        fl =  spm_input('Freedman-Lane',start+1,'b',{'no' 'yes'},[false true],1);
        start = start + 1;
    end
    if fl, nPerms = -nPerms; end
    inference = spm_input('Inference',start+1,'b',{'voxel' 'cluster'},[true false],1);
    if inference, inference = 'voxel'; else, inference = 'cluster'; end
    FWE = spm_input('correction',start+2,'b',{'FWE' 'uncorr'},[true false],1);
    if strcmp(inference,'cluster')
        threshVox = spm_input('p(vox)',start+3,'e','0.001',1);
        threshClust = spm_input('p(clust)',start+4,'e','0.05',1);
    else
        threshVox = spm_input('p(vox)',start+3,'e','0.05',1);
        threshClust = [];
    end
end

function [denoisingOptions,cohort,gmMask] = afxLNMSetup(start,afxRsDir)
    % obtain denoising parameters
    denoisingOptions = spm_input('Denoising strategy?',start,'m',{'WMCSF' 'GSR' 'PCA'},{fullfile('denoisingOptions','smoothed_95_WMCSF.mat') fullfile('denoisingOptions','smoothed_95_GSR.mat') fullfile('denoisingOptions','smoothed_99_PCA.mat')},1);
    denoisingOptions = denoisingOptions{1};
    % obtain normative sample to use
    d = [dir(fullfile(afxRsDir,'data','HCP','subjects','subjects*.mat'));dir(fullfile(afxRsDir,'data','NKI_sample','subjects','subjects*.mat'))];
    subjectFiles = fullfile({d.folder},{d.name});
    descr = {};
    for i = 1:length(subjectFiles)
        subjects = load(subjectFiles{i});
        subjects = subjects.subjects(~[subjects.subjects.exclude]);
        tmp = [subjects.info];
        if ~isfield(tmp,'age'), tmp(1).age = NaN; end
        if ~isfield(tmp,'sex'), tmp(1).sex = NaN; end
        if ~isfield(tmp,'handedness'), tmp(1).handedness = NaN; end
        descr{i} = sprintf(['n=%i, age=%.1f',char(177) ,'%.1f, %.1f%% male, handedness (R/A/L): %i/%i/%i'],length(subjects),mean([tmp.age]),std([tmp.age]),sum([tmp.sex])/length(subjects)*100,sum([tmp.handedness]==1),sum([tmp.handedness]==.5),sum([tmp.handedness]==0));
    end
    cohort = spm_input('Subject cohort',start+1,'m',descr,subjectFiles,1);
    cohort = cohort{1};
    gmMask = spm_input('gray matter threshold?',start+2,'e','.5',1);
end

function afxSaveCellArray(fname,heading,dat)
    fid = fopen(fname,'w');
    fprintf(fid,'%s\r\n',heading);
    for i = 1:length(dat)
        fprintf(fid,'%s\r\n',dat{i});
    end
    fclose(fid);
end
