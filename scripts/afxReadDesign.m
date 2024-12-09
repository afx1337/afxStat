function [Y,X,rowLabels,colLabels,XYZ] = afxReadDesign(designFname)
  if ~exist(designFname,'file'), error(['Could not read file ' designFname]); end
  [pth,~,~] = fileparts(designFname);
  [pth,~,~] = fileparts(pth);
  [~,~,raw] = xlsread(designFname);
  fprintf('Read design file %s with %i rows and %i columns.\n',designFname,size(raw,1),size(raw,2));
  X = cell2mat(raw(2:end,2:end));
  colLabels = raw(1,2:end);
  imgFiles = raw(2:end,1);
  for i = 1:length(imgFiles)
      curImg = imgFiles{i};
      curImg = strrep(curImg,'\',filesep());
      curImg = strrep(curImg,'/',filesep());
      if ~exist(curImg,'file')
          [~,name,ext] = fileparts(curImg);
          tmpName = fullfile(pth,'img',[name ext]);
          if exist(tmpName)
              imgFiles{i} = tmpName;
          else
              error(['Image ' [name ext] ' does not exist.']);
          end
      end
  end
  [Y.dat,XYZ,Y.dim,Y.mat] = afxVolumeRead(imgFiles);
  rowLabels = imgFiles;
end
