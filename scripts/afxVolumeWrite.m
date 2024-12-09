function fname = afxVolumeWrite(fname,y,dim,dt,mat,descrip,slopeScaling)
    % afxVolumeWrite(fname,y,dim,dt,mat,descrip,slopeScaling)
    %
    % dt ... uint8, int16*, int32, float32, float64, int8, uint16, uint32
	
    if ~exist('slopeScaling','var'), slopeScaling = true; end
    if ~exist('descrip','var'), descrip = 'afxVolumeWrite'; end
    Vo = spm_vol();
    if slopeScaling
        Vo = rmfield(Vo,'pinfo'); % enables automatic slope scaling
    else
        Vo(1).pinfo(1,1) = 1; % slope
        Vo(1).pinfo(2,1) = 0; % intercept
    end
    Vo(1).fname = fname;
    Vo(1).dim = dim;
    Vo(1).mat = mat;
    Vo(1).dt = [spm_type(dt) spm_platform('bigend')];
    %Vo(1).n = Vfunc.n;
    Vo(1).descrip = descrip;
    spm_write_vol(Vo,reshape(y,dim));
end
