function [Y A Xtrue Xtrue480 Xfft mask ROIs nd splineInterp] = GenerateFesslerPhantom(nc,nt,SNR,path,seed)
% Syntax:   [Y A Xtrue Xtrue480 Xfft mask ROIs nd splineInterp] = GenerateFesslerPhantom(nc,nt,SNR);
%           [Y A Xtrue Xtrue480 Xfft mask ROIs nd splineInterp] = GenerateFesslerPhantom(nc,nt,SNR,path);
%           [Y A Xtrue Xtrue480 Xfft mask ROIs nd splineInterp] = GenerateFesslerPhantom(nc,nt,SNR,path,seed);
%
%                  Y = ny x 1
%                  A = ny x (nd * nt)
%              Xtrue = nx x ny x nt
%           Xtrue480 = nx x ny x 480
%               Xfft = nx x ny x nt
%               mask = nx x ny
%               ROIs = {nx x ny}

% Set random seed
if (exist('seed','var') && ~isempty(seed))
    rng(seed);
end

% Load dynamic object
if (exist('path','var') && ~isempty(path))
    % Load dynamic object from file
    load(path); % dyn_obj dce
else
    % Generate dynamic object from brainweb phantom
    path = '../phantom - fessler/'; % Path to .fld/.rawb files
    tic;
    [dyn_obj dce] = ir_mri_dce_obj1('chat',0,'brainweb_dir',path,'phase','');
    toc
end

% Save full truth
Xtrue480 = dyn_obj;

% Add phase to dynamic object
nx = size(dyn_obj,1);
ny = size(dyn_obj,2);
[xx yy] = ndgrid((1:nx) / nx,(1:ny) / ny);
phase =	exp((2i * pi) * (xx .* yy));
dyn_obj = dyn_obj .* repmat(single(phase), [1 1 size(dyn_obj,3)]);

% Generate k-space data
mask = (conv2(dce.labels,ones(3),'same') > 0);
[samp1 samp2] = ir_mri_dce_samp1(size(dyn_obj), ...
                                'n_tr_merge',dce.n_tr_merge, ...
                                'Nframe',nt, ...
                                'chat',1);
nd = nnz(mask);
smap = mri_sensemap_sim('nx',nx,'ny',ny,'ncoil',nc,'coil_distance',1.5);
Ytrue = ir_mri_dce_kspace1(dyn_obj,samp2,'Nframe',nt,'smap',smap,'chat',0);
Ytrue = permute(reshape(masker(Ytrue,samp1),[],nt,nc),[1 3 2]);

% Generate "nt truth"
Xtrue = reshape(dyn_obj,[nx ny (size(dyn_obj,3) / nt) nt]);
Xtrue = squeeze(mean(Xtrue,3));

% Generate observations
SNR2sigma = @(SNR,Ytrue) exp(-SNR / 20) * (norm(Ytrue(:)) / sqrt(numel(Ytrue))) / sqrt(2);
Y = Ytrue + SNR2sigma(SNR,Ytrue) * (randn(size(Ytrue)) + 1i * randn(size(Ytrue)));

% Non-iterative reconstruction
kfull = data_share_fill(reshape(permute(Y,[1 3 2]),[],nc),samp1);
xc = permute(ir_fftshift2(ifft2(ir_fftshift2(kfull))),[1 2 4 3]);
Xfft = squeeze(sum(xc .* conj(repmat(smap,[1 1 1 nt])),3)) ./ ...
       repmat(sum(abs(smap).^2,3),[1 1 nt]);
Y = Y(:);

% Spline interpolator
tt = ((1:nt) - 0.5) / nt * dce.duration_s / 60;
splineInterp = @(X) spline(tt,X,dce.ti);

% Get lesion ROIs
ROIs = cell(3,1);
for i = 1:3
    lesioni = 10 + i;
    ROIs{i} = (dce.labels == lesioni);
end     

%--------------------------------------------------------------------------
% Generate system matrix
%--------------------------------------------------------------------------
% Coil maps
tmp = cell(nc,1);
for ic = 1:nc
    diag = masker(smap(:,:,ic),mask);
    tmp{ic} = Gdiag(diag,'mask',mask);
end
Acoil = vertcat(tmp{:});
Acoil = kronI(nt,Acoil);

% Sampling pattern
tmp = cell(nt,1);
for ir = 1:nt
    tmp{ir} = Gdft('mask',mask, ...
                   'samp',samp1(:,:,ir), ...
                   'fftshift',true, ...
                   'ifftshift',true);
    tmp{ir} = kronI(nc,tmp{ir});
end
Aft = block_diag(tmp{:});
A = Aft * Acoil;
%--------------------------------------------------------------------------
