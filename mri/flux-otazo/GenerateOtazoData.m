function [Xtrue X0 ROIs Y A normA] = GenerateOtazoData(p,SNR,seed,inpath)
% Syntax:   [Xtrue X0 ROIs Y A normA] = GenerateOtazoData(p,SNR,seed,inpath);

% Set random seed
rng(seed);

% Load ground truth data
data  = load(inpath);
Xtrue = data.Xtrue;
ROIs  = data.ROIs;

% Undersample data
[Y mask] = subsampleData(data.Ytrue,p,data.ps);

% Add noise
inds      = repmat(mask,[1 1 1 data.nc]);
np        = nnz(inds);
SNR2sigma = @(SNR,Y) exp(-SNR / 20) * (norm(Y(:)) / sqrt(numel(Y))) / sqrt(2);
Y(inds)   = Y(inds) + SNR2sigma(SNR,Y(inds)) * (randn(np,1) + 1i * randn(np,1));

% Generate system matrix
A = Emat_xyt(mask,data.samp);

% Compute norm(A)
%{
m      = data.ny * data.nx * data.nt;
AAtfcn = @(x) reshape(A' * (A * reshape(x,[ny nx nt])),m,1);
opts   = struct('issym',true,'isreal',false,'maxit',5,'disp',1);
normA  = sqrt(eigs(AAtfcn,m,1,'LM',opts)) %#ok
%}
normA  = data.normA; % Use fully-sampled norm(A) as upper bound

%--------------------------------------------------------------------------
% Compute initial iterate
%--------------------------------------------------------------------------
%{
% Method #1: Back-projection
X0 = A' * Y;
%}

% Method #2: FFT reconstruction
center = @(X) fftshift(fftshift(X,1),2);
samp   = flipdim(flipdim(data.samp,1),2);
kfull  = data_share_fill(reshape(Y(abs(Y) ~= 0),[],data.nc),mask);
xc     = permute(center(ifft2(center(kfull))),[1 2 4 3]);
Xfft   = squeeze(sum(xc .* conj(repmat(samp,[1 1 1 data.nt])),3)) ./ ...
         repmat(sum(abs(samp).^2,3),[1 1 data.nt]);
X0     = flipdim(flipdim(Xfft,1),2);

%{
% Method #3: Least-squares
X0 = A' * Y;
tau = 1.5 / norm(A)^2; % numerator in (0,2)
Niters = 20;
for i = 1:Niters
    X0last = X0;
    X0 = X0 - tau * (A' * (A * X0 - Y));
    delta = norm(X0(:) - X0last(:)) / norm(X0last(:));
    fprintf('Iteration %i/%i (delta = %.3g)\n',i,Niters,delta);
end
%}

% Scale initial iterate
X0 = X0 * data.nx;
X0 = X0 * fminsearch(@(K) norm(abs(Xtrue(:)) - abs(K * X0(:))),1);
%--------------------------------------------------------------------------
