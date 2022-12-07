function [Y, A, normA, Xtrue, X0] = generateCardiacPerfData(p,SNR,seed,data)
% Syntax:   [Y, A, normA, Xtrue, X0] = generateCardiacPerfData(p,SNR,seed,data);

% Set random seed
rng(seed);

% Load ground truth data
if ischar(data)
    % Path supplied, so load data
    data = load(data);
end
Xtrue = data.Xtrue;

% Undersample data
[Y mask] = subsampleData(data.Ytrue,p,data.ps);

% Add noise, if requested
if ~isempty(SNR) && isfinite(SNR)
    inds      = repmat(mask,[1 1 1 data.nc]);
    np        = nnz(inds);
    SNR2sigma = @(SNR,Y) exp(-SNR / 20) * (norm(Y(:)) / sqrt(numel(Y))) / sqrt(2);
    Y(inds)   = Y(inds) + SNR2sigma(SNR,Y(inds)) * (randn(np,1) + 1i * randn(np,1));
end

% System matrix
A = Emat_xyt(mask,data.samp);

% Compute norm(A)
%{
dim    = [data.ny, data.nx, data.nt];
m      = prod(dim);
AAtfcn = @(x) reshape(A' * (A * reshape(x,dim)),m,1);
opts   = struct('issym',true,'isreal',false,'maxit',25,'disp',1);
normA  = sqrt(eigs(AAtfcn,m,1,'LM',opts)) %#ok
%}
normA = data.normA; % Use fully-sampled norm(A) as upper bound

%--------------------------------------------------------------------------
% Compute initial iterate
%--------------------------------------------------------------------------
%{
% Method #1: Back-projection
X0 = A' * Y;
%}

% Method #2: FFT reconstruction
samp   = flipdim(flipdim(data.samp,1),2);
kfull  = data_share_fill(reshape(Y(abs(Y) ~= 0),[],data.nc),mask);
center = @(X) fftshift(fftshift(X,1),2);
xc     = permute(center(ifft2(center(kfull))),[1 2 4 3]);
Xfft   = squeeze(sum(xc .* conj(repmat(samp,[1 1 1 data.nt])),3)) ./ ...
         repmat(sum(abs(samp).^2,3),[1 1 data.nt]);
X0     = flipdim(flipdim(Xfft,1),2);

%{
% Method #3: Least-squares
X0 = A' * Y;
tau = 1.5 / normA^2; % Numerator in (0,2)
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
