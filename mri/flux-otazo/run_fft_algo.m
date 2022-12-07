function [NRMSE recon Xtrue ROIs] = run_fft_algo(p,SNR,seed,inpath)
% Syntax:   [NRMSE recon Xtrue ROIs] = run_fft_algo(p,SNR,seed,inpath);

% Generate Otazo data
[Xtrue, X0, ROIs, ~, ~, ~] = GenerateOtazoData(p,SNR,seed,inpath);
[ny nx nt] = size(Xtrue);
nd = ny * nx;
Xtrue = reshape(Xtrue,[nd nt]);
X0 = reshape(X0,[nd nt]);

% Return FFT reconstruction
NRMSEfcn = @(X,Xt) norm(X(:) - Xt(:)) / norm(Xt(:));
NRMSE.X = NRMSEfcn(X0,Xtrue);
recon.X = X0;
