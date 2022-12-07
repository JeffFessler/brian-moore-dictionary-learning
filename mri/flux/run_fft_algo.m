function [Xtrue Xfft mask ROIs nd splineInterp] = run_fft_algo(nc,nt,SNR,seed,Xtrue_full,dce)
% Syntax:   [Xtrue Xfft mask ROIs nd splineInterp] = run_fft_algo(nc,nt,SNR,seed);
%           [Xtrue Xfft mask ROIs nd splineInterp] = run_fft_algo(nc,nt,SNR,seed,Xtrue_full,dce);

% Load Fessler phantom
[Xtrue Xfft mask ROIs nd splineInterp] = GenerateFesslerPhantom(nc,nt,SNR,seed,Xtrue_full,dce);
