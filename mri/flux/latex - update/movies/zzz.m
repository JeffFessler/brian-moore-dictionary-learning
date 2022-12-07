%% Ground truth data
%
% M = [Ground_Truth]
%

% Playback reconstruction as movie
load('ground_truth.mat');
PlayMovie(M,fps,[nan round(c * size(M,1))]);

%% FFT-based reconstruction
%
% Setup:
%   nc = 8, nt = 20, SNR = 50
%
% M = [Ground_Truth   Xhat_fft]
%

% Playback reconstruction as movie
load('FFT_recon.mat');
PlayMovie(M,fps,[nan round(c * size(M,1))]);

%% Hand-tuned OptShrink L + S
%
% Setup:
%   nc = 8, nt = 20, SNR = 50
% OptShrink:
%   r = 1, lambdaS = 0.010, NRMSE = 13.4%
%
% M = [Ground_Truth   Xhat_fft   Xhat_opt   Lhat_opt   Shat_opt]
%

% Playback reconstruction as movie
load('LpS_handopt_OptShrink.mat');
PlayMovie(M,fps,[nan round(c * size(M,1))]);

%% ROI3-Optimized OptShrink/SVT L + S
%
% Setup:
%   nc = 8, nt = 20, SNR = 50
% OptShrink:
%   r = 5, lambdaS = 0.004, NRMSE = 17.4%
% SVT:
%   lambdaL = 0.001, lambdaS = 0.004, NRMSE = 17.9%
%
% M = [Ground_Truth   Xhat_opt   Lhat_opt   Shat_opt]
%     [  Xhat_fft     Xhat_svt   Lhat_svt   Shat_svt]
%

% Playback reconstruction as movie
load('LpS_ROI3_optimized.mat');
PlayMovie(M,fps,[nan round(c * size(M,1))]);

%% X-Optimized OptShrink/SVT L + S
%
% Setup:
%   nc = 8, nt = 20, SNR = 50
% OptShrink:
%   r = 1, lambdaS = 0.720, NRMSE = 13.2%
% SVT:
%   lambdaL = 0.027, lambdaS = 0.720, NRMSE = 13.3%
%
% M = [Ground_Truth   Xhat_opt   Lhat_opt   Shat_opt]
%     [  Xhat_fft     Xhat_svt   Lhat_svt   Shat_svt]
%

% Playback reconstruction as movie
load('LpS_X_optimized.mat');
PlayMovie(M,fps,[nan round(c * size(M,1))]);
