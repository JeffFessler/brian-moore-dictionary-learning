function [Xtrue, Y, M] = loadData(path,idim,T,dt,p,seed,SNR)
% Syntax: [Xtrue, Y, M] = loadData(path,idim,T,dt,p,seed);
%         [Xtrue, Y, M] = loadData(path,idim,T,dt,p,seed,SNR);

% Parse inputs
ADD_NOISE = exist('SNR','var') && ~isempty(SNR);

% Set random seed
rng(seed);

% Load data
Xtrue = loadVideo(path,idim,T,dt);
[ny, nx, nt] = size(Xtrue);

% Add noise, if requested
if ADD_NOISE && isfinite(SNR)
    sigma = 10^(-SNR / 20) * norm(Xtrue(:)) / sqrt(numel(Xtrue));
    Y = Xtrue + sigma * randn(ny,nx,nt);
else
    Y = Xtrue;
end

% Missing data
M = (rand(ny,nx,nt) > p);
Y(~M) = 0;
