function [NRMSE recon Xtrue ROIs] = run_svt_algo(p,SNR,seed,inpath,lambdaL,proxS,lambdaS,algoStr,T)
% Syntax:   [NRMSE recon Xtrue ROIs] = run_svt_algo(p,SNR,seed,inpath,lambdaL,proxS,lambdaS,algoStr,T);

% Generate Otazo data
[Xtrue X0 ROIs Y A normA] = GenerateOtazoData(p,SNR,seed,inpath);
[ny nx nt] = size(Xtrue);
nd = ny * nx;
Xtrue = reshape(Xtrue,[nd nt]);
X0 = reshape(X0,[nd nt]);

% Compute parameter scaling constants
CL = svd(Xtrue); CL = CL(2);
switch lower(proxS)
    case 'soft'
        % Ell-1 regularizer
        CS = max(abs(Xtrue(:))); 
    case 'mixed21'
        % Mixed 2/1 regularizer
        CS = max(sqrt(sum(abs(Xtrue).^2,2)));
    otherwise
        % Unsupported sparse regularizer
        error('Sparse regularizer "%s" is not supported',proxS);
end

% Perform reconstruction
method = struct();
method.normA = normA;
method.proxL = 'SVT';
method.lambdaL = CL * lambdaL;
method.proxS = proxS;
method.lambdaS = CS * lambdaS;
method.X0 = X0;
method.dataSize = '3D'; % {'2D','3D'}
truth = struct('Xtrue',Xtrue);
[NRMSE recon] = RunLSAlgo(algoStr,Y,A,T,ny,nx,nt,method,truth);
