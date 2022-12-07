function [NRMSE recon Xtrue ROIs] = run_optshrink_algo(p,SNR,seed,inpath,r,proxS,lambdaS,algoStr,T)
% Syntax:   [NRMSE recon Xtrue ROIs] = run_optshrink_algo(p,SNR,seed,inpath,r,proxS,lambdaS,algoStr,T);

% Generate Otazo data
[Xtrue X0 ROIs Y A normA] = GenerateOtazoData(p,SNR,seed,inpath);
[ny nx nt] = size(Xtrue);
nd = ny * nx;
Xtrue = reshape(Xtrue,[nd nt]);
X0 = reshape(X0,[nd nt]);

% Make sure computation is valid
if (r > min(nd,nt))
    % Quick return
    NRMSE = [];
    recon = struct('X',nan(size(Xtrue)));
    return;
end

% Compute parameter scaling constants
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
method.proxL = 'OptShrink';
method.r = r;
method.proxS = proxS;
method.lambdaS = CS * lambdaS;
method.X0 = X0;
method.dataSize = '3D'; % {'2D','3D'}
truth = struct('Xtrue',Xtrue);
[NRMSE recon] = RunLSAlgo(algoStr,Y,A,T,ny,nx,nt,method,truth);
