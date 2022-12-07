function [NRMSE recon Xtrue Xfft mask ROIs nd splineInterp] = run_svt_algo(nc,nt,SNR,seed,lambdaL,proxS,lambdaS,algoStr,T,Xtrue_full,dce)
% Syntax:   [NRMSE recon Xtrue Xfft mask ROIs nd splineInterp] = run_svt_algo(nc,nt,SNR,seed,lambdaL,proxS,lambdaS,algoStr,T,Xtrue_full,dce);

% Load Fessler phantom
[Xtrue Xfft mask ROIs nd splineInterp Y A] = GenerateFesslerPhantom(nc,nt,SNR,seed,Xtrue_full,dce);
Xtrue_mat = masker(Xtrue,mask);
nf = size(Xtrue_full,3);

% Compute parameter scaling constants
CL = svd(Xtrue_mat); CL = CL(2);
switch lower(proxS)
    case 'soft'
        % Ell-1 regularizer
        CS = max(abs(Xtrue_mat(:))); 
    case 'mixed21'
        % Mixed 2/1 regularizer
        CS = max(sqrt(sum(abs(Xtrue_mat).^2,2)));
    otherwise
        % Unsupported sparse regularizer
        error('Sparse regularizer "%s" is not supported',proxS);
end

% Perform reconstruction
method = struct();
method.normA = get_normA(nt,nf);
method.proxL = 'SVT';
method.lambdaL = CL * lambdaL;
method.proxS = proxS;
method.lambdaS = CS * lambdaS;
method.X0 = masker(Xfft,mask);
method.dataSize = '2D'; % {'2D','3D'}
truth = struct('Xtrue',Xtrue_mat);
[NRMSE recon] = RunLSAlgo(algoStr,Y,A,T,nd,1,nt,method,truth);
