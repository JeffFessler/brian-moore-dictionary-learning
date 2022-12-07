function [Lhat,Shat,Dhat,Bhat,ite,time,cost,mse,delta,sparsity] = runLADS(A,Y,Xtrue,L0,S0,Lreg,lambdaS,lambdaB,p,nIters)
% Syntax:   [Lhat,Shat,Dhat,Bhat,ite,time,cost,mse,delta,sparsity] = runLADS(A,Y,Xtrue,L0,S0,lambdaL,lambdaS,lambdaB,p,nIters);
%           [Lhat,Shat,Dhat,Bhat,ite,time,cost,mse,delta,sparsity] = runLADS(A,Y,Xtrue,L0,S0,{'svt',lambdaL},lambdaS,lambdaB,p,nIters);
%           [Lhat,Shat,Dhat,Bhat,ite,time,cost,mse,delta,sparsity] = runLADS(A,Y,Xtrue,L0,S0,{'opt',r},lambdaS,lambdaB,p,nIters);

% Parse inputs
if ~iscell(Lreg)
    % Default to SVT
    lambdaL = Lreg;
    Lreg = {'svt',lambdaL};
end

% Run LADS
param.E         = A;
param.d         = Y;
param.L0        = L0;
param.S0        = S0;
param.Lfull     = Xtrue;
param.nite      = nIters;
param.niteDL    = 1;
param.niteLS    = 5;
param.tol       = -1;
param.patchsize = [8, 8, 5];
param.slidedis  = [2, 2];
param.D         = dctmtx(prod(param.patchsize));
if strncmpi(Lreg{1},'svt',3)
    % SVT
    lambdaL = Lreg{2};
    param.lambda_L  = lambdaL;
elseif strncmpi(Lreg{1},'opt',3)
    % OptShrink
    r = Lreg{2};
    param.r  = r;
end
param.lambda_S  = lambdaS;
param.lambda_B  = lambdaB; % vector
param.p         = p;
[Lhat,Shat,Dhat,Bhat,ite,time,cost,mse,delta,sparsity] = lads_ist2(param);
