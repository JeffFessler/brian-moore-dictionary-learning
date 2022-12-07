function [Lhat,Shat,ite,time,cost,mse,delta] = runLpS(A,Y,Xtrue,L0,S0,Lreg,lambdaS,nIters)
% Syntax:   [Lhat,Shat,ite,time,cost,mse,delta] = runLpS(A,Y,Xtrue,L0,S0,lambdaL,lambdaS,nIters);
%           [Lhat,Shat,ite,time,cost,mse,delta] = runLpS(A,Y,Xtrue,L0,S0,{'svt',lambdaL},lambdaS,nIters);
%           [Lhat,Shat,ite,time,cost,mse,delta] = runLpS(A,Y,Xtrue,L0,S0,{'opt',r},lambdaS,nIters);

% Parse inputs
if ~iscell(Lreg)
    % Default to SVT
    lambdaL = Lreg;
    Lreg = {'svt',lambdaL};
end

% Run L + S
param.E     = A;
param.d     = Y;
param.T     = TempFFT(3);
param.L0    = L0;
param.S0    = S0;
param.Lfull = Xtrue;
param.nite  = nIters;
param.tol   = -1;
if strncmpi(Lreg{1},'svt',3)
    % SVT
    lambdaL = Lreg{2};
    param.lambda_L  = lambdaL;
elseif strncmpi(Lreg{1},'opt',3)
    % OptShrink
    r = Lreg{2};
    param.r  = r;
end
param.lambda_S = lambdaS;
[Lhat,Shat,ite,time,cost,mse,delta] = lps_ist2(param);
