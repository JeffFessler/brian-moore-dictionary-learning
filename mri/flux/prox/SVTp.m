function Yhat = SVTp(Y,p,lambda)
% Syntax:   Yhat = SVTp(Y,p,lambda);
%
% Ell-p singular value shrinkage
%

% Perform ell-p singular value shrinkage
[U S V] = svd(Y,'econ');
Yhat = U * diag(pShrink(diag(S),p,lambda)) * V';
