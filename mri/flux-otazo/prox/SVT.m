function Yhat = SVT(Y,lambda)
% Syntax:   Yhat = SVT(Y,lambda);
%
% Singular value thresholding
%

% Perform singular value thresholding
[U S V] = svd(Y,'econ');
Yhat = U * diag(soft(diag(S),lambda)) * V';
