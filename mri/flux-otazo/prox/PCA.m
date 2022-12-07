function Yhat = PCA(Y,r)
% Syntax:   Yhat = PCA(Y,r);
%
% PCA
%

% Perform PCA
[U S V] = svd(Y,'econ');
Yhat = U(:,1:r) * S(1:r,1:r) * V(:,1:r)';
