function Yhat = EuclideanBallProject(Y,M,eps)
% Syntax:   Yhat = EuclideanBallProject(Y,M,eps);
%
% Euclidean ball projection
%

% Orthogonal projection onto Euclidean ball
YmM = Y - M;
nrm = norm(YmM,'fro');
if (nrm > eps)
    Yhat = M + (eps / nrm) * YmM;
else
    Yhat = Y;
end
