function Yhat = mixed21(Y,lambda)
% Syntax:   Yhat = mixed21(Y,lambda);
%
% Mixed 2-1 norm
% 
% Column-form:
% \|Y\|_{2,1} = lambda(1) * \|Y(:,1)\|_2 + ... + lambda(n) * \|Y(:,n)\|_2
% 
% Row-form:
% \|Y\|_{2,1} = lambda(1) * \|Y(1,:)\|_2 + ... + lambda(m) * \|Y(m,:)\|_2
%

%{
% Sparse columns
Yhat = bsxfun(@times,max(0,1 - lambda(:)' ./ sqrt(sum(abs(Y).^2,1))),Y);
Yhat(isnan(Yhat)) = 0;
%}

% Sparse rows
Yhat = bsxfun(@times,max(0,1 - lambda(:) ./ sqrt(sum(abs(Y).^2,2))),Y);
Yhat(isnan(Yhat)) = 0;
