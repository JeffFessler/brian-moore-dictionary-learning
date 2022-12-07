%%

d = 10;
m = 5;
n = 7;

A = randn(d,n);
[D, ~] = qr(randn(d,m),0);
B = randn(m,n);

T1 = norm(A - D * B,'fro')^2 %#ok
T2 = norm(D' * A - B,'fro')^2 + trace((eye(d) - D * D') * (A * A')) %#ok

err = abs(T1 - T2) %#ok

%%

d = 10;
m = 5;
n = 7;

A = randn(d,n);
B = randn(m,n);

[U, ~, V] = svd(A * B','econ');
D = U * V';

[U, ~, V] = svd(A * B');
T = [eye(m); zeros(d - m,m)];
D2 = U * T * V';

norm(D - D2)
