%%

p = 17;
m = 5;
n = 2;
gamma = 0.9;

t = 15;

Y = randn(p,n,t);
D = randn(p,m);
B = randn(m,n,t);

psi = 0;
Qt = zeros(p,p);
Rt = zeros(p,m);
St = zeros(m,m);
for j = 1:t
    Yj = Y(:,:,j);
    Bj = B(:,:,j);
    
    psi = psi + gamma^(t - j) * norm(Yj - D * Bj,'fro')^2;
    
    Qlast = Qt;
    Rlast = Rt;
    Slast = St;
    Qt = gamma * Qt + Yj * Yj';
    Rt = gamma * Rt + Yj * Bj';
    St = gamma * St + Bj * Bj';
end

f1 = @(Q,R,S) trace(Q) - 2 * trace(D' * R) + trace(D' * D * S);
f2 = @(Q,R,S) trace(Q) + trace(D' * (D * S - 2 * R));

psi_pred1 = f1(Qt,Rt,St);
psi_pred2 = f2(Qt,Rt,St);
psi_pred3 = gamma * f1(Qlast,Rlast,Slast) + norm(Y(:,:,t) - D * B(:,:,t),'fro')^2;
psi_pred4 = gamma * f2(Qlast,Rlast,Slast) + norm(Y(:,:,t) - D * B(:,:,t),'fro')^2;

ERR1 = norm(psi - psi_pred1) %#ok
ERR2 = norm(psi - psi_pred2) %#ok
ERR3 = norm(psi - psi_pred3) %#ok
ERR4 = norm(psi - psi_pred4) %#ok

%%

dt = 2;
T = 3;
nt = 10;
gamma = 0.50;
SNR = inf;

X1 = zeros(1,1,nt);
X2 = zeros(1,1,nt);
X3 = zeros(1,1,nt);
for t = 1:dt:(nt - T + 1)
    %Xnew = permute(t:(t + T - 1),[1, 3, 2]);
    Xnew = corrupt(repmat(t,[1, 1, T]),SNR);
    X1 = updateRecon(X1,Xnew,0,t,dt);
    X2 = updateRecon(X2,Xnew,gamma,t,dt);
    X3 = updateRecon(X3,Xnew,1,t,dt);
end
x1 = squeeze(X1);
x2 = squeeze(X2);
x3 = squeeze(X3);

x1_x2_x3 = [x1, x2, x3] %#ok

figure();
phndl = zeros(1,3);
phndl(1) = plot(x1,'b-o'); hold on;
phndl(2) = plot(x2,'r-x');
phndl(3) = plot(x3,'g-*');
legend(phndl,'latest','gamma','avg');
axis tight; padAxis();
