%%

ps = 1 / 24;
SNR = inf;
seed = 1;
inpath = 'otazo_full.mat';

[Y, A, normA, Xtrue, X0] = generateCardiacPerfData(ps,SNR,seed,inpath);

nc = size(Y,4);
x = nan(1,nc);
for i = 1:nc
    x(i) = nnz(Y(1,:,1,i) ~= 0);
end

figure;

subplot(1,2,1);
imshow(double(abs(Y(:,:,1,1)) ~= 0));

subplot(1,2,2);
PMFplot(PMFgen(x),'width',0);

%%

% Knobs
f1 = 7;
f2 = 13;
gap = 1;
val = 1;

data13 = load('/Users/Brian/Desktop/data13.mat');
data18 = load('/Users/Brian/Desktop/data18.mat');
dataXX = load('otazo_full.mat');

accel13 = round(1 / data13.params.ps) %#ok
accel18 = round(1 / data18.params.ps) %#ok

nrmse13 = data13.stats.nrmse(end) %#ok
nrmse18 = data18.stats.nrmse(end) %#ok

Xtrue  = dataXX.Xtrue;
Xhat13 = reshape(data13.Lhat + data13.Shat,size(Xtrue));
Xhat18 = reshape(data18.Lhat + data18.Shat,size(Xtrue));

nrmse13b = norm(Xhat13(:) - Xtrue(:)) / norm(Xtrue(:)) %#ok
nrmse18b = norm(Xhat18(:) - Xtrue(:)) / norm(Xtrue(:)) %#ok

L = max(abs(Xtrue(:)));
Xtruef  = min(abs(Xtrue)  / L,1);
Xhat13f = min(abs(Xhat13) / L,1);
Xhat18f = min(abs(Xhat18) / L,1);

C = {Xtruef(:,:,f1), Xhat13f(:,:,f1), Xhat18f(:,:,f1); ...
     Xtruef(:,:,f2), Xhat13f(:,:,f2), Xhat18f(:,:,f2)};

 
M = cell2mov(C,gap,val);
imshow(M);
