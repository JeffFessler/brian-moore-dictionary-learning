%%

dt = 2;
T  = 5;
nt = 50;

nt = nt - mod(nt - T,dt);

tt = 1:dt:(nt + 1 - T);

ni = numel(tt);

L = nan(ni,T);
for i = 1:ni
    t = tt(i);
    L(i,:) = t:(t + T - 1);
end

counts = histc(L(:),1:nt);

plot(1:nt,counts,'bo');
axis([1, nt, 0, max(counts)]); padAxis();

%%

gamma = 0.5;
dt = 1;
T = 5;

Xhat = zeros(3,2,10);
Y = randn(3,2,5);

for t = 1:5
    Xhat = updateRecon(Xhat,Y,gamma,t,dt);
end

actual = Xhat(:,:,5)

g = permute(gamma.^(0:4),[1, 3, 2]);

expected = sum(bsxfun(@times,g,Y),3) / sum(g)
