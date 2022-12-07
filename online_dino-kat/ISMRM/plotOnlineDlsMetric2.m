function plotOnlineDlsMetric2(metric,nItersi,nIters,dt,nReps)
% Syntax:   plotOnlineDlsMetric2(metric,nItersi,nIters,dt,nReps);

% Parse inputs
metric = reshape(metric(:),[],nReps);
nt = (size(metric,1) - nItersi) / nIters;

% Plot online metric
tt = 1:nIters;
phndl = zeros(1,nReps);
pstr = cell(1,nReps);
cm = linspecer(nReps);
for j = 1:nReps
    phndl(j) = plot(1:nItersi,metric(1:nItersi,j),'-','Color',cm(j,:));
    pstr{j} = sprintf('pass %d',j);
    hold on;
    for k = 1:nt
        x0 = k * dt;
        y0 = nItersi + (k - 1) * nIters;
        plot(x0 + tt,metric(y0 + tt,j),'-','Color',cm(j,:));
    end
end
%legend(phndl,pstr{:});
