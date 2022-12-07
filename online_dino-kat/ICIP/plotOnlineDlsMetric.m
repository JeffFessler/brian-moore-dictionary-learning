function plotOnlineDlsMetric(metric,nItersi,nIters,dt,varargin)
% Syntax:   plotOnlineDlsMetric(metric,nItersi,nIters,dt,...);

% Plot online metric
plot(1:nItersi,metric(1:nItersi),varargin{:}); hold on;
nt = (numel(metric) - nItersi) / nIters;
tt = 1:nIters;
for k = 1:nt
    x0 = k * dt;
    y0 = nItersi + (k - 1) * nIters;
    plot(x0 + tt,metric(y0 + tt),varargin{:});
end
