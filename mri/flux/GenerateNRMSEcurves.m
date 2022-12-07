function GenerateNRMSEcurves(fftpath,svtpath,optpath,outpath)
% Syntax:   GenerateNRMSEcurves(fftpath,svtpath,optpath,outpath);

% Output variables list
vars = {};

% Load FFT data
if ~isempty(fftpath)
    ltimer = tic;
    fftspec = struct('dims',1:4,'avg',4,'min',[],'data',5:9,'len',9);
    [NRMSE_fft labels_fft] = loadNRMSEdata(fftpath,fftspec); %#ok
    vars((end + 1):(end + 2)) = {'NRMSE_fft','labels_fft'};
    fprintf('FFT NRMSE curves computed [Time = %.2fs]\n',toc(ltimer));
end

% Load SVT data
if ~isempty(svtpath)
    ltimer = tic;
    %svtspec = struct('dims',1:6,'avg',4,'min',5:6,'data',7:11,'len',11);
    svtspec = struct('dims',1:6,'avg',4,'min',[],'data',7:11,'len',11);
    [NRMSE_svt labels_svt idxopt_svt] = loadNRMSEdata(svtpath,svtspec); %#ok
    vars((end + 1):(end + 3)) = {'NRMSE_svt','labels_svt','idxopt_svt'};
    fprintf('SVT NRMSE curves computed [Time = %.2fs]\n',toc(ltimer));
end

% Load OptShrink data
if ~isempty(optpath)
    ltimer = tic;
    %optspec = struct('dims',1:6,'avg',4,'min',5:6,'data',7:11,'len',11);
    optspec = struct('dims',1:6,'avg',4,'min',[],'data',7:11,'len',11);
    [NRMSE_opt labels_opt idxopt_opt] = loadNRMSEdata(optpath,optspec); %#ok
    vars((end + 1):(end + 3)) = {'NRMSE_opt','labels_opt','idxopt_opt'};
    fprintf('OptShrink NRMSE curves computed [Time = %.2fs]\n',toc(ltimer));
end

% Make output directory, if necessary
[outdir,~,~] = fileparts(outpath);
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% Save data to file
wtimer = tic;
save(outpath,vars{:});
fprintf('NRMSE curves written to "%s" [Time = %.2fs]\n',outpath,toc(wtimer));

end

function [NRMSE labels idxopt] = loadNRMSEdata(path,spec)
% Syntax:   [labels NRMSE idxopt] = loadNRMSEdata(path,spec);

% Load data
paths = findMatchingFiles(path);
data = zeros(0,spec.len);
for i = 1:length(paths)
    datai = load(paths{i});
    data = [data; datai.data]; %#ok
end
labs = datai.labels;

% Organize data
Nidxs = length(spec.dims);
inds_list = cell(1,Nidxs);
for i = 1:Nidxs
    [labels.(labs{spec.dims(i)}),~,inds_list{i}] = unique(data(:,spec.dims(i)));
end
sz = cellfun(@(f)length(labels.(f)),fields(labels));
sz = sz(:)'; % row vector
inds = sub2ind(sz,inds_list{:});
Ndata = length(spec.data);
NRMSE = cell(Ndata,1);
idxopt = cell(Ndata,1);
for i = 1:Ndata
    % Assimilate data
    NRMSE{i} = nan(sz);
    NRMSE{i}(inds) = data(:,spec.data(i));
    
    % Average over requested dimensions
    for j = 1:length(spec.avg)
        NRMSE{i} = nanmean(NRMSE{i},spec.avg(j));
    end
    
    % Minimize over requested dimensions
    [NRMSE{i} idxopt{i}] = multi_min(NRMSE{i},spec.min);
    %[NRMSE{i} idxopt{i}] = ndmin(NRMSE{i},spec.min);
end

end
