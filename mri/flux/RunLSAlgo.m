function [NRMSE recon] = RunLSAlgo(algoStr,varargin)
%
% Syntax:   [NRMSE recon] = RunLSAlgo(algoStr,Y,A,T,nx,ny,nt,method,truth);
%           
% Inputs:   algoStr = {'L+S PGM','L+S ADMM','L&S ALM','L&S ADMM'}
%           
%           method.normA        % System matrix norm
%           method.proxL        % {'OptShrink','SVT','SVTp','PCA'}
%           method.r            % {'OptShrink','PCA'} only
%           method.p            % 'SVTp' only
%           method.lambdaL      % {'SVT','SVTp'}
%           method.proxS        % {'soft','mixed21'}
%           method.lambdaS      %  All
%           method.nb           % [nby nbx] locally low-rank block dims
%           method.mask         % Logical mask (Fessler data)
%           method.X0           % Initial iterate
%           method.dataSize     % {'2D','3D'}
%           
%           truth.Xtrue         % Optional
%           truth.Ltrue         % Optional
%           truth.Strue         % Optional
%           
% Outputs:  NRMSE.X
%           NRMSE.L
%           NRMSE.S
%           
%           recon.X
%           recon.L
%           recon.S
%
    
    % Run specified algorithm
    switch algoStr
        case 'L+S PGM'
            % Peform L + S via PGM
            [NRMSE recon] = Run_LpS_PGM(varargin{:});
        case 'L+S ADMM'
            % Peform L + S via ADMM
            [NRMSE recon] = Run_LpS_MRM(varargin{:});
        case 'L&S ALM'
            % Peform L & S via ALM
            [NRMSE recon] = Run_LaS_ALM(varargin{:});
        case 'L&S ADMM'
            % Peform L & S via ADMM
            [NRMSE recon] = Run_LaS_ADMM(varargin{:});
    end
end

function [NRMSE recon] = Run_LpS_PGM(Y,A,T,nx,ny,nt,method,truth)
% Perform L + S via PGM
    
    % Parse inputs
    nd = nx * ny;
    
    % Construct algo arguments
    switch upper(method.dataSize)
        case '3D'
           % 3D data
            model.A = @(X) A * reshape(X,[ny nx nt]);
            model.At = @(Y) reshape(A' * Y,[nd nt]);
            %model.T = @(X) T * reshape(X,[ny nx nt]);
            %model.Tt = @(Y) reshape(T' * Y,[nd nt]);
            model.T = @(X) T * X;
            model.Tt = @(Y) T' * Y;
        case '2D'
           % 2D data
            model.A = @(X) A * X(:);
            model.At = @(Y) reshape(A' * Y(:),[nd nt]);
            model.T = @(X) T * X;
            model.Tt = @(Y) T' * Y;
    end
    if isfield(method,'nb')
        model.nb = method.nb;
    end
    if isfield(method,'mask')
        model.mask = method.mask;
    end
    if isfield(method,'X0')
        model.M0 = method.X0;
    end
    
    % Construct optimization parameters
    opt.methodL = method.proxL;
    opt.methodS = method.proxS;
    opt.tau = 0.75 / abs(method.normA)^2;
    opt.tol = -1; % 1e-3
    opt.Nmin = 5;
    opt.Nmax = 50;
    
    % Set parameters
    model.lambdaS = method.lambdaS / opt.tau;
    switch lower(opt.methodL)
        case 'svt'
            model.lambdaL = method.lambdaL / opt.tau;
        case 'svtp'
            model.p = method.p;
            model.lambdaL = method.lambdaL / opt.tau;
        case {'optshrink','pca'}
            model.r = method.r;
    end
    
    % Run algo
    [L S NRMSE] = LpS_MRI(Y,model,opt,truth);
    recon = struct('X',L + S,'L',L,'S',S);
end

function [NRMSE recon] = Run_LpS_MRM(Y,A,T,nx,ny,nt,method,truth)
% Perform L + S via ADMM
    
    % Parse inputs
    nd = nx * ny;
    
    % Construct optimization parameters
    opt.rho = 1;
    opt.epsilon = -1; % 1e-3
    opt.Nmin = 5;
    opt.Nmax = 50;
    
    % Define low-rank regularizer
    if ~isfield(method,'mask')
        method.mask = [];
    end
    if isfield(method,'nb')
        % Locally low-rank model
        method.proxL = [method.proxL '_LLR'];
        switch lower(method.proxL)
            case 'svt_llr'
                params = {method.mask method.nb (opt.rho * method.lambdaL)};
            case 'svtp_llr'
                params = {method.mask method.nb method.p (opt.rho * method.lambdaL)};
            case {'optshrink_llr','pca_llr'}
                params = {method.mask method.nb method.r};
        end
    else
        % Standard (full) matrix model
        switch lower(method.proxL)
            case 'svt'
                params = {(opt.rho * method.lambdaL)};
            case 'svtp'
                params = {method.p (opt.rho * method.lambdaL)};
            case {'optshrink','pca'}
                params = {method.r};
        end
    end
    
    % Construct X model
    X11.prox = method.proxL;
    X11.params = params;
    X21.prox = method.proxS;
    X21.params = {(opt.rho * method.lambdaS)};
    switch upper(method.dataSize)
        case '3D'
           % 3D data
            %X21.U = @(X) T * reshape(X,[ny nx nt]);
            %X21.Ut = @(Y) reshape(T' * Y,[nd nt]);
            X21.U = @(X) T * X;
            X21.Ut = @(Y) T' * Y;
        case '2D'
           % 2D data
            X21.U = @(X) T * X;
            X21.Ut = @(Y) T' * Y;
    end
    modelX = {{X11} {X21}};
    
    % Construct algo model
    switch upper(method.dataSize)
        case '3D'
           % 3D data
            model.A = @(X) A * reshape(X,[ny nx nt]);
            model.At = @(Y) reshape(A' * Y,[nd nt]);
        case '2D'
           % 2D data
            model.A = @(X) A * X(:);
            model.At = @(Y) reshape(A' * Y(:),[nd nt]);
    end
    model.X = modelX;
    if isfield(method,'X0')
        model.X0 = method.X0;
    end
    
    % Construct GD parameters
    GDopt.k = 2;
    GDopt.normA = method.normA;
    GDopt.tau = 0.95;
    
    % Run algorithm
    if ~isempty(truth) && isfield(truth,'Xtrue')
        truth = struct('Xtrue',truth.Xtrue);
    end
    [X NRMSE] = MSM_MRI(Y,model,opt,GDopt,truth);
    %[X NRMSE] = MSMn_MRI(Y,model,opt,GDopt,truth);
    NRMSE = struct('X',NRMSE.X,'L',NRMSE.Xi{1},'S',NRMSE.Xi{2});
    recon = struct('X',X{1} + X{2},'L',X{1},'S',X{2});
end

function [NRMSE recon] = Run_LaS_ALM(Y,A,T,nx,ny,nt,method,truth)
% Perform L & S via ALM
    
    % Parse inputs
    nd = nx * ny;
    
    % Construct algo model
    switch upper(method.dataSize)
        case '3D'
           % 3D data
            model.A = @(X) A * reshape(X,[ny nx nt]);
            model.At = @(Y) reshape(A' * Y,[nd nt]);
            %model.T = @(X) T * reshape(X,[ny nx nt]);
            %model.Tt = @(Y) reshape(T' * Y,[nd nt]);
            model.T = @(X) T * X;
            model.Tt = @(Y) T' * Y;
        case '2D'
           % 2D data
            model.A = @(X) A * X(:);
            model.At = @(Y) reshape(A' * Y(:),[nd nt]);
            model.T = @(X) T * X;
            model.Tt = @(Y) T' * Y;
    end
    if isfield(method,'nb')
        model.nb = method.nb;
    end
    if isfield(method,'mask')
        model.mask = method.mask;
    end
    if isfield(method,'X0')
        model.X0 = method.X0;
    end
    
    % Construct optimization parameters
    opt.methodL = method.proxL;
    opt.methodS = method.proxS;
    opt.beta0 = 1e6;
    opt.betaInc = 25;
    opt.Ninner = 5;
    opt.Nouter = 10;
    
    % Construct GD parameters
    GDopt.k = 2;
    GDopt.normA = method.normA;
    GDopt.tau = 0.95;
    
    % Set optimization parameters
    model.lambdaS = method.lambdaS;
    switch lower(opt.methodL)
        case 'svt'
            model.lambdaL = method.lambdaL;
        case 'svtp'
            model.p = method.p;
            model.lambdaL = method.lambdaL;
        case {'optshrink','pca'}
            model.r = method.r;
    end
    
    % Run algorithm
    if ~isempty(truth) && isfield(truth,'Xtrue')
        Xtrue = truth.Xtrue;
    else
        Xtrue = [];
    end
    [X NRMSE] = LaS_MRI(Y,model,opt,GDopt,Xtrue);
    NRMSE = struct('X',NRMSE(:),'L',[],'S',[]);
    recon = struct('X',X,'L',[],'S',[]);
end

function [NRMSE recon] = Run_LaS_ADMM(Y,A,T,nx,ny,nt,method,truth)
% Perform L & S via ADMM
    
    % Parse inputs
    nd = nx * ny;
    
    % Construct optimization parameters
    opt.rho = 1;
    opt.epsilon = -1; % 1e-3
    opt.Nmin = 5;
    opt.Nmax = 50;
    
    % Define low-rank regularizer
    if ~isfield(method,'mask')
        method.mask = [];
    end
    if isfield(method,'nb')
        % Locally low-rank model
        method.proxL = [method.proxL '_LLR'];
        switch lower(method.proxL)
            case 'svt_llr'
                params = {method.mask method.nb (opt.rho * method.lambdaL)};
            case 'svtp_llr'
                params = {method.mask method.nb method.p (opt.rho * method.lambdaL)};
            case {'optshrink_llr','pca_llr'}
                params = {method.mask method.nb method.r};
        end
    else
        % Standard (full) matrix model
        switch lower(method.proxL)
            case 'svt'
                params = {(opt.rho * method.lambdaL)};
            case 'svtp'
                params = {method.p (opt.rho * method.lambdaL)};
            case {'optshrink','pca'}
                params = {method.r};
        end
    end
    
    % Construct X model
    X11.prox = method.proxL;
    X11.params = params;
    X12.prox = method.proxS;
    X12.params = {(opt.rho * method.lambdaS)};
    switch upper(method.dataSize)
        case '3D'
            % 3D data
            %X12.U = @(X) T * reshape(X,[ny nx nt]);
            %X12.Ut = @(Y) reshape(T' * Y,[nd nt]);
            X12.U = @(X) T * X;
            X12.Ut = @(Y) T' * Y;
        case '2D'
            % 2D data
            X12.U = @(X) T * X;
            X12.Ut = @(Y) T' * Y;
    end
    modelX{1} = {X11 X12};
    
    % Construct algo model
    switch upper(method.dataSize)
        case '3D'
            % 3D data
            model.A = @(X) A * reshape(X,[ny nx nt]);
            model.At = @(Y) reshape(A' * Y,[nd nt]);
        case '2D'
            % 2D data
            model.A = @(X) A * X(:);
            model.At = @(Y) reshape(A' * Y(:),[nd nt]);
    end
    model.X = modelX;
    if isfield(method,'X0')
        model.X0 = method.X0;
    end
    
    % Construct GD parameters
    GDopt.k = 2;
    GDopt.normA = method.normA;
    GDopt.tau = 0.95;
    
    % Run algorithm
    if ~isempty(truth) && isfield(truth,'Xtrue')
        truth = struct('Xtrue',truth.Xtrue);
    end
    [X NRMSE] = MSM_MRI(Y,model,opt,GDopt,truth);
    %[X NRMSE] = MSMn_MRI(Y,model,opt,GDopt,truth);
    NRMSE = struct('X',NRMSE.X(:),'L',[],'S',[]);
    recon = struct('X',X{1},'L',[],'S',[]);
end
