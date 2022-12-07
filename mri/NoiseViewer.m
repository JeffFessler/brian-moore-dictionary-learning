function NoiseViewer(ISTn,SLRNn,dims)
% load('NoisyCardiacPerfComparison.mat');
% NoiseViewer(ISTn,SLRNn);
% NoiseViewer(ISTn,SLRNn,[33 96 33 96]);

% Plot formatting
dh = 0.03;
bottom = 3 * dh;
margin = 0.03;

% Get sizes (assumed same for ISTn and SLRNn)
[Nparam1 Nparam2] = size(ISTn.sigma);
[nx ny nt] = size(ISTn.L{1});
if (nargin < 3)
    dims = [1 nx 1 ny];
end
frame = 1;

%--------------------------------------------------------------------------
% Create figure window
%--------------------------------------------------------------------------
% Create figure
figure('name',sprintf('L + S Algorithms: ROI = [%i %i %i %i]',dims(1),dims(2),dims(3),dims(4)));

% Axes
h = tight_subplot(2,3,margin,[bottom margin],[margin margin]);

% Sliders
ISTslider1 = uicontrol('Style','slider','Min',1,'Max',Nparam1,'Value',ceil(Nparam1 / 2),'SliderStep',[1 1] / (Nparam1 - 1),'Units','Normalized','Position',[0 dh .25 dh],'Callback',@UpdateIST);
ISTslider2 = uicontrol('Style','slider','Min',1,'Max',Nparam2,'Value',ceil(Nparam2 / 2),'SliderStep',[1 1] / (Nparam2 - 1),'Units','Normalized','Position',[0 0 .25 dh],'Callback',@UpdateIST);
SLRNslider1 = uicontrol('Style','slider','Min',1,'Max',Nparam1,'Value',ceil(Nparam1 / 2),'SliderStep',[1 1] / (Nparam1 - 1),'Units','Normalized','Position',[.35 dh .25 dh],'Callback',@UpdateSLRN);
SLRNslider2 = uicontrol('Style','slider','Min',1,'Max',Nparam2,'Value',ceil(Nparam2 / 2),'SliderStep',[1 1] / (Nparam2 - 1),'Units','Normalized','Position',[.35 0 .25 dh],'Callback',@UpdateSLRN);
frameslider = uicontrol('Style','slider','Min',1,'Max',nt,'Value',1,'SliderStep',[1 1] / (nt - 1),'Units','Normalized','Position',[.7 0 0.3 dh],'Callback',@UpdateAll);

% Text boxes
ISTtext1 = uicontrol('Style','text','Units','Normalized','Position',[.25 dh .1 dh],'String','','BackgroundColor',[.8 .8 .8]);
ISTtext2 = uicontrol('Style','text','Units','Normalized','Position',[.25 0 .1 dh],'String','','BackgroundColor',[.8 .8 .8]);
SLRNtext1 = uicontrol('Style','text','Units','Normalized','Position',[.6 dh .1 dh],'String','','BackgroundColor',[.8 .8 .8]);
SLRNtext2 = uicontrol('Style','text','Units','Normalized','Position',[.6 0 .1 dh],'String','','BackgroundColor',[.8 .8 .8]);
frametext = uicontrol('Style','text','Units','Normalized','Position',[.7 dh .3 dh],'String','','BackgroundColor',[.8 .8 .8]);

% Draw everything
UpdateAll();
%--------------------------------------------------------------------------

%
% Nested functions
%

% Update everything
function UpdateAll(varargin)
    % Update frame
    frame = round(get(frameslider,'Value'));
    set(frametext,'String',sprintf('frame %i/%i',frame,nt));
    
    % Update ISTn images
    UpdateIST();
    
    % Update SLRNn images
    UpdateSLRN();
end

% Update ISTn images
function UpdateIST(varargin)
    % Get param values
    val1 = round(get(ISTslider1,'Value'));
    val2 = round(get(ISTslider2,'Value'));
    
    % Update slider labels
    set(ISTtext1,'String',sprintf('sigma = %.4f',ISTn.sigma(val1,val2)));
    set(ISTtext2,'String',sprintf('P = %.4f',ISTn.P(val1,val2)));
    
    % Get images
    L = ISTn.L{val1,val2};
    S = ISTn.S{val1,val2};
    M = L + S;
    
    % Format images
    [M L S] = FormatImages(M,L,S);
    
    % Extract ROI and frame
    M = M(dims(1):dims(2),dims(3):dims(4),frame);
    L = L(dims(1):dims(2),dims(3):dims(4),frame);
    S = S(dims(1):dims(2),dims(3):dims(4),frame);
    
    % Update images
    axes(h(1));
    imshow(M,[0 1]);
    title('ITS   [L + S]');
    axes(h(2));
    imshow(L,[0 1]);
    title('ITS   [L]');
    axes(h(3));
    imshow(S,[0 1]);
    title('ITS   [S]');
    
    % Draw everything
    drawnow;
end

% Update SLRNn images
function UpdateSLRN(varargin)
    % Get param values
    val1 = round(get(SLRNslider1,'Value'));
    val2 = round(get(SLRNslider2,'Value'));
    
    % Update slider labels
    set(SLRNtext1,'String',sprintf('sigma = %i',SLRNn.sigma(val1,val2)));
    set(SLRNtext2,'String',sprintf('P = %.4f',SLRNn.P(val1,val2)));
    
    % Get images
    L = SLRNn.L{val1,val2};
    S = SLRNn.S{val1,val2};
    M = L + S;
    
    % Format images
    [M L S] = FormatImages(M,L,S);
    
    % Extract ROI and frame
    M = M(dims(1):dims(2),dims(3):dims(4),frame);
    L = L(dims(1):dims(2),dims(3):dims(4),frame);
    S = S(dims(1):dims(2),dims(3):dims(4),frame);
    
    % Update images
    axes(h(4));
    imshow(M,[0 1]);
    title('SLRNn   [L + S]');
    axes(h(5));
    imshow(L,[0 1]);
    title('SLRNn   [L]');
    axes(h(6));
    imshow(S,[0 1]);
    title('SLRNn   [S]');
    
    % Draw everything
    drawnow;
end

end

function [M L S] = FormatImages(M,L,S)

M = abs(M);
L = abs(L);
S = abs(S);

m = max([max(M(:)) max(L(:)) max(S(:))]);

M = M / m;
L = L / m;
S = S / m;

end
