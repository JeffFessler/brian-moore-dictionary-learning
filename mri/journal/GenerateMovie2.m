function GenerateMovie2(X,outpath,fps,watermark,res,lims)
%--------------------------------------------------------------------------
% Syntax:       GenerateMovie2(X,outpath,fps,watermark);
%               GenerateMovie2(X,outpath,fps,watermark,res);
%               GenerateMovie2(X,outpath,fps,watermark,res,lims);
%               
% Inputs:       X is an (ny x nx x nt) cube containing nt video frames
%               
%               outpath is the desired output filename/path (without
%               extension)
%               
%               fps is the desired frames per second for the movie
%               
%               watermark is a string (or cell array of strings) to display
%               as a watermark through the image
%               
%               [OPTIONAL] res is the desired vertical resolution, in
%               pixels. The default value (or when res = []) is ny
%               
%               [OPTIONAL] lims = [minlim maxlim] contains the desired
%               limits for image display. The default value (or when lims 
%               = []) is lims = [min(X(:)) max(X(:))]
%               
% Description:  This function generates an .avi video of the input data X
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         August 25, 2014
%--------------------------------------------------------------------------

% HACK ALERT: Must use "software" OpenGL for getframe() to work...
opengl software;

% Watermark knobs
fontName = 'Helvetica';
fontWeight = 'normal';
fontSize = 0.1;
color = [1 1 1];
angle = 45;
alpha = 0.85;

% Video formats
%format = 'Motion JPEG AVI';
format = 'Uncompressed AVI';
%format = 'Archival';
%format = 'Motion JPEG 2000';
%format = 'MPEG-4';

% Parse inputs
if ischar(watermark)
    watermark = {watermark};
end
if ~exist('res','var') || isempty(res)
    % Default resolution
    res = size(X,1);
end
[ny nx] = size(imresize(X(:,:,1),[res nan],'nearest'));
nt = size(X,3);
if ~exist('lims','var') || isempty(lims)
    % Default display limits
    lims = [min(X(:)) max(X(:))];
end

% Create video object (could also use avifile() function)
vidObj = VideoWriter(outpath,format);
vidObj.FrameRate = fps;

% Open the video file
vidObj.open();

% Set up a nice figure for recording
scrsz = get(0,'ScreenSize');
xyc = 0.5 * scrsz(3:4);
dim = [nx ny];
f = figure('Name','', ...
           'MenuBar','None', ...
           'NumberTitle','off', ...
           'Position',[(xyc - 0.5 * dim) dim], ...
           'Renderer','OpenGL', ...
           'Resize','off');
ax(1) = axes('Position',[0 0 1 1]);
ax(2) = axes('Position',[0 0 1 1]);
ax(3) = axes('Position',[0 0 1 1]);
axis(ax,'off');
uistack(ax(1),'bottom');
uistack(ax(3),'top');

% Add watermark
text(0.5,0.5,watermark,'Parent',ax(2), ...
                       'Rotation',angle, ...
                       'HorizontalAlignment','center', ...
                       'VerticalAlignment','middle', ...
                       'FontUnits','normalized', ...
                       'FontSize',fontSize, ...
                       'FontWeight',fontWeight, ...
                       'FontName',fontName, ...
                       'Color',color);
axis(ax(2),[0 1 0 1]);

% Record the movie
for k = 1:nt
    % Format the kth frame
    Xk = imresize(X(:,:,k),[ny nx],'nearest');
    Xk(Xk < lims(1)) = lims(1);
    Xk(Xk > lims(2)) = lims(2);
    Xk = Xk - lims(1);
    Xk = im2uint8(Xk / lims(2));
    
    % Display the kth frame
    colormap(gray(256));
    image([1 nx],[1 ny],Xk,'Parent',ax(1));
    image([1 nx],[1 ny],Xk,'AlphaData',alpha,'Parent',ax(3));
    axis(ax([1 3]),[1 nx 1 ny]);
    
    % Format axes for printing
    set(ax,'Position',[0 0 1 1]);
    axis(ax,'off');
    
    % Write frame to file
    drawnow; pause(0.01);
    frame = getframe(f);
    vidObj.writeVideo(frame);
end
close(f);

% Close the .avi file
vidObj.close();

% Return to OpenGL hardware 
opengl hardware;
