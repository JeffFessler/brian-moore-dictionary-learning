function PlayMovie(M,fps,dim)
%--------------------------------------------------------------------------
% Syntax:       PlayMovie(M);
%               PlayMovie(M,fps,width);
%               PlayMovie(M,fps,[width NaN]);
%               PlayMovie(M,fps,[NaN height]);
%               PlayMovie(M,fps,[width height]);
%               
% Inputs:       M is an (ny x nx x nt) cube containing nt video frames,
%               each with a height of ny pixels and a width of nx pixels
%               
%               [OPTIONAL] fps is the desired frames per second during
%               playback. The default value is fps = 24
%               
%               [OPTIONAL] [width height] are the desired frame dimensions
%               during playback, in pixels. When either value is NaN or is
%               omitted, the appropriate aspect ratio-preserving value is
%               computed internally
%               
% Description:  This function plays back the data in the cube M as a movie
%               
% Controls:     {Backspace,Enter,Delete}    Start/stop movie playback
%               {Left,Up} Arrow             Go back one frame
%               {Right,Down} Arrow          Advance one frame
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         February 20, 2015
%--------------------------------------------------------------------------

% Parse inputs
[ny nx nt] = size(M);
if ~exist('fps','var')
    % Default frame rate
    fps = 24;
end
if ~exist('dim','var') || all(isnan(dim(:)))
    % Default display dimensions
    dim = [nx ny];
elseif isnan(dim(1))
    % Compute aspect-preserving width
    dim(1) = round((nx / ny) * dim(2));
elseif isnan(dim(2)) || (length(dim) == 1)
    % Compute aspect-preserving height
    dim(2) = round((ny / nx) * dim(1));
end

% Scale to [0 255]
M = M - min(M(:));
M = round((255 / max(M(:))) * double(M));

% Initialize figure
scrsz = get(0,'Screensize');
fig = figure('MenuBar','None', ...
             'NumberTitle','off', ...
             'DockControl','off', ...
             'name','', ...
             'Position',[0.5 * (scrsz(3:4) - dim) dim], ...
             'KeyPressFcn',@(s,e)handleKeyPress(e), ...
             'Visible','off');
ax = axes('Position',[0 0 1 1]);
imh = image([1 dim(1)] - 0.5,[1 dim(2)] - 0.5,zeros(dim),'Parent',ax);
colormap(gray(256)); % Grayscale
axis(ax,'off');

% Initialize frame timer
timerobj = timer('ExecutionMode','FixedRate', ...
                 'StartDelay',0, ...
                 'Period',round(1000 / fps) / 1000, ...
                 'TasksToExecute',Inf, ...
                 'TimerFcn',@(s,e)displayFrame());

% Play movie
idx = 0;
lock = false;
set(fig,'Visible','on');
start(timerobj);

%--------------------------------------------------------------------------
% Nested functions
%--------------------------------------------------------------------------

% Display the next frame
function displayFrame
    % Update frame number
    idx = idx + 1;
    
    % Update display
    if ~lock && ishandle(fig)
        % Set lock
        lock = true;
        
        % Display next frame
        set(imh,'CData',M(:,:,idx));
        set(fig,'name',sprintf('Frame %i/%i',idx,nt));
        
        % Release lock
        lock = false;
    end
    
    % Check for stopping condition
    if ~ishandle(fig) || (idx >= nt)
        % Stop timer
        stop(timerobj);
        delete(timerobj);
    end
end

% Handle figure keypress
function handleKeyPress(e)
    % Parse keypress
    keyChar = e.Character;
    if isempty(keyChar)
        % Quick return
        return;
    end
    
    % Get animation state
    isRunning = strcmpi(get(timerobj,'Running'),'on');
    
    % Handle keypress
    switch double(keyChar)
        case {8 13 127}
            % Backspace/enter/delete key
            if isRunning == true
                % Pause animation
                stop(timerobj);
            else
                % Start animation
                start(timerobj);
            end
        case {28 30}
            % Left/up key
            if isRunning
                % Pause animation
                stop(timerobj);
            end
            
            % Previous frame
            idx = max(0,idx - 2);
            displayFrame();
        case {29 31}
            % Right/down key
            if isRunning
                % Pause animation
                stop(timerobj);
            end
            
            % Next frame
            displayFrame();
    end
end

end
