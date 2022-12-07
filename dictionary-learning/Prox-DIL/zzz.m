%% First order Lipschitz constant estimation

% B step size
if SECOND_ORDER || (it == 1)
    % Second-order
    LB = norm(D)^2;
else
    % First order
    if accel
        LB = norm(GBbar(:) - GBbarlast(:)) / norm(Bbar(:) - Bbarlast(:));
        GBbarlast = GBbar;
        Bbarlast = Bbar;
    else
        LB = norm(GB(:) - GBlast(:)) / norm(B(:) - Blast(:));
        GBlast = GB;
    end
end
tauB = tau / LB;

% B update
if accel
    % Accelerated proximal gradient
    Bbar  = B + ((tlast - 1) / t) * (B - Blast);
    Blast = B;
    GBbar = D' * (D * Bbar - Y);
    B     = shrink(Bbar - tauB * GBbar,tauB);
else
    % Proximal gradient
    Blast = B;
    GB    = D' * (D * B - Y);
    B     = shrink(B - tauB * GB,tauB);
end

% D step size
if SECOND_ORDER || (it == 1)
    % Second-order
    LD = norm(B)^2;
else
    % First order
    if accel
        LD = norm(GDbar(:) - GDbarlast(:)) / norm(Dbar(:) - Dbarlast(:));
        GDbarlast = GDbar;
        Dbarlast = Dbar;
    else
        LD = norm(GD(:) - GDlast(:)) / norm(D(:) - Dlast(:));
        GDlast = GD;
    end
end
tauD = tau / LD;

% D update
if accel
    % Accelerated proximal gradient
    Dbar   = D + ((tlast - 1) / t) * (D - Dlast);
    Dlast  = D;
    GDbar  = (Dbar * B - Y) * B';
    D = projectAtoms(Dbar - tauD * GDbar,dr,ddim,d0);
else
    % Proximal gradient
    Dlast = D;
    GD    = (D * B - Y) * B';
    D     = projectAtoms(D - tauD * GD,dr,ddim,d0);
end
