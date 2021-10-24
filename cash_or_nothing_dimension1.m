function u = cash_or_nothing_dimension1(x, hx, K, c, sigx, r, T, dt)
    [Nt, Nx] = deal(round(T/dt), length(x));
    
    % Create Tri-diagonal matrix
    ax = (r*x(3:Nx).*hx(3:Nx)-(sigx*x(3:Nx)).^2)./(hx(2:Nx-1).*(hx(2:Nx-1)+hx(3:Nx)));
    bx = 1/dt+((sigx*x(2:Nx)).^2-r*x(2:Nx).*(hx(2:Nx)-hx(1:Nx-1)))./(hx(2:Nx).*hx(1:Nx-1));
    gx = (-r*x(2:Nx).*hx(1:Nx-1)-(sigx*x(2:Nx)).^2)./(hx(2:Nx).*(hx(1:Nx-1)+hx(2:Nx)));
    bx(Nx-1) = bx(Nx-1) + gx(Nx-1);
    A = diag(ax,-1) + diag(bx,0) + diag(gx(1:Nx-2), 1);

    Ax = A + eye(Nx - 1)*r;
    u = zeros(Nx, 1);
    
    for i = 1 : Nx
        if x(i) >= K
            u(i) = c;
        end
    end
    
    for n = 1 : Nt
        f = u(2:Nx)/dt; 
        u(2:Nx) = Ax\f;
    end
end