function u = cash_or_nothing_dimension1(x0, x, hx, K, c, sigx, r, T, dt)
    Nt = round(T/dt);
    Nx = length(x);
    A = tridiagonal_matrix(Nx, x, hx, sigx, r, dt);

    Ax = A + eye(Nx - 1) * r;
    u = zeros(Nx, 1);
    for i = 1 : Nx
        if x(i) >= K
            u(i) = c;
        end
    end
    
    for n = 1 : Nt
        f = u(2 : Nx)/dt; 
        u(2 : Nx) = Ax \ f;
    end
end