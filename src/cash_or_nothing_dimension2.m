function[u] = cash_or_nothing_dimension2(x0, x, hx, K, c, sigx, rho, r, T, dt)
    Nt = round(T/dt);
    Nx = length(x);

    A = tridiagonal_matrix(Nx, x, hx, sigx, r, dt);
    [y0, y, Ny, hy, sigy] = deal(x0, x, Nx, hx, sigx);
    B = tridiagonal_matrix(Ny, y, hy, sigy, r, dt);

    Ax = A + eye(Nx - 1) * 0.5 * r; 
    Ay = B + eye(Ny - 1) * 0.5 * r;    
    u = zeros(Nx + 1, Ny + 1);
    
    for i = 1 : Nx
        for j = 1 : Ny
            if x(i) >= K && y(j) >= K
                u(i, j) = c;
            end
        end
    end

    u(2 : Nx, Ny + 1) = u(2 : Nx, Ny); 
    u(Nx + 1, 2 : Ny + 1) = u(Nx, 2 : Ny + 1);
    [v, fx, fy] = deal(u, zeros(Nx - 1, 1), zeros(Ny - 1, 1));
    
    for n = 1 : Nt
        for j = 2 : Ny
            for i = 2 : Nx
                fx(i-1) = u(i,j)/dt+0.5*rho*sigx*sigy*x(i)*y(j)...
                *(u(i+1,j+1)+u(i-1,j-1)-u(i-1,j+1)-u(i+1,j-1))...
                /(hx(i-1)*hy(j)+hx(i)*hy(j)+hx(i)*hy(j-1)+hx(i-1)*hy(j-1));
            end
            v(2 : Nx, j) = Ax \ fx;
        end
        
        v(2 : Nx, Ny + 1) = v(2 : Nx, Ny); 
        v(Nx + 1, 2 : Ny + 1) = v(Nx, 2: Ny + 1);

        for i = 2 : Nx
            for j = 2 : Ny
                fy(j-1) = v(i,j)/dt+0.5*rho*sigx*sigy*x(i)*y(j)*...
                (v(i+1,j+1)+v(i-1,j-1)-v(i-1,j+1)-v(i+1,j-1))...
                /(hy(j-1)*hx(i)+hy(j)*hx(i)+hy(j)*hx(i-1)+hy(j-1)*hx(i-1));
            end
            u(i, 2 : Ny) = Ay \ fy;
        end
        
        u(2 : Nx, Ny + 1) = u(2 : Nx, Ny); 
        u(Nx + 1, 2 : Ny + 1) = u(Nx, 2 : Ny + 1);
    end
end