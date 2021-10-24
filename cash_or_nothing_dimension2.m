function u = cash_or_nothing_dimension2(x, hx, K, c, sigx, rho, r, T, dt)
    [Nt, Nx] = deal(round(T/dt), length(x));
    
    % Create Tri-diagonal matrix
    ax = (r*x(3:Nx).*hx(3:Nx)-(sigx*x(3:Nx)).^2)./(hx(2:Nx-1).*(hx(2:Nx-1)+hx(3:Nx)));
    bx = 1/dt+((sigx*x(2:Nx)).^2-r*x(2:Nx).*(hx(2:Nx)-hx(1:Nx-1)))./(hx(2:Nx).*hx(1:Nx-1));
    gx = (-r*x(2:Nx).*hx(1:Nx-1)-(sigx*x(2:Nx)).^2)./(hx(2:Nx).*(hx(1:Nx-1)+hx(2:Nx)));
    bx(Nx-1) = bx(Nx-1)+gx(Nx-1);
    A = diag(ax, -1) + diag(bx, 0) + diag(gx(1:Nx-2), 1);

    y = x; Ny = Nx; hy = hx; sigy = sigx;
    ay = (r*y(3:Ny).*hy(3:Ny)-(sigy*y(3:Ny)).^2)./(hy(2:Ny-1).*(hy(2:Ny-1)+hy(3:Ny)));
    by = 1/dt+((sigy*y(2:Ny)).^2-r*y(2:Ny).*(hy(2:Ny)-hy(1:Ny-1)))...
    ./(hy(2:Ny).*hy(1:Ny-1));
    gy = (-r*y(2:Ny).*hy(1:Ny-1)-(sigy*y(2:Ny)).^2)./(hy(2:Ny).*(hy(1:Ny-1)+hy(2:Ny)));
    by(Ny-1) = by(Ny-1) + gy(Ny-1);
    B = diag(ay, -1) + diag(by, 0) + diag(gy(1:Ny-2), 1);
    
    Ax = A + eye(Nx-1)*0.5*r; 
    Ay = B + eye(Ny-1)*0.5*r;    
    u = zeros(Nx+1, Ny+1);
    
    for i = 1 : Nx
        for j = 1 : Ny
            if x(i) >= K && y(j) >= K
                u(i,j) = c;
            end
        end
    end

    u(2:Nx, Ny+1) = u(2:Nx,Ny); 
    u(Nx+1, 2:Ny+1) = u(Nx,2:Ny+1);
    v=u;    fx=zeros(Nx-1,1);   fy=zeros(Ny-1,1);
    for n=1:Nt
        for j=2:Ny
            for i=2:Nx
                fx(i-1)=u(i,j)/dt+0.5*rho*sigx*sigy*x(i)*y(j)...
                *(u(i+1,j+1)+u(i-1,j-1)-u(i-1,j+1)-u(i+1,j-1))...
                /(hx(i-1)*hy(j)+hx(i)*hy(j)+hx(i)*hy(j-1)+hx(i-1)*hy(j-1));
            end
            v(2:Nx,j)=Ax\fx;
        end
        v(2:Nx,Ny+1)=v(2:Nx,Ny); v(Nx+1,2:Ny+1)=v(Nx,2:Ny+1);

        for i=2:Nx
            for j=2:Ny
                fy(j-1)=v(i,j)/dt+0.5*rho*sigx*sigy*x(i)*y(j)*...
                (v(i+1,j+1)+v(i-1,j-1)-v(i-1,j+1)-v(i+1,j-1))...
                /(hy(j-1)*hx(i)+hy(j)*hx(i)+hy(j)*hx(i-1)+hy(j-1)*hx(i-1));
            end
            u(i,2:Ny)=Ay\fy;
        end
        u(2:Nx,Ny+1)=u(2:Nx,Ny); u(Nx+1,2:Ny+1)=u(Nx,2:Ny+1);
    end
end