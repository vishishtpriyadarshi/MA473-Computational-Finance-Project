function[u] = cash_or_nothing_dimension3(x0, x, hx, K, c, sigx, rhoxy, rhoyz, rhoxz, r, T, dt)
    Nt=round(T/dt);
    Nx = length(x);
    
    A = matrix_A(Nx, x, hx, sigx, r, dt);
    
    y0=x0; y=x; Ny=Nx; hy=hx; sigy=sigx;
    B = matrix_B(Ny, y, hy, sigy, r, dt);
       
    z0=x0; z=x; Nz=Nx; hz=hx; sigz=sigx;
    C = matrix_C(Nz, z, hz, sigz, r, dt);
    
    Ax=A+eye(Nx-1)*r/3; Ay=B+eye(Ny-1)*r/3; Az=C+eye(Nz-1)*r/3;
    u=zeros(Nx+1,Ny+1,Nz+1);
    
    for i=1:Nx
        for j=1:Ny
            for k=1:Nz
                if x(i)>=K && y(j)>=K && z(k)>=K
                    u(i,j,k)=c;
                end
            end
        end
    end

    u(2:Nx,2:Ny,Nz+1)=u(2:Nx,2:Ny,Nz);
    u(Nx+1,2:Ny,2:Nz+1)=u(Nx,2:Ny,2:Nz+1);
    u(2:Nx+1,Ny+1,2:Nz+1)=u(2:Nx+1,Ny,2:Nz+1);
    v=u; fx=zeros(Nx-1,1); fy=zeros(Ny-1,1); fz=zeros(Nz-1,1);

    for n=1:Nt
        for j=2:Ny
            for k=2:Nz
                for i=2:Nx
                    fx(i-1)=(1/3)*(rhoxy*sigx*sigy*x(i)*y(j)...
                    *(u(i+1,j+1,k)-u(i+1,j-1,k)-u(i-1,j+1,k)+u(i-1,j-1,k))...
                    /(hx(i-1)*hy(j)+hx(i)*hy(j)+hx(i)*hy(j-1)+hx(i-1)*hy(j-1))...
                    +rhoxz*sigx*sigz*x(i)*z(k)*(u(i+1,j,k+1)-u(i+1,j,k-1)-u(i-1,j,k+1)+u(i-1,j,k-1))...
                    /(hx(i-1)*hz(k)+hx(i)*hz(k)+hx(i)*hz(k-1)+hx(i-1)*hz(k-1))...
                    +rhoyz*sigy*sigz*y(j)*z(k)*(u(i,j+1,k+1)-u(i,j+1,k-1)-u(i,j-1,k+1)+u(i,j-1,k-1))...
                    /(hy(j-1)*hz(k)+hy(j)*hz(k)+hy(j)*hz(k-1)+hy(j-1)*hz(k-1)))+u(i,j,k)/dt;
                end
                v(2:Nx,j,k)=Ax\fx;
            end
        end
        v(2:Nx,2:Ny,Nz+1)=v(2:Nx,2:Ny,Nz);
        v(Nx+1,2:Ny,2:Nz+1)=v(Nx,2:Ny,2:Nz+1);
        v(2:Nx+1,Ny+1,2:Nz+1)=v(2:Nx+1,Ny,2:Nz+1);

        for k=2:Nz
            for i=2:Nx
                for j=2:Ny
                    fy(j-1)=(1/3)*(rhoxy*sigx*sigy*x(i)*y(j)...
                    *(v(i+1,j+1,k)-v(i+1,j-1,k)-v(i-1,j+1,k)+v(i-1,j-1,k))...
                    /(hx(i-1)*hy(j)+hx(i)*hy(j)+hx(i)*hy(j-1)+hx(i-1)*hy(j-1))...
                    +rhoxz*sigx*sigz*x(i)*z(k)*(v(i+1,j,k+1)-v(i+1,j,k-1)-v(i-1,j,k+1)+v(i-1,j,k-1))...
                    /(hx(i-1)*hz(k)+hx(i)*hz(k)+hx(i)*hz(k-1)+hx(i-1)*hz(k-1))...
                    +rhoyz*sigy*sigz*y(j)*z(k)*(v(i,j+1,k+1)-v(i,j+1,k-1)-v(i,j-1,k+1)+v(i,j-1,k-1))...
                    /(hy(j-1)*hz(k)+hy(j)*hz(k)+hy(j)*hz(k-1)+hy(j-1)*hz(k-1)))+v(i,j,k)/dt;
                end
            u(i,2:Ny,k)=Ay\fy;
            end
        end
        u(2:Nx,2:Ny,Nz+1)=u(2:Nx,2:Ny,Nz);
        u(Nx+1,2:Ny,2:Nz+1)=u(Nx,2:Ny,2:Nz+1);
        u(2:Nx+1,Ny+1,2:Nz+1)=u(2:Nx+1,Ny,2:Nz+1);

        for j=2:Ny
            for i=2:Nx
                for k=2:Nz
                    fz(k-1)=(1/3)*(rhoxy*sigx*sigy*x(i)*y(j)...
                    *(u(i+1,j+1,k)-u(i+1,j-1,k)-u(i-1,j+1,k)+u(i-1,j-1,k))...
                    /(hx(i-1)*hy(j)+hx(i)*hy(j)+hx(i)*hy(j-1)+hx(i-1)*hy(j-1))...
                    +rhoxz*sigx*sigz*x(i)*z(k)*(u(i+1,j,k+1)-u(i+1,j,k-1)- u(i-1,j,k+1)+u(i-1,j,k-1))...
                    /(hx(i-1)*hz(k)+hx(i)*hz(k)+hx(i)*hz(k-1)+hx(i-1)*hz(k-1))...
                    +rhoyz*sigy*sigz*y(j)*z(k)*(u(i,j+1,k+1)-u(i,j+1,k-1)-u(i,j-1,k+1)+u(i,j-1,k-1))...
                    /(hy(j-1)*hz(k)+hy(j)*hz(k)+hy(j)*hz(k-1)+hy(j-1)*hz(k-1)))+u(i,j,k)/dt;
                end
                v(i,j,2:Nz)=Az\fz;
            end
        end
        v(2:Nx,2:Ny,Nz+1)=v(2:Nx,2:Ny,Nz);
        v(Nx+1,2:Ny,2:Nz+1)=v(Nx,2:Ny,2:Nz+1);
        v(2:Nx+1,Ny+1,2:Nz+1)=v(2:Nx+1,Ny,2:Nz+1);
        u=v;
    end
    threeD_price=interp3(x,y,z,u(1:Nx,1:Ny,1:Nz),x0,y0,z0,'linear');
    disp(threeD_price);
end