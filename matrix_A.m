function[A] = matrix_A(Nx, x, hx, sigx, r, dt)
    ax=(r*x(3:Nx).*hx(3:Nx)-(sigx*x(3:Nx)).^2)./(hx(2:Nx-1).*(hx(2:Nx-1)+hx(3:Nx)));
    bx=1/dt+((sigx*x(2:Nx)).^2-r*x(2:Nx).*(hx(2:Nx)-hx(1:Nx-1)))./(hx(2:Nx).*hx(1:Nx-1));
    gx=(-r*x(2:Nx).*hx(1:Nx-1)-(sigx*x(2:Nx)).^2)./(hx(2:Nx).*(hx(1:Nx-1)+hx(2:Nx)));
    bx(Nx-1)=bx(Nx-1)+gx(Nx-1);
    A = diag(ax,-1)+diag(bx,0)+diag(gx(1:Nx-2),1);
end