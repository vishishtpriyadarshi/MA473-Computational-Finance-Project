function[B] = matrix_B(Ny, y, hy, sigy, r, dt)
    ay=(r*y(3:Ny).*hy(3:Ny)-(sigy*y(3:Ny)).^2)./(hy(2:Ny-1).*(hy(2:Ny-1)+hy(3:Ny)));
    by=1/dt+((sigy*y(2:Ny)).^2-r*y(2:Ny).*(hy(2:Ny)-hy(1:Ny-1)))...
    ./(hy(2:Ny).*hy(1:Ny-1));
    gy=(-r*y(2:Ny).*hy(1:Ny-1)-(sigy*y(2:Ny)).^2)./(hy(2:Ny).*(hy(1:Ny-1)+hy(2:Ny)));
    by(Ny-1)=by(Ny-1)+gy(Ny-1);
    B=diag(ay,-1)+diag(by,0)+diag(gy(1:Ny-2),1);
end