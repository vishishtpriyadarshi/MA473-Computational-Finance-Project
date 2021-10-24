function[C] = matrix_C(Nz, z, hz, sigz, r, dt)
    az=(r*z(3:Nz).*hz(3:Nz)-(sigz*z(3:Nz)).^2)./(hz(2:Nz-1).*(hz(2:Nz-1)+hz(3:Nz)));
    bz=1/dt+((sigz*z(2:Nz)).^2-r*z(2:Nz).*(hz(2:Nz)-hz(1:Nz-1)))./(hz(2:Nz).*hz(1:Nz-1));
    gz=(-r*z(2:Nz).*hz(1:Nz-1)-(sigz*z(2:Nz)).^2)./(hz(2:Nz).*(hz(1:Nz-1)+hz(2:Nz)));
    bz(Nz-1)=bz(Nz-1)+gz(Nz-1);
    C=diag(az,-1)+diag(bz,0)+diag(gz(1:Nz-2),1);
end