% Numerical solution for cash-or-nothing
function[u] = cash_or_nothing_dimension1(x0, x, hx, K, c, sigx, r, T, dt)
    % x0=100;
    % L=300;
    % x=[0 1.5:4:77.5 80.5:3:119.5 122.5:4:L-1.5 L];
    % hx=diff(x); hx=[hx,hx(end)];
    % K=100; c=100; sigx=0.3; r=0.03;
    % T=1; dt=0.5/365; 
    Nt=round(T/dt);
    Nx = length(x);
    A = matrix_A(Nx, x, hx, sigx, r, dt);

    Ax=A+eye(Nx-1)*r;
    u=zeros(Nx,1);
    for i=1:Nx
        if x(i)>=K
            u(i)=c;
        end
    end
    for n=1:Nt
        f=u(2:Nx)/dt; u(2:Nx)= Ax\f;
    end
    oneD_price=interp1(x,u,x0,'linear');
    disp(oneD_price);
end