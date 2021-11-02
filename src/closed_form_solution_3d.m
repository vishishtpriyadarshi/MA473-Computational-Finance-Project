function price = closed_form_solution_3d(x0, y0, z0, K, c, r, T, sigma_x, sigma_y, sigma_z, rho_xy, rho_xz, rho_yz)
    corr = [1, rho_xy, rho_xz; rho_xy , 1, rho_yz; rho_xz rho_yz, 1];
    price = zeros(length(x0), length(x0), length(x0));

    for i = 1:length(x0)
        for j = 1:length(x0)
            for k = 1:length(x0)
                dx=(log(x0(i)/K)+(r-sigma_x^2/2)*T)/(sigma_x*sqrt(T));
                dy=(log(y0(j)/K)+(r-sigma_y^2/2)*T)/(sigma_y*sqrt(T));
                dz=(log(z0(k)/K)+(r-sigma_z^2/2)*T)/(sigma_z*sqrt(T));
                price(i, j, k) = c*exp(-r*T)*mvncdf([dx,dy,dz],[0,0,0],corr);
            end
        end
    end    
end