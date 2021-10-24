function price = closed_form_solution_2d(x0, y0, K, c, sigma_x, sigma_y, r, T, rho)
    corr = [1 rho; rho 1];
    price = zeros(length(x0));
    
    for i = 1 : length(x0)
        for j = 1 : length(y0)
            dx = (log(x0(i)/K) + (r - sigma_x^2/2)*T) / (sigma_x*sqrt(T));
            dy = (log(y0(j)/K) + (r - sigma_y^2/2)*T) / (sigma_y*sqrt(T));
            price(i, j) = c * exp(-r*T) * mvncdf([dx,dy], [0,0], corr);
        end
    end
end
