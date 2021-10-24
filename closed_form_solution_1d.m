function price = closed_form_solution_1d(x0, K, c, sigma, r, T)
    d = (log(x0/K) + (r - 0.5*sigma^2) * T) / (sigma*sqrt(T));
    normcdf = 1/2 * erfc(-d/sqrt(2));
    price = c*exp(-r*T) * normcdf; 
end
