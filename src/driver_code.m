prompt = "Enter the Number of Assets (1, 2 or 3): ";
dimension = input(prompt);

[x0, L] = deal(100, 300);
[numerical_solution, error_arr] = deal([]);

for i = 1 : 1
    fprintf("Processing Grid - %d...\n", i);
    
    % Non-uniform grid
    x = get_grid(i, L);
    hx = diff(x); hx = [hx, hx(end)];

    % Specified values
    [K, c, sigx, r] = deal(100, 100, 0.3, 0.03);
    [T, dt] = deal(1, 0.5/365);
    L = 300;

    if dimension == 1
        u_computed = cash_or_nothing_dimension1(x0, x, hx, K, c, sigx, r, T, dt);
        u_exact = closed_form_solution_1d(x', K, c, sigx, r, T);
        
        req_price = interp1(x, u_computed, x0,'linear');
        if i == 1
            exact_sol = closed_form_solution_1d(x0, K, c, sigx, r, T);
            fprintf("Exact solution, i.e, u(%d, T) = ", x0); disp(exact_sol);
        end
    elseif dimension == 2
        rho = 0.5;
        u_computed = cash_or_nothing_dimension2(x0, x, hx, K, c, sigx, rho, r, T, dt);
        u_exact = closed_form_solution_2d(x', x', K, c, sigx, sigx, r, T, rho);
        
        req_price = interp2(x, x, u_computed(1:length(x), 1:length(x)), x0, x0, 'linear');
        if i == 1
            exact_sol = closed_form_solution_2d(x0, x0, K, c, sigx, sigx, r, T, rho);
            fprintf("Exact solution, i.e, u(%d, %d, T) = ", x0, x0); disp(exact_sol);
        end
    else
        rhoxy = 0.5;
        rhoyz = 0.5;
        rhoxz = 0.5;
        u_computed = cash_or_nothing_dimension3(x0, x, hx, K, c, sigx, rhoxy, rhoyz, rhoxz, r, T, dt);
        u_exact = closed_form_solution_3d(x', x', x', K, c, r, T, sigx, sigx, sigx, rhoxy, rhoxz, rhoyz);

        req_price = interp3(x, x, x, u_computed(1:length(x), 1:length(x), 1:length(x)), x0, x0, x0, 'linear');
        if i == 1
            exact_sol = closed_form_solution_3d(x0, x0, x0, K, c, r, T, sigx, sigx, sigx, rhoxy, rhoxz, rhoyz);
            fprintf("Exact solution, i.e, u(%d, %d, %d, T) = ", x0, x0, x0); disp(exact_sol);
        end
    end

    error = compute_error(x, u_computed, u_exact, dimension);
    numerical_solution(end + 1) = req_price;
    error_arr(end + 1) = error;
    
    if i == 1
        plot_graphs(x, x, x, u_computed, "Computed Solution", dimension);
        plot_graphs(x, x, x, u_exact, "Exact Solution", dimension);
    end
end


T = [[1; 2; 3], numerical_solution', error_arr'];
T = array2table(T, 'VariableNames', {'Grid';'Numerical solution';'Error'});
disp(T);


function x = get_grid(choice, L)
    if choice == 1
        x = [0 1.5:4:77.5 80.5:3:119.5 122.5:4:L-1.5 L];
    elseif choice == 2
        x = [0 1:3:79 81:2:121 124:3:L L];
    else
        x = [0:0.5:80.5 81.5:1:120.5 122.5:2:L L];
    end
end
