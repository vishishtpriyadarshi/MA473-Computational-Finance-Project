function error = compute_error(x, u_computed, u_exact, dim)
    error = 0;
    fprintf("\nComputing error...\n\n");
    
    start_idx = find(x > 80); start_idx = start_idx(1);
    end_idx = find(x < 120); end_idx = end_idx(end);
    N = end_idx - start_idx + 1;
        
    if dim == 1
        for i = start_idx : end_idx
            error = error + ((u_computed(i) - u_exact(i)) / u_exact(i))^2;
        end
        
        error = sqrt(error/N);
    elseif dim == 2
        for i = start_idx : end_idx
            for j = start_idx : end_idx
                error = error + ((u_computed(i, j) - u_exact(i, j)) / u_exact(i, j))^2;
            end
        end
        
        error = sqrt(error/(N^2));
    else
        for i = start_idx : end_idx
            for j = start_idx : end_idx
                for k = start_idx : end_idx
                    error = error + ((u_computed(i, j, k) - u_exact(i, j, k)) / u_exact(i, j, k))^2;
                end
            end
        end
        
        error = sqrt(error/(N^3));
    end
end
