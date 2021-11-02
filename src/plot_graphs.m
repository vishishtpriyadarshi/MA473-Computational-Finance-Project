function plot_graphs(x, y, z, u, plot_title, dimension)
    figure();
    
    if dimension == 1
        plot(x, u);
        zlabel("u(x, y)");
    elseif dimension == 2
        [x, y] = meshgrid(x, y);
        u = u(1:length(x), 1:length(x));
        surf(x, y, u);
        zlabel("u(x, y)");
    else
        [x1, y1] = meshgrid(x, y);
        u = u(1:length(x1), 1:length(x1), 1:length(x1));
        surf(x1, y1, u);
        xlabel("x"); ylabel("y"); zlabel("u(x, y, z)");
        
        [y2, z2] = meshgrid(y, z);
        surf(y2, z2, u);
        xlabel("y"); ylabel("z"); zlabel("u(x, y, z)");
        
        [x3, z3] = meshgrid(x, z);
        surf(x3, z3, u);
        zlabel("u(x, y, z)");
    end
    
    title(plot_title);
    xlabel("x");ylabel("y");
    colorbar;
end
