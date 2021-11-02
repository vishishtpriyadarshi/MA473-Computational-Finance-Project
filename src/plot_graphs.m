function plot_graphs(x, y, u, plot_title, dimension)
    figure();
    
    if dimension == 1
        plot(x, u);
        zlabel("u(x, y)");
    elseif dimension == 2
        [x, y] = meshgrid(x, y);
        u = u(1:length(x), 1:length(x));
        surf(x, y, u);
    end
    
    title(plot_title);
    xlabel("x"); ylabel("y"); zlabel("u(x, y)");
    colorbar;
end
