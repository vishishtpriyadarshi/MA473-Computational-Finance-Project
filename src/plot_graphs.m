function plot_graphs(x, y, z, u, plot_title, dimension)
    figure();
    
    if dimension == 1
        plot(x, u);
        zlabel("u(x, y)");
    elseif dimension == 2
        [x, y] = meshgrid(x, y);
        u = u(1:length(x), 1:length(x));
        surf(x, y, u);
    else
        [y1, z1] = meshgrid(y, z);
        u1 = squeeze(u(ceil(length(x)/2), 1:length(x), 1:length(x)));
        surf(y1, z1, u1); title(plot_title); xlabel("y"); ylabel("z"); zlabel("u(y, z)");
        
        figure();
        [x2, z2] = meshgrid(x, z);
        u2 = squeeze(u(1:length(x), ceil(length(x)/2), 1:length(x)));
        surf(x2, z2, u2); title(plot_title); xlabel("x"); ylabel("z"); zlabel("u(x, z)");
        
        figure();
        [x3, y3] = meshgrid(x, y);
        u3 = squeeze(u(1:length(x), 1:length(x), ceil(length(x)/2)));
        surf(x3, y3, u3);
    end
    
    title(plot_title);
    xlabel("x"); ylabel("y"); zlabel("u(x, y)");
    colorbar;
end
