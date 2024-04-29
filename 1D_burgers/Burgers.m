%% Burgers eq. u_t - (eps*u_x)_x = 0

clear all
close all
prob = 1;
plot_freq = 2000000;
if prob == 1
    xl = -1;
    xr = 1;
    grids = 40*[1 2 4 8 16]+2;
    % grids = 201;
    % grids = 11;
    c = 2;
    epsilons = 0.1;
    T = 0.4;
    T0 = 0;
elseif prob == 2
    xl = 0;
    xr = 2*pi;
    grids = 201;
    c = 2;
    epsilons = [1 0.1 0.001, 0];
    % epsilons = 0;
    T = 2;
    T0 = 0;
elseif prob == 3
    xl = 0;
    xr = 2*pi;
    grids = 40*[1 2 4 8 16] +1;
    c = 2;
    epsilons = 0.5;
    T = 2;
    T0 = 0;
end


u0 = @(x,t, eps) c - tanh((x+0.5)./(2*eps));
u0sin = @(x,t) sin(x);
u_exact = @(x,t, eps) c - tanh((x+0.5-c.*t)./(2*eps));


for i = 1:length(grids)
    for k = 1:length(epsilons)
        t = 0;
        N = grids(i);
        h = (xr-xl)/(N-1);
        if prob == 3
            eps = epsilons*h;
        else
            eps = epsilons(k);
        end
        x = linspace(xl, xr, N)';
        M = mass_matrix_assembler(x);
        M = sparse(M);
        S = stiffness_matrix_assembler(x);
        S = sparse(S);
        A = advection_matrix_assembler(x);
        A = sparse(A);

        F = @(u, t) M\(-A*(0.5*u.^2) - eps*S*u);
        if prob == 1
            u = u0(x,T0, eps);
        else
            u = u0sin(x);
        end
        plot(u)       
        dt = 0.8*(((xr-xl)/(N-1))^2)/3;
        dt = T/(ceil(T/dt));
        timesteps = T/dt;
        plotnum = 1;
        while t < T
            u = rungekutta_4(u, t, dt, F);
            if prob == 1

                u(end) = u_exact(xr,t+dt, eps);
                u(1) = u_exact(xl,t+dt, eps);
            else
                u(end) = 0;
                u(1) = 0;
            end
            if mod(plotnum, plot_freq) == 0
                plot(x,u, '-x')

                if prob == 1
                    hold on
                    plot(x, u_exact(x,t,eps), '-o')
                    hold off
                    legend('sim', 'analytic')
                end

                pause(0.1)
            end
            t = t+dt;
            plotnum = plotnum +1;
        end
        if prob == 1 || prob == 3
            err(i) = norm(u-u_exact(x,t,eps),2)/(grids(i)-1);
            % err(i) =  norm(u-u_exact(x,t,eps),Inf);
            u_sol{i} = u;
        else
            err(k) = norm(u-u_exact(x,t,eps),2);
            u_sol{k} = u;
        end
    end
end

if prob == 1
    save('A1', 'err', 'u_sol', 'grids', 'xr','xl', 'epsilons')
elseif prob == 2
    save('A2', 'u_sol', 'xr', 'xl', 'epsilons', 'grids');
elseif prob == 3
    save('A3', 'u_sol', 'xr', 'xl', 'epsilons', 'grids');
end
