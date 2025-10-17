function [t,x] = rk_solver(f, x0, u, tfinal, delta_t, BT)
    % f: function handle for the system of ODEs f(x, u)
    % x0: the initial state of the system
    % u: function handle for the input u(t)
    % tfinal: the final time
    % delta_t: the time step
    % BT: Butcher tableau (struct with fields a, b, c for Runge-Kutta coefficients)

    xk=x0;

    % Extract Butcher tableau coefficients
    a = BT.a;
    b = BT.b;
    c = BT.c;

    % Number of stages and number of state variables
    n_stage = length(c);
    n_x = length(x0);
    n_K = length(a(:, 1));

    % Number of time steps
    n_step = length(0:delta_t:tfinal);
 
    % Preallocate state array
    x = zeros(n_x, n_step);
    K = zeros(n_K, n_x);
    t = zeros(n_step, 1);
    

    for i = 1:n_step-1
        if i == 1
            x(:, i) = xk;
        end
        for n = 1:n_stage
            K(n,:) = f(x(:, i) + (delta_t*a(n, :)*K)', u(t(i) + c(n)*delta_t));
        end
        %x(:, i + 1) = x(:, i) + delta_t * sum(b .* K, 2);
        x(:, i + 1) = x(:, i) + (delta_t.*b*K)';
        t(i + 1) = t(i) + delta_t;
    end
end