function [t, x] = RKsolver(f, u, x0, T_final, dT, ButcherArray)

n_step = length(0:dT:T_final);
a = ButcherArray.a;
b = ButcherArray.b;
c = ButcherArray.c;

n_stage = length(c);
n_x = length(x0);
n_K = length(a(:, 1));

K = zeros(n_K, n_x);
x = zeros(n_x, n_step);
t = zeros(n_step, 1);

for k = 1:n_step-1
    if k == 1
        x(:, k) = x0;
    end
    for n = 1:n_stage
        K(n, :) = f(x(:, k) + (dT*a(n, :)*K)', u(t(k) + c(n)*dT));
    end
    x(:, k + 1) = x(:, k) + (dT.*b*K)';
    t(k + 1) = t(k) + dT;
end