function [K, x, t] = IRKsolver(f, u, x0, ButcherArray, T_final, dT, tol)

n_step = T_final/dT;
a = ButcherArray.a;
b = ButcherArray.b;
c = ButcherArray.c;

n_stage = length(c);
n_x = length(x0);
n_K = length(a(:, 1));

K = sym('K', [n_K, n_x], 'real');
x = sym('x', [n_x, 1], 'real');
r = sym('r', [n_stage, n_x], 'real');
syms u_s t real

if size(r) == [2 1]
    r = f(x + dT*a*K, u(t + c*dT)) - K;
else
    for i = 1:n_stage
        r(i, :) = f(x + dT*a(i, :)'*K(i, :), u(t + c(i)*dT)) - K(i, :)';
    end
end

r = reshape(r, [numel(r), 1]);
K_vec = reshape(K, [numel(K), 1]);
dr = jacobian(r, K_vec);

matlabFunction(r, dr, 'file', 'rFileIEuler', 'vars', {u_s, x, K});

t = (0:dT:T_final);
x = [x0, zeros(n_x, n_step - 1)];
K = ones(n_K, n_x);

% Loop for the Implicit Euler
for k = 1:n_step
    
    iter = true;
    while iter
        [r, dr] = rFileIEuler(u, x(:, k), K);
        K = K - reshape(dr\r, [size(K)]);

        norm(r);
        if norm(r) < tol
            iter = false;
        end
    end
    x(:, k + 1) = x(:, k) + (dT.*b*K)';
end

