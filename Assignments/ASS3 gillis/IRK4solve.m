function [K, x, t] = IRK4solve(f, x0, u, tfinal, delta_t, BT, Tol)

A = BT.A;
b = BT.b;
c = BT.c;

% Number of stages and number of state variables
n_step = tfinal/delta_t;
s = length(c); % Stages (2 for IRK4 required)
o = 2 * s; % 4th-order method
nx = length(x0); % Number of states
nK = length(A(:, 1));

% Defining r(xNext, x, u): rIRK4_file
x = sym('x', [nx, 1],'real');
% Define K and r as nxs-matrices
K_mat = sym('K', [nK, nx], 'real');
r_mat = sym('r', [s, nx],'real');
syms t real;
syms u_s real


if size(r_mat) == [2 1]
    r_mat = f(x + delta_t*A*K_mat, u(t + c*delta_t)) - K_mat;
else
    for i = 1:s
        r_mat(i,:) = f(x + delta_t*A(i,:)'*K_mat(i,:), u(t + c(i)*delta_t)) - K_mat(i,:)';
    end
end

% Reshape K and r to vectors to enable computation of Jacobian
%r = reshape(r_mat, [nx * s, 1]);
r = reshape(r_mat, [numel(r_mat), 1]);
K = reshape(K_mat, [numel(K_mat), 1]);

dr = jacobian(r, K);
matlabFunction(r, dr, 'file', 'rIRK4_file', 'vars', {u_s, x, K_mat});

t = (0:delta_t:tfinal);
x = [x0, zeros(nx, n_step - 1)];
K = ones(nK, nx);


% Loop for the IRK4 method
for k = 1:n_step
    % Newton iteration
    iter = true;
    %K = ones(nK, nx); % (for lambda)
    %K = repelem(x(:, k), s);  % Initialize K with the current state (Just for vanderpol)
    while iter
        [r, dr] = rIRK4_file(u,x(:,k),K);  % Compute r and dr (vanderpol))
        
        
        % Update K using Newton's method
        K = K - reshape(dr\r, [size(K)]);
        norm(r);    
        % Check convergence
        if norm(r) < Tol
            iter = false;
        end
    end
    
    x(:, k + 1) = x(:, k) + (delta_t.*b*K)';
end