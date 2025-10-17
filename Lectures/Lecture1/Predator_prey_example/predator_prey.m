%% Parameters
alpha = 1;
beta = 1;
delta = 1;
gamma = 1;
x1_0 = 2; % initial prey
x2_0 = 2; % initial predator

%% Derivative function
derivative = @(t, X) [X(1)*(alpha - beta*X(2));
                       X(2)*(-delta + gamma*X(1))];

%% Time setup
Nt = 1000;
tmax = 30;
t = linspace(0, tmax, Nt);
X0 = [x1_0, x2_0];

%% Solve using ode45
[t_sol, X_sol] = ode45(derivative, t, X0);
x1 = X_sol(:,1);
x2 = X_sol(:,2);

%% Plot ode45 results
figure;
plot(t_sol, x1, 'b', 'DisplayName','Prey'); hold on;
plot(t_sol, x2, 'r', 'DisplayName','Predator');
grid on;
xlabel('Time t [days]');
ylabel('Population');
title('ode45 method');
legend;

%% Euler method
Euler = @(func, X0, t) euler_solver(func, X0, t);

Xe = Euler(derivative, X0, t);

figure;
plot(t, Xe(:,1), 'b', 'DisplayName','Prey'); hold on;
plot(t, Xe(:,2), 'r', 'DisplayName','Predator');
grid on;
xlabel('Time t [s]');
ylabel('Population');
ylim([0 6]);
title('Euler method');
legend;

%% Varying beta
beta_values = [0.5 1 1.5 2];

% Prey population
figure;
hold on;
for b = beta_values
    deriv = @(t,X) [X(1)*(alpha - b*X(2)); X(2)*(-delta + gamma*X(1))];
    [t_sol, X_sol] = ode45(deriv, t, X0);
    plot(t_sol, X_sol(:,1), 'DisplayName', sprintf('\\beta=%.1f',b));
end
xlabel('Time');
ylabel('Prey population');
title('Prey population for varying \beta');
grid on;
legend;

% Predator population
figure;
hold on;
for b = beta_values
    deriv = @(t,X) [X(1)*(alpha - b*X(2)); X(2)*(-delta + gamma*X(1))];
    [t_sol, X_sol] = ode45(deriv, t, X0);
    plot(t_sol, X_sol(:,2), 'DisplayName', sprintf('\\beta=%.1f',b));
end
xlabel('Time');
ylabel('Predator population');
title('Predator population for varying \beta');
grid on;
legend;

%% Phase plane for varying initial prey
x1_initial_values = [1 2 3 5 8];
figure;
hold on;
for x1_0i = x1_initial_values
    X0i = [x1_0i, x2_0];
    [~, X_sol] = ode45(derivative, t, X0i);
    plot(X_sol(:,1), X_sol(:,2), 'DisplayName', sprintf('x1(0)=%d', x1_0i));
end
xlabel('Prey population');
ylabel('Predator population');
title('Predator vs Prey (phase plane)');
grid on;
legend;
axis equal;


%% Euler solver function
function X = euler_solver(func, X0, t)
    dt = t(2) - t(1);
    nt = length(t);
    n_vars = length(X0);
    X = zeros(nt, n_vars);
    X(1,:) = X0;
    for i = 1:nt-1
        X(i+1,:) = X(i,:) + func(t(i), X(i,:)')'*dt;
    end
end
