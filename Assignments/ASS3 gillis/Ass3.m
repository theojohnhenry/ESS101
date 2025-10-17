%% Assignment 3 - ESS101 
% Gillis Magnusson
% Leonardo Barberán
close all
clear all
clc

%% 2a) Accuracy and stability
syms t x u real
lambda = -2;
dt = 1e-1;
t_final = 2;
x0 = 1;

% Defining the butcher tableaus as structs
% One for explicit euler
ee.a = 0;
ee.b = 1;
ee.c = 0;

% One for RK2
rk2.a = [0 0;1/2 0];
rk2.b = [0 1];
rk2.c = [0 ; 1/2];

% One for RK4
rk4.a = [0 0 0 0;0.5 0 0 0;0 0.5 0 0;0 0 1 0];
rk4.b = [1/6 1/3 1/3 1/6];
rk4.c = [0;1/2;1/2;1];

f(x, u) = lambda*x;
f_GRK = matlabFunction(f);
u(t) = 0;
matlabFunction(u);

%t_exact = 0:dt:t_final;
t_exact_2a = 0:dt:t_final;
f_exact_2a = x0.*exp(lambda.*t_exact_2a);

%creating a matrix of answers for each method
[t_ee_2a,x_ee_2a]=rk_solver(f_GRK,x0,u,t_final,dt,ee);
[t_rk2_2a,x_rk2_2a]=rk_solver(f_GRK,x0,u,t_final,dt,rk2);
[t_rk4_2a,x_rk4_2a]=rk_solver(f_GRK,x0,u,t_final,dt,rk4);

% figure(1)
% plot(x_ee,t_ee)
% hold on
% figure(2)
% plot(x_rk2,t_rk2)
figure(3)
plot(t_rk4_2a, x_rk4_2a)

%% 2 b)
% Plot the accuracy (i.e. the error relative to the true solution of (14)) of your integrators vs time step ∆t.
% Compare your results with the theoretical order of accuracy of the various schemes.
% Hint: a plot similar to Figure 6.5 in the Lecture Notes is useful

% using logarithmic error
dt = logspace(-2,-4,10);
n_step = length(dt);

% creating an zeros-vector to begin with
error_ee = zeros(1, n_step);
error_rk2 = zeros(1, n_step);
error_rk4 = zeros(1, n_step);

% loopingggg
for i = 1:n_step
    [t_ee_2b, x_ee_2b] =   rk_solver(f_GRK, x0, u, t_final, dt(i), ee);
    [t_rk2_2b, x_rk2_2b] = rk_solver(f_GRK, x0, u, t_final, dt(i), rk2);
    [t_rk4_2b, x_rk4_2b] = rk_solver(f_GRK, x0, u, t_final, dt(i), rk4);
    f2a_exact = x0.*exp(lambda.*(0:dt(i):t_final));
    error_EE(i) = norm(f2a_exact - x_ee_2b);
    error_RK2(i) = norm(f2a_exact - x_rk2_2b);
    error_RK4(i) = norm(f2a_exact - x_rk4_2b);
end
figure1 = figure;
subplot(2, 1, 1)
plot(t_exact_2a, f_exact_2a)
hold on
plot(t_ee_2a,x_ee_2a, 'r+-')
plot(t_rk2_2a,x_rk2_2a, 'bo', 'MarkerSize', 6)
plot(t_rk4_2a,x_rk4_2a, 'k*', 'MarkerSize', 6)
xlim([0, 2])
xlabel('Time [s]')
ylabel('x-value')
legend('Exact value', 'Explicit Euler', 'RK2', 'RK4')
title('Comparison between exact function value and EE/RK2/RK4')

subplot(2, 1, 2) 
loglog(dt, error_ee, '.-')
hold on
loglog(dt, error_rk2, '.-')
loglog(dt, error_rk4, '.-')
set(gca, 'XDir', 'reverse')
legend('Error EE', 'Error RK2', 'Error RK4')
xlabel('\Delta t')
ylabel('Error')
title('Accuracy of EE/RK2/RK4')

%% 2 c)
dt = 1e-1;
order = [1, 2, 4];
lambda = zeros(size(order));

for i = 1:length(order)
    S = 0;
    while abs(S) <=1 
        lambda(i) = lambda(i) - 0.001;
        S = 0;
        for k = 0:order(i)
            S = S + (lambda(i)*dt)^k/factorial(k);
        end 
    end
end
fprintf('Explicit Euler becomes unstable at lambda = %.3f, RK2 at %.3f and RK4 at %.3f.\n', lambda(1), lambda(2), lambda(3))
%% 3a) Van-Der-Pol Oscillator
xy0 = [1; 0];
t_final = 25;
u = 5;

% RK4 tableau
rk4.a = [0 0 0 0;0.5 0 0 0;0 0.5 0 0;0 0 1 0];
rk4.b = [1/6 1/3 1/3 1/6];
rk4.c = [0;1/2;1/2;1];

dt = 0.09;
t_time = 0:dt:t_final;

% Initialize X matrix to store results
X = zeros(2, length(t_time));  % 2 rows for x and y, columns for time steps
X(:, 1) = xy0;      % Set the initial conditions


xin = xy0;
x = sym('x', [2, 1], 'real');
t = sym('t', [1, 1], 'real');
syms U
U(t) = 0;
matlabFunction(U);

% for ODE45 
options = odeset();
[t_ode_3a, x_ode_3a] = ode45(@(t, x) vpd_system(t, x, u), t_time, xy0);

% Create a figure with subplots
figure2=figure;
% First subplot for x(t)
subplot(2, 1, 1);  % 2 rows, 1 column, position 1
plot(t_ode_3a, x_ode_3a(:,1), 'b', 'DisplayName', 'x(t)');
xlabel('Time');
ylabel('x(t)');
title('Van-der-Pol Oscillator: x(t) (ODE45)');
grid on;

% Second subplot for y(t)
subplot(2, 1, 2);  % 2 rows, 1 column, position 2
plot(t_ode_3a, x_ode_3a(:,2), 'r', 'DisplayName', 'y(t)');
xlabel('Time');
ylabel('y(t)');
title('Van-der-Pol Oscillator: y(t) (ODE45)');
grid on;



%% b) 
[t_rk4_3b, x_rk4_3b] = rk_solver(@(x, U) vpd_system(t_time, x, u), xy0, U, t_final, dt, rk4);
options = odeset('AbsTol', 1e-8, 'RelTol', 1e-8);

[t_ode_3b, x_ode_3b] = ode45(@(t, x) vpd_system(t, x, u), t_time, xy0);



% For RK4 
figure3=figure;

% First subplot for x(t)
subplot(2, 1, 1);  % 2 rows, 1 column, position 1
plot(t_rk4_3b, x_rk4_3b(1,:), 'b', 'DisplayName', 'x(t)');
hold on
plot(t_ode_3b, x_ode_3b(:,1), 'g', 'DisplayName', 'x(t)');
legend('RK4','ODE45');
xlabel('Time');
ylabel('x(t)');
title('Van-der-Pol Oscillator: x(t) using RK4 scheme');
grid on;

% Second subplot for y(t)
subplot(2, 1, 2);  % 2 rows, 1 column, position 2
plot(t_rk4_3b, x_rk4_3b(2,:), 'r', 'DisplayName', 'y(t)');
hold on
plot(t_ode_3b, x_ode_3b(:,2), 'm', 'DisplayName', 'y(t)');
legend('RK4','ODE45');
xlabel('Time');
ylabel('y(t)');
title('Van-der-Pol Oscillator: y(t) using RK4 scheme');
grid on;

%% 4a) IRK

% can be seen in IRK4solve.m

%% Task 4 b)
t_final = 1;
dt = 1e-2;
lambda = -2;
TOL = 1e-3;

syms x u f t
f(x, u) = lambda*x;
f = matlabFunction(f);
x0 = 1;
u(t) = 0;
matlabFunction(u);

% Defining the Butcher array for IRK4
irk4.A = [1/4, 1/4 - sqrt(3)/6;
     1/4 + sqrt(3)/6, 1/4];
irk4.b = [1/2, 1/2];
irk4.c = [1/2 - sqrt(3)/6, 1/2 + sqrt(3)/6]';

[K_irk4_4b, x_irk4_4b, t_irk4_4b] = IRK4solve(f, x0, u, t_final, dt, irk4, TOL);

figure4 = figure;
title('IRK4 v. RK4 on test system')
plot(t_irk4_4b, x_irk4_4b)
hold on
plot(t_rk4_2a, x_rk4_2a, 'k*', 'MarkerSize', 3)
xlim([0 1])
xlabel('Time [s]')
ylabel('x-value')
legend('IRK4', 'RK4')

%% c)

xy0 = [1;0];
t_final = 25;
dt=1e-2;
x = sym('x', [2, 1], 'real');
t = sym('t', [1, 1], 'real');
syms U
U(t) = 0;
matlabFunction(U);
u = 5;
t_time = 0:dt:t_final;
[K_irk4_4c, x_irk4_4c, t_irk4_4c] = IRK4solve(@(x, U) vpd_system(t_time, x, u), xy0, U, t_final, dt, irk4, TOL);

figure5 = figure;
title('IRK4 v. RK4 on test system')
plot(t_irk4_4c, x_irk4_4c)
hold on
plot(t_rk4_3b, x_rk4_3b, 'k*', 'MarkerSize', 3)
xlim([0 1])
xlabel('Time [s]')
ylabel('x-value')
legend('IRK4', 'RK4')




