close all;
clear;
clc;

%% Task 2 a)
syms t x u real

lambda = -2;
dT = 1e-1;
T_final = 2;
x0 = 1;
f(x, u) = lambda*x;
f_2a = matlabFunction(f);
u(t) = 0;
matlabFunction(u);

T_exact = 0:dT:T_final;
f_exact2a = x0.*exp(lambda.*T_exact);


EE.c = 0; EE.a = 0; EE.b = 1;
RK2.c = [0, 0.5]';
RK2.a = [0, 0; 0.5, 0];
RK2.b = [0, 1];
RK4.c = [0, 0.5, 0.5, 1]'; 
RK4.a = [0, 0, 0, 0;0.5, 0, 0, 0;0, 0.5, 0, 0;0, 0, 1, 0];
RK4.b = [1/6, 1/3, 1/3, 1/6];

[t2aEE, x2aEE] = RKsolver(f_2a, u, x0, T_final, dT, EE);
[t2aRK2, x2aRK2] = RKsolver(f_2a, u, x0, T_final, dT, RK2);
[t2aRK4, x2aRK4] = RKsolver(f_2a, u, x0, T_final, dT, RK4);

%% Task 2 b)

%%%% BE AWARE THAT THIS PART TAKES QUITE A WHILE %%%%%
dt = logspace(-2, -4, 10);
n_step = length(dt);
eEE = zeros(1, n_step);
eRK2 = zeros(1, n_step);
eRK4 = zeros(1, n_step);

for i = 1:n_step
    [t2bEE, x2bEE] = RKsolver(f_2a, u, x0, T_final, dt(i), EE);
    [t2bRK2, x2bRK2] = RKsolver(f_2a, u, x0, T_final, dt(i), RK2);
    [t2bRK4, x2bRK4] = RKsolver(f_2a, u, x0, T_final, dt(i), RK4);
    f_exact = x0.*exp(lambda.*(0:dt(i):T_final));
    eEE(i) = norm(f_exact - x2bEE);
    eRK2(i) = norm(f_exact - x2bRK2);
    eRK4(i) = norm(f_exact - x2bRK4);
end

figure1 = figure;
subplot(2, 1, 1)
plot(T_exact, f_exact2a)
hold on
plot(t2aEE, x2aEE, 'r+-')
plot(t2aRK2, x2aRK2, 'bo', 'MarkerSize', 6)
plot(t2aRK4, x2aRK4, 'k*', 'MarkerSize', 6)
xlim([0, 2])
xlabel('Time [s]')
ylabel('x-value')
legend('Exact value', 'Explicit Euler', 'RK2', 'RK4')
title('Comparison between exact function value and EE/RK2/RK4')

subplot(2, 1, 2) 
loglog(dt, eEE, '.-')
hold on
loglog(dt, eRK2, '.-')
loglog(dt, eRK4, '.-')
set(gca, 'XDir', 'reverse')
legend('Error EE', 'Error RK2', 'Error RK4')
xlabel('\Delta t')
ylabel('Error')
title('Accuracy of EE/RK2/RK4')


%% Task 2 c)
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

%% Task 3 a)
x = sym('x', [2, 1], 'real');
t = sym('t', [1, 1], 'real');
u_val = 5;

x0 = [1; 0];
tspan = [0, 25];

options = odeset();
[t3a, x3a] = ode45(@(t, x)VdPFun(u_val, t, x), tspan, x0, options);

figure2 = figure;
subplot(2, 1, 1)
plot(t3a, x3a(:, 1))
xlabel('Time [s]')
ylabel('x_1-value')

subplot(2, 1, 2)
plot(t3a, x3a(:, 2), 'Color', "#D95319")
xlabel('Time [s]')
ylabel('x_2-value')

sgtitle('Evaluation of VdP system using ODE45')


%% Task 3 b)
dT = 0.08;
tstop = 25;
[t3bRK, x3bRK] = RKsolver(@(x, u)VdPFun(u_val, tspan, x), u, x0, tstop, dT, RK4);

options = odeset('AbsTol', 1e-8, 'RelTol', 1e-8);
[t3bODE, x3bODE] = ode45(@(t, x)VdPFun(u_val, t, x), tspan, x0, options);

figure3 = figure;
subplot(2, 1, 1)
plot(t3bRK, x3bRK(1, :), 'k*', 'MarkerSize', 3)
hold on
plot(t3bODE, x3bODE(:, 1))
xlim([0, 25])
xlabel('Time [s]')
ylabel('x_1-value')
legend('RK4', 'MATLAB ODE45')

subplot(2, 1, 2)
plot(t3bRK, x3bRK(2, :), 'k*', 'MarkerSize', 3)
hold on
plot(t3bODE, x3bODE(:, 2))
xlim([0, 25])
xlabel('Time [s]')
ylabel('x_2-value')
legend('RK4', 'MATLAB ODE45')

sgtitle('Comparison between own RK4 method and ODE45')



%% Task 4 a)

% Can be found as IRKsolver.m


%% Task 4 b)
t_final = 1;
dT = 1e-2;
lambda = -2;

syms x u_test
f(x, u_test) = lambda*x;
f = matlabFunction(f);
x0 = 1;

IRK4.c = [(1/2) - (sqrt(3)/6), (1/2) + (sqrt(3)/6)]';
IRK4.a = [1/4, (1/4) - (sqrt(3)/6);
          (1/4) + (sqrt(3)/6), 1/4];
IRK4.b = [1/2, 1/2];

[K4b, x4b, t4b] = IRKsolver(f, u, x0, IRK4, t_final, dT, 1e-3);

figure4 = figure;
title('IRK4 v. RK4 on test system')
plot(t4b, x4b)
hold on
plot(t2aRK4, x2aRK4, 'k*', 'MarkerSize', 3)
xlim([0 1])
xlabel('Time [s]')
ylabel('x-value')
legend('IRK4', 'RK4')


%% Task 4 c)
x0 = [1, 0]';
t_final = 25;
dT = 1e-2;
tol = 1e-3;

[K4c, x4c, t4c] = IRKsolver(@(x, u)VdPFun(u_val, tspan, x), u, x0, IRK4, t_final, dT, tol);


figure5 = figure;
subplot(2, 1, 1)
plot(t4c, x4c(1, :))
hold on
plot(t3bRK, x3bRK(1, :))
xlim([0, 25])
xlabel('Time [s]')
ylabel('x_1-value')
legend('IRK4', 'RK4')

subplot(2, 1, 2)
plot(t4c, x4c(2, :))
hold on
plot(t3bRK, x3bRK(2, :))
xlim([0, 25])
xlabel('Time [s]')
ylabel('x_2-value')
legend('IRK4', 'RK4')

sgtitle('IRK4 v. RK4 on VdP system')





















