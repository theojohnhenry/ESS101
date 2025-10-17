clear all
close all
clc
%% 1.3 a)
a0 = -2.923;
b0 = 7.18;
c0 = 2.8;
mu = 0;
sigma_a = 3.8;
N = 100;
rng(1);
x = 50.*rand(N, 1);
e_a = randn(N, 1).*sigma_a^2;
y_a = a0 + b0.*x + e_a;
LSest_a = polyfit(x, y_a, 1);
xFit_a = linspace(min(x), max(x), N);
yFit_a = polyval(LSest_a, xFit_a);
figure(1);
plot(x, y_a, '.', 'MarkerSize', 10)
hold on
plot(xFit_a, yFit_a, 'r', 'LineWidth', 1)
xlabel('Random x')
ylabel('y(x)')
legend('Generated data', 'Best fit', 'Location', 'se')
%% 1.3 b)
sigma_b = 12.8;
e_b = randn(N, 1).*sigma_b^2;
y_b = a0 + b0.*x + c0.*x.^2 + e_b;
LSest_b = polyfit(x, y_b, 2);
xFit_b = linspace(min(x), max(x), N);
yFit_b = polyval(LSest_b, xFit_b);
figure(2);
plot(x, y_b, '.', 'MarkerSize', 10)
hold on
plot(xFit_b, yFit_b, 'r', 'LineWidth', 1)
xlabel('Random x')
ylabel('y(x)')
legend('Generated data', 'Best fit', 'Location', 'se')
LSest_b_lin = polyfit(x, y_b, 1);
yFit_b_lin = polyval(LSest_b_lin, xFit_b);
for i = 1:N
eps_lin(i) = (y_b(i) - (LSest_b_lin(1)*x(i) + LSest_b_lin(2)))^2;
eps(i) = (y_b(i) - (LSest_b(3) + LSest_b(2)*x(i) + LSest_b(1)*x(i)^2))^2;
end
eps_lin = sum(eps_lin)/N;
eps = sum(eps)/N;
figure(3)
plot(x, y_b, '.', 'MarkerSize', 10)
hold on
plot(xFit_b, yFit_b_lin, 'r', 'LineWidth', 1)
xlabel('Random x')
ylabel('y(x)')
legend('Generated data', 'Best fit', 'Location', 'se')
%% 2.3 - 12a
inp = load('input.mat');
inp = inp.u;
outp = load('output.mat');
outp = outp.y;
uest = inp(1:length(inp)/2);
uval = inp((length(inp)/2 + 1):end);
yest = outp(1:length(outp)/2);
yval = outp((length(outp)/2 + 1):end);
H = zeros(length(yest), 3);
H(1,:) = [0 0 0];
for i = 2:length(yest)
H(i, 1) = -yest(i - 1);
end
for i = 3:length(yest)
H(i, 2) = -yest(i - 2);
end
H(:, 3) = uest;
th = (H'*H)\(H.'*yest);
ans_12a = [th(1), th(2), th(3)]; % [a1_hat, a2_hat, b0_hat]
Ha = H;
clear th
clear H
%% 2.3 - 12b
H = zeros(length(yest), 4);
H(1,:) = [0 0 0 0];
for i = 2:length(yest)
H(i, 1) = -yest(i - 1);
end
for i = 3:length(yest)
H(i, 2) = -yest(i - 2);
end
for i = 2:length(yest)
H(i, 4) = uest(i - 1);
end
H(:, 3) = uest;
th = (H'*H)\H'*yest;
ans_12b = [th(1), th(2), th(3), th(4)]; % [a1_hat, a2_hat, b0_hat, b1_hat]
Hb = H;
clear th
clear H
%% 2.3 - 12c
H = zeros(length(yest), 4);
H(1,:) = [0 0 0 0];
for i = 2:length(yest)
H(i, 1) = -yest(i - 1);
end
for i = 3:length(yest)
H(i, 2) = -yest(i - 2);
end
for i = 4:length(yest)
H(i, 3) = -yest(i - 3);
end
for i = 2:length(uest)
H(i, 4) = uest(i - 1);
end
th = (H'*H)\H'*yest;
ans_12c = [th(1), th(2), th(3), th(4)]; % [a1_hat, a2_hat, a3_hat, b1_hat]
% Inverting the values because thats how it works aparently
% ans_12a = -ans_12a;
% ans_12b = -ans_12b;
% ans_12c = -ans_12c;
ahats = [ans_12a ans_12b ans_12c]
%% 2.3 b) - Prediction
N = length(yval);
ypred = [zeros(N, 1), zeros(N, 1), zeros(N, 1)]; % [a, b, c]
ypred(1, 1) = ans_12a(3)*uval(1);
ypred(2, 1) = ans_12a(3)*uval(2) - ans_12a(1)*yval(1);
ypred(1, 2) = ans_12b(3)*uval(1);
ypred(2, 2) = ans_12b(3)*uval(2) + ans_12b(4)*uval(1) - ans_12b(1)*yval(1);
ypred(1, 3) = 0;
ypred(2, 3) = ans_12c(4)*uval(1) - ans_12c(1)*yval(1);
ypred(3, 3) = ans_12c(4)*uval(2) - ans_12c(1)*yval(2) - ans_12c(2)*yval(1);
for i = 3:N
ypred(i, 1) = ans_12a(3)*uval(i) - ans_12a(1)*yval(i - 1) - ans_12a(2)*yval(i - 2);
ypred(i, 2) = ans_12b(3)*uval(i) + ans_12b(4)*uval(i - 1) - ans_12b(1)*yval(i - 1) - ans_12b(2)*yval(i - 2);
end
for i = 4:N
ypred(i, 3) = ans_12c(4)*uval(i - 1) - ans_12c(1)*yval(i - 1) - ans_12c(2)*yval(i - 2) - ans_12c(3)*yval(i - 3);
end
predRMSE = rmse(yval, ypred)
%% 2.3 b) - Simulation
ysim = [zeros(N, 1), zeros(N, 1), zeros(N, 1)]; % [a, b, c]
ysim(1, 1) = ans_12a(3)*uval(1);
ysim(2, 1) = ans_12a(3)*uval(2) - ans_12a(1)*ysim(1,1);
ysim(1, 2) = ans_12b(3)*uval(1);
ysim(2, 2) = ans_12b(3)*uval(2) + ans_12b(4)*uval(1) - ans_12b(1)*ysim(1, 2);
ysim(1, 3) = 0;
ysim(2, 3) = ans_12c(4)*uval(1) - ans_12c(1)*ysim(1,3);
ysim(3, 3) = ans_12c(4)*uval(2) - ans_12c(1)*ysim(2,3) - ans_12c(2)*ysim(1,3);
for i = 3:N
ysim(i, 1) = ans_12a(3)*uval(i) - ans_12a(1)*ysim(i-1, 1) - ans_12a(2)*ysim(i- 2, 1);
ysim(i, 2) = ans_12b(3)*uval(i) + ans_12b(4)*uval(i - 1) - ans_12b(1)*ysim(i-1, 2) - ans_12b(2)*ysim(i-2, 2);
end
for i = 4:N
ysim(i, 3) = ans_12c(4)*uval(i - 1) - ans_12c(1)*ysim(i-1, 3) - ans_12c(2)*ysim(i-2, 3) - ans_12c(3)*ysim(i-3, 3);
end
% Compare with real data and compute least-squares error
simERROR = yval - ysim;
simRMSE = rms(simERROR)
