clc; clear; close all;

N = 1000;                      % length of time series
rng(1);                       % for reproducibility

% external input: interest rate shock
u = 0.01*randn(N,1);          % small random variations

% true ARX system (hidden from "you")
a1 = -0.6; a2 = 0.2; b0 = 0.5;
y = zeros(N,1);
y(1:2) = [1000; 1002];        % initial index values

for t = 3:N
    y(t) = -a1*y(t-1) - a2*y(t-2) + b0*u(t) + 2*randn();  % system with noise
end

plot(y); xlabel('Time'); ylabel('Index value');
title('Synthetic stock index data');

phi1 = y(2:end-1);   % y(t-1)
phi2 = y(1:end-2);   % y(t-2)
phi3 = u(3:end);     % u(t)
PHI = [phi1 phi2 phi3];
Y   = y(3:end);

thetahat = (PHI' * PHI) \ (PHI' * Y);
disp('Estimated parameters:');
disp(thetahat');



y_pred = PHI * thetahat;

figure;
plot(3:N, Y, 'b', 'DisplayName','True index'); hold on
plot(3:N, y_pred, 'r--', 'DisplayName','1-step prediction');
legend; xlabel('Time'); ylabel('Index');
title('ARX prediction vs true stock index');

y_sim = zeros(N,1);
y_sim(1:2) = y(1:2);   % initialize

for t = 3:N
    phi_t = [y_sim(t-1); y_sim(t-2); u(t)];
    y_sim(t) = thetahat' * phi_t;
end

figure;
plot(y, 'b', 'DisplayName','True index'); hold on
plot(y_sim, 'r--', 'DisplayName','Simulated ARX model');
legend; xlabel('Time'); ylabel('Index');
title('ARX simulation vs true stock index');% Forecast horizon
N_forecast = N + 50;  

% Allocate forecast array
y_fore = zeros(N_forecast,1);  

% Copy known true data into the beginning
y_fore(1:N) = y;  

% Simulate beyond known horizon using only model
for t = N+1:N_forecast
    phi_t = [y_fore(t-1); y_fore(t-2); u(mod(t-1,N)+1)]; % re-use inputs cyclically
    y_fore(t) = thetahat' * phi_t;
end

% Plot: compare true vs forecast
figure;
plot(1:N, y, 'b', 'DisplayName','True index'); hold on
plot(1:N_forecast, y_fore, 'r--', 'DisplayName','Forecast');
xline(N,'k--','Forecast starts here');
xlabel('Time'); ylabel('Index');
legend; grid on;
title('ARX Forecast beyond known data');