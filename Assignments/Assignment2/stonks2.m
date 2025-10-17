clc; clear; close all;

N = 200;                  % number of time steps
rng(2);                   % reproducibility

% external input: interest rate shock
u = 0.01*randn(N,1);      

% smoother ARX parameters
a1 = -0.95; a2 = -0.04; b0 = 0.2;
y = zeros(N,1);
y(1:2) = [1000; 1001];    % initial index values

for t = 3:N
    y(t) = -a1*y(t-1) - a2*y(t-2) + b0*u(t) + 0.5*randn();
end

plot(y); xlabel('Time'); ylabel('Index value');
title('Synthetic stock index (ARX dynamics)');

phi1 = y(2:end-1);   
phi2 = y(1:end-2);   
phi3 = u(3:end);     
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
y_sim(1:2) = y(1:2);

for t = 3:N
    phi_t = [y_sim(t-1); y_sim(t-2); u(t)];
    y_sim(t) = thetahat' * phi_t;
end

figure;
plot(y, 'b', 'DisplayName','True index'); hold on
plot(y_sim, 'r--', 'DisplayName','Simulated ARX model');
legend; xlabel('Time'); ylabel('Index');
title('ARX simulation vs true stock index');