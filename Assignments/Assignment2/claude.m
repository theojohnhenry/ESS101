clc, clear, clf
load("output.mat"), load("input.mat");

% 75/25 split of training data and validation
N = length(y);
N_est = 3000; %variable parameter for split
N_val = length(y) - N_est ;

u_est = u(1:N_est); y_est = y(1:N_est);
u_val = u(N_est+1:N); y_val = y(N_est+1:N);

span_est = linspace(1,N_est,N_est); %1-3000
span_val = linspace(N_est+1,N,N_val); %3001-4000

scatter(span_est,y_est); %est data in blue
hold on
scatter(span_val, y_val); %val data in orange
hold on

% Model C: Need 3 past y values and u(t-1)
% Change from your original: PHI should be negative y terms, and only u(t-1)
PHI = [-y_est(3:N_est-1), -y_est(2:N_est-2), -y_est(1:N_est-3), u_est(3:N_est-1)];
y_target = y_est(4:N_est);  % Changed from 3:N_est to 4:N_est

thetahat = (PHI' * PHI) \ (PHI' * y_target)

%% 1-step pred. RMSE
% For validation: same structure
PHI_Val = [-y_val(3:N_val-1), -y_val(2:N_val-2), -y_val(1:N_val-3), u_val(3:N_val-1)];
y_val_target = y_val(4:N_val);  % Changed to match
y_val_pred = PHI_Val * thetahat;

rmse_pred = sqrt(mean((y_val_target - y_val_pred).^2))