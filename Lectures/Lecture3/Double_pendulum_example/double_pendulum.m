% Double Pendulum Simulation in MATLAB
clear; clc; close all;

% Pendulum parameters
L1 = 1; L2 = 1;   % lengths (m)
m1 = 1; m2 = 1;   % masses (kg)
g = 9.81;         % gravitational acceleration (m/s^2)

% Initial conditions: [theta1, dtheta1, theta2, dtheta2]
y0 = [3*pi/7; 0; 3*pi/4; 0];

% Time setup
tmax = 30; dt = 0.01;
tspan = 0:dt:tmax;

% Solve ODE
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
[time, y] = ode45(@(t,y) deriv(t,y,L1,L2,m1,m2,g), tspan, y0, opts);

theta1 = y(:,1);
theta2 = y(:,3);

% Cartesian coordinates
x1 = L1*sin(theta1);
y1 = -L1*cos(theta1);
x2 = x1 + L2*sin(theta2);
y2 = y1 - L2*cos(theta2);

% Energy drift check
E0 = calc_E(y0,L1,L2,m1,m2,g);
E = arrayfun(@(k) calc_E(y(k,:),L1,L2,m1,m2,g), 1:length(y));
EDRIFT = 0.05;
if max(abs(E - E0)) > EDRIFT
    error('Maximum energy drift exceeded.');
end

% Plot settings
r = 0.05; % bob radius
fps = 10;
di = round(1/(fps*dt));
trail_secs = 1;
max_trail = round(trail_secs/dt);

figure('Color','w');
for i = 1:di:length(time)
    clf; hold on;

    % Draw rods
    plot([0 x1(i) x2(i)], [0 y1(i) y2(i)], 'k-', 'LineWidth', 2);

    % Draw anchor and bobs
    viscircles([0 0], r/2,'Color','k','LineWidth',1,'EnhanceVisibility',false);
    viscircles([x1(i) y1(i)], r,'Color','b','LineWidth',1,'EnhanceVisibility',false);
    viscircles([x2(i) y2(i)], r,'Color','r','LineWidth',1,'EnhanceVisibility',false);

    % Trail fading effect
    ns = 20;
    s = floor(max_trail/ns);
    for j = 1:ns
        imin = i - (ns-j)*s;
        if imin < 1, continue; end
        imax = min(imin+s, length(x2));
        alpha = (j/ns)^2; % fading
        plot(x2(imin:imax), y2(imin:imax), 'r-', 'LineWidth', 2, ...
             'Color', [1 0 0 alpha]);
    end

    % Formatting
    axis equal off;
    xlim([-L1-L2-r, L1+L2+r]);
    ylim([-L1-L2-r, L1+L2+r]);

    drawnow;
end


% ================= Helper Functions ================= %
function dydt = deriv(~, y, L1, L2, m1, m2, g)
    % Unpack
    theta1 = y(1); z1 = y(2);
    theta2 = y(3); z2 = y(4);

    c = cos(theta1-theta2);
    s = sin(theta1-theta2);

    theta1dot = z1;
    z1dot = (m2*g*sin(theta2)*c - m2*s*(L1*z1^2*c + L2*z2^2) - ...
            (m1+m2)*g*sin(theta1)) / L1 / (m1 + m2*s^2);
    theta2dot = z2;
    z2dot = ((m1+m2)*(L1*z1^2*s - g*sin(theta2) + g*sin(theta1)*c) + ...
            m2*L2*z2^2*s*c) / L2 / (m1 + m2*s^2);

    dydt = [theta1dot; z1dot; theta2dot; z2dot];
end

function E = calc_E(y, L1, L2, m1, m2, g)
    th1 = y(1); th1d = y(2);
    th2 = y(3); th2d = y(4);

    V = -(m1+m2)*L1*g*cos(th1) - m2*L2*g*cos(th2);
    T = 0.5*m1*(L1*th1d)^2 + 0.5*m2*((L1*th1d)^2 + (L2*th2d)^2 + ...
        2*L1*L2*th1d*th2d*cos(th1-th2));
    E = T + V;
end
