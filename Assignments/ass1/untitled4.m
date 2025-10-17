syms L m1 m2 z real
syms p1 p2 [3 1] real  % components: p11, p12, p13, etc.
syms qdot [6 1] real

e = p1 - p2;

% Define T, V, C directly without jacobian, since you want compact form
T = (1/2) * (m1*(qdot(1:3).' * qdot(1:3)) + m2*(qdot(4:6).' * qdot(4:6)));
V = m1*g*[0 0 1]*p1 + m2*g*[0 0 1]*p2;
C = (1/2)*(e.'*e - L^2);

Lag = T - V - z*C