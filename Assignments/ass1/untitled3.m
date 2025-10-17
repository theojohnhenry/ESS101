clc,clear
g=9.82 

syms L m1 m2 z real
syms x1 x2 y1 y2 z1 z2 real
syms x1dot y1dot z1dot x2dot y2dot z2dot real

p1 = [x1;y1;z1]
p2 = [x2;y2;z2]

e = p1 - p2

q = [p1;p2]
qdot = [x1dot; y1dot; z1dot; x2dot; y2dot; z2dot];

J1 = jacobian(p1, q)
J2 = jacobian(p2, q)

W = m1 * J1' * J1 + m2 * J2' * J2 

T = (1/2) * qdot' * W * qdot
V = m1*g*[0 0 1]*p1 + m2*g*[0 0 1]*p2
C = (1/2) * (e'* e - L^2)


Lag = T - V - z'*C