clc, clear

g = 9.82

syms theta phi real
syms x y z real
syms xdot ydot zdot 
syms thetadot phidot
syms L m1 m2 real
syms xbis ybis zbis thetabis phibis


p1 = [x; y; z]
p3 = L * [sin(phi)*cos(theta);
    sin(phi)*sin(theta);
    -cos(phi)]

p2 = p1 + p3
q = [p1; theta; phi]


J = jacobian(p2,q) % dp2/dq

qdot = [xdot; ydot; zdot; thetadot; phidot] %hastighet

pdot = J * qdot

T2 = (1/2) * m2 * pdot' * pdot
V2 = m2*g*[0 0 1]*p2

W = m1 * J' * J + m2 * J' * J

