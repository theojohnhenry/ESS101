clc, clear

g = 9.82

syms theta phi real
syms x y z real
syms xdot ydot zdot real
syms thetadot phidot real
syms L m1 m2 real
syms xbis ybis zbis thetabis phibis real
syms ux uy uz real



p1 = [x; y; z]
p3 = L * [sin(phi)*cos(theta);
    sin(phi)*sin(theta);
    -cos(phi)]

p2 = p1 + p3
q = [p1; theta; phi]


J1 = jacobian(p1,q) 
J2 = jacobian(p2,q) 

qdot = [xdot; ydot; zdot; thetadot; phidot] %hastighet
qbis = [xbis; ybis; zbis; thetabis; phibis]

W = m1 * J1' * J1 + m2 * J2' * J2 
Wq = W * qdot
coriolis = jacobian(Wq, q)

T = (1/2) * qdot' * W * qdot
V = m1*g*[0 0 1]*p1 + m2*g*[0 0 1]*p2
Lag = T - V

nablaqT = jacobian(T,q)
nablaqL = jacobian(V,q)

dLagdq = jacobian(Lag, q)'
dLagdqdot = jacobian(Lag, qdot)'

dLagdqdot_dt = jacobian(dLagdqdot,[q; qdot]) * [qdot; qbis]

EulaG = dLagdqdot_dt - dLagdq

U = [ux; uy; uz]

Q = J1.' * U

dVdq = jacobian(V, q)

M = W
b = Q + nablaqT - nablaqL - coriolis

%We now have the expression Mq'' = b(q,q',q'')
