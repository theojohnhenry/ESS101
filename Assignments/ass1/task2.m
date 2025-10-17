clc, clear 

syms p1 p2 [3 1] real
syms p1prim p2prim [3 1] real
syms p1bis p2bis [3 1] real
syms m1 m2 z g L ux uy uz real

Q = [ux;uy;uz;0;0;0]
q = [p1; p2]
qprim = [p1prim;p2prim] 
qbis = [p1bis;p2bis]
% leftside
diagElems = [m1,m1,m1,m2,m2,m2]
M = diag(diagElems)

e = p1 - p2;
C = (1/2)*(e.'*e - L^2);
gradc = jacobian(C, [p1; p2])'
%rightside
dp1dq = jacobian(p1,q)
dp2dq = jacobian(p2,q)

W = m1 * dp1dq' * dp1dq + m2 * dp2dq' * dp2dq

T = (1/2) * qprim' * W * qprim
V = m1*g*[0 0 1]*p1 + m2*g*[0 0 1]*p2

Lag = T - V

H1 = Q - jacobian((M*qprim),qprim)*qprim + jacobian(T, q)' + jacobian(V,q)'
h2 = jacobian(C,q)*qprim
H2 = -jacobian(h2,q)*qprim

finalA = [M gradc
        gradc' 0] 

finalB = [qbis; z]
 
finalC = [H1;H2]

finalB = finalA\finalC