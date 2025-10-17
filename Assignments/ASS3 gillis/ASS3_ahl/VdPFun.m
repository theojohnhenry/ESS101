function F = VdPFun(u_val, t, x)

F = [x(2)
     u_val*(1-x(1)^2)*x(2)-x(1)];
end