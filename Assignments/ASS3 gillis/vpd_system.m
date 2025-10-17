function vdp_dx = vpd_system(t,x,u)
vdp_dx = [x(2); 
    u*(1 - x(1)^2) * x(2) - x(1)];