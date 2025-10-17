clc,clf,clear
%% butcher tableaus
rk1_a = zeros(1);
rk1_b = [1];
rk1_c = [0];

rk2_a = zeros(2); rk2_a(2,1) = 1/2;
rk2_b = [0 1];
rk2_c = [0; 1/2];

rk4_a = zeros(4); rk4_a(2,1) = 1/2; rk4_a(3,2) = 1/2; rk4_a(4,3) = 1;
rk4_b = [1/6 1/3 1/3 1/6];
rk4_c = [0;1/2;1/2;1];

%%
syms x

lambda = -2;
tf = 2;
dt = 10^-1;
x0 = 1;

f=lambda*x; %xdot, system dynamics
mf_f = matlabFunction(f);

N = round(tf/dt); %antal steg


xlin = 0:dt:tf-dt;
%scatter(xlin,mf_f(xlin),'blue');
hold on;

x_rk1 = [x0];
x_rk2 = [x0];
x_rk4 = [x0];


for k=1:N-1
   xnext1 = rk_step(mf_f, x_rk1(k),dt,rk1_a, rk1_b, rk1_c);
   x_rk1 = [x_rk1, xnext1];
   
   xnext2 = rk_step(mf_f, x_rk2(k),dt,rk2_a, rk2_b, rk2_c);
   x_rk2 = [x_rk2, xnext2];

   xnext4 = rk_step(mf_f, x_rk4(k),dt,rk4_a, rk4_b, rk4_c);
   x_rk4 = [x_rk4, xnext4];

end

scatter(xlin, x_rk1,'red')
hold on;
scatter(xlin, x_rk2,'green')
hold on;
scatter(xlin, x_rk4,'blue')
hold on;

plot(xlin, exp(lambda*xlin), 'yellow')
