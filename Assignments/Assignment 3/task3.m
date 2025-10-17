clc, clear, clf
%our system

u = 5;
tf = 25;
x0 = 1;
y0 = 0;
%----
deltat = 0.1; % STEP LENGTH

% xprim = f(t,X)
f = @(t,X)[X(2);u*(1-X(1)^2)*X(2)-X(1)];

%ode of vanderpol

options = odeset();
[ode_t, ode_X] = ode45(f, [0 tf], [x0;y0], options);
scatter(ode_t, ode_X(:,1))
xlabel('t')
ylabel('x ') 

%i observe a cyclical pattern ?

%our ode for rk4

Xk = [x0;y0]; %seeding with start values
i=1;
rk_t = zeros(tf/)
rk_X = zeros(0:deltat:tf)

for k=0:deltat:tf

    k1 = f(k,Xk);
    k2 = f(k + (deltat/2), Xk + (deltat/2) * k1);
    k3 = f(k + (deltat/2), Xk + (deltat/2) * k2);
    k4 = f(k + deltat, Xk + deltat*k3);

    Xnext = Xk + deltat * ( ...
        (1/6 * k1) ...
        + (1/3 * k2) ...
        + (1/3 * k3) ...
        + (1/6 * k4) ...
    );
    Xk = Xnext;

    % store it 
    i=i+1;

end
