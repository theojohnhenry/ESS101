function x_next = rk_step(f,x,dt,a,b,c)
    %tar in dynamics
    %nuvarnade state
    %timestep
    %a b c konstanter i butcher
    K = []
    s = length(b);
    %make k vector
    for i=1:1:s
        sum = 0;
        for j=1:1:s
            sum = sum + a(i,j)*K(j);
        end
        Ki = f(x+dt*sum);
        K = [K, Ki];
    end

    jonny=0;
    for i=1:1:s
        jonny=jonny+b(i)*K(i);
    end

    x_next = x+dt*jonny;
end