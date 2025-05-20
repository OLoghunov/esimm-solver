sigma = 10;
rho = 28;
beta = 8/3;
params = [sigma rho beta];

x = 1;
y = 1;
z = 1;
init_conditions = [x, y, z];

t_start = 0.1;
t_step = 1e-2;
t_end = 50;
t_span = t_start:t_step:t_end;

ode = @(x)lorenz(x);
boot = @(ode, y, h)rk4(ode, y, h);

sol_num = zeros(length(init_conditions), length(t_span));
j = 1;
Y = init_conditions';

for t = t_span
    Y = boot(ode, Y, t_step);
    sol_num(:,j) = Y;
    j = j + 1;
end

plot3(sol_num(1,:), sol_num(2,:), sol_num(3,:))

function dXdt = lorenz(X)
    x = X(1);
    y = X(2);
    z = X(3);

    sigma = 10;
    rho = 28;
    beta = 8/3;

    dxdt = sigma * (y - x);
    dydt = x * (rho - z) - y;
    dzdt = x * y - beta * z;

    dXdt = [dxdt; dydt; dzdt];
end

function res = rk4(ode, y, h)
    k1 = h * ode(y);
    k2 = h * ode(y + k1/2);
    k3 = h * ode(y + k2/2);
    k4 = h * ode(y + k3);
    res = y + (k1 + 2*k2 + 2*k3 + k4) / 6;
end







