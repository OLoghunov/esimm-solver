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
multistep = @(ode, y, y_1, y_2, y_3, h)ab4(ode, y, y_1, y_2, y_3, h);

sol_num = zeros(length(init_conditions), length(t_span));
j = 1;
Y = init_conditions';

% for t = t_span
%     Y = boot(ode, Y, t_step);
%     sol_num(:,j) = Y;
%     j = j + 1;
% end

% multistep boot
sol_num(:,1) = init_conditions';
sol_num(:,2) = boot(ode, Y, t_step);
sol_num(:,3) = boot(ode, y_1, t_step);
sol_num(:,4) = boot(ode, y_3, t_step);
j = 5;

% multistep calcs
for t = t_span
    Y = multistep(ode, Y, sol_num(:,j-2), sol_num(:,j-3), sol_num(:,j-4), t_step);
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

function res = rk2(ode, y, h)
    k1 = ode(y);
    k2 = ode(y + k1*h/2);
    res = y + h*k2;       
end

function res = rk4(ode, y, h)
    k1 = h * ode(y);
    k2 = h * ode(y + k1/2);
    k3 = h * ode(y + k2/2);
    k4 = h * ode(y + k3);
    res = y + (k1 + 2*k2 + 2*k3 + k4) / 6;
end

function res = ab4(ode, y, y_1, y_2, y_3, h)

    k1 = ode(y);
    k2 = ode(y_1);
    k3 = ode(y_2);
    k4 = ode(y_3);

    res = y + h * (55/24*k1 - 59/24*k2 + 37/24*k3 - 9/24*k4);
end







