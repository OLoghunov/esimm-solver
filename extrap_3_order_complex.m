% Конфиг
X0 = [1.1; 1; 1];

k11 = 8/7;
k12 = -1/7;

alpha1 = 8/7 + 1.075i;
alpha2 = -1/7 - 12.075i;

Tmax = 30;

h_span = logspace(-4, -0.5, 20);

j = 1;
E = zeros(1, length(h_span));

elapsed_extr_3 = [];
max_error_extr_3 = [];

elapsed_extr_3_complex = [];
max_error_extr_3_complex = [];

elapsed_extr_4 = [];
max_error_extr_4 = [];

elapsed_cd = [];
max_error_cd = [];

% Цикл по шагам
for h = h_span
    
    t = 0:h:Tmax;
    x_ord_3 = zeros(length(X0), length(t));
    x2 = zeros(length(X0), 2 * length(t) - 1);

    x_ord_3(:, 1) = X0;
    x2(:, 1) = X0;

    x_ref = x_ord_3;
    x_cd = x_ord_3;
    x_ord_4 = x_ord_3;
    x_ord_3_complex = x_ord_3;

    tic;
    % Разгон методом более высокого порядка (RK7) для 3 порядка
    for i = 2:2
        x_ord_3(:, i) = rk7(x_ord_3(:, i - 1), h);
    end
    
    % Экстраполяция с h
    for i = 3:length(t)
        T11 = cd(x_ord_3(:, i-1), h);
        T12 = cd(x_ord_3(:, i-2), 2 * h);

        T21 = k11 * T11 + k12 * T12;
        x_ord_3(:, i) = T21;
    end

    elapsed_extr_3 = [elapsed_extr_3; toc];

    tic;
    % Разгон методом более высокого порядка (RK7) для 4 порядка
    for i = 2:3
        x_ord_4(:, i) = rk7(x_ord_4(:, i - 1), h);
    end
    
    % Экстраполяция с h
    for i = 4:length(t)
        T11 = cd(x_ord_4(:, i-1), h);
        T12 = cd(x_ord_4(:, i-2), 2 * h);
        T13 = cd(x_ord_4(:, i-3), 3 * h);

        T21 = 8/7 * T11 - 1/7 * T12;
        T22 = 27/26 * T11 - 1/26 * T13;

        T31 = 189/85 * T21 - 104/85 * T22;
        x_ord_4(:, i) = T31;
    end

    elapsed_extr_4 = [elapsed_extr_4; toc];

    tic;
    % Разгон методом более высокого порядка (RK7) для 3 порядка с
    % комплексными коэффициентами
    for i = 2:2
        x_ord_3_complex(:, i) = rk7(x_ord_3_complex(:, i - 1), h);
    end
    
    % Экстраполяция с h
    for i = 3:length(t)
        T11 = cd(x_ord_3_complex(:, i-1), h);
        T12 = cd(x_ord_3_complex(:, i-2), 2 * h);

        T21 = alpha1 * T11 + alpha2 * T12;
        x_ord_3_complex(:, i) = real(T21);
    end

    elapsed_extr_3_complex = [elapsed_extr_3_complex; toc];

    for i = 2:4
        x2(:, i) = rk7(x2(:, i - 1), h / 2);
    end

    % Экстраполяция с h/2 для 3 порядка (для Order plot)
    for i = 5:2:2 * length(t) - 1
        T11 = cd(x2(:, i-1), h / 2);
        T12 = cd(x2(:, i-2), 2 * h / 2);

        T21 = k11 * T11 + k12 * T12;
        x2(:, i) = real(T21);

        T11 = cd(x2(:, i), h / 2);
        T12 = cd(x2(:, i-1), 2 * h / 2);

        T21 = k11 * T11 + k12 * T12;
        x2(:, i+1) = real(T21);
    end

    % Моделирование через RK7 для референса
    for i = 2:length(t)
        x_ref(:,i) = rk7(x_ref(:,i-1), h);
    end

    tic;
    for i = 2:length(t)
        x_cd(:,i) = cd(x_cd(:,i-1), h);
    end
    elapsed_cd = [elapsed_cd; toc];
    
    x2 = x2(:, 2 * (1:length(t)) - 1); % берем только нечетные точки

    % e(h)/e(h/2) для текущего h
    e1 = vecnorm(x_ord_3 - x_ref);
    e2 = vecnorm(x2 - x_ref);
    E(j) =  max(e1) / max(e2);
    j = j + 1;

    max_error_extr_3 = [max_error_extr_3; max(vecnorm(x_ord_3 - x_ref))];
    max_error_extr_4 = [max_error_extr_4; max(vecnorm(x_ord_4 - x_ref))];
    max_error_extr_3_complex = ...
        [max_error_extr_3_complex; max(vecnorm(x_ord_3_complex - x_ref))];
    max_error_cd = [max_error_cd; max(vecnorm(x_cd - x_ref))];
end

% График порядка
figure(2);
semilogx(h_span, log2(E));
xlabel('Step $h$', Interpreter='latex');
ylabel('Method order', Interpreter='latex');
grid on;
title('Order Plot');

% График эффективности
figure(3);
loglog(elapsed_cd, max_error_cd, 'o-', 'LineWidth', 1.5);
hold on
loglog(elapsed_extr_3, max_error_extr_3, 'o-', 'LineWidth', 1.5);
loglog(elapsed_extr_4, max_error_extr_4, 'o-', 'LineWidth', 1.5);
loglog(elapsed_extr_3_complex, max_error_extr_3_complex, 'o-', 'LineWidth', 1.5);
hold off
xlabel('Computational Time (s)', 'Interpreter', 'latex');
ylabel('Max Error', 'Interpreter', 'latex');
title('Efficiency Plot');
legend('cd', 'cd-extr-ord-3', 'cd-extr-ord-4', 'cd-extr-ord-3-complex')
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Метод CD с системой Нозе-Гувера
function res = cd(X, h)
    x = X(1) + h/2*X(2);
    y = X(2) + h/2*(-x + X(2)*X(3));
    z = X(3) + h*(1-y^2);
    
    y = (y - h/2*x)/(1-h/2*z);
    x = x + h/2*(y);

    res = [x; y; z];
end

% Метод RK7 (референсный)
function y = rk7(y, h)
    func = @noze_hoover;
    c_1_11 = 41.0 / 840.0;
    c6 = 34.0 / 105.0;
    c_7_8= 9.0 / 35.0;
    c_9_10 = 9.0 / 280.0;
    
    a2 = 2.0 / 27.0;
    
    b31 = 1.0 / 36.0;
    b32 = 3.0 / 36.0;
    b41 = 1.0 / 24.0;
    b43 = 3.0 / 24.0;
    b51 = 20.0 / 48.0;
    b53 = -75.0 / 48.0;
    b54 = 75.0 / 48.0;
    b61 = 1.0 / 20.0;
    b64 = 5.0 / 20.0;
    b65 = 4.0 / 20.0;
    b71 = -25.0 / 108.0;
    b74 =  125.0 / 108.0;
    b75 = -260.0 / 108.0;
    b76 =  250.0 / 108.0;
    b81 = 31.0/300.0;
    b85 = 61.0/225.0;
    b86 = -2.0/9.0;
    b87 = 13.0/900.0;
    b91 = 2.0;
    b94 = -53.0/6.0;
    b95 = 704.0 / 45.0;
    b96 = -107.0 / 9.0;
    b97 = 67.0 / 90.0;
    b98 = 3.0;
    b10_1 = -91.0 / 108.0;
    b10_4 = 23.0 / 108.0;
    b10_5 = -976.0 / 135.0;
    b10_6 = 311.0 / 54.0;
    b10_7 = -19.0 / 60.0;
    b10_8 = 17.0 / 6.0;
    b10_9 = -1.0 / 12.0;
    b11_1 = 2383.0 / 4100.0;
    b11_4 = -341.0 / 164.0;
    b11_5 = 4496.0 / 1025.0;
    b11_6 = -301.0 / 82.0;
    b11_7 = 2133.0 / 4100.0;
    b11_8 = 45.0 / 82.0;
    b11_9 = 45.0 / 164.0;
    b11_10 = 18.0 / 41.0;

    h2_7 = a2 * h;
    
    k1 = func(y);
    k2 = func(y + h2_7 * k1);
    k3 = func(y + h * (b31*k1 + b32*k2));
    k4 = func(y + h * (b41*k1 + b43*k3));
    k5 = func(y + h * (b51*k1 + b53*k3 + b54*k4));
    k6 = func(y + h * (b61*k1 + b64*k4 + b65*k5));
    k7 = func(y + h * (b71*k1 + b74*k4 + b75*k5 + b76*k6));
    k8 = func(y + h * (b81*k1 + b85*k5 + b86*k6 + b87*k7));
    k9 = func(y + h * (b91*k1 + b94*k4 + b95*k5 + b96*k6 + b97*k7 + b98*k8));
    k10 = func(y + h * (b10_1*k1 + b10_4*k4 + b10_5*k5 + b10_6*k6...
                                              + b10_7*k7 + b10_8*k8 + b10_9*k9));
    k11 = func(y + h * (b11_1*k1 + b11_4*k4 + b11_5*k5 + b11_6*k6...
                               + b11_7*k7 + b11_8*k8 + b11_9*k9 + b11_10 * k10));
    y = y +  h * (c_1_11 * (k1 + k11)  + c6 * k6 + c_7_8 * (k7 + k8)...
                                               + c_9_10 * (k9 + k10));
end

% Система Нозе-Гувера
function dx = noze_hoover(X)
    dx = X;
    dx(1) = X(2);
    dx(2) = -X(1) + X(2) * X(3);
    dx(3) = 1 - X(2)^2;
end
