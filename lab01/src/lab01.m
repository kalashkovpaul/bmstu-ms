%% Read
file = fopen('selection.txt', 'r');
X = fscanf(file, '%f,');
fclose(file);
%% Calculation
M_max = max(X);
M_min = min(X);
R = M_max - M_min;
n = length(X);
mu = sum(X) / n;
Ssquare = sum((X - mu).^2) / (n - 1);
m = floor(log2(n)) + 2;
step = R / m;
sorted_X = sort(X);
intervals = cell(1, m);
i = 1;
for cur = (M_min + step):step:M_max
    last = cur - step;
    intervals(i) = {X((last < X) & (X < cur))};
    i = i + 1;
end
intervals{m} = [intervals{m}; X(X == M_max)];
density = 1:m;
for i = 1:m
    density(i) = length(intervals{i}) / n / step;
end
%% probability density function
figure('Name', 'Функция плотности распределения вероятностей');
histogram(X, m, 'BinEdges', M_min:step:M_max, 'Normalization', 'pdf');
hold on;
x = (M_min - step):0.1:(M_max+step);
f = pdf('Normal', x, mu, Ssquare);
plot(x,f,'LineWidth',1.5, 'color', 'red')
hold off;
%% probability function
figure('Name', 'Функция распределения вероятностей')
x = (M_min - step):0.1:(M_max+step);
[F_empirical, x_empirical] = ecdf(X);
F = normcdf((x - mu) / sqrt(Ssquare));
hold on;
plot(x_empirical, F_empirical,'LineWidth', 1.5);
plot(x, F,'LineWidth', 1.5, 'color', 'red');
hold off;