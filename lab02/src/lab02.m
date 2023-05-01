% Read
file = fopen('selection.txt', 'r');
X = fscanf(file, '%f,');
fclose(file);
gamma = 0.9;

% 1-2
[muhat, muci] = my_normfit_mu(X, 1 - gamma);
[s2hat, s2ci] = my_normfit_s2(X, 1 - gamma);

% 3
process_mu(X, gamma, muhat);
process_s2(X, gamma, s2hat);


function [muhat, muci] = normfit_mu(X, alpha)
    [muhat, ~, muci, ~] = normfit(X, alpha);
end

function [s2hat, s2ci] = normfit_s2(X, alpha)
    [~, sigmahat, ~, sigmaci] = normfit(X, alpha);
    s2hat = sigmahat ^ 2;
    s2ci = sigmaci .^ 2;
end

function [muhat, muci] = my_normfit_mu(X, alpha)
    muhat = mean(X);
    s = std(X);
    gamma = 1 - alpha;
    n = length(X);
    mu_bottom = muhat + s * tinv((1 - gamma) / 2, n - 1) / sqrt(n);
    mu_top = muhat + s * tinv((1 + gamma) / 2, n - 1) / sqrt(n);
    muci = [mu_bottom, mu_top];
end

function [s2hat, s2ci] = my_normfit_s2(X, alpha)
    s2hat = var(X);
    gamma = 1 - alpha;
    n = length(X);
    s2_top = (n - 1) * s2hat / chi2inv((1 - gamma) / 2, n - 1);
    s2_bottom = (n - 1) * s2hat / chi2inv((1 + gamma) / 2, n - 1);
    s2ci = [s2_bottom, s2_top];
end

function process_parameter(X, gamma, est, fit, line_legend, est_legend, top_legend, bottom_legend)
    N = length(X);
    figure;
    hold on;
    grid on;
    plot([1, N], [est, est]);
    ests = [];
    cis_bottom = [];
    cis_top = [];
    for n = 1:N
        [est, cis] = fit(X(1:n), 1 - gamma);
        ests = [ests, est];
        cis_bottom = [cis_bottom, cis(1)];
        cis_top = [cis_top, cis(2)];
    end
    plot(20:N, ests(20:N));
    plot(20:N, cis_bottom(20:N));
    plot(20:N, cis_top(20:N));
    l = legend(line_legend, est_legend, top_legend, bottom_legend);
    set(l, 'Interpreter', 'latex', 'fontsize', 18);
    hold off;
end

function process_mu(X, gamma, muhat)
    process_parameter(X, gamma, muhat, @my_normfit_mu, '$\hat\mu(\vec x_N)$', '$\hat\mu(\vec x_n)$', '$\underline\mu(\vec x_n)$', '$\overline\mu(\vec x_n)$');
end

function process_s2(X, gamma, S2)
    process_parameter(X, gamma, S2, @my_normfit_s2, '$\hat\sigma^2(\vec x_N)$', '$\hat\sigma^2(\vec x_n)$', '$\underline\sigma^2(\vec x_n)$', '$\overline\sigma^2(\vec x_n)$');
end