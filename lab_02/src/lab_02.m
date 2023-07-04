clear all;
X = [1.52,1.26,2.17,1.75,-0.19,2.24,2.76,1.52,1.89,3.10,2.61,1.18,1.83,...
    1.85,3.39,2.31,2.99,1.61,2.57,1.81,1.73,1.89,-0.00,2.27,1.61,2.57,...
    2.54,1.67,1.49,0.12,-0.04,1.36,2.04,2.04,-0.05,0.67,1.32,0.78,0.89,...
    2.73,1.51,1.48,1.67,2.18,1.70,4.20,1.81,2.66,1.72,0.77,3.16,1.86,...
    3.66,4.30,0.98,3.00,0.99,1.72,2.71,2.47,2.56,1.99,0.23,0.66,...
    2.47,2.71,2.28,2.59,3.30,2.08,0.90,0.49,2.38,0.71,0.10,1.50,...
    0.21,0.44,3.94,1.50,1.70,-0.73,1.76,2.71,1.95,-0.71,1.32,...
    3.95,2.64,-0.04,3.24,1.67,2.31,0.18,0.79,3.26,3.44,2.64,0.89,...
    2.47,4.02,2.12,0.61,2.59,1.44,1.82,2.94,3.03,1.97,2.30,0.80,0.52,...
    1.21,2.13,2.82,1.56,2.84,3.54,0.86,0.42]; 
gamma = 0.9;
% 1-2
[muhat, muci] = my_normfit_mu(X, 1 - gamma);
[s2hat, s2ci] = my_normfit_s2(X, 1 - gamma);

% 3
fprintf("n = %f\n", length(X));
fprintf("mu = %f\n", muhat);
fprintf("S^2 = %f\n", s2hat);
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