function [tau, phi, ang_alpha, hist, X] = alg1_jpta(cfg, bbar, method)
%ALG1_JPTA  算法 1：tau / phi / angle(alpha) 交替优化的迭代设计。
%
%   [tau, phi, ang_alpha, hist] = alg1_jpta(cfg, bbar, method)
%
%   method : 'ls'     -> 式(8a) 线搜索（FFT 网格，精确）
%            'lemma1' -> Lemma 1（式(9)(10) 的 wLS 近似，快）
%
%   每次迭代：
%     1) 固定 angle(alpha)，按 (8a) 或 Lemma 1 求 tau_n，再按 (8b) 求 phi_m
%     2) 把 tau 整体平移到可行区间中央（tau_bar），并同步 angle(alpha)
%     3) 固定 tau、phi，按 (11) 重算 angle(alpha)
%   第 2 步的 tau_bar 公式带 min/max 包裹，不能简化成"减去均值"，
%   否则 kappa 小时 Fobj 会偏低；居中后必须同步 angle(alpha) += 2*pi*f_k*tau_bar，
%   否则整体相位错位、Fobj 直接掉。

if nargin < 3 || isempty(method)
    method = cfg.method;
end

M = cfg.M; N = cfg.N;
ang_alpha = zeros(cfg.K,1);          % 算法 1 的初值
tau = zeros(N,1);
phi = zeros(M,1);
hist = zeros(cfg.maxiter,1);

for it = 1:cfg.maxiter
    % ---- 1) 固定 angle(alpha)，求 tau 与 phi --------------------------
    switch method
        case 'ls'
            tau = solve_tau_ls(cfg, bbar, ang_alpha);
        case 'lemma1'
            tau = solve_tau_wls(cfg, bbar, ang_alpha);
        otherwise
            error('alg1_jpta: 未知 method=%s', method);
    end
    phi = solve_phi(cfg, bbar, ang_alpha, tau);

    % ---- 2) 把 tau 推到可行区间中央 ----------------------------------
    tau_min = min(tau);
    tau_max = max(tau);
    half    = cfg.kappa/(2*cfg.W);
    tau_bar = min(mean(tau), half + tau_min);
    tau_bar = max(tau_bar, tau_max - half);
    tau = tau - tau_bar;
    ang_alpha = ang_alpha + 2*pi*cfg.f*tau_bar;

    % ---- 3) 固定 tau、phi，重算 angle(alpha) --------------------------
    ang_alpha = update_alpha(cfg, bbar, tau, phi);

    X = jpta_forward(cfg, tau, phi);
    hist(it) = fobj(cfg, bbar, X);
end
end
