function [tau, phi, ang_alpha] = alg2_heur1(cfg)
%ALG2_HEUR1  算法 2：行为 1 的启发式解（闭式，不迭代）。
%
%   [tau, phi, ang_alpha] = alg2_heur1(cfg)
%
%   论文原文（Algorithm 2, Heuristic Solution for Behavior 1）：
%       tau_n = sum_{m in M_n} m*[sin(theta0-Dtheta/2)*fmin
%                                - sin(theta0+Dtheta/2)*fmax] / (2*W*f0*|M_n|)
%       tau <- tau - mean(tau);  tau <- clip(tau, +-kappa/(2W))
%       phi_m = pi*(m-1)*sin(theta0) + 2*pi*f0*tau_n,  wrapped to (-pi,pi]
%       alpha_k 由式 (11) 给出
%
%   物理含义：行为 1 的期望相位 psi_m(f) 在带内是直线，其斜率正比于
%   (m-1)*[sin(theta0-Dtheta/2)*fmin - sin(theta0+Dtheta/2)*fmax]，除以
%   2*W*f0 即得该天线所需的等效时延。同一 TTD 下的天线取斜率均值。
%
%   注意两处记号：
%     * 分子对 m 求和后除以 |M_n|，等价于取组内 m 的均值；
%     * 论文写 m 从 1 到 M，因此 (m-1) 出现在 phi 里而 m 出现在 tau 里，
%       这是原文写法，照抄（不统一成 m-1，否则 Fobj 会偏离图 5(a)）。
%
%   输出 tau 的裁剪区间与式 (4) 一致（对称），见 solve_tau_ls.m 的说明。

M = cfg.M;
N = cfg.N;

[~, bbar] = gen_beam(cfg, 1);

% fmax = f_[(K-1)/2]、fmin = f_[(1-K)/2] 的论文写法，按实际带边取
% （与附录 A 的 kappa_req 推导一致：带边 95 / 104.995 GHz）。
fmax = cfg.f(end);
fmin = cfg.f(1);

th_lo = cfg.theta0 - cfg.dtheta/2;
th_hi = cfg.theta0 + cfg.dtheta/2;
c     = (sin(th_lo)*fmin - sin(th_hi)*fmax)/(2*cfg.W*cfg.f0);

tau = zeros(N,1);
for n = 1:N
    rows   = find(cfg.gidx == n);      % 1-based 天线号，即论文的 m
    tau(n) = sum(rows)*c/numel(rows);  % sum_{m in M_n} m / |M_n| * c
end

% 居中 + 裁剪（顺序与算法 2 一致）
tau = tau - mean(tau);
half = cfg.kappa/(2*cfg.W);
tau  = max(min(tau, half), -half);

% 移相器：对齐 k = 0 处的波束中心
phi = zeros(M,1);
for n = 1:N
    rows     = find(cfg.gidx == n);
    ph       = pi*(rows-1)*sin(cfg.theta0) + 2*pi*cfg.f0*tau(n);
    phi(rows) = mod(ph + pi, 2*pi) - pi;
end

ang_alpha = update_alpha(cfg, bbar, tau, phi);
end
