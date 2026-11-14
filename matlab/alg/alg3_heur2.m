function [tau, phi, ang_alpha] = alg3_heur2(cfg)
%ALG3_HEUR2  算法 3：行为 2 的启发式解（闭式，不迭代）。
%
%   [tau, phi, ang_alpha] = alg3_heur2(cfg)
%
%   论文原文（Algorithm 3, Heuristic Solution for Behavior 2）：
%       [eb]_m = ( e^{j*pi*m*sin(theta2)} + e^{j*pi*m*sin(theta1)} ) / sqrt(2M)
%       tau_n  = -3/(2*pi*W) * angle( sum_{m in M_n} [eb]_m^*
%                            * e^{j*pi*m*sin(theta2)*f_[K/3]/f0} )
%
%   ⚠ 伪代码里那个弯钩是"取相位"算子 angle(.)，不是乘号也不是除法。
%   漏掉它会让 tau 变成复数/被全部裁到边界（实测 Fobj 从 0.56 掉到 0.40）。
%   带上 angle 后 |tau_n| <= 3/(2W)，span <= 3/W，正好对上 Remark 2 的
%   "behavior 2 最大时延 3/W"，也能复现图 6(b) 的 kappa~2 饱和。
%       tau <- tau - mean(tau);  tau <- clip(tau, +-kappa/(2W))
%       phi_m = angle([eb]_m) + 2*pi*f0*tau_n,  wrapped to (-pi,pi]
%       alpha_k 由式 (11) 给出
%
%   两处按勘误处理：
%     1) [eb]_m 的两项，原文都印成 sin(theta2)，第二项应为 sin(theta1)
%        （见已发布 Errata，/errata/#J17）。否则 eb 退化为单波束，
%        行为 2 的阶跃结构丢掉。
%     2) 系数 -3/(2*pi*W)：3 对应"斜率跨 3/W"，与 Remark 2 的
%        "behavior 2 最大时延 3/W" 自洽；正文那句 kappa = 2 是同类笔误。
%
%   记号：eb 用 m（不是 m-1），phi 也用 m，均照原文。若改成 (m-1)，
%   kappa=6.4 时的 Fobj 从 ~0.56 掉到 0.2 以下，可用来判别。

M = cfg.M;
N = cfg.N;

[~, bbar] = gen_beam(cfg, 2);

m  = (1:M).';
eb = (exp(1j*pi*m*sin(cfg.theta2)) + exp(1j*pi*m*sin(cfg.theta1)))/sqrt(2*M);

% f_[K/3]：k = floor(K/3)
kq  = floor(cfg.K/3);
fq  = cfg.f0 + cfg.W*kq/cfg.K;

tau = zeros(N,1);
for n = 1:N
    rows   = find(cfg.gidx == n);
    tau(n) = -3/(2*pi*cfg.W) * angle(sum(conj(eb(rows)) ...
                 .* exp(1j*pi*rows*sin(cfg.theta2)*fq/cfg.f0)));
end

tau  = tau - mean(tau);
half = cfg.kappa/(2*cfg.W);
tau  = max(min(tau, half), -half);

phi = zeros(M,1);
for n = 1:N
    rows      = find(cfg.gidx == n);
    ph        = angle(eb(rows)) + 2*pi*cfg.f0*tau(n);
    phi(rows) = mod(ph + pi, 2*pi) - pi;
end

ang_alpha = update_alpha(cfg, bbar, tau, phi);
end
