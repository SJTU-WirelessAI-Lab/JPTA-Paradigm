function tau = solve_tau_wls(cfg, bbar, ang_alpha)
%SOLVE_TAU_WLS  Lemma 1（式(9)(10)）：相位解缠 + 2 参数加权最小二乘 + 取模 + 裁剪。
%
%   对每根天线 m：
%     1) 沿 k 递增做相位解缠 U( angle([bbar_k]_m) - ang_alpha_k )，得 y_k = -U(...)
%     2) 以 2*pi*f_k*tau - phi_m = y_k 为模型做加权 (w_k*|[bbar_k]_m|) 两参数 LS
%     3) 式(10a) 取模到 [-K/(2W), K/(2W))，再式(10b) 裁剪到 [-kappa/(2W), kappa/(2W)]
%   顺序不能反：先取模再裁剪。
%
%   符号说明：式(9) 方括号内是 [2*pi*f_k*tau - phi_m + U(angle(bbar)-ang_alpha)]，
%   令 y_k = -U(...) 即得上述模型；这与"bbar_k'TPd_k 应等于 e^{-j*ang_alpha_k}"
%   的推导一致，不要"顺手化简"成 +ang_alpha，否则 Fobj 会掉到 0.3 以下。

f  = cfg.f;
w  = cfg.omega;
% 频率去中心化 + 归一化。2*pi*f_k*tau - phi_m 中斜率与截距在 95-105 GHz 上
% 高度相关（f 的相对变化只有 10%），而 2*pi*f 的绝对值又高达 6.6e11，
% 直接解正规方程会病态到 RCOND ~ 1e-21。令
%     x_k = 2*pi*(f_k - mean(f))/s,  s = max|2*pi*(f_k - mean(f))|
% 则截距自动吸收 2*pi*mean(f)*tau，tau 的估计不变，条件数回到 O(1)。
Xc = 2*pi*(f - mean(f));
s  = max(abs(Xc));
Xn = Xc/s;
tau = zeros(cfg.N,1);

for n = 1:cfg.N
    rows = find(cfg.gidx == n);
    a1 = 0; a2 = 0; a3 = 0; a4 = 0; a5 = 0;
    for m = rows(:).'
        y  = -unwrap(angle(bbar(m,:).') - ang_alpha);     % K x 1
        wk = w.*abs(bbar(m,:).');                         % K x 1 加权
        a1 = a1 + sum(wk);
        a2 = a2 + sum(wk.*Xn);
        a3 = a3 + sum(wk.*Xn.^2);
        a4 = a4 + sum(wk.*y);
        a5 = a5 + sum(wk.*Xn.*y);
    end
    % 正规方程：[a3 -a2; -a2 a1] * [tau*s; phi] = [a5; -a4]
    sol    = [a3, -a2; -a2, a1] \ [a5; -a4];
    tau(n) = sol(1)/s;
end

% 式(10a)：取模到 [-K/(2W), K/(2W))
tau = mod(tau + cfg.K/(2*cfg.W), cfg.K/cfg.W) - cfg.K/(2*cfg.W);
% 式(10b)：裁剪到可行区间
tau = max(min(tau, cfg.kappa/(2*cfg.W)), -cfg.kappa/(2*cfg.W));
end
