function phi = solve_phi(cfg, bbar, ang_alpha, tau)
%SOLVE_PHI  式(8b)：给定 tau 与 digital 相位后的移相器相位。
%
%   phi_m = angle( sum_k omega_k * conj(e^{j*ang_alpha_k}) * [bbar_k]_m
%                  * e^{j*2*pi*f_k*tau_n} ),   m in M_n
%
%   逐 TTD 计算：同组天线共用 tau_n，但 bbar 逐天线不同。

C = bbar.*(cfg.omega.*exp(-1j*ang_alpha)).';   % M x K：bbar_k,m * w_k*e^{-j*ang_alpha}
phi = zeros(cfg.M,1);

for n = 1:cfg.N
    rows = find(cfg.gidx == n);
    phi(rows) = angle(C(rows,:)*exp(1j*2*pi*cfg.f*tau(n)));
end
end
