function X = jpta_forward(cfg, tau, phi)
%JPTA_FORWARD  式(1)：每个子载波上的单位范数模拟波束 T*P*d_k，M x K。
%
%   x_k = (1/sqrt(M)) * diag(exp(j*phi)) * P * [exp(-j*2*pi*f_k*tau_n)] * alpha_k*s_k
%
%   T = (1/sqrt(M))*diag(exp(j*phi))，d_k = exp(-j*2*pi*f_k*tau)。
%   这样 ||T*P*d_k|| = ||P*d_k||/sqrt(M) = sqrt(M)/sqrt(M) = 1，
%   与论文"||TPd_k||=1"一致。
%
%   这是整个复现里唯一的"物理"实现，所有算法与所有图都调它。

D = exp(-1j*2*pi*(tau(:))*(cfg.f(:).'));   % N x K，TTD 的频响
X = (cfg.P*D).*exp(1j*phi(:));             % M x K，移相器
X = X/sqrt(cfg.M);                         % 式(1) 的 1/sqrt(M)
end
